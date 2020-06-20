using PyCall
using DataFrames
using JSONTables
using Plots
using DFTK
using JSON3


const DEBUG = false


debugsuffix() = DEBUG ? "_DEBUG" : ""

function all_elements()
    kgrids = open(JSON3.read, joinpath(@__DIR__, "kgrids.json"))
    string.(keys(kgrids))
end

function loadjson(file::AbstractString)
    open(file) do fp
        DataFrame(jsontable(read(fp, String)))
    end
end

function savejson(file::AbstractString, df::DataFrame)
    open(file, "w") do fp
        write(fp, arraytable(df))
    end
end

function run_dftk(ase_atoms)
    # Pretty much follows the values from SupplMat-Abinit-HGHsc.pdf

    # Try also core=:semicore
    getpsp(el) = load_psp(el.symbol, functional="pbe", core=:fullcore)

    magnetic_moments = load_magnetic_moments(ase_atoms)
    if iszero(sum(pair -> sum(pair[2]), magnetic_moments))
        empty!(magnetic_moments)
    end

    lattice = load_lattice(ase_atoms)
    atoms   = [ElementPsp(el.symbol, psp=getpsp(el)) => position
               for (el, position) in load_atoms(ase_atoms)]
    model   = model_PBE(lattice, atoms, smearing=Smearing.FermiDirac(),
                        temperature=0.01 * DFTK.units.eV,
                        magnetic_moments=magnetic_moments)
    model.spin_polarization != :collinear && empty!(magnetic_moments)

    if DEBUG
        basis = PlaneWaveBasis(model, 10, kgrid=[3, 3, 3])
    else
        symbol = atoms[1][1].symbol
        kgrid  = open(JSON3.read, joinpath(@__DIR__, "kgrids.json"))[symbol]
        basis  = PlaneWaveBasis(model, 250 * DFTK.units.Ry, kgrid=kgrid)
    end
    ρspin  = guess_spin_density(basis, magnetic_moments)
    scfres = self_consistent_field(basis, ρspin=ρspin, tol=1e-8, mixing=HybridMixing())

    scfres.energies.total / DFTK.units.eV
end


function compute_volumes(symbol)
    # Follows https://wiki.fysik.dtu.dk/ase/tutorials/deltacodesdft/deltacodesdft.html
    dcdft = pyimport("ase.collections").dcdft

    # Table of volume and energy per atom
    df = DataFrame(code=String[], symbol=String[], volume=Float64[], energy=Float64[])
    for cellfactor in (94:2:108) ./ 100
        println("#\n# Working on $symbol $(100cellfactor)%\n#")

        atoms = dcdft.__getitem__(symbol)
        atoms.set_cell(atoms.cell * cbrt(cellfactor), scale_atoms=true)
        push!(df, ("DFTK", symbol,
                   atoms.get_volume() / length(atoms),
                   run_dftk(atoms) / length(atoms)))

        println()
        println()
    end
    df
end

function fit_eos(df)
    # For each Symbol fit an eos and return a new DataFrame
    EquationOfState = pyimport("ase.eos").EquationOfState

    # See https://juliadata.github.io/DataFrames.jl/stable/man/split_apply_combine/
    combine(groupby(df, [:symbol, :code])) do gdf
        eos = EquationOfState(gdf.volume, gdf.energy, "birchmurnaghan")
        eos.fit(warn=true)
        e0, B, Bp, v0 = eos.eos_parameters
        (energy=e0, B=B, Bp=Bp, volume=v0)
    end
end


function reference_eos()
    dcdft = pyimport("ase.collections").dcdft
    kJ = pyimport("ase.units").kJ
    df = DataFrame(symbol=String[], code=String[], energy=Union{Float64,Missing}[],
                   B=Float64[], Bp=Float64[], volume=Float64[])
    # df uses the units of ASE:
    # energy  Total energy in eV
    # B       Bulk modulus in eV/Å^3
    # Bp      Pressure derivative of bulk modulus
    # volume  Volume per atom in Å^3

    for (symbol, data) in dcdft.data
        for code in ("wien2k", "exp")
            if "$(code)_B" in keys(data)
                # Note unit conversion: Bulk modulus in dcdft.data stored in GPa,
                # but we need eV / Å^3 for consistency with ASE
                push!(df, (symbol, "$(code)",
                           get(data, "$(code)_energy", missing),
                           data["$(code)_B"] * 1e-24 * kJ,
                           get(data, "$(code)_Bp",     missing),
                           get(data, "$(code)_volume", missing)))
            end
        end
    end

    df
end

function amend_deltas(df)
    function calc_delta(compare, symbol, volume, B, Bp)
        delta = pyimport("ase.utils.deltacodesdft").delta
        dfcompare = filter(:code => isequal(compare), df)
        compare = filter(:symbol => isequal(symbol), dfcompare)
        isempty(compare) && return missing
        ctuple = (only(compare.volume), only(compare.B), only(compare.Bp))
        1000 * delta(ctuple..., volume, B, Bp)
    end

    transform(df,
        [:symbol, :volume, :B, :Bp] => (
            (symbol, volume, B, Bp) -> calc_delta.("wien2k", symbol, volume, B, Bp)
        ) => "Δ wien2k",
        [:symbol, :volume, :B, :Bp] => (
            (symbol, volume, B, Bp) -> calc_delta.("exp", symbol, volume, B, Bp)
        ) => "Δ exp"
    )
end

function volumedata(elements=all_elements())
    datafile = "volumedata$(debugsuffix()).json"
    !isfile(datafile) && savejson(datafile, DataFrame())
    dfvol = loadjson(datafile)

    for symbol in elements
        !isempty(dfvol) && symbol in dfvol.symbol && continue
        res = compute_volumes(symbol)

        # Update dfvol and make a checkpoint save
        dfvol = loadjson(datafile)
        append!(dfvol, res)
        savejson(datafile, dfvol)
    end
    dfvol
end


function compute_deltas(elements=all_elements())
    dfeos = append!(reference_eos(), fit_eos(volumedata(elements)))
    dfdelta = amend_deltas(dfeos)
    savejson("delta$(debugsuffix()).json", dfdelta)
    dfdelta
end


function plot_case(df, symbol)
    birchmurnaghan = pyimport("ase.eos").birchmurnaghan

    df = filter(:symbol => isequal(symbol), df)
    vmin = minimum(df.volume)
    vmax = maximum(df.volume)
    v = collect((0.94 * vmin):0.05:(1.06 * vmax))

    p = plot()
    for (i, code) in enumerate(("wien2k", "DFTK", "exp"))
        dfcode = filter(:code => isequal(code), df)
        isempty(dfcode) && continue
        ctuple = (only(dfcode.volume), only(dfcode.B), only(dfcode.Bp))
        println(ctuple)
        plot!(p, v, birchmurnaghan(v, 0.0, ctuple...), label=code, color=i)
    end

    dfdftk = filter(:code => isequal("DFTK"), df)
    if !isempty(dfdftk)
        eDFTK = only(dfdftk.energy)
        computed = filter(:symbol => isequal(symbol), volumedata())
        println(eDFTK, " ", computed.energy .- eDFTK)
        scatter!(p, computed.volume, computed.energy .- eDFTK)
    end

    xlabel!(p, "volume [Å^3]")
    ylabel!(p, "energy [eV/atom]")
    p
end


function main()
    DEBUG && @warn "Running in debug mode. Results not converged."
    elements = isempty(ARGS) ? all_elements() : ARGS
    compute_deltas(elements)
end
(abspath(PROGRAM_FILE) == @__FILE__) && main()
