# Delta-Test for DFTK

This Julia script performs the famous [Delta test](https://molmod.ugent.be/deltacodesdft)
for the [density-functional toolkit (DFTK)](https://dftk.org).

## Running the code
Running the code requires an installation of
[Julia 1.4.0](https://julialang.org/downloads/#current_stable_release),
of [DFTK](https://docs.dftk.org/dev/guide/installation/)
and of [ASE](https://wiki.fysik.dtu.dk/ase/).
With this setup can generate the data for the deltatest by executing:
```bash
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
julia run_deltatest.jl  # Generate data
```
