# glowify_func
Generate interlacing cuts in a grid that span a surface described by a function. 

## Requirements:

Currently, this code only runs on [julia v0.6](https://julialang.org/downloads/oldreleases.html). It won't run properly on v0.7+ for now. 

Once you have julia v0.6, run
```
julia> Pkg.clone(" https://github.com/prcamp/PointsNCurves.git")
```
from the julia command prompt (alternatively you can clone the package directly into your julia package directory).

To open the jupyter notebook, install `IJulia` using 
```
julia> Pkg.add("IJulia")
```
after which you should be able to launch a jupyter session with
```
julia> using IJulia
julia> notebook()
```

See the jupyter notebook for further details on using these functions.
