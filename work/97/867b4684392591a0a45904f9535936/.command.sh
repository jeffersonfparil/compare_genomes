#!/usr/bin/env bash
cd /data-weedomics-3
echo 'using Pkg; Pkg.add(["Plots", "DataFrames", "CSV", "ProgressMeter", "JLD2"])' > install_julia_packages.jl
julia install_julia_packages.jl
rm install_julia_packages.jl
