#!/usr/bin/env bash

mkdir -p jlexamples
./ipynb-to-jl.sh
# rm -f jlexamples/*checkpoint*
julia -p auto --project="." -e 'include("./run_notebooks.jl")'
rm -rf jlexamples/
