#!/bin/env bash

mkdir -p jlexamples
./ipynb-to-jl.sh "full_fledged_schema_examples"
./ipynb-to-jl.sh "full_fledged_schema_examples_new"
rm jlexamples/full_fledged_schema_examples*checkpoint*
julia -p auto --project="." -e 'include("./run_notebooks.jl")'
rm -rf jlexamples/
