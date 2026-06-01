#!/bin/env bash

subdir=$1
find ./examples/$subdir -type f -name "*.ipynb" | while read f; do
    fname=$(basename "$f")
    dname=$(dirname "$f" | cut -c12-) # Get directory name; remove ./examples/
    mkdir -p "./jlexamples/$dname"
    jq -j '.cells
           | map( select(.cell_type == "code") | .source + ["\n\n"] )
           | .[][]' "$f" > ./jlexamples/"$dname"/"$fname".jl;
done
