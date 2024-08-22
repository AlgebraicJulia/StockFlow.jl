#!/bin/env bash

subdir=$1
find ./examples/$subdir -type f -name "*.ipynb" | while read f; do
    fname=$(basename "$f")
    jq -j '.cells
           | map( select(.cell_type == "code") | .source + ["\n\n"] )
           | .[][]' "$f" > ./jlexamples/"$subdir"."$fname".jl;
done
