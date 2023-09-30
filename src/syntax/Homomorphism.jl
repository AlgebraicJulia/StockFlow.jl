module Homomorphism

using ...StockFlow
using ..Syntax
using ..Syntax: NothingFunction
using MLStyle
using Catlab.CategoricalAlgebra

export @hom

macro hom(sf1, sf2, block)
  escaped_block = Expr(:quote, block)
  quote 
      hom($(esc(sf1)), $(esc(sf2)), $escaped_block);
  end
end

function hom(sf1, sf2, block)
  Base.remove_linenums!(block)
  stocks::Vector{Pair{Symbol,Symbol}} = []
  params::Vector{Pair{Symbol,Symbol}} = []
  dyvars::Vector{Pair{Symbol,Symbol}} = []
  flows::Vector{Pair{Symbol,Symbol}} = []
  sums::Vector{Pair{Symbol,Symbol}} = []
  current_phase = (_, _) -> ()
  for statement in block.args
      @match statement begin
          QuoteNode(:stocks) => begin
            current_phase = s -> push!(stocks, s)
          end
          QuoteNode(:parameters) => begin
            current_phase = p -> push!(params, p)
          end
          QuoteNode(:dynamic_variables) => begin
            current_phase = d -> push!(dyvars, d)
          end
          QuoteNode(:flows) => begin
            current_phase = f -> push!(flows, f)
          end
          QuoteNode(:sums) => begin
            current_phase = s -> push!(sums, s)
          end
          QuoteNode(kw) =>
            error("Unknown block type for homomorphism syntax: " * String(kw))
          :($A => $B) => current_phase(A => B)
          _ => error("Unknown symbol format.  Must be A => B.")
      end
  end
  explicit_mappings =  Dict(:S => stocks, :F => flows, :V => dyvars, :P => params, :SV => sums)
  apply_hom(explicit_mappings, sf1, sf2)
end

function names_to_index(sf, func)
  Dict(name => i for (i, name) in enumerate(func(sf)))
end

function all_names_to_index(sf)
  Dict(:S => names_to_index(sf, snames),
  :F => names_to_index(sf, fnames),
  :V => names_to_index(sf, vnames),
  :P => names_to_index(sf, pnames),
  :SV => names_to_index(sf, svnames)
  )
end


function apply_hom(hom::Dict{Symbol, Vector{Pair{Symbol,Symbol}}}, sf1, sf2)
  srcnames = all_names_to_index(sf1)
  tgtnames = all_names_to_index(sf2)
  println(srcnames)
  println(tgtnames)

  homstocks = [srcnames[:S][name1] => tgtnames[:S][name2] for (name1, name2) in hom[:S]]
  homflows = [srcnames[:F][name1] => tgtnames[:F][name2] for (name1, name2) in hom[:F]]
  homparams = [srcnames[:P][name1] => tgtnames[:P][name2] for (name1, name2) in hom[:P]]
  homdyvars = [srcnames[:V][name1] => tgtnames[:V][name2] for (name1, name2) in hom[:V]]
  homsums = [srcnames[:SV][name1] => tgtnames[:SV][name2] for (name1, name2) in hom[:SV]]


  nec_maps = Dict{Symbol, Vector{Int64}}(:S => map(first, sort!(homstocks, by=x -> x[2])), 
  :F => map(first, sort!(homflows, by=x -> x[2])), 
  :SV => map(first, sort!(homsums, by=x -> x[2])), 
  :P => map(first, sort!(homparams, by=x -> x[2])), 
  :V => map(first, sort!(homdyvars, by=x -> x[2]))
  )
  println(nec_maps)

  links = infer_links(sf1, sf2, nec_maps)


  ACSetTransformation(sf1, sf2, ; filter(kv -> !isempty(kv[2]), nec_maps)..., filter(kv -> !isempty(kv[2]), links)..., :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)


end


end