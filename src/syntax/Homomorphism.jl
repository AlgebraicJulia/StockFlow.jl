module Homomorphism

using ...StockFlow
using ..Syntax
import ..Syntax: infer_links, substitute_symbols, DSLArgument, NothingFunction, invert_vector, complete_mappings
using MLStyle
using Catlab.CategoricalAlgebra

export @hom

macro hom(sf1, sf2, block)
  escaped_block = Expr(:quote, block)
  quote 
      hom($(esc(sf1)), $(esc(sf2)), $escaped_block);
  end
end

function interpret_homomorphism_syntax(line :: Expr)::Vector{DSLArgument}
  @match line begin
    :($s => $t) => [DSLArgument(s,t)]
    :($(shead...), $s => $t) => [[[DSLArgument(ss, t) for ss in shead] ; DSLArgument(s, t)]]
    _ => error("Unknown symbol format.  Must be A => B.")
  end
end


function hom(sf1, sf2, block)
  Base.remove_linenums!(block)
  stocks::Vector{DSLArgument} = []
  params::Vector{DSLArgument} = []
  dyvars::Vector{DSLArgument} = []
  flows::Vector{DSLArgument} = []
  sums::Vector{DSLArgument} = []
  current_phase = (_, _) -> ()
  for statement in block.args
      @match statement begin
          QuoteNode(:stocks) => begin
            current_phase = s -> append!(stocks, interpret_homomorphism_syntax(s))
          end
          QuoteNode(:parameters) => begin
            current_phase = p -> append!(params, interpret_homomorphism_syntax(p))
          end
          QuoteNode(:dynamic_variables) => begin
            current_phase = d -> append!(dyvars, interpret_homomorphism_syntax(d))
          end
          QuoteNode(:flows) => begin
            current_phase = f -> append!(flows, interpret_homomorphism_syntax(f))
          end
          QuoteNode(:sums) => begin
            current_phase = s -> append!(sums, interpret_homomorphism_syntax(s))
          end
          QuoteNode(kw) =>
            error("Unknown block type for homomorphism syntax: " * String(kw))
         _ => current_phase(statement)
      end
  end

  # # TODO: Assert that tempstratdefault not in any names
  # push!(stocks, DSLArgument(:_, -1, Set{Symbol}()))
  # push!(params, DSLArgument(:_, -1, Set{Symbol}()))
  # push!(dyvars, DSLArgument(:_, -1, Set{Symbol}()))
  # push!(flows, DSLArgument(:_, -1, Set{Symbol}()))
  # push!(sums, DSLArgument(:_, -1, Set{Symbol}()))
  sf1_snames = invert_vector(snames(sf1), Symbol)
  sf1_svnames = invert_vector(svnames(sf1), Symbol)
  sf1_vnames = invert_vector(vnames(sf1), Symbol)
  sf1_fnames = invert_vector(fnames(sf1), Symbol)
  sf1_pnames = invert_vector(pnames(sf1), Symbol)

  push!(sf1_snames, :_ => -1)
  push!(sf1_svnames, :_ => -1)
  push!(sf1_vnames, :_ => -1)
  push!(sf1_fnames, :_ => -1)
  push!(sf1_pnames, :_ => -1)



  master_vector = [substitute_symbols(sf1_snames, invert_vector(snames(sf2), Symbol), stocks),
    substitute_symbols(sf1_svnames, invert_vector(svnames(sf2), Symbol), sums),
    substitute_symbols(sf1_vnames, invert_vector(vnames(sf2), Symbol), dyvars),
    substitute_symbols(sf1_fnames, invert_vector(fnames(sf2), Symbol), flows),
    substitute_symbols(sf1_pnames, invert_vector(pnames(sf2), Symbol), params),
  ]::Vector{Dict{Int, Int}}

  # map(x -> push!(x, -1 => :_), master_vector)

  new_master_dict = complete_mappings(sf1, sf2, master_vector)

  links = infer_links(sf1, sf2, new_master_dict)

  no_attribute_type = map(sf2, Name=NothingFunction, 
    Op=NothingFunction, Position=NothingFunction)

  transformation = ACSetTransformation(sf1, no_attribute_type, ; filter(kv -> !isempty(kv[2]), new_master_dict)..., filter(kv -> !isempty(kv[2]), links)..., :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)
  @assert is_natural(transformation)
  return transformation

  
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


# function apply_hom(hom::Dict{Symbol, Vector{Pair{Symbol,Symbol}}}, sf1, sf2)
#   srcnames = all_names_to_index(sf1)
#   tgtnames = all_names_to_index(sf2)

#   homstocks = [srcnames[:S][name1] => tgtnames[:S][name2] for (name1, name2) in hom[:S]]

#   homflows = [srcnames[:F][name1] => tgtnames[:F][name2] for (name1, name2) in hom[:F]]
#   homparams = [srcnames[:P][name1] => tgtnames[:P][name2] for (name1, name2) in hom[:P]]
#   homdyvars = [srcnames[:V][name1] => tgtnames[:V][name2] for (name1, name2) in hom[:V]]
#   homsums = [srcnames[:SV][name1] => tgtnames[:SV][name2] for (name1, name2) in hom[:SV]]

#   # nec_maps = Dict{Symbol, Vector{Int64}}(
#   #   :S => map()
#   # )

#   nec_maps = Dict{Symbol, Vector{Int64}}(:S => map(last, sort!(homstocks, by=x -> x[1])), 
#   :F => map(last, sort!(homflows, by=x -> x[1])), 
#   :SV => map(last, sort!(homsums, by=x -> x[1])), 
#   :P => map(last, sort!(homparams, by=x -> x[1])), 
#   :V => map(last, sort!(homdyvars, by=x -> x[1]))
#   )

#   links = infer_links(sf1, sf2, nec_maps)




# end


end