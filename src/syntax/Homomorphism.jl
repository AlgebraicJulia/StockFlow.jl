module Homomorphism

using ...StockFlow
import ...StockFlow: state_dict

using ..Syntax

import ..Syntax: infer_links, substitute_symbols, DSLArgument, NothingFunction, 
    complete_mappings
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

  sf1_snames = Dict{Symbol, Int}(state_dict(snames(sf1)))
  sf1_svnames = Dict{Symbol, Int}(state_dict(svnames(sf1)))
  sf1_vnames = Dict{Symbol, Int}(state_dict(vnames(sf1)))
  sf1_fnames = Dict{Symbol, Int}(state_dict(fnames(sf1)))
  sf1_pnames = Dict{Symbol, Int}(state_dict(pnames(sf1)))

  @assert (:_ ∉ keys(sf1_snames) && :_ ∉ keys(sf1_svnames) 
    && :_ ∉ keys(sf1_vnames) && :_ ∉ keys(sf1_fnames) && 
    :_ ∉ keys(sf1_pnames))  "A stockflow contains :_ !  \
    Please change temp_strat_default to a different symbol or \
    rename offending object."

  push!(sf1_snames, :_ => -1)
  push!(sf1_svnames, :_ => -1)
  push!(sf1_vnames, :_ => -1)
  push!(sf1_fnames, :_ => -1)
  push!(sf1_pnames, :_ => -1)



  master_vector = [substitute_symbols(sf1_snames, Dict{Symbol, Int}(state_dict(snames(sf2))), stocks),
    substitute_symbols(sf1_svnames, Dict{Symbol, Int}(state_dict(svnames(sf2))), sums),
    substitute_symbols(sf1_vnames, Dict{Symbol, Int}(state_dict(vnames(sf2))), dyvars),
    substitute_symbols(sf1_fnames, Dict{Symbol, Int}(state_dict(fnames(sf2))), flows),
    substitute_symbols(sf1_pnames, Dict{Symbol, Int}(state_dict(pnames(sf2))), params),
  ]::Vector{Dict{Int, Int}}


  new_master_dict = complete_mappings(sf1, sf2, master_vector)

  links = infer_links(sf1, sf2, new_master_dict)

  no_attribute_type = map(sf2, Name=NothingFunction, 
    Op=NothingFunction, Position=NothingFunction)

  transformation = ACSetTransformation(sf1, no_attribute_type, ; filter(kv -> !isempty(kv[2]), new_master_dict)..., filter(kv -> !isempty(kv[2]), links)..., :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)
  @assert is_natural(transformation)
  return transformation

  
end

function all_names_to_index(sf)
  Dict(:S => state_dict(snames(sf))
  :F => state_dict(fnames(sf))
  :V => state_dict(vnames(sf))
  :P => state_dict(pnames(sf))
  :SV => state_dict(svnames(sf))
  )
end




end