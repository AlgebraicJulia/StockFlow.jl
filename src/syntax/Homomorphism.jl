module Homomorphism

export @hom


using ...StockFlow
using MLStyle
using ..Syntax

using Catlab.CategoricalAlgebra

import ..Syntax: STRICT_MAPPINGS, STRICT_MATCHES, USE_ISSUB, ISSUB_DEFAULT, infer_links, substitute_symbols




"""
I mean, may as well.  Comes for free after doing the stratification work.  Just need to parse the notation somewhat differently.
Gonna be a little annoying to deal with all these changes in a pull request, sorry about that.

@homomorphism (src, tgt) begin
    :stocks
    A => A′

    :flows
    f1′ => f2′
   ⋮ 
end


"""
macro hom(sf, block)
    @assert sf.head == :tuple && length(sf.args) == 2
    @assert all(x -> x ∈ names(Main), sf.args)

    sfsrc, sftgt = map(x -> getfield(Main, x), sf.args)
    @assert all(x -> typeof(x) <: AbstractStockAndFlowStructureF, [sfsrc, sftgt]) # The infer links call is particular to stockflows 

    Base.remove_linenums!(block)

       
    src_snames = Dict(S => i for (i, S) in enumerate(snames(sfsrc)))
    src_svnames = Dict(S => i for (i, S) in enumerate(svnames(sfsrc)))
    src_vnames = Dict(S => i for (i, S) in enumerate(vnames(sfsrc)))
    src_fnames = Dict(S => i for (i, S) in enumerate(fnames(sfsrc)))
    src_pnames = Dict(S => i for (i, S) in enumerate(pnames(sfsrc)))

    tgt_snames = Dict(S => i for (i, S) in enumerate(snames(sftgt)))
    tgt_svnames = Dict(S => i for (i, S) in enumerate(svnames(sftgt)))
    tgt_vnames = Dict(S => i for (i, S) in enumerate(vnames(sftgt)))
    tgt_fnames = Dict(S => i for (i, S) in enumerate(fnames(sftgt)))
    tgt_pnames = Dict(S => i for (i, S) in enumerate(pnames(sftgt)))

    src_stock_mappings_dict::Dict{Int, Int} = Dict()
    src_flow_mappings_dict::Dict{Int, Int} = Dict()
    src_dyvar_mappings_dict::Dict{Int, Int} = Dict()
    src_param_mappings_dict::Dict{Int, Int} = Dict()
    src_sum_mappings_dict::Dict{Int, Int} = Dict()
    

    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> read_homomorphism_line_and_update_dictionaries!(s, src_snames, tgt_snames, src_stock_mappings_dict)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_pnames, tgt_pnames, src_param_mappings_dict)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> read_homomorphism_line_and_update_dictionaries!(v, src_vnames, tgt_vnames, src_dyvar_mappings_dict)
            end            
            QuoteNode(:flows) => begin
                current_phase = f -> read_homomorphism_line_and_update_dictionaries!(f, src_fnames, tgt_fnames, src_flow_mappings_dict)
            end                    
            QuoteNode(:sums) => begin
                current_phase = sv -> read_homomorphism_line_and_update_dictionaries!(sv, src_svnames, tgt_svnames, src_sum_mappings_dict)
            end                    
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    if !STRICT_MAPPINGS
        one_tgt_stock = length(snames(sftgt)) == 1 ? 1 : 0 
        one_tgt_flow = length(fnames(sftgt)) == 1 ? 1 : 0
        one_tgt_dyvar = length(vnames(sftgt)) == 1 ? 1 : 0
        one_tgt_param = length(pnames(sftgt)) == 1 ? 1 : 0
        one_tgt_sum = length(svnames(sftgt)) == 1 ? 1 : 0
    else
        one_tgt_stock = one_tgt_flow = one_tgt_dyvar = one_tgt_param = one_tgt_sum = 0
    end

    src_stock_mappings::Vector{Int} = [get(src_stock_mappings_dict, i, one_tgt_stock)  for i in 1:ns(sfsrc)]
    src_flow_mappings::Vector{Int} =  [get(src_flow_mappings_dict, i, one_tgt_flow)  for i in 1:nf(sfsrc)]
    src_dyvar_mappings::Vector{Int} = [get(src_dyvar_mappings_dict, i, one_tgt_dyvar)  for i in 1:nvb(sfsrc)]
    src_param_mappings::Vector{Int} = [get(src_param_mappings_dict, i, one_tgt_param)  for i in 1:np(sfsrc)]
    src_sum_mappings::Vector{Int} =  [get(src_sum_mappings_dict, i, one_tgt_sum)  for i in 1:nsv(sfsrc)]




    all_mappings = [src_stock_mappings..., src_flow_mappings..., src_dyvar_mappings..., src_param_mappings..., src_sum_mappings...]

    @assert(all(x -> x != 0, all_mappings))


    nothing_function = x -> nothing
    no_attribute_tgt = map(sftgt, Name=name->nothing, Op=op->nothing, Position=pos->nothing)

    src_necmaps = Dict(:S => src_stock_mappings, :F => src_flow_mappings, :V => src_dyvar_mappings, :P => src_param_mappings, :SV => src_sum_mappings)    
    src_inferred_links = infer_links(sfsrc, sftgt, src_necmaps)
    src_to_tgt = ACSetTransformation(sfsrc, no_attribute_tgt; src_necmaps..., src_inferred_links..., Op = nothing_function, Position = nothing_function, Name = nothing_function)

    @assert is_natural(src_to_tgt) # unfortunately, at this point, all the attributes will be nothing.

    return src_to_tgt

end





function read_homomorphism_line_and_update_dictionaries!(line::Expr, src_names::Dict{Symbol, Int}, tgt_names::Dict{Symbol, Int}, src_mappings::Dict{Int, Int})
    current_src_symbol_dict = interpret_homomorphism_notation(line) 

    current_src_dict = substitute_symbols(src_names, tgt_names, current_src_symbol_dict ; use_substr_prefix=USE_ISSUB, issubstr_prefix=ISSUB_DEFAULT)
    
    if STRICT_MATCHES
        @assert (all(x -> x ∉ keys(src_mappings), keys(current_src_dict))) # check that we're not overwriting a value which has already been assigned
        merge!(src_mappings, current_src_dict) # accumulate dictionary keys

    else
        mergewith!((x, y) -> x, src_mappings, current_src_dict)
    end
end



function interpret_homomorphism_notation(mapping_pair::Expr)
    @match mapping_pair begin
        :($s => $t) => return (Dict(s => t))
        :($(shead...), $s => $t) => return push!(Dict(ss => t for ss in shead), s => t)
        _ => error("Unknown line format found in homomorphism notation.") 
    end

end



end
