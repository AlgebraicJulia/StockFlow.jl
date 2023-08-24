module Homomorphism

using ...StockFlow
using MLStyle
using ..Syntax

import Base.get
import Catlab.CategoricalAlgebra.CSets.ACSetTransformation
import Catlab.CategoricalAlgebra.Limits.pullback
import Catlab.CategoricalAlgebra.FreeDiagrams.apex


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
macro homomorphism(sf, block) begin
    @assert sf.head == :tuple && length(sf.args) == 2
    @assert all(x -> x ∈ names(Main), sf.args)

    sfsrc, sftgt = map(x -> getfield(Main, x), sf.args)
    @assert all(x -> typeof(x) <: AbstractStockAndFlowStructureF, [sfsrc, sftgt]) # The infer links call is particular to stockflows 

    Base.remove_linenums!(block)

       
    src_snames = Dict(S => i for (i, S) in enumerate(snames(src)))
    src_svnames = Dict(S => i for (i, S) in enumerate(svnames(src)))
    src_vnames = Dict(S => i for (i, S) in enumerate(vnames(src)))
    src_fnames = Dict(S => i for (i, S) in enumerate(fnames(src)))
    src_pnames = Dict(S => i for (i, S) in enumerate(pnames(src)))

    tgt_snames = Dict(S => i for (i, S) in enumerate(snames(tgt)))
    tgt_svnames = Dict(S => i for (i, S) in enumerate(svnames(tgt)))
    tgt_vnames = Dict(S => i for (i, S) in enumerate(vnames(tgt)))
    tgt_fnames = Dict(S => i for (i, S) in enumerate(fnames(tgt)))
    tgt_pnames = Dict(S => i for (i, S) in enumerate(pnames(tgt)))


  current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_snames, tgt_snames, src_stock_mappings_dict)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_pnames, tgt_pnames, src_param_mappings_dict)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_vnames, tgt_vnames, src_dyvar_mappings_dict)
            end            
            QuoteNode(:flows) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_fnames, tgt_fnames, src_flow_mappings_dict)
            end                    
            QuoteNode(:sums) => begin
                current_phase = p -> read_homomorphism_line_and_update_dictionaries!(p, src_svnames, tgt_svnames, src_sum_mappings_dict)
            end                    
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    if !STRICT_MAPPINGS
        one_type_stock = length(snames(tgt)) == 1 ? 1 : 0 
        one_type_flow = length(fnames(tgt)) == 1 ? 1 : 0
        one_type_dyvar = length(vnames(tgt)) == 1 ? 1 : 0
        one_type_param = length(pnames(tgt)) == 1 ? 1 : 0
        one_type_sum = length(svnames(tgt)) == 1 ? 1 : 0
    else
        one_type_stock = one_type_flow = one_type_dyvar = one_type_param = one_type_sum = 0
    end

    src_stock_mappings::Vector{Int} = [get(src_stock_mappings_dict, i, one_tgt_stock)  for i in 1:ns(src)]
    src_flow_mappings::Vector{Int} =  [get(src_flow_mappings_dict, i, one_tgt_flow)  for i in 1:nf(src)]
    src_dyvar_mappings::Vector{Int} = [get(src_dyvar_mappings_dict, i, one_tgt_dyvar)  for i in 1:nvb(src)]
    src_param_mappings::Vector{Int} = [get(src_param_mappings_dict, i, one_tgt_param)  for i in 1:np(src)]
    src_sum_mappings::Vector{Int} =  [get(src_sum_mappings_dict, i, one_tgt_sum)  for i in 1:nsv(src)]




    all_mappings = [src_stock_mappings..., src_flow_mappings..., src_dyvar_mappings..., src_param_mappings..., src_sum_mappings...]

    @assert(all(x -> x != 0, all_mappings))


    nothing_function = x -> nothing
    no_attribute_tgt = map(sftgt, Name=name->nothing, Op=op->nothing, Position=pos->nothing)

    src_necmaps = Dict(:S => src_stock_mappings, :F => src_flow_mappings, :V => src_dyvar_mappings, :P => src_param_mappings, :SV => src_sum_mappings)    
    src_inferred_links = infer_links(src, type, src_necmaps)
    src_to_tgt = ACSetTransformation(src, no_attribute_tgt; src_necmaps..., src_inferred_links..., Op = nothing_function, Position = nothing_function, Name = nothing_function)

    return src_to_tgt


end

end





function read_homomorphism_line_and_update_dictionaries!(line::Expr, src_names::Dict{Symbol, Int}, tgt_names::Dict{Symbol, Int}, src_mappings::Dict{Int, Int})
    current_src_symbol_dict = interpret_homomorphism_notation(line) 

    current_src_dict = substitute_symbols_homomorphism(src_names, tgt_names, current_src_symbol_dict ; use_substr_prefix=USE_ISSUB, issubstr_prefix=ISSUB_DEFAULT)
    
    if STRICT_MATCHES, current_tgt_symbol_dict
        @assert (all(x -> x ∉ keys(src_mappings), keys(current_src_dict))) # check that we're not overwriting a value which has already been assigned
        merge!(src_mappings, current_src_dict) # accumulate dictionary keys

    else
        mergewith!((x, y) -> x, src_mappings, current_src_dict)
    end
end





# LMAO
function interpret_homomorphism_notation(mapping_pair::Expr)
    @match mapping_pair begin
        :($s => $t) => return (Dict(s => t))
        :($(shead...), $s => $t) => return push!(Dict(ss => t for ss in shead), s => t)
        _ => error("Unknown line format found in homomorphism notation.") 
    end

end









function substitute_symbols_homomorphism(s, t, ds ; use_substr_prefix=true, issubstr_prefix="Ξ")

    if !use_substr_prefix # this bit isn't necessary, as it's covered by the else block, but it's way simpler, and there may be cases where we don't want to check for substrings
        new_src_dict = Dict(s[strata_symbol] => t[type_symbol] for (strata_symbol, type_symbol) in ds)
        # new_aggregate_dict = Dict(a[aggregate_symbol] => t[type_symbol] for (aggregate_symbol, type_symbol) in da)
        return new_src_dict
    else

        @assert(allequal(values(da)))


        t_original_value::Symbol = only(Set(values(ds))) 
        t_val_string = string(t_original_value)


        if startswith(t_val_string, issubstr_prefix)
            t_match_string = chopprefix(t_val_string, issubstr_prefix)
            tgt_index = only(filter(((key, value),) ->  occursin(t_match_string, string(key)), t)).second # grab the only t value with t_match_string as a substring (there could be multiple, in which case throw error)
        else
            tgt_index = t[t_original_value]
        end



        new_src_dict = Dict()

        for src_key::Symbol in keys(ds)
            src_key_string = string(src_key)

            if startswith(src_key_string, issubstr_prefix)
                src_match_string = chopprefix(src_key_string, issubstr_prefix)
                push!(new_src_dict, [s[key] => tgt_index for (key,_) in filter(((key, value),) -> occursin(src_match_string, string(key)), s)]...) 
            else
                push!(new_src_dict, s[src_key] => tgt_index)
            end
        end

        return new_src_dict

    
        

    end
end

end
