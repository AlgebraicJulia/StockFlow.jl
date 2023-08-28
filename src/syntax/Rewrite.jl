module Rewrite
export sfrewrite




struct RewriteModifications
    added::Vector{Dict{Symbol, Vector{Dict}}} # :S => [Dict(:sname => :pop), ...], :LS = [(:lss => 1, :lsv => 1), ...], ...
    removed::Vector{Dict{Symbol, Vector{Dict}}} 
    swapped::Vector{Dict{Symbol, Symbol}} # to be parsed at the end and be added to added/removed, if hasn't been already.
end


function sfrewrite(sf, block)
    Base.remove_linenums!(block)

    sf_snames = Dict(S => i for (i, S) in enumerate(snames(sf)))
    sf_svnames = Dict(S => i for (i, S) in enumerate(svnames(sf)))
    sf_vnames = Dict(S => i for (i, S) in enumerate(vnames(sf)))
    sf_fnames = Dict(S => i for (i, S) in enumerate(fnames(sf)))
    sf_pnames = Dict(S => i for (i, S) in enumerate(pnames(sf)))
    
    name_vector = [collect(keys(names))... for names ∈ [sf_snames, sf_svnames, sf_vnames, sf_fnames, sf_pnames]]
    println(name_vector)
    @assert allunique(name_vector), "Not all names are unique!"

    L = StockAndFlowF()
    I = StockAndFlowF()
    R = StockAndFlowF()

    # So the way this is gonna work:
    # prefix an item with ! to indicate it's being removed from I
    # prefix an item with ⊕ to indicate it's being added to R
    # under swapping:
    # - use X => Y to indicate all instances of X in sf are to be swapped with Y in R
    # - use v1 = A + B => C + B to indicate this particular link is being swapped
    # (this makes it easier to swap links)
    # when all links to an object have been removed, it can be inferred that the object itself can be removed, but this isn't currently done.

    # Once we gather all the information:
    # - Create stockflow L from all links attached to parts indicated to be removed and swapped, without them removed
    # - Create stockflow I from all links attached to parts to be removed or swapped, with the removed and swapped parts removed
    # - Create stockflow R from all links attached to parts to be added or swapped to.
    # Define hom(I, L), hom(I, R)
    # Rewrite
    


    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:swapping) => begin
                current_phase = sw -> read_rewrite_line_and_update_dictionaries!(sw, name_vector, sf) # this one will be different than the others
            end
            QuoteNode(:stocks) => begin
                current_phase = s -> read_rewrite_line_and_update_dictionaries!(s, strata_snames, type_snames, aggregate_snames, strata_stock_mappings_dict, aggregate_stock_mappings_dict)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> read_rewrite_line_and_update_dictionaries!(p, strata_pnames, type_pnames, aggregate_pnames, strata_param_mappings_dict, aggregate_param_mappings_dict)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> read_rewrite_line_and_update_dictionaries!(v, strata_vnames, type_vnames, aggregate_vnames, strata_dyvar_mappings_dict, aggregate_dyvar_mappings_dict)
            end            
            QuoteNode(:flows) => begin
                current_phase = f -> read_rewrite_line_and_update_dictionaries!(f, strata_fnames, type_fnames, aggregate_fnames, strata_flow_mappings_dict, aggregate_flow_mappings_dict)
            end                    
                  
            QuoteNode(:sums) => begin
                current_phase = sv -> read_rewrite_line_and_update_dictionaries!(sv, strata_svnames, type_svnames, aggregate_svnames, strata_sum_mappings_dict, aggregate_sum_mappings_dict)
            end                    
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end
end


end