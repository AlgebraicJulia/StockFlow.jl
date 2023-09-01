module Stratification
export sfstratify

using ...StockFlow
using ..Syntax
using MLStyle
import Base: get
using Catlab.CategoricalAlgebra
import ..Syntax: infer_links, substitute_symbols, DSLArgument, NothingFunction, invert_vector


"""
    interpret_stratification_notation(mapping_pair::Expr)::Tuple{Vector{DSLArgument}, Vector{DSLArgument}}
Take an expression of the form a1, ..., => t <= s1, ..., where every element is a symbol, and return a 2-tuple of form ((a1 => t, a2 => t, ...), (s1 => t, ...))
"""
function interpret_stratification_notation(mapping_pair::Expr)::Tuple{Vector{DSLArgument}, Vector{DSLArgument}}
    @match mapping_pair begin


        :($s => $t <= $a) => return ([DSLArgument(s,t)], [DSLArgument(a,t)])
        :($s => $t <= $a, $(atail...)) => ([DSLArgument(s,t)], [DSLArgument(a,t) ; [DSLArgument(as,t) for as in atail] ])#return (Dict(unwrap_key_expression(s, t)), push!(Dict(unwrap_key_expression(as, t) for as in atail), unwrap_key_expression(a, t)))
        :($(shead...), $s => $t <= $a) => ([[DSLArgument(ss, t) for ss in shead] ; DSLArgument(s, t)], [DSLArgument(a, t)])#return (push!(Dict(unwrap_key_expression(ss, t) for ss in shead), unwrap_key_expression(s, t)), Dict(unwrap_key_expression(a, t)))

        if mapping_pair.head == :tuple end => begin
            middle_index = findfirst(x -> typeof(x) == Expr && length(x.args) == 3, mapping_pair.args) # still isn't specific enough
            if isnothing(middle_index)
                error("Malformed line $mapping_pair, could not find center.")
            end
            @match mapping_pair.args[middle_index] begin
                :($stail => $t <= $ahead) => begin
                sdict = [[DSLArgument(ss, t) for ss in mapping_pair.args[1:middle_index-1]] ; DSLArgument(stail, t)]
                adict = [DSLArgument(ahead, t) ; [DSLArgument(as, t) for as in mapping_pair.args[middle_index+1:end]]]
                    return (sdict, adict)
                end
                _ => "Unknown format found for match; middle three values formatted incorrectly."
            end
        end
        _ => error("Unknown line format found in stratification notation.") 
    end
end




function read_stratification_line_and_update_dictionaries!(line::Expr, strata_names::Dict{Symbol, Int}, type_names::Dict{Symbol, Int}, aggregate_names::Dict{Symbol, Int}, strata_mappings::Dict{Int, Int}, aggregate_mappings::Dict{Int, Int} ; strict_matches = false, use_flags = true)
    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(line)

    current_strata_dict = substitute_symbols(strata_names, type_names, current_strata_symbol_dict ; use_flags=use_flags)
    current_aggregate_dict = substitute_symbols(aggregate_names, type_names, current_aggregate_symbol_dict ; use_flags=use_flags)
    
    if strict_matches
        @assert (all(x -> x ∉ keys(strata_mappings), keys(current_strata_dict))) "Attempt to overwrite a mapping in strata!"
        # check that we're not overwriting a value which has already been assigned
        merge!(strata_mappings, current_strata_dict) # accumulate dictionary keys


        @assert (all(x -> x ∉ keys(aggregate_mappings), keys(current_aggregate_dict))) "Attempt to overwrite a mapping in aggregate!"
        merge!(aggregate_mappings, current_aggregate_dict)

    else 
        mergewith!((x, y) -> x, strata_mappings, current_strata_dict) # alternatively, can use: only ∘ first
        mergewith!((x, y) -> x, aggregate_mappings, current_aggregate_dict)
    end

end

"""

"""
function print_unmapped(mappings::Vector{Pair{Vector{Int}, Vector{Symbol}}}, name="STOCKFLOW")
    for (ints, dicts) in mappings
        for (i, val) in enumerate(ints)
            if val == 0
                println("UNMAPPED IN $(name):")
                println(dicts[i])
            end
        end
    end
end

"""
Iterates over each line in a stratification quoteblock and updates the appropriate dictionaries
"""
function iterate_over_stratification_lines!(block, strata_names, type_names, aggregate_names, strata_mappings, aggregate_mappings; strict_matches=false, use_flags=true)
    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> read_stratification_line_and_update_dictionaries!(s, strata_names[1], type_names[1], aggregate_names[1], strata_mappings[1], aggregate_mappings[1]; strict_matches=strict_matches, use_flags=use_flags)
            end
            QuoteNode(:sums) => begin
                current_phase = sv -> read_stratification_line_and_update_dictionaries!(sv, strata_names[2], type_names[2], aggregate_names[2], strata_mappings[2], aggregate_mappings[2]; strict_matches=strict_matches, use_flags=use_flags)
            end        
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> read_stratification_line_and_update_dictionaries!(v, strata_names[3], type_names[3], aggregate_names[3], strata_mappings[3], aggregate_mappings[3]; strict_matches=strict_matches, use_flags=use_flags)
            end       
            QuoteNode(:flows) => begin
                current_phase = f -> read_stratification_line_and_update_dictionaries!(f, strata_names[4], type_names[4], aggregate_names[4], strata_mappings[4], aggregate_mappings[4]; strict_matches=strict_matches, use_flags=use_flags)
            end                    
            QuoteNode(:parameters) => begin
                current_phase = p -> read_stratification_line_and_update_dictionaries!(p, strata_names[5], type_names[5], aggregate_names[5], strata_mappings[5], aggregate_mappings[5]; strict_matches=strict_matches, use_flags=use_flags)
            end
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end
end


"""
    sfstratify(strata, type, aggregate, block ; kwargs)

    Ok, so the general idea here is:
    1. Grab all names from strata, type and aggregate, and create dictionaries which map them to their indices
    2. iterate over each line in the block
        2a. Split each line into a dictionary which maps all strata to that type and all aggregate to that type
        2b. Convert from two Symbol => Symbol dictionaries to two Int => Int dictionaries, using the dictionaries from step 1
            2bα. If applicable, for symbols with temp_strat_default as a prefix, find all matching symbols in the symbol dictionaries, and map all those
        2c. Accumulate respective dictionaries (optionally, only allow first match vs throw an error (strict_matches = false vs true))
    3. Create an array of 0s for stocks, flows, parameters, dyvars and sums for strata and aggregate.   Insert into arrays all values from the two Int => Int dictionaries
        3a. If strict_mappings = false, if there only exists one option in type to map to, and it hasn't been explicitly specified, add it.  If strict_mappings = true and it hasn't been specified, throw an error.
    4. Do a once-over of arrays and ensure there aren't any zeroes (unmapped values) remaining (helps with debugging when you screw up stratifying)
    5. Deal with attributes (create a copy of type sf with attributes mapped to nothing)
    6. Infer LS, LSV, etc., if possible.
    7. Construct strata -> type and aggregate -> type ACSetTransformations
    8. Return pullback (with flattened attributes)
"""
function sfstratify(strata::AbstractStockAndFlowStructureF, type::AbstractStockAndFlowStructureF, aggregate::AbstractStockAndFlowStructureF, block::Expr ; strict_mappings = false, strict_matches = false, temp_strat_default = :_, use_temp_strat_default = true, use_flags = true, return_homs = false)

    Base.remove_linenums!(block)

    # STEP 1
    
    # invert_vector: Vector{K} -> Dict{K, Int} where int is original index and all K (symbols, in this case) are unique.
    strata_snames::Dict{Symbol, Int} = invert_vector(snames(strata)) 
    strata_svnames::Dict{Symbol, Int} = invert_vector(svnames(strata))
    strata_vnames::Dict{Symbol, Int} = invert_vector(vnames(strata))
    strata_fnames::Dict{Symbol, Int} = invert_vector(fnames(strata))
    strata_pnames::Dict{Symbol, Int} = invert_vector(pnames(strata))

    type_snames::Dict{Symbol, Int} = invert_vector(snames(type)) 
    type_svnames::Dict{Symbol, Int} = invert_vector(svnames(type))
    type_vnames::Dict{Symbol, Int} = invert_vector(vnames(type))
    type_fnames::Dict{Symbol, Int} = invert_vector(fnames(type))
    type_pnames::Dict{Symbol, Int} = invert_vector(pnames(type))

    aggregate_snames::Dict{Symbol, Int} = invert_vector(snames(aggregate))
    aggregate_svnames::Dict{Symbol, Int} = invert_vector(svnames(aggregate))
    aggregate_vnames::Dict{Symbol, Int} = invert_vector(vnames(aggregate))
    aggregate_fnames::Dict{Symbol, Int} = invert_vector(fnames(aggregate))
    aggregate_pnames::Dict{Symbol, Int} = invert_vector(pnames(aggregate))


    strata_stock_mappings_dict::Dict{Int, Int} = Dict()
    strata_flow_mappings_dict::Dict{Int, Int} = Dict()
    strata_dyvar_mappings_dict::Dict{Int, Int} = Dict()
    strata_param_mappings_dict::Dict{Int, Int} = Dict()
    strata_sum_mappings_dict::Dict{Int, Int} = Dict()
    
    aggregate_stock_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_flow_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_dyvar_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_param_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_sum_mappings_dict::Dict{Int, Int} = Dict()
    

    strata_all_name_mappings::Vector{Dict{Symbol, Int}} = [strata_snames, strata_svnames, strata_vnames, strata_fnames, strata_pnames]
    type_all_name_mappings::Vector{Dict{Symbol, Int}} = [type_snames, type_svnames, type_vnames, type_fnames, type_pnames]
    aggregate_all_name_mappings::Vector{Dict{Symbol, Int}} = [aggregate_snames, aggregate_svnames, aggregate_vnames, aggregate_fnames, aggregate_pnames]

    strata_all_index_mappings = [strata_stock_mappings_dict, strata_sum_mappings_dict, strata_dyvar_mappings_dict, strata_flow_mappings_dict, strata_param_mappings_dict]
    aggregate_all_index_mappings = [aggregate_stock_mappings_dict, aggregate_sum_mappings_dict, aggregate_dyvar_mappings_dict, aggregate_flow_mappings_dict, aggregate_param_mappings_dict]


    if use_temp_strat_default

        strata_all_names::Vector{Vector{Symbol}} = [snames(strata), svnames(strata), vnames(strata), fnames(strata), pnames(strata)]
        aggregate_all_names::Vector{Vector{Symbol}} = [snames(aggregate), svnames(aggregate), vnames(aggregate), fnames(aggregate), pnames(aggregate)]

        @assert all(x -> temp_strat_default ∉ keys(x), strata_all_names) "Strata contains $temp_strat_default !  Please change temp_strat_default to a different symbol or rename offending object."
        @assert all(x -> temp_strat_default ∉ keys(x), aggregate_all_names) "Aggregate contains $temp_strat_default !  Please change temp_strat_default to a different symbol or rename offending object."

        map(x -> (push!(x, (temp_strat_default => -1))), strata_all_name_mappings)
        map(x -> (push!(x, (temp_strat_default => -1))), aggregate_all_name_mappings)
    end

    # STEP 2
    iterate_over_stratification_lines!(block, strata_all_name_mappings, type_all_name_mappings, aggregate_all_name_mappings, strata_all_index_mappings, aggregate_all_index_mappings ; strict_matches=strict_matches, use_flags=use_flags)


    # get the default value, if it has been assigned.  Use 0 if it hasn't.
    default_index_strata_stock = get(strata_stock_mappings_dict, -1, 0)
    default_index_strata_flow = get(strata_flow_mappings_dict, -1, 0)
    default_index_strata_dyvar = get(strata_dyvar_mappings_dict, -1, 0)
    default_index_strata_param = get(strata_param_mappings_dict, -1, 0)
    default_index_strata_sum = get(strata_sum_mappings_dict, -1, 0)

    default_index_aggregate_stock = get(aggregate_stock_mappings_dict, -1, 0)
    default_index_aggregate_flow = get(aggregate_flow_mappings_dict, -1, 0)
    default_index_aggregate_dyvar = get(aggregate_dyvar_mappings_dict, -1, 0)
    default_index_aggregate_param = get(aggregate_param_mappings_dict, -1, 0)
    default_index_aggregate_sum = get(aggregate_sum_mappings_dict, -1, 0)


    # STEP 3
    if !strict_mappings
        one_type_stock = length(snames(type)) == 1 ? 1 : 0 # if there is only one stock, it needs to have index 1
        one_type_flow = length(fnames(type)) == 1 ? 1 : 0
        one_type_dyvar = length(vnames(type)) == 1 ? 1 : 0
        one_type_param = length(pnames(type)) == 1 ? 1 : 0
        one_type_sum = length(svnames(type)) == 1 ? 1 : 0
    else
        one_type_stock = one_type_flow = one_type_dyvar = one_type_param = one_type_sum = 0
    end

    # Convert back to vectors.  If you find a zero, check if there's a default and use that.  If there isn't a default, check if there's only one option and use that.
    # Otherwise, there's an unassigned value which can't be inferred.
    # Taking max because it's less verbose than ternary and accomplishes the same thing:
    # - If both default_index and one_type are mapped, they must be mapped to the same thing, because one_type being mapped implies there's only one option.
    # - If only one_type is mapped, then it will be positive, and default_infex will be 0
    # - If only default_index is mapped, it will be positive and one_type will be 0

    strata_stock_mappings::Vector{Int} = [get(strata_stock_mappings_dict, i, max(default_index_strata_stock, one_type_stock))  for i in 1:ns(strata)]
    strata_flow_mappings::Vector{Int} =  [get(strata_flow_mappings_dict, i, max(default_index_strata_flow, one_type_flow))  for i in 1:nf(strata)]
    strata_dyvar_mappings::Vector{Int} = [get(strata_dyvar_mappings_dict, i, max(default_index_strata_dyvar, one_type_dyvar))  for i in 1:nvb(strata)]
    strata_param_mappings::Vector{Int} = [get(strata_param_mappings_dict, i,  max(default_index_strata_param, one_type_param))  for i in 1:np(strata)]
    strata_sum_mappings::Vector{Int} =  [get(strata_sum_mappings_dict, i,  max(default_index_strata_sum, one_type_sum))  for i in 1:nsv(strata)]

    aggregate_stock_mappings::Vector{Int} = [get(aggregate_stock_mappings_dict, i, max(default_index_aggregate_stock, one_type_stock))  for i in 1:ns(aggregate)]
    aggregate_flow_mappings::Vector{Int} =  [get(aggregate_flow_mappings_dict, i, max(default_index_aggregate_flow, one_type_flow))  for i in 1:nf(aggregate)]
    aggregate_dyvar_mappings::Vector{Int} = [get(aggregate_dyvar_mappings_dict, i, max(default_index_aggregate_dyvar, one_type_dyvar))  for i in 1:nvb(aggregate)]
    aggregate_param_mappings::Vector{Int} = [get(aggregate_param_mappings_dict, i,  max(default_index_aggregate_param, one_type_param))  for i in 1:np(aggregate)]
    aggregate_sum_mappings::Vector{Int} =  [get(aggregate_sum_mappings_dict, i,  max(default_index_aggregate_sum, one_type_sum))  for i in 1:nsv(aggregate)]


    # This bit is a bit verbose, but makes debugging when making a stratification easier.  Tells you exactly which ones you forgot to map.

    all_mappings = [strata_stock_mappings..., strata_flow_mappings..., strata_dyvar_mappings..., strata_param_mappings..., strata_sum_mappings..., aggregate_stock_mappings..., aggregate_flow_mappings..., aggregate_dyvar_mappings..., aggregate_param_mappings..., aggregate_sum_mappings...]
    
    strata_mappings::Vector{Pair{Vector{Int}, Vector{Symbol}}} = [strata_stock_mappings => snames(strata), strata_flow_mappings => fnames(strata), strata_dyvar_mappings => vnames(strata), strata_param_mappings => pnames(strata), strata_sum_mappings => svnames(strata)]
    aggregate_mappings::Vector{Pair{Vector{Int}, Vector{Symbol}}} = [aggregate_stock_mappings => snames(aggregate), aggregate_flow_mappings => fnames(aggregate), aggregate_dyvar_mappings => vnames(aggregate), aggregate_param_mappings => pnames(aggregate), aggregate_sum_mappings => svnames(aggregate)]

    # STEP 4

    #unmapped: 
    if !(all(x -> x != 0, all_mappings))
        print_unmapped(strata_mappings, "STRATA")
        print_unmapped(aggregate_mappings, "AGGREGATE")
        error("There is an unmapped value!")
    end
    

    # STEP 5
    # NothingFunction(x...) = nothing;
    no_attribute_type = map(type, Name=NothingFunction, Op=NothingFunction, Position=NothingFunction)

    # STEP 6/7 
    # This is where we pull out the magic to infer links.
    #
    #  A <- C -> B
    #  ||        ||
    #  v         v
    #  A'<- C'-> B'
    #
    # implies
    #
    #  A <- C -> B
    #  ||  ||    ||
    #  v    v    v
    #  A'<- C'-> B'
    #
    
    strata_necmaps = Dict(:S => strata_stock_mappings, :F => strata_flow_mappings, :V => strata_dyvar_mappings, :P => strata_param_mappings, :SV => strata_sum_mappings)    
    strata_inferred_links = infer_links(strata, type, strata_necmaps)
    strata_to_type = ACSetTransformation(strata, no_attribute_type; strata_necmaps..., strata_inferred_links..., Op = NothingFunction, Position = NothingFunction, Name = NothingFunction)

    
    aggregate_necmaps = Dict(:S => aggregate_stock_mappings, :F => aggregate_flow_mappings, :V => aggregate_dyvar_mappings, :P => aggregate_param_mappings, :SV => aggregate_sum_mappings)
    aggregate_inferred_links = infer_links(aggregate, type, aggregate_necmaps)
    aggregate_to_type = ACSetTransformation(aggregate, no_attribute_type; aggregate_necmaps..., aggregate_inferred_links..., Op = NothingFunction, Position = NothingFunction, Name =NothingFunction)

    

    # STEP 8
    pullback_model = pullback(strata_to_type, aggregate_to_type) |> apex |> rebuildStratifiedModelByFlattenSymbols;

    if return_homs
        return pullback_model, strata_to_type, aggregate_to_type
    else
        return pullback_model
    end
    
end





end