module Stratify
export @stratify

using ...StockFlow
using ..Syntax
using MLStyle
import Base.get
using Catlab.CategoricalAlgebra
import ..Syntax: STRICT_MAPPINGS, STRICT_MATCHES, USE_ISSUB, ISSUB_DEFAULT, infer_links, substitute_symbols


RETURN_HOMS = false


"""
Take an expression of the form a1, ..., => t <= s1, ..., where every element is a symbol, and return a 2-tuple of dictionaries of form ((a1 => t, a2 => t, ...), (s1 => t, ...))
"""
function interpret_stratification_notation(mapping_pair::Expr)
    @match mapping_pair begin


        :($s => $t <= $a) => return (Dict(s => t), Dict(a => t))
        :($s => $t <= $a, $(atail...)) => return (Dict(s => t), push!(Dict(as => t for as in atail), a => t))
        :($(shead...), $s => $t <= $a) => return (push!(Dict(ss => t for ss in shead), s => t), Dict(a => t))

        if mapping_pair.head == :tuple end => begin
            middle_index = findfirst(x -> typeof(x) == Expr, mapping_pair.args)
            @match mapping_pair.args[middle_index] begin
                :($stail => $t <= $ahead) => begin
                    sdict = push!(Dict(ss => t for ss in mapping_pair.args[1:middle_index-1]), stail => t)
                    adict =  push!(Dict(ss => t for as in mapping_pair.args[middle_index:end]), ahead => t)
                    return (sdict, adict)
                end
                _ => "Unknown format found for match; middle three values formatted incorrectly."
            end
        end
        _ => error("Unknown line format found in stratification notation.") 
    end
end




function read_stratification_line_and_update_dictionaries!(line::Expr, strata_names::Dict{Symbol, Int}, type_names::Dict{Symbol, Int}, aggregate_names::Dict{Symbol, Int}, strata_mappings::Dict{Int, Int}, aggregate_mappings::Dict{Int, Int})
    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(line)

    current_strata_dict = substitute_symbols(strata_names, type_names, current_strata_symbol_dict ; use_substr_prefix=USE_ISSUB, issubstr_prefix=ISSUB_DEFAULT)
    current_aggregate_dict = substitute_symbols(aggregate_names, type_names, current_aggregate_symbol_dict ; use_substr_prefix=USE_ISSUB, issubstr_prefix=ISSUB_DEFAULT)
    # current_strata_dict, current_aggregate_dict = substitute_symbols_stratification(strata_names, type_names, aggregate_names, current_strata_symbol_dict, current_aggregate_symbol_dict ; use_substr_prefix=USE_ISSUB, issubstr_prefix=ISSUB_DEFAULT)
    
    if STRICT_MATCHES
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
    @stratify (strata, type, aggregate) begin ... end

    Ok, so the general idea here is:
    1. Grab all names from strata, type and aggregate, and create dictionaries which map them to their indices
    2. iterate over each line in the block
        2a. Split each line into a dictionary which maps all strata to that type and all aggregate to that type
        2b. Convert from two Symbol => Symbol dictionaries to two Int => Int dictionaries, using the dictionaries from step 1
            2bα. Map symbols without ISSUB_DEFAULT to Int
            2bβ. If applicable, for symbols with ISSUB_DEFAULT as a prefix, find all matching symbols in the symbol dictionaries, and map all those
        2c. Accumulate respective dictionaries (optionally, only allow first match vs throw an error (STRICT_MATCHES = false vs true))
    3. Create an array of 0s for stocks, flows, parameters, dyvars and sums for strata and aggregate.   Insert into arrays all values from the two Int => Int dictionaries
        3a. If STRICT_MAPPINGS = false, if there only exists one option in type to map to, and it hasn't been explicitly specified, add it.  If STRICT_MAPPINGS = true and it hasn't been specified, throw an error.
    4. Do a once-over of arrays and ensure there aren't any zeroes (unmapped values) remaining
    5. Deal with attributes (create a copy of type sf with attributes mapped to nothing)
    6. Infer LS, LSV, etc.
    7. Construct strata -> type and aggregate -> type ACSetTransformations (Maybe not ACSet, because we don't care about attributes)
    8. Return pullback (with flattened attributes)
"""
macro stratify(sf, block) # Trying to be very vigilant about catching errors.


    @assert sf.head == :tuple && length(sf.args) == 3 "Incorrect arguments!  Expected tuple of length three, got: $(sf)"
    @assert all(x -> x ∈ names(Main), sf.args) "Stockflows not in Main namespace!  Did you forget to define them?  If doing it at runtime, need an @eval call on stratify!"
    
    # Maybe want it to be not just in Main?

    strata, type, aggregate = map(x -> getfield(Main, x), sf.args) # Attempting to be clever and get around an eval call
    # note, because of this, you're probably gonna want to run this in a jupyter notebook, or use an @eval call.  It won't work at normal runtime.



    @assert all(x -> typeof(x) <: AbstractStockAndFlowStructureF, [strata, type, aggregate]) # Unsure if we want to be more or less strict with the type check
    # they all need parameter fields right now.  Could probably generalize.

    Base.remove_linenums!(block)

    # STEP 1
       
    strata_snames = Dict(S => i for (i, S) in enumerate(snames(strata)))
    strata_svnames = Dict(S => i for (i, S) in enumerate(svnames(strata)))
    strata_vnames = Dict(S => i for (i, S) in enumerate(vnames(strata)))
    strata_fnames = Dict(S => i for (i, S) in enumerate(fnames(strata)))
    strata_pnames = Dict(S => i for (i, S) in enumerate(pnames(strata)))

    type_snames = Dict(S => i for (i, S) in enumerate(snames(type)))
    type_svnames = Dict(S => i for (i, S) in enumerate(svnames(type)))
    type_vnames = Dict(S => i for (i, S) in enumerate(vnames(type)))
    type_fnames = Dict(S => i for (i, S) in enumerate(fnames(type)))
    type_pnames = Dict(S => i for (i, S) in enumerate(pnames(type)))

    aggregate_snames = Dict(S => i for (i, S) in enumerate(snames(aggregate)))
    aggregate_svnames = Dict(S => i for (i, S) in enumerate(svnames(aggregate)))
    aggregate_vnames = Dict(S => i for (i, S) in enumerate(vnames(aggregate)))
    aggregate_fnames = Dict(S => i for (i, S) in enumerate(fnames(aggregate)))
    aggregate_pnames = Dict(S => i for (i, S) in enumerate(pnames(aggregate)))


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
    

    strata_all_names = [snames(strata), svnames(strata), vnames(strata), fnames(strata), pnames(strata)]
    type_all_names = [snames(type), svnames(type), vnames(type), fnames(type), pnames(type)]
    aggregate_all_names = [snames(aggregate), svnames(aggregate), vnames(aggregate), fnames(aggregate), pnames(aggregate)]

    @assert all(map(x -> allunique(x), strata_all_names)) "Not all names in strata model are unique!"
    @assert all(map(x -> allunique(x), type_all_names)) "Not all names in type model are unique!"
    @assert all(map(x -> allunique(x), aggregate_all_names)) "Not all names in aggregate model are unique!"
    # ensure names in each stockflow are unique within themselves.
    # there's the possibility this would be called after mapping to nothing, so if you get an error here, check that the names aren't all nothing


    all_names =  [strata_all_names..., type_all_names..., aggregate_all_names...]

    if USE_ISSUB
        @assert all(ks -> all(k -> !startswith(string(k), ISSUB_DEFAULT), values(ks)), all_names) "A name starts with the prefix used to identify substrings!"
        # ensure that no name in stockflows start with prefix used to identify substrings
    end


    # STEP 2
    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> read_stratification_line_and_update_dictionaries!(s, strata_snames, type_snames, aggregate_snames, strata_stock_mappings_dict, aggregate_stock_mappings_dict)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> read_stratification_line_and_update_dictionaries!(p, strata_pnames, type_pnames, aggregate_pnames, strata_param_mappings_dict, aggregate_param_mappings_dict)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> read_stratification_line_and_update_dictionaries!(v, strata_vnames, type_vnames, aggregate_vnames, strata_dyvar_mappings_dict, aggregate_dyvar_mappings_dict)
            end            
            QuoteNode(:flows) => begin
                current_phase = f -> read_stratification_line_and_update_dictionaries!(f, strata_fnames, type_fnames, aggregate_fnames, strata_flow_mappings_dict, aggregate_flow_mappings_dict)
            end                    
                  
            QuoteNode(:sums) => begin
                current_phase = sv -> read_stratification_line_and_update_dictionaries!(sv, strata_svnames, type_svnames, aggregate_svnames, strata_sum_mappings_dict, aggregate_sum_mappings_dict)
            end                    
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    # STEP 3
    if !STRICT_MAPPINGS
        one_type_stock = length(snames(type)) == 1 ? 1 : 0 # if there is only one stock, it needs to have index 1
        one_type_flow = length(fnames(type)) == 1 ? 1 : 0
        one_type_dyvar = length(vnames(type)) == 1 ? 1 : 0
        one_type_param = length(pnames(type)) == 1 ? 1 : 0
        one_type_sum = length(svnames(type)) == 1 ? 1 : 0
    else
        one_type_stock = one_type_flow = one_type_dyvar = one_type_param = one_type_sum = 0
    end

    strata_stock_mappings::Vector{Int} = [get(strata_stock_mappings_dict, i, one_type_stock)  for i in 1:ns(strata)] # if key isn't in dictionary, returns one_type_stock
    strata_flow_mappings::Vector{Int} =  [get(strata_flow_mappings_dict, i, one_type_flow)  for i in 1:nf(strata)]
    strata_dyvar_mappings::Vector{Int} = [get(strata_dyvar_mappings_dict, i, one_type_dyvar)  for i in 1:nvb(strata)]
    strata_param_mappings::Vector{Int} = [get(strata_param_mappings_dict, i, one_type_param)  for i in 1:np(strata)]
    strata_sum_mappings::Vector{Int} =  [get(strata_sum_mappings_dict, i, one_type_sum)  for i in 1:nsv(strata)]

    aggregate_stock_mappings::Vector{Int} = [get(aggregate_stock_mappings_dict, i, one_type_stock)  for i in 1:ns(aggregate)]
    aggregate_flow_mappings::Vector{Int} =  [get(aggregate_flow_mappings_dict, i, one_type_flow)  for i in 1:nf(aggregate)]
    aggregate_dyvar_mappings::Vector{Int} = [get(aggregate_dyvar_mappings_dict, i, one_type_dyvar)  for i in 1:nvb(aggregate)]
    aggregate_param_mappings::Vector{Int} = [get(aggregate_param_mappings_dict, i, one_type_param)  for i in 1:np(aggregate)]
    aggregate_sum_mappings::Vector{Int} =  [get(aggregate_sum_mappings_dict, i, one_type_sum)  for i in 1:nsv(aggregate)]
    
    all_mappings = [strata_stock_mappings..., strata_flow_mappings..., strata_dyvar_mappings..., strata_param_mappings..., strata_sum_mappings..., aggregate_stock_mappings..., aggregate_flow_mappings..., aggregate_dyvar_mappings..., aggregate_param_mappings..., aggregate_sum_mappings...]
    
    strata_mappings = [strata_stock_mappings => snames(strata), strata_flow_mappings => fnames(strata), strata_dyvar_mappings => vnames(strata), strata_param_mappings => pnames(strata), strata_sum_mappings => svnames(strata)]
    aggregate_mappings = [aggregate_stock_mappings => snames(aggregate), aggregate_flow_mappings => fnames(aggregate), aggregate_dyvar_mappings => vnames(aggregate), aggregate_param_mappings => pnames(aggregate), aggregate_sum_mappings => svnames(aggregate)]

    # STEP 4

    #unmapped: 
    if !(all(x -> x != 0, all_mappings))
        for (ints, dicts) in strata_mappings
            for (i, val) in enumerate(ints)
                if val == 0
                    println("UNMAPPED IN STRATA:")
                    println(dicts[i])
                end
            end
        end
        for (ints, dicts) in aggregate_mappings
            for (i, val) in enumerate(ints)
                if val == 0
                    println("UNMAPPED IN AGGREGATE:")
                    println(dicts[i])
                end
            end
        end
        error("There is an unmapped value!")
    end
    

    # STEP 5
    nothing_function = x -> nothing
    no_attribute_type = map(type, Name=name->nothing, Op=op->nothing, Position=pos->nothing)

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
    strata_to_type = ACSetTransformation(strata, no_attribute_type; strata_necmaps..., strata_inferred_links..., Op = nothing_function, Position = nothing_function, Name = nothing_function)

    
    aggregate_necmaps = Dict(:S => aggregate_stock_mappings, :F => aggregate_flow_mappings, :V => aggregate_dyvar_mappings, :P => aggregate_param_mappings, :SV => aggregate_sum_mappings)
    aggregate_inferred_links = infer_links(aggregate, type, aggregate_necmaps)
    aggregate_to_type = ACSetTransformation(aggregate, no_attribute_type; aggregate_necmaps..., aggregate_inferred_links..., Op = nothing_function, Position = nothing_function, Name = nothing_function)

    

    # STEP 8
    pullback_model = pullback(strata_to_type, aggregate_to_type) |> apex |> rebuildStratifiedModelByFlattenSymbols;

    if RETURN_HOMS
        return pullback_model, strata_to_type, aggregate_to_type
    else
        return pullback_model
    end
    
end

end