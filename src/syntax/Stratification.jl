module Stratification
export sfstratify, @stratify, @n_stratify

using ...StockFlow
using ..Syntax
using MLStyle
import Base: get
using Catlab.CategoricalAlgebra
import ..Syntax: infer_links, substitute_symbols, DSLArgument, NothingFunction, invert_vector



struct SFNames

    sf::AbstractStockAndFlowF 

    snames::Vector{Symbol} 
    svnames::Vector{Symbol}
    vnames::Vector{Symbol}
    fnames::Vector{Symbol}
    pnames::Vector{Symbol}

    # name -> index
    s::Dict{Symbol, Int}
    sv::Dict{Symbol, Int}
    v::Dict{Symbol, Int}
    f::Dict{Symbol, Int}
    p::Dict{Symbol, Int}

    # index -> new index
    ms::Dict{Int, Int}
    msv::Dict{Int, Int}
    mv::Dict{Int, Int}
    mf::Dict{Int, Int}
    mp::Dict{Int, Int}

    # index -> new index, where the first index is the actual index of the vector, the second is the int at that location
    mvs::Vector{Int}
    mvsv::Vector{Int}
    mvv::Vector{Int}
    mvf::Vector{Int}
    mvp::Vector{Int}


    SFNames(sfarg::AbstractStockAndFlowF) = (new(sfarg,
     snames(sfarg), svnames(sfarg), vnames(sfarg), fnames(sfarg), pnames(sfarg),
     invert_vector(snames(sfarg)), invert_vector(svnames(sfarg)), invert_vector(vnames(sfarg)), invert_vector(fnames(sfarg)), invert_vector(pnames(sfarg)),
     Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}(), Dict{Int, Int}(),
     Vector{Int}(),Vector{Int}(), Vector{Int}(), Vector{Int}(), Vector{Int}()))
end

function get_mappings(sfn::SFNames)::NTuple{5, Dict{Int, Int}}
    return sfn.ms, sfn.msv, sfn.mv, sfn.mf, sfn.mp
end

function get_mapped_vectors(sfn::SFNames)::NTuple{5, Vector{Int}}
    return sfn.mvs, sfn.mvsv, sfn.mvv, sfn.mvf, sfn.mvp
end

function get_mappings_infer_links_format(sfn::SFNames)::Dict{Symbol, Vector{Int}}
    Dict(:S => sfn.mvs, :SV => sfn.mvsv, :V => sfn.mvv, :F => sfn.mvf, :P => sfn.mvp)
end

function all_unique_names(sfn::SFNames)::Bool # Unnecessary, this is checked in invert_vector
    return allunique(sfn.snames) && allunique(sfn.svnames) && allunique(vnames) && allunique(fnames) && allunique(pnames)
end

function no_temp_strat_default_in_names(sfn::SFNames, temp_strat_default)::Bool
    return temp_strat_default ∉ keys(sfn.s) && temp_strat_default ∉ keys(sfn.sv) && temp_strat_default ∉ keys(sfn.v) && temp_strat_default ∉ keys(sfn.f) && temp_strat_default ∉ keys(sfn.p)
end

function add_temp_strat_default!(sfn::SFNames, temp_strat_default)
    push!(sfn.s, temp_strat_default => -1)
    push!(sfn.sv, temp_strat_default => -1)
    push!(sfn.v, temp_strat_default => -1)
    push!(sfn.f, temp_strat_default => -1)
    push!(sfn.p, temp_strat_default => -1)
end

function is_all_mapped(sfn::SFNames)::Bool
    return all(vec -> 0 ∉ vec, get_mapped_vectors(sfn))
end

function get_names(sfn::SFNames)::NTuple{5, Vector{Symbol}}
    return sfn.snames, sfn.svnames, sfn.vnames, sfn.fnames, sfn.pnames
end




"""
    interpret_stratification_standard_notation(mapping_pair::Expr)::Tuple{Vector{DSLArgument}, Vector{DSLArgument}}
Take an expression of the form a1, ..., => t <= s1, ..., where every element is a symbol, and return a 2-tuple of form ((a1 => t, a2 => t, ...), (s1 => t, ...))
"""
function interpret_stratification_standard_notation(mapping_pair::Expr)::Vector{Vector{DSLArgument}}
    @match mapping_pair begin


        :($s => $t <= $a) => return [[DSLArgument(s,t)], [DSLArgument(a,t)]]
        :($s => $t <= $a, $(atail...)) => [[DSLArgument(s,t)], [DSLArgument(a,t) ; [DSLArgument(as,t) for as in atail] ]]
        :($(shead...), $s => $t <= $a) => [[[DSLArgument(ss, t) for ss in shead] ; DSLArgument(s, t)], [DSLArgument(a, t)]]

        if mapping_pair.head == :tuple end => begin
            middle_index = findfirst(x -> typeof(x) == Expr && length(x.args) == 3, mapping_pair.args) # still isn't specific enough
            if isnothing(middle_index)
                error("Malformed line $mapping_pair, could not find center.")
            end
            @match mapping_pair.args[middle_index] begin
                :($stail => $t <= $ahead) => begin
                sdict = [[DSLArgument(ss, t) for ss in mapping_pair.args[1:middle_index-1]] ; DSLArgument(stail, t)]
                adict = [DSLArgument(ahead, t) ; [DSLArgument(as, t) for as in mapping_pair.args[middle_index+1:end]]]
                    return [sdict, adict]
                end
                _ => "Unknown format found for match; middle three values formatted incorrectly."
            end
        end
        _ => error("Unknown line format found in stratification notation.") 
    end
end



function interpret_stratification_generalized_notation(mapping_pair::Expr)::Vector{Vector{DSLArgument}}
    # TODO: add a crap ton of assert statements.
    # We're assuming it's gonna take the form [(a1, ..., an), ..., (k1, ..., km)] => t1

    @assert length(mapping_pair.args) == 3 && typeof(mapping_pair.args[3]) == Symbol && length(mapping_pair.args[2]) && typeof(mapping_pair.args[2].args) == Vector

    # TODO: Include assert that length(mapping_pair.args[2].args) is the same as the number of stockflows
    # Maybe include an assert that each element of the vector contains only symbols

    other = mapping_pair.args[2].args # needs to be a vector of tuples of symbols
    type = mapping_pair.args[3] # needs to be a symbol
    return [[DSLArgument(sym, type) for sym ∈ tup.args] for tup ∈ other]
end





"""
Gets mapping information from each line and updates dictionaries.  If a symbol already has a mapping and another is found, keep the first, or throw an error if strict_matches = true.
"""
function read_stratification_line_and_update_dictionaries!(line::Expr, other_names::Vector{Dict{Symbol, Int}}, type_names::Dict{Symbol, Int}, other_mappings::Vector{Dict{Int, Int}} ; use_standard_stratification_syntax = true, strict_matches = false, use_flags = true)
    if use_standard_stratification_syntax
        interpret_stratification_notation_function = interpret_stratification_standard_notation
    else
        interpret_stratification_notation_function = interpret_stratification_generalized_notation
    end
    
    current_symbol_dict::Vector{Vector{DSLArgument}} = interpret_stratification_notation_function(line)

    current_mapping_dict::Vector{Dict{Int, Int}} = ((x, y) -> substitute_symbols(x,type_names, y; use_flags=use_flags)).(other_names, current_symbol_dict)

    ((cumulative_dict, new_dict) -> mergewith!((cv, nv) -> cv, cumulative_dict, new_dict)).(other_mappings, current_mapping_dict)

end

"""
Print all symbols such that the corresponding int is 0, representing an unmapped object.
"""
function print_unmapped(SFNames, name="STOCKFLOW")
    for (indices, names) ∈ zip(SFNames.get_mapped_vectors, SFNames.get_names)
        for (i, val) ∈ enumerate(indices)
            if val == 0
                println("UNMAPPED IN $(name):")
                println(names[i])
            end
        end
    end
end

"""
Iterates over each line in a stratification syntax block and updates the appropriate dictionaries.
"""
function iterate_over_stratification_lines!(block, other_names::Vector{SFNames}, type_names::SFNames; use_standard_stratification_syntax=true, strict_matches=false, use_flags=true)
    
    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> read_stratification_line_and_update_dictionaries!(s, getfield.(other_names, :s), type_names.s, getfield.(other_names, :ms); use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)
            end
            QuoteNode(:sums) => begin
                current_phase = sv -> read_stratification_line_and_update_dictionaries!(sv, getfield.(other_names, :sv), type_names.sv, getfield.(other_names, :msv); use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)
            end        
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> read_stratification_line_and_update_dictionaries!(v, getfield.(other_names, :v), type_names.v, getfield.(other_names, :mv); use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)
            end       
            QuoteNode(:flows) => begin
                current_phase = f -> read_stratification_line_and_update_dictionaries!(f,getfield.(other_names, :f), type_names.f, getfield.(other_names, :mf); use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)
            end                    
            QuoteNode(:parameters) => begin
                current_phase = p -> read_stratification_line_and_update_dictionaries!(p, getfield.(other_names, :p), type_names.p, getfield.(other_names, :mp); use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)
            end
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end
end

"""
Apply default mappings, infer mapping if there's only a single option, and convert from Dict{Int, Int} to Vector{Int}
"""
function complete_mappings!(sfm::SFNames, sftype::AbstractStockAndFlowF; strict_mappings = false)
    # get the default value, if it has been assigned.  Use 0 if it hasn't.
    all_index_mappings = get_mappings(sfm)

    default_index_stock = get(all_index_mappings[1], -1, 0)
    default_index_sum = get(all_index_mappings[2], -1, 0)
    default_index_dyvar = get(all_index_mappings[3], -1, 0)
    default_index_flow = get(all_index_mappings[4], -1, 0)
    default_index_param = get(all_index_mappings[5], -1, 0)



    # STEP 3
    if !strict_mappings
        one_type_stock = length(snames(sftype)) == 1 ? 1 : 0 # if there is only one stock, it needs to have index 1
        one_type_flow = length(fnames(sftype)) == 1 ? 1 : 0
        one_type_dyvar = length(vnames(sftype)) == 1 ? 1 : 0
        one_type_param = length(pnames(sftype)) == 1 ? 1 : 0
        one_type_sum = length(svnames(sftype)) == 1 ? 1 : 0
    else
        one_type_stock = one_type_flow = one_type_dyvar = one_type_param = one_type_sum = 0
    end

    # Convert back to vectors.  If you find a zero, check if there's a default and use that.  If there isn't a default, check if there's only one option and use that.
    # Otherwise, there's an unassigned value which can't be inferred.
    # Taking max because it's less verbose than ternary and accomplishes the same thing:
    # - If both default_index and one_type are mapped, they must be mapped to the same thing, because one_type being mapped implies there's only one option.
    # - If only one_type is mapped, then it will be positive, and default_infex will be 0
    # - If only default_index is mapped, it will be positive and one_type will be 0

    stock_mappings::Vector{Int} = [get(all_index_mappings[1], i, max(default_index_stock, one_type_stock))  for i ∈ eachindex(sfm.snames)]
    sum_mappings::Vector{Int} =  [get(all_index_mappings[2], i,  max(default_index_sum, one_type_sum))  for i ∈ eachindex(sfm.svnames)]
    dyvar_mappings::Vector{Int} = [get(all_index_mappings[3], i, max(default_index_dyvar, one_type_dyvar))  for i ∈ eachindex(sfm.vnames)]
    flow_mappings::Vector{Int} =  [get(all_index_mappings[4], i, max(default_index_flow, one_type_flow))  for i ∈ eachindex(sfm.fnames)]
    param_mappings::Vector{Int} = [get(all_index_mappings[5], i,  max(default_index_param, one_type_param))  for i ∈ eachindex(sfm.pnames)]

    append!(sfm.mvs, stock_mappings)
    append!(sfm.mvsv, sum_mappings)
    append!(sfm.mvv, dyvar_mappings)
    append!(sfm.mvf, flow_mappings)
    append!(sfm.mvp, param_mappings)

end


"""
    sfstratify(strata, type, aggregate, block ; kwargs)

    1. Grab all names from strata, type and aggregate, and create dictionaries which map them to their indices
    2. Iterate over each line in the block
        2a. Split each line into a dictionary which maps all strata to that type and all aggregate to that type
        2b. Convert from two Symbol => Symbol dictionaries to two Int => Int dictionaries, using the dictionaries from step 1
            2bα. If applicable, for symbols with ~ as a prefix, find all symbols with matching substrings in the symbol dictionaries, and map all those
        2c. Accumulate respective dictionaries (optionally, only allow first match vs throw an error (strict_matches = false vs true))
    3. Create an array of 0s for stocks, flows, parameters, dyvars and sums for strata and aggregate.   Insert into arrays all values from the two Int => Int dictionaries
        3a. If strict_mappings = false, if there only exists one option in type to map to, and it hasn't been explicitly specified, add it.  If strict_mappings = true and it hasn't been specified, throw an error.
    4. Do a once-over of arrays and ensure there aren't any zeroes (unmapped values) remaining (helps with debugging when you screw up stratifying)
    5. Deal with attributes (create a copy of type sf with attributes mapped to nothing)
    6. Infer LS, LSV, etc.
    7. Construct strata -> type and aggregate -> type ACSetTransformations
    8. Return pullback (with flattened attributes)
"""
function sfstratify(others::Vector{K}, type::K, block::Expr ; use_standard_stratification_syntax = true, strict_mappings = false, strict_matches = false, temp_strat_default = :_, use_temp_strat_default = true, use_flags = true, return_homs = false) where {K <: AbstractStockAndFlowStructureF}

    Base.remove_linenums!(block)

    # STEP 1


    other_names::Vector{SFNames} = [SFNames(sf) for sf ∈ others]
    type_names::SFNames = SFNames(type) # has some unnecessary fields.

    if use_temp_strat_default
        # Applies function to every element in vector.
        @assert all((sfn -> no_temp_strat_default_in_names(sfn, temp_strat_default)).(other_names)) && no_temp_strat_default_in_names(type_names, temp_strat_default) "A stockflow contains $(temp_strat_default) !  Please change temp_strat_default to a different symbol or rename offending object."
        (sfn -> add_temp_strat_default!(sfn, temp_strat_default)).(other_names)
    end


    # STEP 2
    iterate_over_stratification_lines!(block, other_names, type_names ; use_standard_stratification_syntax=use_standard_stratification_syntax, strict_matches=strict_matches, use_flags=use_flags)


    (sfn -> complete_mappings!(sfn, type ; strict_mappings=strict_mappings)).(other_names)

    # STEP 4

    # This bit makes debugging when making a stratification easier.  Tells you exactly which ones you forgot to map.

    #unmapped: 
    if !(all(is_all_mapped.(other_names)))
        for i ∈ eachindex(other_names)
            print_unmapped(other_names[i], "STOCKFLOW $i")
        end
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
    
    generate_all_mappings_function = m -> Dict(infer_links(m.sf, type, get_mappings_infer_links_format(m))..., get_mappings_infer_links_format(m)..., :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)
    all_mappings = generate_all_mappings_function.(other_names)
   
    all_transformations = [ACSetTransformation(sfn.sf, no_attribute_type ; mappings...) for (sfn, mappings) ∈ zip(other_names, all_mappings)]

    # STEP 8

    pullback_model = pullback(all_transformations...) |> apex |> rebuildStratifiedModelByFlattenSymbols;
  


    if return_homs
        return pullback_model, all_transformations
    else
        return pullback_model
    end
    
end

"""
    stratify(strata, type, aggregate, block)
Take three stockflows and a block describing where the first and third map on to the second, and get a new stratified stockflow.
Left side are strata objects, middle are type, right are aggregate.  Each strata and aggregate object is associated with one type object.
The resultant stockflow contains objects which are the product of strata and aggregate objects which map to the same type object.
Use _ to match all objects in that category, ~ as a prefix to match all objects with the following string as a substring.  Objects always go with their first match.
If the type model has a single object in a category, the mapping to it is automatically assumed.  In the below example, we wouldn't need to specify :stocks or :sums.

```julia

@stratify WeightModel l_type ageWeightModel begin
    :stocks
    _ => pop <= _
    
    :flows
    ~Death => f_death <= ~Death
    ~id => f_aging <= ~aging
    ~Becoming => f_fstOrder <= ~id
    _ => f_birth <= f_NB

    
    :dynamic_variables
    v_NewBorn => v_birth <= v_NB
    ~Death => v_death <= ~Death
    ~id  => v_aging <= v_agingCA, v_agingAS
    v_BecomingOverWeight, v_BecomingObese => v_fstOrder <= v_idC, v_idA, v_idS
    
    :parameters
    μ => μ <= μ
    δw, δo => δ <= δC, δA, δS
    rw, ro => rFstOrder <= r
    rage => rage <= rageCA, rageAS
    
    :sums
    N => N <= N
    
end 
```
"""
macro stratify(strata, type, aggregate, block)
    escaped_block = Expr(:quote, block)
    quote
        sfstratify([$(esc(strata)), $(esc(aggregate))], $(esc(type)), $(esc(escaped_block)))
    end
end

macro n_stratify(args...)
    if length(args) < 2
        return :(MethodError("Too few arguments provided!  Please provide some number of stockflows, then the type stock flow, then a quote block."))
    else
        # TODO: Potentially Expr(:quote, args[end])
        if length(args) == 2
            return Expr(:call, :sfnstratify, :(Vector{AbstractStockAndFlowF}()), esc(args[1]), esc(args[2]))
        else 
            return Expr(:call, :sfnstratify, Expr(:vect, esc.(args[1:end-2])...), esc(args[end-1]), esc(args[end]))
        end
    end
end


end