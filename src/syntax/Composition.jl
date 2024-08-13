module Composition
export sfcompose, @compose

using ...StockFlow
using ..Syntax
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

import ..Syntax: create_foot
import Catlab.Programs.RelationalPrograms: UntypedUnnamedRelationDiagram

using MLStyle


"""
Construct a uwd to compose your open stockflows
"""
function create_uwd(;
    Box::Vector{Symbol} = Vector{Symbol}(), # stockflows
    Port::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}(), # stockflow => foot number, for each foot on stockflow
    OuterPort::Vector{Int} = Vector{Int}(),  # unique feet number (1:n)
    Junction::Vector{Symbol} = Vector{Symbol}() # A symbol for each (unique) foot
    )

    uwd = UntypedUnnamedRelationDiagram{Symbol, Symbol}(0)
    add_parts!(uwd, :Box, length(Box), name=Box)
    add_parts!(uwd, :Junction, length(Junction), variable=Junction)
    add_parts!(uwd, :Port, length(Port), box=map(first, Port), junction=map(last, Port))
    add_parts!(uwd, :OuterPort, length(OuterPort), outer_junction=OuterPort)
    return uwd
end



function interpret_composition_notation(statement::Expr, create_foot_function)
    sf_symbols = nothing
    foot_expressions = nothing

    @match statement begin
        :($sfs ^ $foot) => begin
            sf_symbols = sfs
            foot_expressions = foot
        end
        :($sfs ^ $stock => $sum) => begin
            sf_symbols = sfs
            foot_expressions = Expr(:call, :(=>), stock, sum)
        end

        :($sfs ^ $stock1 => $sum1, $(foot...)) => begin
            sf_symbols = sfs
            foot_expressions = :($stock1 => $sum1, $(foot...))
        end
        _ => error("Unknown composition syntax: $statement")
    end

    sf_foot = create_foot_function(foot_expressions)
    sf_vector = nothing

    @match sf_symbols begin
        ::Symbol => begin sf_vector = [sf_symbols] end
        ::Expr => begin sf_vector = collect(sf_symbols.args) end
        _ => error("Unknown composition syntax in $K arguments: $sf_symbols")
    end

    return sf_vector, sf_foot
end



"""
sirv = sfcompose(sir, svi, quote
    (sr, sv)
    (sr, sv) ^ S => N, I => N
end)

Cannot use () => () as a foot, 
the length of the first tuple must be the same as the number of stock flows
given as argument, and every foot can only be used once.
"""
function sfcompose(sfs::Vector, block::Expr, main_type, foot_type, create_foot_function; RETURN_UWD = false) # (sf1, sf2, ..., block)

    Base.remove_linenums!(block)
    if length(sfs) == 0
        return main_type()
    end

    if length(block.args) == 0 # equivalent to coproduct
        # If sf_names is empty, there can't be any arguments after it.
        sf_names = [gensym() for _ in 1:length(sfs)]
    else
        sf_names = block.args[1].args
    end
  
    # If the tuple of aliases is smaller than 
    if length(sf_names) < length(sfs)
        len_diff = length(sfs) - length(sf_names)
        append!(sf_names, [gensym() for _ in 1:len_diff])
    end


    @assert allunique(sf_names) "Not all choices of names for stock flows \
    are unique!"

    empty_foot =  foot_type()


    # symbol representation of sf => (sf itself, sf's feet)
    # Every sf has empty foot as first foot to get around being unable to create OpenStockAndFlowF without feet
    sf_map::Dict{Symbol, Tuple{main_type, Vector{foot_type}}} = 
        Dict(sf_names[i] => (sfs[i], [empty_foot]) for i ∈ eachindex(sf_names)) # map the symbols to their corresponding stockflows

    # all feet
    feet_index_dict::Dict{foot_type, Int} = Dict(empty_foot => 1)
    for statement in block.args[2:end]
        stockflows, foot = interpret_composition_notation(statement, create_foot_function)
        # adding new foot to list
        @assert (foot ∉ keys(feet_index_dict)) "Foot has already been used,\
         or you are using an empty foot!"
        push!(feet_index_dict, foot => length(feet_index_dict) + 1)
        for stockflow in stockflows
            # adding this foot to each stock flow to its left
            push!(sf_map[stockflow][2], foot)
        end
    end

    Box::Vector{Symbol} = sf_names


    Port = Vector{Tuple{Int, Int}}()

    # for each (alias, (sf, sf's feet)) 
    #   for each foot in sf's feet
    #       grab the index of a foot
    #       grab the index of alias
    #       push (alias index, foot index) to port




    for (k, v) ∈ sf_map # TODO: Just find a better way to do this.
        for foot ∈ v[2]
            push!(Port, (findfirst(x -> x == k, sf_names), feet_index_dict[foot]))
        end
    end

    Junction::Vector{Symbol} = [gensym() for _ ∈ 1:length(feet_index_dict)]
    OuterPort::Vector{Int} = collect(1:length(feet_index_dict))

    uwd = create_uwd(Box=Box, Port=Port, Junction=Junction, OuterPort=OuterPort)

    open_stockflows::AbstractDict = Dict(sf_key => Open(sf_val[1], sf_val[2]...) for (sf_key, sf_val) ∈ sf_map)

    if RETURN_UWD
        return apex(oapply(uwd, open_stockflows)), uwd
    else
        return apex(oapply(uwd, open_stockflows))
    end

end


"""
Compose models.  Works with Stockflow and Causal Loop, but not together.

Cannot use () => () as a foot.

```julia
# Stockflow
sirv = @compose sir svi begin
    (sr, sv)
    (sr, sv) ^ S => N, I => N
end

# Causal Loop with polarities
ABCD = @compose ABC BCD begin
    (ABC, BCD)
    (ABC, BCD) ^ B => -C
end

# Causal Loop without polarities
ABCD2 = @compose ABC2 BCD2 begin
    (ABC, BCD)
    (ABC, BCD) ^ B => C
end
```

"""
macro compose(args...)
    if length(args) == 0
        return :(MethodError("No arguments provided!  Please provide some \
        number of stockflows, then a quote block."))
    end
    escaped_block = Expr(:quote, args[end])
    sfs = esc.(args[1:end-1])
    quote
        if length($sfs) == 0
            :(MethodError("Could not infer type given no arguments.  Please provide at \
            least one stockflow or causal loop diagram"))
        else
            local model_type = typeof($(sfs[1]))
            if !(all(x -> typeof(x) == model_type, $(sfs)))
                :(MethodError("Arguments to composition have inconsistent types."))
            end

            if model_type <: AbstractStockAndFlowF
                sfcompose([$(sfs...)], $escaped_block, StockAndFlowF, StockAndFlow0, create_foot)
            elseif model_type <: CausalLoopPol
                sfcompose([$(sfs...)], $escaped_block, CausalLoopPol, CausalLoopPol, x -> to_clp(cl_macro(x)))
            elseif model_type <: CausalLoopPM
                sfcompose([$(sfs...)], $escaped_block, CausalLoopPM, CausalLoopPM, cl_macro)
            elseif model_type <: CausalLoop
                sfcompose([$(sfs...)], $escaped_block, CausalLoop, CausalLoop, (x -> (cl_macro(x, true))))
            else
                :(MethodError("Invalid type $(model_type) for composition syntax."))
            end
        end
    end
end



end