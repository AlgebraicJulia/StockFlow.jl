module Composition
export sfcompose, @compose

using ...StockFlow
using ..Syntax
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

import ..Syntax: create_foot
import Catlab.Programs.RelationalPrograms: UntypedUnnamedRelationDiagram



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

"""
Parse expression of form A ^ B => C, extract sf A and foot B => C
"""
function interpret_center_of_composition_statement(center::Expr)::Tuple{Symbol, Expr} # sf, foot defintion
    @assert length(center.args) == 3 && center.args[1] == :(=>) && typeof(center.args[2]) == Expr "Invalid argument: expected A ^ B => C, A ^ () => C or A ^ B => (), got $center"
        # third argument can be symbol or (), the latter of which is an Expr
        center_caret_statement = center.args[2]
        @assert length(center_caret_statement.args) == 3 && center_caret_statement.args[1] == :^ && typeof(center_caret_statement.args[2]) == Symbol "Invalid center argument: expected A ^ B or A ^ (), got $center"
        # third argument here, too, can be symbol or ()
        return (center_caret_statement.args[2], Expr(:call, :(=>), center_caret_statement.args[3], center.args[3]))
end

"""
Go line by line and associate stockflows and feet
"""
function interpret_composition_notation(mapping_pair::Expr)::Tuple{Vector{Symbol}, StockAndFlow0}

    if mapping_pair.head == :call # (A ^ B => C) case (incl where B or C are ())
        sf, foot_def = interpret_center_of_composition_statement(mapping_pair)
        return [sf], create_foot(foot_def)
    end

    expr_args = mapping_pair.args
    stockflows = collect(Base.Iterators.takewhile(x -> typeof(x) == Symbol, expr_args))
    center_index = length(stockflows) + 1
    @assert center_index <= length(expr_args) "A tuple is an invalid expression for composition syntax.  Expected argument of form sf1, sf2, ... ^ stock1 => sum1, stock2 => sum2, ..."
    center = expr_args[center_index]

    foot_temp = Vector{Expr}()

    sf, foot_def = interpret_center_of_composition_statement(center)
    push!(foot_temp, foot_def)
    push!(stockflows, sf)
    append!(foot_temp, expr_args[center_index+1:end])

    return (stockflows, create_foot(Expr(:tuple, foot_temp...)))
end


"""
sirv = sfcompose(sir, svi, quote
    (sr, sv)
    sr, sv ^ S => N, I => N
end)

Cannot use () => () as a foot, 
the length of the first tuple must be the same as the number of stock flows
given as argument, and every foot can only be used once.
"""
function sfcompose(sfs::Vector{K}, block::Expr; RETURN_UWD = false) where {K <: AbstractStockAndFlowF} # (sf1, sf2, ..., block)

    Base.remove_linenums!(block)
    if length(sfs) == 0
        return StockAndFlowF()
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


    # symbol representation of sf => (sf itself, sf's feet)
    sf_map::Dict{Symbol, Tuple{AbstractStockAndFlowF, Vector{StockAndFlow0}}} = 
        Dict(sf_names[i] => (sfs[i], Vector{StockAndFlow0}()) for i ∈ eachindex(sf_names)) # map the symbols to their corresponding stockflows

    # all feet
    feet_index_dict::Dict{StockAndFlow0, Int} = Dict{StockAndFlow0, Int}()
    for statement in block.args[2:end]
        stockflows, foot = interpret_composition_notation(statement)
        # adding new foot to list
        @assert (foot ∉ keys(feet_index_dict)) "Foot has already been used!"
        push!(feet_index_dict, foot => length(feet_index_dict) + 1)
        for stockflow in stockflows
            # adding this foot to each stock flow to its left
            push!(sf_map[stockflow][2], foot)
        end
    end


    sf_without_feet = Vector{Symbol}()



    Port = Vector{Tuple{Int, Int}}()

    # for each (alias, (sf, sf's feet)) 
    #   for each foot in sf's feet
    #       grab the index of a foot
    #       grab the index of alias
    #       push (alias index, foot index) to port




    for (k, v) ∈ sf_map 
        if length(v[2]) == 0
            pushfirst!(sf_without_feet, k)
        end
        for foot ∈ v[2]
            push!(Port, (findfirst(x -> x == k, sf_names), feet_index_dict[foot]))
        end
    end

    Box::Vector{Symbol} = setdiff(sf_names, sf_without_feet)


    Junction::Vector{Symbol} = [gensym() for _ ∈ 1:length(feet_index_dict)]
    OuterPort::Vector{Int} = collect(1:length(feet_index_dict))

    uwd = create_uwd(Box=Box, Port=Port, Junction=Junction, OuterPort=OuterPort)
    stockflows_with_feet = Dict(filter(((k,v),) -> length(v[2]) != 0, sf_map))
    if length(stockflows_with_feet) == 0
        composed_sf = StockAndFlowF()
    else
        open_stockflows::AbstractDict = Dict(sf_key => Open(sf_val[1], sf_val[2]...) for (sf_key, sf_val) ∈ stockflows_with_feet)
    
        composed_sf = apex(oapply(uwd, open_stockflows))
    end

    for sf in sf_without_feet
        composed_sf = apex(coproduct(composed_sf, sf_map[sf][1]))
    end

    if RETURN_UWD
        return composed_sf, uwd
    else
        return composed_sf
    end

end


"""
Compose models.
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
            sfcompose(Vector{StockAndFlowF}(), $escaped_block)
        else
            sfcompose([$(sfs...)], $escaped_block)
        end
    end
end


end