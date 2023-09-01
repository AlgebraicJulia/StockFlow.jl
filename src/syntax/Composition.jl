module Composition
export sfcompose

using ...StockFlow
using ..Syntax
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams

import ..Syntax: create_foot
import Catlab.Programs.RelationalPrograms: UntypedUnnamedRelationDiagram


RETURN_UWD = false

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
the length of the first tuple must be the same as the number of stock flows given as argument,
and every foot can only be used once.
"""
function sfcompose(args...) #(sf1, sf2, ..., block)

    @assert length(args) > 0 "Didn't get any arguments!"

    block = args[end]
    @assert typeof(block) == Expr "Didn't get an expression block for last argument!"
    
    sfs = args[1:end-1]


    @assert all(sf -> typeof(sf) <: AbstractStockAndFlowF, sfs) "Not all arguments before the block are stock flows!"


    Base.remove_linenums!(block)
    
    sf_names::Vector{Symbol} = block.args[1].args # first line are the names you want to use for the ordered arguments.
    # That is, first line needs to be a tuple, with the first argument being what you'll call the first stockflow

    if length(sfs) == 0 # Composing 0 stock flows should give you an empty stock flow
        return StockAndFlowF()
   end

   @assert length(sf_names) == length(sfs) "The number of symbols on the first line is not the same as the number of stock flow arguments provided.  Stockflow #: $(length(sfs)) Symbol #: $(length(sf_names))"



    @assert allunique(sf_names) "Not all choices of names for stock flows are unique!"


    empty_foot =  (@foot () => ())


    # symbol representation of sf => (sf itself, sf's feet)
    # Every sf has empty foot as first foot to get around being unable to create OpenStockAndFlowF without feet
    sf_map::Dict{Symbol, Tuple{AbstractStockAndFlowF, Vector{StockAndFlow0}}} = Dict(sf_names[i] => (sfs[i], [empty_foot]) for i ∈ eachindex(sf_names)) # map the symbols to their corresponding stockflows
    
    # all feet
    feet_index_dict::Dict{StockAndFlow0, Int} = Dict(empty_foot => 1)
    for statement in block.args[2:end]
        stockflows, foot = interpret_composition_notation(statement)
        # adding new foot to list
        @assert (foot ∉ keys(feet_index_dict)) "Foot has already been used, or you are using an empty foot!"
        push!(feet_index_dict, foot => length(feet_index_dict) + 1)
        for stockflow in stockflows
            # adding this foot to each stock flow to its left
            push!(sf_map[stockflow][2], foot)
        end
    end

    Box::Vector{Symbol} = sf_names


    Port = Vector{Tuple{Int, Int}}()

    for (k, v) ∈ sf_map # TODO: Just find a better way to do this.
        for foot ∈ v[2]
            push!(Port, (findfirst(x -> x == k, sf_names), feet_index_dict[foot]))
        end
    end

    Junction::Vector{Symbol} = [gensym() for _ ∈ 1:length(feet_index_dict)]
    OuterPort::Vector{Int} = collect(1:length(feet_index_dict))

    uwd = create_uwd(Box=Box, Port=Port, Junction=Junction, OuterPort=OuterPort)

    # I'd prefer this to be a vector, but oapply didn't like that
    # I'd also prefer that I don't include the empty foot, but Open doesn't want to accept stockflows with no feet.
    # open_stockflows::AbstractDict = Dict(sf_key => Open(sf_val, foot_dict[sf_val]...,) for (sf_key, sf_val) ∈ sf_map)

    open_stockflows::AbstractDict = Dict(sf_key => Open(sf_val[1], sf_val[2]...) for (sf_key, sf_val) ∈ sf_map)

    if RETURN_UWD # UWD might be a bit screwed up from the empty foot being first.
        return apex(oapply(uwd, open_stockflows)), uwd
    else
        return apex(oapply(uwd, open_stockflows))
    end

end

end