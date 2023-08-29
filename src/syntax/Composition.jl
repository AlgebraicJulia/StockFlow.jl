module Composition
export sfcompose

using ...StockFlow
import StockFlow: OpenStockAndFlowF
using ..Syntax
using MLStyle
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams


import ..Syntax: create_foot, match_foot_format, STRICT_MAPPINGS, STRICT_MATCHES, USE_ISSUB, ISSUB_DEFAULT, infer_links, substitute_symbols, iterate_stockflow_quoteblocks, DSLArgument, NothingFunction


import Catlab.Programs.RelationalPrograms: UntypedUnnamedRelationDiagram


RETURN_UWD = false

function create_uwd(;Box::Vector{Symbol} = Vector{Symbol}(), Port::Vector{Tuple{Int, Int}} = Vector{Tuple{Int, Int}}(), OuterPort::Vector{Int} = Vector{Int}(), Junction::Vector{Symbol} = Vector{Symbol}())
    uwd = UntypedUnnamedRelationDiagram{Symbol, Symbol}(0)
    add_parts!(uwd, :Box, length(Box), name=Box)
    add_parts!(uwd, :Junction, length(Junction), variable=Junction)
    add_parts!(uwd, :Port, length(Port), box=map(first, Port), junction=map(last, Port))
    add_parts!(uwd, :OuterPort, length(OuterPort), outer_junction=OuterPort)
    return uwd
end




function interpret_composition_notation(mapping_pair::Expr)::Tuple{Vector{Symbol}, StockAndFlow0}
    stockflows = Vector{Symbol}()
    foot_temp = Vector{Expr}()
    @match mapping_pair.args begin
        [:(=>), :($sf ^ $stock), sum] => begin
            # println(stock)
            # println(sum)
            push!(stockflows, sf)
            push!(foot_temp, Expr(:call, :(=>), stock, sum))
        end
        Many[
            sf::Symbol && Do(push!(stockflows, sf)) ||
            Expr(:call, :(=>), Expr(:call, :^, sf, stock), sum) && Do(push!(stockflows, sf), push!(foot_temp, Expr(:call, :(=>), stock, sum))) ||
            foot::Expr && Do(push!(foot_temp, foot)) 
        ] => true
    end

    return stockflows, create_foot(Expr(:tuple, foot_temp...))
end





function sfcompose(args...) #(sf1, sf2, ..., block)

    sfs = args[1:end-1]
    block = args[end]

    # @assert sf.head == :tuple
    # @assert all(x -> x ∈ names(Main), sf.args)
    # sfs = map(x -> getfield(Main, x), sf.args)

    # sf_map = Dict(sym => getfield(Main, sym) for sym ∈ sf.args)

    # @assert all(x -> typeof(x) <: AbstractStockAndFlowStructureF, sfs)
    # sf_tuple = sf.args

    Base.remove_linenums!(block)
    
    # @assert allunique(sf_tuple)
    sf_names::Vector{Symbol} = block.args[1].args
    # println(sf_names)
    sf_map::Dict{Symbol, AbstractStockAndFlowF} = Dict(sf_names[i] => sfs[i] for i ∈ eachindex(sf_names))
    

    foot_dict = Dict{AbstractStockAndFlowStructureF, Vector{StockAndFlow0}}(x => [] for x in sfs) # intialize sf => foot dict
    feet_index_dict::Dict{StockAndFlow0, Int} = Dict()
    for statement in block.args[2:end]
        stockflows, foot = interpret_composition_notation(statement)
        push!(feet_index_dict, foot => length(feet_index_dict) + 1)
        for stockflow in stockflows
            push!(foot_dict[sf_map[stockflow]], foot)
        end
    end

    Box::Vector{Symbol} = sf_names
    Port = [(i,feet_index_dict[j]) for (i,key) ∈ enumerate(sf_names) for j ∈ foot_dict[sf_map[key]]]
    Junction = [gensym() for _ ∈ 1:length(feet_index_dict)]
    OuterPort = collect(1:length(feet_index_dict))

    uwd = create_uwd(Box=Box, Port=Port, Junction=Junction, OuterPort=OuterPort)
    # open_stockflows::AbstractDict = Dict(sf => Open(sf, feet...) for (sf, feet) ∈ foot_dict)

    open_stockflows::AbstractDict = Dict(sf_key => Open(sf_val, foot_dict[sf_val]...) for (sf_key, sf_val) ∈ sf_map)

    # create_uwd(Box=sf_vector, Junction=[gensym() for _ ∈ 1:length(foot_vector)], Port=flattened_foot_int_indexed_vector, OuterPort=collect(1:length(foot_vector)))


    if RETURN_UWD
        return apex(oapply(uwd, open_stockflows)), uwd
    else
        return apex(oapply(uwd, open_stockflows))
    end


    

end


    # # I really fumble at the 10 yard line here
    # # this whole bit needs to be reworked.

    # flattened_foot_vector = Vector{Tuple{AbstractStockAndFlowF, StockAndFlow0}}()

    # # open_stockflows = Dict{Symbol, OpenStockAndFlowF}()
    # foot_set = Set{StockAndFlow0}()
    # for (stockflow, feet) ∈ foot_dict
    #     union!(foot_set, feet)

    #     # println(foot_dict[stockflow])
    #     # println("DSAD")
    #     # append!(foot_dict[stockflow], feet)
    #     # println(feet)

    #     # push!(open_stockflows(Open(stockflow, feet...)))
    #     append!(flattened_foot_vector, [(stockflow, foot) for foot in feet])
    #     # println([(stockflow, foot) for foot in feet])

    # end
    # # println(flattened_foot_vector)


    # open_stockflows::AbstractDict = Dict(sf => Open(sf_map[sf], feet...) for (sf, feet) ∈ foot_dict)

    # # open_stockflows::Vector{OpenStockAndFlowF} = [Open(sf_map[sf], feet...) for (sf, feet) ∈ foot_dict]

    # # println(open_stockflows[1].feet[1] == @foot A => ())


    # sf_vector::Vector{Symbol} = collect(sf_names)
    # sf_dict = Dict(stockflow => i for (i, stockflow) ∈ enumerate(sf_vector)) # terrible


    # # sf_tuple_dict = Dict{Symbol, Int}(stockflow => i for (i, stockflow) ∈ sf_vector)

    # foot_vector::Vector{StockAndFlow0} = collect(foot_set)


    # foot_dict2 = Dict(foot => i for (i, foot) ∈ enumerate(foot_vector)) # also terrible



    # # ok the hash for two feet that have the same values is the same, this might be an issue.
    # # In fact hash for any stock flow with same arguments is same
    # flattened_foot_int_indexed_vector::Vector{Tuple{Int, Int}} = [(sf_dict[stockflow_symbol], foot_dict2[foot]) for (stockflow_symbol, foot) ∈ flattened_foot_vector]



    # # foot_index_dict = Dict{StockAndFlow0, Int}(foot => i for (i, foot) ∈ enumerate(foot_vector))

    # # port_vector = [(]

    # # println(flattened_foot_vector)

    # # println(sf_vector, [gensym() for _ ∈ 1:length(foot_vector)], flattened_foot_int_indexed_vector, collect(1:length(flattened_foot_int_indexed_vector)) )

    # uwd = create_uwd(Box=sf_vector, Junction=[gensym() for _ ∈ 1:length(foot_vector)], Port=flattened_foot_int_indexed_vector, OuterPort=collect(1:length(foot_vector)))

    # # println(uwd)


end