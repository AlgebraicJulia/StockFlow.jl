module StockFlowUnitsSyntax

using ...StockFlow.StockFlowUnits

import ...StockFlow: StockAndFlowF, state_dict, ns, np

using ..Syntax
import ..Syntax: StockAndFlowBlock, parse_dyvar!, parse_flow!, parse_sum!, stock_and_flow_syntax_to_arguments, get

using MLStyle
using MLStyle.Modules.AST

using Catlab.CategoricalAlgebra


export @stock_and_flow_U


function parse_stock_units!(stocks, s::Union{Expr, Symbol})
    s_dict = @capture ($Stock:$units) s
    push!(stocks, (s_dict[:Stock], s_dict[:units]))
end

function parse_param_units!(params, p::Union{Expr, Symbol})
    p_dict = @capture ($Param:$units) p
    push!(params, (p_dict[:Param], p_dict[:units]))
end












function parse_stock_and_flow_units_syntax(statements::Vector{Any})
    stocks::Vector{Tuple{Symbol, Union{Symbol, Expr}}} = []
    params::Vector{Tuple{Symbol, Union{Symbol, Expr}}} = []
    dyvars::Vector{Tuple{Symbol,Expr}} = []
    flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
    sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []
    current_phase = (_, _) -> ()
    for statement in statements
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> parse_stock_units!(stocks, s)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> parse_param_units!(params, p)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = d -> parse_dyvar!(dyvars, d)
            end
            QuoteNode(:flows) => begin
                current_phase = f -> parse_flow!(flows, f)
            end
            QuoteNode(:sums) => begin
                current_phase = s -> parse_sum!(sums, s)
            end
            QuoteNode(kw) =>
                error("Unknown block type for Stock and Flow syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    s = StockAndFlowBlock(map(((s, u),) -> s,  stocks), map(((p, u),) -> p,  params), dyvars, flows, sums)
    sunits = map(((s, u),) -> u,  stocks)
    punits = map(((p, u),) -> u,  params)
    return s, sunits, punits
end


macro stock_and_flow_U(block)
    Base.remove_linenums!(block)
    block_args = block.args
    return quote
      local syntax_lines, sunits, punits = parse_stock_and_flow_units_syntax($block_args)
      local saff_args = stock_and_flow_syntax_to_arguments(syntax_lines)
      local sff = StockAndFlowF(
          saff_args.stocks,
          saff_args.params,
          map(kv -> kv.first => get(kv.second), saff_args.dyvars),
          saff_args.flows,
          saff_args.sums,
      )
      FtoU(sff, sunits, punits)
      #add_units_sff!(sff, sunits, punits)
      
      

      
    end
end


function add_units_sff!(sf, sunits, punits)




    cunits = collect(Set(sunits) ∪ Set(punits))
    cunit_dict = state_dict(cunits)

      

    cunit_exp_dict = Dict{Union{Symbol, Expr}, Dict{Symbol, Float64}}()
    unit_set = Set{Symbol}()
    
    for cu in cunits
        cunit_exp_dict[cu] = extract_exponents(cu)
        union!(unit_set, keys(cunit_exp_dict[cu]))
    end

    
    units = collect(unit_set)
    
    unit_dict = state_dict(units)
    

    
    add_units!(sf, length(units), uname=([Symbol(x) for x in units]))
    add_cunits!(sf, length(cunits), cuname=([Symbol(x) for x in cunits]))
    
    for cu in cunits
        for (u, exp) in pairs(cunit_exp_dict[cu])
            add_UCUlink!(sf, unit_dict[u], cunit_dict[cu], exp=exp)
        end
    end
    


    scu = collect(map(x -> cunit_dict[x], sunits))
    pcu = collect(map(x -> cunit_dict[x], punits))
    
    #set_subpart!(sf, :scu, scu)
    #set_subpart!(sf, :pcu, pcu)

    for i in 1:ns(sf)
        sf.subparts.scu[i] = scu[i]
    end
    for i in 1:np(sf)
        sf.subparts.pcu[i] = pcu[i]
    end
    #set_scu!(sf, scu)
    #set_pcu!(sf, pcu)
        
        
end
        







end
