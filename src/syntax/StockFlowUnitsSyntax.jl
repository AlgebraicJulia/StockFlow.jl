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
    @match s begin
        stock::Symbol => begin
            push!(stocks, (stock, :NONE))
        end
        :($stock:$unit) => begin
            push!(stocks, (stock, unit))
        end
        _ => return error("Unknown argument $s in stock syntax.")
    end
    # s_dict = @capture ($Stock:$units) s
    # push!(stocks, (s_dict[:Stock], s_dict[:units]))
end

function parse_param_units!(params, p::Union{Expr, Symbol})
    @match p begin
        param::Symbol => begin
            push!(params, (param, :NONE))
        end
        :($param:$unit) => begin
            push!(params, (param, unit))
        end
        _ => return error("Unknown argument $p in param syntax.")
    end
    # p_dict = @capture ($Param:$units) p
    # push!(params, (p_dict[:Param], p_dict[:units]))
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




end
