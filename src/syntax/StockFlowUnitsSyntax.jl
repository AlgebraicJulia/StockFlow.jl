
using ..StockFlow.StockFlowUnits

using MLStyle.Modules.AST

export @stock_and_flow_U, @foot_U


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
    end
end


function match_foot_format_U(footblock::Union{Expr, Symbol})
    @match footblock begin
        ::Symbol => ((), (), (), (), footblock)

        :(()              => ())              => ((), (), (), (), ())
        :($(s :: Symbol)  => ())              => (s, (), (), s => :NONE, ())
        :(()              => $(sv :: Symbol)) => ((), sv, (), (), ())
        :($(s :: Symbol)  => $(sv :: Symbol)) => (s, sv, s => sv, s => :NONE, ())

        :(()              => () : $cu)              => ((), (), (), :NONE => cu, ())
        :($(s :: Symbol)  => () : $cu)              => (s, (), (), s => cu, ())
        :(()              => $(sv :: Symbol) : $cu) => ((), sv, (), :NONE => cu, ())
        :($(s :: Symbol)  => $(sv :: Symbol): $cu) => (s, sv, s => sv, s => cu, ())

        :($(s :: Symbol)  => sv)              => error("Non-symbolic second argument of foot: $sv")
        :($s              => $(sv :: Symbol)) => error("Non-symbolic first argument of foot: $s")
        :($s              => $sv)             => error("Foot definition requires symbolic names. Received: $s, $sv")
        Expr(:call, name, args...)            => error("Received: $name called with $args. Expected foot definition of form: A => B.")
        _                                     => error("Invalid foot definition.")
    end
end

function create_foot_U(block::Expr)
    @match block.head begin

        :tuple => begin
            if isempty(block.args) # case for create_foot(:())
                error("Cannot create foot with no arguments.")
            end
            foot_s = Vector{Symbol}()
            foot_sv = Vector{Symbol}()
            foot_ssv = Vector{Pair{Symbol, Symbol}}()
            foot_cu = Vector{Pair{Symbol, Union{Symbol, Expr}}}()
            foot_u = Vector{Symbol}()
            for (s, sv, ssv, cu, u) ∈ map(match_foot_format_U, block.args)
                if s != () push!(foot_s, s) end
                if sv != () push!(foot_sv, sv) end
                if ssv != () push!(foot_ssv, ssv) end
                if cu != () push!(foot_cu, cu) end
                if u != () push!(foot_u, u) end
            end
            return footU(unique(foot_s), unique(foot_sv), foot_ssv, foot_cu, unique(foot_u))
        end
        :call => footU(match_foot_format_U(block)...)
        _ => error("Invalid expression type $(block.head).  Expecting tuple or call.")
    end
end




macro foot_U(block::Expr)
    Base.remove_linenums!(block)
    return create_foot_U(block)
end



