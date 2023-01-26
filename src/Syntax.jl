module Syntax

using StockFlow

struct StockAndFlowSyntax
    stocks::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}}
    params::Array{Symbol}
    dyvars::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}}
    flows::Array{Pair{Symbol,Symbol}}
    sums::Array{Symbol}
end

function stocks_phase(stocks, stock)
    stock_name = stock.args[1]
    stock_syms = Tuple(Symbol(x) for x in stock.args[2].args)
    push!(stocks, (stock_name, stock_syms))
end

function params_phase(params, param)
    push!(params, param)
end

function dyvars_phase(dyvars, dyvar)
    dyvar_name = dyvar.args[1]
    dyvar_def = Tuple(d for d in dyvar.args[2].args)
    push!(dyvars, (dyvar_name, dyvar_def))
end

function flows_phase(flows, flow)
    push!(flows, (flow.args[2] => flow.args[3]))
end

function sums_phase(sums, sym)
    push!(sums, sym)
end

function collect_elements(statements)
    Base.remove_linenums!(statements)
    stocks::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}} = []
    params::Array{Symbol} = []
    dyvars::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}} = []
    flows::Array{Pair{Symbol,Symbol}} = []
    sums::Array{Symbol} = []
    current_phase = (_, _) -> ()
    for statement in statements
        if typeof(statement) <: QuoteNode
            kw = statement.value
            if kw == :stocks
                current_phase = s -> stocks_phase(stocks, s)
            elseif kw == :parameters
                current_phase = p -> params_phase(params, p)
            elseif kw == :dynamic_variables
                current_phase = d -> dyvars_phase(dyvars, d)
            elseif kw == :flows
                current_phase = f -> flows_phase(flows, f)
            elseif kw == :sums
                current_phase = s -> sums_phase(sums, s)
            else
                throw("Unknown keyword for Stock and Flow syntax: " * String(kw))
            end
        else
            current_phase(statement)
        end
    end

    return StockAndFlowSyntax(stocks, params, dyvars, flows, sums)
end

macro stock_and_flow(block)
    Base.remove_linenums!(block)
    syntax = collect_elements(block.args)
    stocks = Tuple(stock_name => outputs for (stock_name, outputs) in syntax.stocks)
    params = syntax.params
    dyvars = Tuple(dyvar_name => ((dyvar_param1, dyvar_param2) => dyvar_func) for (dyvar_name, (dyvar_func, dyvar_param1, dyvar_param2)) in syntax.dyvars)
    flows  = Tuple(f for f in syntax.flows)
    sums   = syntax.sums
    return StockAndFlowF(stocks, params, dyvars, flows, sums)
end

# New syntax
SIR = @stock_and_flow begin
    :stocks
    S = (F_NONE, inf, N)
    I = (inf, rec, N)
    R = (rec, F_NONE, N)

    :parameters
    c
    beta
    tRec

    :dynamic_variables
    v_prevalence = I / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I / tRec

    :flows
    inf => v_newInfections
    rec => v_newRecovery

    :sums
    N
end

# Current:
SIR_curr = StockAndFlowF((:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),# stocks
    (:c, :beta, :tRec),# parameters
    (:v_prevalence => ((:I, :N) => :/), :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*), :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
        :v_newInfections => ((:S, :v_perSIncidenceRate) => :*), :v_newRecovery => ((:I, :tRec) => :/)),# dynamical variables
    (:inf => :v_newInfections, :rec => :v_newRecovery),# flows
    (:N))# sum dynamical variables


end
