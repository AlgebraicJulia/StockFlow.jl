module Syntax

using StockFlow
using MLStyle

function replace_in_last(exprs, symtoreplace, targetsym)
    idx = lastindex(exprs)
    expr = last(exprs)
    try
        updatedexpr = @match expr begin
            Expr(:(=), sym, expr) => begin
                if sym == symtoreplace
                    :($targetsym = $expr)
                end
            end
        end
        exprs[idx] = updatedexpr
    catch _
    end
end

function infix_expression_to_binops(expression, finalsym=nothing)
    exprs = []
    function loop(e)
        @match e begin
            ::Symbol =>
                e
            Expr(:call, f, a) => begin
                asym = loop(a)
                varname = gensym("")
                push!(exprs, :($varname = $f($asym)))
                varname
            end
            Expr(:call, f, a, b) => begin
                asym = loop(a)
                bsym = loop(b)
                varname = gensym("")
                push!(exprs, :($varname = $f($asym, $bsym)))
                varname
            end
            Expr(:call, f, args...) => begin
                argsyms = map(loop, args)
                lastsym = gensym("")
                a = popfirst!(argsyms)
                b = popfirst!(argsyms)
                symexpr = :($lastsym = $f($a, $b))
                push!(exprs, symexpr)
                for argsym in argsyms
                    currsym = gensym("")
                    push!(exprs, :($currsym = $f($lastsym, $argsym)))
                    lastsym = currsym
                end
                lastsym
            end
            Expr(en, _, _, _) || Expr(en, _, _) => begin
                throw("Unhandled expression cannot be converted into form f(a, b) " * String(en))
            end
        end
    end
    lastsym = loop(expression)
    if finalsym !== nothing
        replace_in_last(exprs, lastsym, finalsym)
        lastsym = finalsym
    end
    exprs, lastsym
end

function extract_flow_name_and_equation(flow)
    @match flow begin
        :($flow_name($expr)) =>
            (flow_name, expr)
        Expr(en, _, _, _) || Expr(en, _, _) =>
            throw("Unhandled expression in flow definition " * String(en))
    end
end

function parse_flow(flow_definition)
    @match flow_definition begin
        :(TODO => $flow => $stock_out) || :(☁ => $flow => $stock_out) =>
            (:F_NONE, flow, stock_out)
        :($stock_in => $flow => TODO) || :($stock_in => $flow => ☁) =>
            (stock_in, flow, :F_NONE)
        :($stock_in => $flow => $stock_out) =>
            (stock_in, flow, stock_out)
        Expr(en, _, _, _) || Expr(en, _, _) =>
            throw("Unhandled expression in flow definition " * String(en))
    end
end
struct StockAndFlowArguments
    stocks::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}}
    params::Array{Symbol}
    dyvars::Array{Tuple{Symbol,Tuple{Symbol,Symbol,Symbol}}}
    flows::Array{Pair{Symbol,Symbol}}
    sums::Array{Symbol}
end

struct StockAndFlowSyntax
    stocks::Array{Symbol}
    params::Array{Symbol}
    dyvars::Array{Tuple{Symbol,Expr}}
    flows::Array{Tuple{Symbol,Expr,Symbol}}
    sums::Array{Tuple{Symbol,Expr}}
end

function stocks_phase(stocks, stock)
    push!(stocks, stock)
end

function params_phase(params, param)
    push!(params, param)
end

function dyvars_phase(dyvars, dyvar)
    @match dyvar begin
        :($dyvar_name = $dyvar_def) =>
            push!(dyvars, (dyvar_name, dyvar_def))
        Expr(c, _, _) || Expr(c, _, _, _) =>
            throw("Unhandled expression in dynamic variable definition " * String(c))
    end
end

function flows_phase(flows, flow)
    push!(flows, parse_flow(flow))
end

function sums_phase(sums, sum)
    @match sum begin
        :($sum_name = $equation) =>
            push!(sums, (sum_name, equation))
        Expr(c, _, _) || Expr(c, _, _, _) =>
            throw("Unhandled expression in sum defintion " * String(c))
    end
end

function collect_elements(statements)
    stocks::Array{Symbol} = []
    params::Array{Symbol} = []
    dyvars::Array{Tuple{Symbol,Expr}} = []
    flows::Array{Tuple{Symbol,Expr,Symbol}} = []
    sums::Array{Tuple{Symbol,Expr}} = []
    current_phase = (_, _) -> ()
    for statement in statements
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> stocks_phase(stocks, s)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> params_phase(params, p)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = d -> dyvars_phase(dyvars, d)
            end
            QuoteNode(:flows) => begin
                current_phase = f -> flows_phase(flows, f)
            end
            QuoteNode(:sums) => begin
                current_phase = s -> sums_phase(sums, s)
            end
            QuoteNode(kw) =>
                throw("Unknown block type for Stock and Flow syntax: " * String(kw))
            _ =>
                current_phase(statement)
        end
    end

    s = StockAndFlowSyntax(stocks, params, dyvars, flows, sums)
    println(s)
    return s
end

function assemble_stock_and_flow_f(syntax_elements)
end

macro stock_and_flow(block)
    Base.remove_linenums!(block)
    syntax = collect_elements(block.args)
    #    stocks = Tuple(stock_name => outputs for (stock_name, outputs) in syntax.stocks)
    #    params = syntax.params
    #    dyvars = Tuple(dyvar_name => ((dyvar_param1, dyvar_param2) => dyvar_func) for (dyvar_name, (dyvar_func, dyvar_param1, dyvar_param2)) in syntax.dyvars)
    #    flows = Tuple(f for f in syntax.flows)
    #    sums = syntax.sums
    #    return StockAndFlowF(stocks, params, dyvars, flows, sums)
end

test_infix_expr = :(a + b * c - d - e - f(g) + h(j, k) + l + m * n * o * p / q / r + s - t - u * v - w / x * y + z)
Base.remove_linenums!(test_infix_expr)

SIR = @stock_and_flow begin
    :stocks
    S
    I
    R

    :parameters
    c
    beta
    tRec
    omega
    alpha

    :dynamic_variables
    v_prevalence = I / N
    v_forceOfInfection = c * v_prevalence * beta

    :flows
    S => inf(S * v_forceOfInfection) => I
    ☁ => births(totalPopulation * alpha) => S
    S => deathsS(S * omega) => ☁
    I => rec(I / tRec) => R
    I => deathsI(I * omega) => ☁
    R => deathsR(R * omega) => ☁


    :sums
    totalPopulation = (S, I, R)
end

# New syntax
SIR_2 = @stock_and_flow begin
    :stocks
    S
    I
    R

    :parameters
    c
    beta
    tRec

    :dynamic_variables
    v_prevalence = I / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS

    :flows
    S => inf(S * v_perSIncidenceRate) => I
    I => rec(I / tRec) => R

    :sums
    N = [S, I, R]
end

# Current:
SIR_curr = StockAndFlowF((:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),# stocks
    (:c, :beta, :tRec),# parameters
    (:v_prevalence => ((:I, :N) => :/), :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*), :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
        :v_newInfections => ((:S, :v_perSIncidenceRate) => :*), :v_newRecovery => ((:I, :tRec) => :/)),# dynamical variables
    (:inf => :v_newInfections, :rec => :v_newRecovery),# flows
    (:N))# sum dynamical variables

end
