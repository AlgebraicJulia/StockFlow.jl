module Syntax

using StockFlow
using MLStyle


struct StockAndFlowArguments
    stocks::Vector{Pair{Symbol,Tuple{Union{Symbol,Vector{Symbol}},Union{Symbol,Vector{Symbol}},Union{Symbol,Vector{Symbol}}}}}
    params::Vector{Symbol}
    dyvars::Vector{Pair{Symbol,Pair{Tuple{Symbol,Symbol},Symbol}}}
    flows::Vector{Pair{Symbol,Symbol}}
    sums::Vector{Symbol}
end

struct StockAndFlowSyntax
    stocks::Vector{Symbol}
    params::Vector{Symbol}
    dyvars::Vector{Tuple{Symbol,Expr}}
    flows::Vector{Tuple{Symbol,Expr,Symbol}}
    sums::Vector{Tuple{Symbol,Vector{Symbol}}}
end

function set_final_binop_varname!(exprs::Vector{Tuple{Symbol,Expr}}, symtoreplace::Symbol, targetsym::Symbol)
    idx = lastindex(exprs)
    (s, expr) = last(exprs)
    if s == symtoreplace
      exprs[idx] = (targetsym, expr)
    end
end

function is_binop(e::Expr)
    @match e begin
        Expr(:call, f, a, b) => true
        _ => false
    end
end
function infix_expression_to_binops(expression::Expr, finalsym=nothing::Union{Nothing,Symbol})
    exprs::Vector{Tuple{Symbol,Expr}} = []
    function loop(e)
        @match e begin
            ::Symbol =>
                e
            Expr(:call, f, a, b) => begin
                asym = loop(a)
                bsym = loop(b)
                varname = gensym("")
                push!(exprs, (varname, :($f($asym, $bsym))))
                varname
            end
            Expr(:call, f, args...) => begin
                argsyms = map(loop, args)
                lastsym = gensym("")
                a = popfirst!(argsyms)
                b = popfirst!(argsyms)
                symexpr = :($f($a, $b))
                push!(exprs, (lastsym, symexpr))
                for argsym in argsyms
                    currsym = gensym("")
                    push!(exprs, (currsym, :($f($lastsym, $argsym))))
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
        set_final_binop_varname!(exprs, lastsym, finalsym)
        lastsym = finalsym
    end
    exprs, lastsym
end

function extract_flow_name_and_equation(flow::Expr)
    @match flow begin
        :($flow_name($expr)) =>
            (flow_name, expr)
        :($flow_name($expr, name=$_)) =>
            (flow_name, expr)
        Expr(en, _, _, _) || Expr(en, _, _) =>
            throw("Unhandled expression in flow name definition " * String(en))
    end
end

function parse_flow_io(flow_definition::Expr)
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

function parse_stock!(stocks::Vector{Symbol}, stock::Symbol)
    push!(stocks, stock)
end

function parse_param!(params::Vector{Symbol}, param::Symbol)
    push!(params, param)
end

function parse_dyvar!(dyvars::Vector{Tuple{Symbol,Expr}}, dyvar::Expr)
    @match dyvar begin
        :($dyvar_name = $dyvar_def) =>
            push!(dyvars, (dyvar_name, dyvar_def))
        Expr(c, _, _) || Expr(c, _, _, _) =>
            throw("Unhandled expression in dynamic variable definition " * String(c))
    end
end

function parse_flow!(flows::Vector{Tuple{Symbol,Expr,Symbol}}, flow::Expr)
    push!(flows, parse_flow_io(flow))
end

function parse_sum!(sums::Vector{Tuple{Symbol,Vector{Symbol}}}, sum::Expr)
    @match sum begin
        :($sum_name = $equation) =>
            push!(sums, (sum_name, equation.args))
        Expr(c, _, _) || Expr(c, _, _, _) =>
            throw("Unhandled expression in sum defintion " * String(c))
    end
end

function parse_stock_and_flow_syntax(statements::Vector{Any})
    stocks::Vector{Symbol} = []
    params::Vector{Symbol} = []
    dyvars::Vector{Tuple{Symbol,Expr}} = []
    flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
    sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []
    current_phase = (_, _) -> ()
    for statement in statements
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> parse_stock!(stocks, s)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> parse_param!(params, p)
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
                throw("Unknown block type for Stock and Flow syntax: " * String(kw))
            _ =>
                current_phase(statement)
        end
    end

    s = StockAndFlowSyntax(stocks, params, dyvars, flows, sums)
    return s
end


function fnone_or_tuple(arrows::Vector{Symbol})
    if isempty(arrows)
        :F_NONE
    elseif length(arrows) == 1
        arrows[1]
    else
        arrows
    end
end

function assemble_stock_defintions(stocks::Vector{Symbol}, flows::Vector{Tuple{Symbol,Expr,Symbol}}, sum_variables::Vector{Tuple{Symbol,Vector{Symbol}}})
    formatted_stocks = []
    for stock in stocks
        input_arrows::Vector{Symbol} = []
        output_arrows::Vector{Symbol} = []
        sum_arrows::Vector{Symbol} = []
        for (start_object, flow, end_object) in flows
            (flow_name, _) = extract_flow_name_and_equation(flow)
            if start_object == stock
                push!(output_arrows, flow_name)
            end
            if end_object == stock
                push!(input_arrows, flow_name)
            end
        end
        for (sum_variable_name, inputs) in sum_variables
            if stock in inputs
                push!(sum_arrows, sum_variable_name)
            end
        end
        push!(formatted_stocks, (stock => (fnone_or_tuple(input_arrows), fnone_or_tuple(output_arrows), fnone_or_tuple(sum_arrows))))
    end
    return formatted_stocks
end

function unfold_dyvars_to_binops(dyvars::Vector{Tuple{Symbol,Expr}})
    syms::Vector{Pair{Symbol,Pair{Tuple{Symbol,Symbol},Symbol}}} = []
    for (dyvar_name, dyvar_definition) in dyvars
        if is_binop(dyvar_definition)
            @match dyvar_definition begin
                Expr(:call, op, a, b) => begin
                    push!(syms, (dyvar_name => ((a, b) => op)))
                end
                Expr(c, _, _) || Expr(c, _, _, _) =>
                    throw("Unhandled expression in dynamic variable definition " * String(c))
            end
        else
            (binops, _) = infix_expression_to_binops(dyvar_definition, dyvar_name)
            binops_syms = unfold_dyvars_to_binops(binops)
            syms = vcat(syms, binops_syms)
        end
    end
    return syms
end

function disassemble_flow_equation(flow_expression, dyvars)
    @match flow_expression begin
        :($flow_name($expr)) => begin
          if typeof(expr) <: Symbol
              # TODO: search for the 'op' in dyvars to see if it's defined as a dynamic variable.
              #       if so, we can assemble it into 'flow_name => dyvar_name' in stock and flow
              #       Elsewise, err.
              return ([], flow_name => expr)
          else
              (additional_dyvars, var_name) = infix_expression_to_binops(expr)
              return (additional_dyvars, flow_name => var_name)
          end
        end
        :($flow_name($expr, name=$sym)) => begin
            (additional_dyvars, var_name) = infix_expression_to_binops(expr, sym)
            return (additional_dyvars, flow_name => sym)
        end
        Expr(c, _, _) || Expr(c, _, _, _) => begin
            throw("Unhandled expression in flow equation definition " * String(c))
        end
    end
end

function assemble_flows(flows::Vector{Tuple{Symbol,Expr,Symbol}}, dyvars)
    flow_definitions = []
    updated_dyvars = []
    for (start_object, flow, end_object) in flows
        (flow_name, equation) = extract_flow_name_and_equation(flow)
        (additional_dyvars, var_name) = infix_expression_to_binops(equation)
        push!(flow_definitions, (flow_name => var_name))
        updated_dyvars = vcat(updated_dyvars, unfold_dyvars_to_binops(additional_dyvars))
    end
    return (flow_definitions, updated_dyvars)
end

function assemble_stock_and_flow_f(syntax_elements::StockAndFlowSyntax)
    stocks = assemble_stock_defintions(syntax_elements.stocks, syntax_elements.flows, syntax_elements.sums)
    params = syntax_elements.params
    binop_dyvars = unfold_dyvars_to_binops(syntax_elements.dyvars)
    (flows, flow_dyvars) = assemble_flows(syntax_elements.flows, binop_dyvars)
    all_dyvars = vcat(binop_dyvars, flow_dyvars)
    sums = [s for (s, _) in syntax_elements.sums]
    return StockAndFlowArguments(stocks, params, all_dyvars, flows, sums)
end

macro stock_and_flow(block)
    Base.remove_linenums!(block)
    syntax = parse_stock_and_flow_syntax(block.args)
    args = assemble_stock_and_flow_f(syntax)
    return StockAndFlowF(args.stocks, args.params, args.dyvars, args.flows, args.sums)
end

test_infix_expr = :(a + b * c - d - e - f(g) + h(j, k) + l + m * n * o * p / q / r + s - t - u * v - w / x * y + z)
Base.remove_linenums!(test_infix_expr)
# New syntax
SIR = @stock_and_flow begin
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
    S => inf(S * v_perSIncidenceRate, name=v_newInfections) => I
    I => rec(I / tRec, name=v_newRecovery) => R

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


# TODO
# SIR = @stock_and_flow begin
#     :stocks
#     S
#     I
#     R
#
#     :parameters
#     c
#     beta
#     tRec
#     omega
#     alpha
#
#     :dynamic_variables
#     v_prevalence = I / N
#     v_forceOfInfection = c * v_prevalence * beta
#
#     :flows
#     S => inf(S * v_forceOfInfection) => I
#     ☁ => births(totalPopulation * alpha) => S
#     S => deathsS(S * omega) => ☁
#     I => rec(I / tRec) => R
#     I => deathsI(I * omega) => ☁
#     R => deathsR(R * omega) => ☁
#
#
#     :sums
#     totalPopulation = (S, I, R)
# end

end
