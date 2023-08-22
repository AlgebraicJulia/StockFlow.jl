""" Alternative syntax for use in the definition of stock and flow models.

### Examples
```julia
# An S-I-R model of infection
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
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I / tRec

    :flows
    S => inf(v_newInfections) => I
    I => rec(v_newRecovery) => R

    :sums
    N = [S, I, R]
end

# Generates:
# SIR = StockAndFlowF(
#     # stocks
#     (:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),
#     # parameters
#     (:c, :beta, :tRec),
#     # dynamical variables
#     (   :v_prevalence => ((:I, :N) => :/),
#         :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*),
#         :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
#         :v_newInfections => ((:S, :v_perSIncidenceRate) => :*),
#         :v_newRecovery => ((:I, :tRec) => :/),
#     ),
#     # flows
#     (:inf => :v_newInfections, :rec => :v_newRecovery),
#     # sum dynamical variables
#     (:N),
# )

# The same model as before, but with the dynamic variables inferred
SIR_2 = @stock_and_flow begin
    :stocks
    S
    I
    R

    :parameters
    c
    beta
    tRec

    # We can leave out dynamic variables and let them be inferred from flows entirely!

    :flows
    S => inf(S * beta * (c * (I / N))) => I
    I => rec(I / tRec) => R

    :sums
    N = [S, I, R]
end

# Another possible S-I-R model definition
SIR_3 = @stock_and_flow begin
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
    v_prevalence = I / totalPopulation
    v_forceOfInfection = c * v_prevalence * beta

    :flows
    S => inf(S * v_forceOfInfection) => I
    ☁ => births(totalPopulation * alpha) => S
    S => deathsS(S * omega) => ☁
    I => rec(I / tRec) => R
    I => deathsI(I * omega) => ☁
    R => deathsR(R * omega) => ☁


    :sums
    totalPopulation = [S, I, R]
end
```
"""
module Syntax
export @stock_and_flow, @foot, @feet, @stratify

using ..StockFlow
using MLStyle
import Base.get

"""
    stock_and_flow(block :: Expr)

Compiles stock and flow syntax of the line-based block form
```julia
  :stocks
    symbol_1
    symbol_2
    ...
    symbol_n

  :parameters
    param_1
    param_2
    ...
    param_n

  :dynamic_variables
    dyvar_1 = symbol_h * param_g ... - symbol_x / param_y
    ...
    dyvar_n = symbol_k * param_j - dyvar_a ... - symbol_p / param_q

  :flows
    symbol_r => flow_name_1(dyvar_k) => symbol_q
    symbol_z => flow_name_2(dyvar_g * param_v) => symbol_p
    ☁       => flow_name_3(symbol_c + dyvar_b) => symbol_r
    symbol_j => flow_name_4(param_l + symbol_m) => TODO
    ...
    symbol_y => flow_name_n(dyvar_f) => ☁
```
into a StockAndFlowF data type for use with the StockFlow.jl modelling system.
"""
macro stock_and_flow(block)
    Base.remove_linenums!(block)
    syntax_lines = parse_stock_and_flow_syntax(block.args)
    saff_args = stock_and_flow_syntax_to_arguments(syntax_lines)
    return StockAndFlowF(
        saff_args.stocks,
        saff_args.params,
        map(kv -> kv.first => get(kv.second), saff_args.dyvars),
        saff_args.flows,
        saff_args.sums,
    )
end

"""
    StockAndFlowBlock

Contains the elements that make up the Stock and Flow block syntax.

### Fields
- `stocks` -- Each stock is defined by a single valid Julia variable name on a line.
- `params` -- Each parameter is defined by a single valid Julia variable name on a line.
- `dyvars` -- Each dynamic variable is defined by a valid Julia assignment statement of the
              form `dyvar = expr`
- `flows`  -- Each flow is defined by a valid Julia pair expression of the form
              `variable_name => flow_name(flow_expression) => variable_name`
- `sums`   -- Each sum is defined by a valid Julia assignment statement of the form
              `sum_name = [a, b, c, ...]`
"""
struct StockAndFlowBlock
    stocks::Vector{Symbol}
    params::Vector{Symbol}
    dyvars::Vector{Tuple{Symbol,Expr}}
    flows::Vector{Tuple{Symbol,Expr,Symbol}}
    sums::Vector{Tuple{Symbol,Vector{Symbol}}}
end

"""
Contains the five parameters required to instantiate a StockAndFlowF data type,
representing a Stock and Flow model.

This can be used to directly instantiate StockAndFlowF.
"""
abstract type DyvarExprT end
struct Binop{P1,P2} <: DyvarExprT
    binop::Pair{Tuple{P1,P2},Symbol}
end
struct Ref <: DyvarExprT
    ref::Pair{Symbol,Symbol}
end
get(r::Ref) = r.ref
get(b::Binop) = b.binop
struct StockAndFlowArguments
    stocks::Vector{
        Pair{
            Symbol,
            Tuple{
                Union{Symbol,Vector{Symbol}},
                Union{Symbol,Vector{Symbol}},
                Union{Symbol,Vector{Symbol}},
            },
        },
    }
    params::Vector{Symbol}
    dyvars::Vector{Pair{Symbol, DyvarExprT}}
    flows::Vector{Pair{Symbol,Symbol}}
    sums::Vector{Symbol}
end

"""
    parse_stock_and_flow_syntax(statements :: Vector{Any})

Given a vector of Julia expressions, attempt to interpret them using the block syntax
defined for Stock and Flow diagrams.

### Input
- `statements` -- A series of Julia expressions, each a line of code in
                  a block of statements.

### Output
A StockAndFlowSyntax data type which contains the syntax pieces required to define
a Stock and Flow model: stocks, parameters, dynamic variables, flows, and sums.
"""
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
                error("Unknown block type for Stock and Flow syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    s = StockAndFlowBlock(stocks, params, dyvars, flows, sums)
    return s
end

"""
    stock_and_flow_syntax_to_arguments(syntax_elements::StockAndFlowSyntax)

Convert the Stock and Flow Syntax elements to parameters suitable for StockAndFlowF
data type instantiation

### Input
- `syntax_elements` -- The output from `parse_stock_and_flow_syntax`,
                       containing the definition of a stock and flow model

### Output
Parameters for instantiation of a StockAndFlowF data type.
"""
function stock_and_flow_syntax_to_arguments(syntax_elements::StockAndFlowBlock)
    stocks = assemble_stock_definitions(
        syntax_elements.stocks,
        syntax_elements.flows,
        syntax_elements.sums,
    )
    params = syntax_elements.params
    dyvars = dyvar_exprs_to_symbolic_repr(syntax_elements.dyvars)
    dyvar_names = [dyvar_name for (dyvar_name, _dyvar_def) in dyvars]
    (flows, flow_dyvars) = create_flow_definitions(syntax_elements.flows, dyvar_names)
    sums = sum_variables(syntax_elements.sums)
    return StockAndFlowArguments(stocks, params, vcat(dyvars, flow_dyvars), flows, sums)
end


"""
    parse_stock!(stocks :: Vector{Symbol}, stock :: Symbol)

Add a stock to the list of stocks.

Stocks in the Stock and Flow syntax are just a single symbol on a line, so no work
is currently done by this function.

### Input
- `stocks` -- A list of stocks already parsed from the block
- `stock` -- A stock symbol

### Output
None. This mutates the given stocks vector.
"""
function parse_stock!(stocks::Vector{Symbol}, stock::Symbol)
    push!(stocks, stock)
end

"""
    parse_param!(params :: Vector{Symbol}, param :: Symbol)

Add a param to the list of params.

Params in the Stock and Flow syntax are just a single symbol on a line,
so no work is currently done by this function.

### Input
- `params` -- A list of params already parsed from the block
- `param` -- A param symbol

### Output
None. This mutates the given params vector.
"""
function parse_param!(params::Vector{Symbol}, param::Symbol)
    push!(params, param)
end

"""
   is_recursive_dyvar(dyvar_name :: Symbol, dyvar_def :: Expr)

Check that the dyvar_name is not used somewhere in the dyvar_def

### Input
- `dyvar_name` -- A dyvar name as a Julia Symbol
- `dyvar_expr` -- A Julia expression

### Output
True if the dyvar name is used in the expression; false elsewise.
"""
function is_recursive_dyvar(dyvar_name, dyvar_def)
    @match dyvar_def begin
        ::Symbol => dyvar_def == dyvar_name
        :($f()) => f == dyvar_name
        Expr(:call, args...) => true in map(arg -> is_recursive_dyvar(dyvar_name, arg), args)
    end
end
"""
    parse_dyvar!(dyvars :: Vector{Tuple{Symbol, Expr}}, dyvar :: Expr)

Extract the dynamic variable name and defining expression from a Julia expression of form
`dyvar = a + b * c ...`, and add it to the vector of already parsed dynamic variables.

### Input
- `dyvars` -- A list of dynamic variables and their defining Julia expressions
              already parsed from the block
- `dyvar` -- A dynamic variable definition of the form `dyvar = defining_expression`

### Output
None. This mutates the given dyvars vector.
"""
function parse_dyvar!(dyvars::Vector{Tuple{Symbol,Expr}}, dyvar::Expr)
    @match dyvar begin
        :($dyvar_name = $dyvar_def) =>
            if !is_recursive_dyvar(dyvar_name, dyvar_def)
                push!(dyvars, (dyvar_name, dyvar_def))
            else
                error("Recursive dyvar detected in Symbol: " * String(dyvar_name))
            end
        Expr(c, _, _) || Expr(c, _, _, _) =>
            error("Unhandled expression in dynamic variable definition " * String(c))
    end
end

"""
    parse_flow_io(flow_definition :: Expr)

Given a flow definition of the form `SYMBOL => flow_name(flow_equation) => SYMBOL`,
return a 3-tuple of its constituent parts: the start symbol, the end symbol,
and the flow equations's definition as an expression.

### Input
- `flow_definition` -- A flow definition of the form
                       `SYMBOL => flow_name(flow_equation) => SYMBOL`,
                       where SYMBOL can be an arbitrary name or special cases of ☁ or TODO,
                       which corresponds to a flow from nowhere.

### Output
A 3-tuple (input, expression, output) of a an input and output symbol (either of which
may be :F_NONE for a flow from nowhere) and the flow equation as a julia expression.

### Examples
```julia-repl
julia> Syntax.parse_flow_io(:(TODO => birthRate(a * b * c) => S))
(:F_NONE, :(birthRate(a * b * c)), :S)
```
"""
function parse_flow(flow_definition::Expr)
    @match flow_definition begin
        :(TODO => $flow => $stock_out) || :(☁ => $flow => $stock_out) =>
            (:F_NONE, flow, stock_out)
        :($stock_in => $flow => TODO) || :($stock_in => $flow => ☁) =>
            (stock_in, flow, :F_NONE)
        :($stock_in => $flow => $stock_out) => (stock_in, flow, stock_out)
        Expr(en, _, _, _) || Expr(en, _, _) =>
            error("Unhandled expression in flow definition " * String(en))
    end
end

"""
    parse_flow!(flows :: Vector{Tuple{Symbol, Expr, Symbol}}, flow :: Expr)

Extract the flow input, output, and defining expression from a Julia expression of
the form `IN_SYMBOL => FLOW_NAME(FLOW_EQUATION) => OUT_SYMBOL`,
and add it to the vector of already parsed flows.

### Input
- `flows` -- A list of parsed flows
- `flow` -- A flow definition of the form `in => flow_name(flow_equation) => out`

### Output
None. This mutates the given flows vector.
"""
function parse_flow!(flows::Vector{Tuple{Symbol,Expr,Symbol}}, flow::Expr)
    parsed_flow = parse_flow(flow)
    push!(flows, parsed_flow)
end

"""
    parse_sum!(sums :: Vector{Tuple{Symbol, Vector{Symbol}}}, sum :: Expr)

Extract the sum name and the stocks that flow into it from a stock expression of the form
`SUM_NAME = [STOCK_1, STOCK_2, STOCK_3, ...]`

### Input
- `sums` -- A list of parsed sum names and their incoming stocks
- `sum`  -- A sum definition as a Julia expression of the form
            `sum_name = [stock_1, stock_2, stock_3, ...]`

### Output
None. This mutates the given sums vector.
"""
function parse_sum!(sums::Vector{Tuple{Symbol,Vector{Symbol}}}, sum::Expr)
    @match sum begin
        :($sum_name = $equation) => push!(sums, (sum_name, equation.args))
        Expr(c, _, _) || Expr(c, _, _, _) =>
            error("Unhandled expression in sum defintion " * String(c))
    end
end

"""
    assemble_stock_definitions( stocks::Vector{Symbol}
                              , flows::Vector{Tuple{Symbol,Expr,Symbol}}
                              , sum_variables::Vector{Tuple{Symbol,Vector{Symbol}}}
                              )

Convert the raw syntax of a Stock and Flow block definition into a series of stock
definitions suitable for input into the StockAndFlowF data type, which is
(:STOCK_NAME => ((INPUT_ARROWS,...), (OUTPUT_ARROWS,...), (SUM_ARROWS,...))).
The input, output, and sum arrows are calculated from the sum and flow definitions.

### Input
- `stocks` -- A list of stock names as symbols
- `flows` -- A list of flow definitions, which may use any of the stock names as
             inputs or outputs.
- `sum_variables` -- A list of sum definitions.

### Output
The raw syntax definitions of a Stock and Flow block rearranged into a
StockAndFlowF stock definition vector.
"""
function assemble_stock_definitions(
    stocks::Vector{Symbol},
    flows::Vector{Tuple{Symbol,Expr,Symbol}},
    sum_variables::Vector{Tuple{Symbol,Vector{Symbol}}},
)
    formatted_stocks = []
    for stock in stocks
        input_arrows::Vector{Symbol} = []
        output_arrows::Vector{Symbol} = []
        sum_arrows::Vector{Symbol} = []
        for (start_object, flow, end_object) in flows
            (flow_name, _) = extract_function_name_and_args_expr(flow)
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
        push!(
            formatted_stocks,
            (
                stock => (
                    fnone_value_or_vector(input_arrows),
                    fnone_value_or_vector(output_arrows),
                    default_or_value_or_vector(sum_arrows; default=:SV_NONE),
                )
            ),
        )
    end
    return formatted_stocks
end

"""
    generate_dyvar_name(name :: Symbol)

Generate a dynamic variable name using the given symbol as a prefix. All dynamic variables are prefixed with 'v_'

### Input
- `name` -- A symbol base for the dyvar

### Output
`name` as a string prefixed with `v_`, for input to `gensym`

### Example
```julia-repl
julia> sym = generate_dyvar_name(:my_var)
"v_my_var"
julia> gensym(sym)
Symbol("##v_my_var#718")
"""
function generate_dyvar_name(name::Symbol)
    name_str = String(name)
    if cmp("v_", name_str) == -1
        name_str
    else
        "v_" * name_str
    end
end

"""
    dyvar_exprs_to_symbolic_repr(dyvars::Vector{Tuple{Symbol,Expr}})

Converts a series of dynamic variable definitions of the form `dyvar = dyvar_expression`
into a form suitable for input into the StockAndFlowF data type:
(:dyvar => ((:arg1, :arg2) => :function_name))

### Input
- `dyvars` -- A vector of pairs of dynamic variable names (as symbols)
              and Julia expressions defining them.

### Output
A vector of dynamic variable definitions suitable for input to StockAndFlowF.
"""
function dyvar_exprs_to_symbolic_repr(dyvars::Vector{Tuple{Symbol,Expr}})
    syms::Vector{Pair{Symbol,DyvarExprT}} = []
    for (dyvar_name, dyvar_definition) in dyvars
        if is_binop_or_unary(dyvar_definition)
            @match dyvar_definition begin
                Expr(:call, op, a) => push!(syms, (dyvar_name => Ref(a => op)))
                Expr(:call, op, a, b) => begin
                    push!(syms, (dyvar_name => Binop((a, b) => op)))
                end
                Expr(c, _, _) || Expr(c, _, _, _) => error(
                    "Unhandled expression in dynamic variable definition " * String(c),
                )
            end
        else
            (binops, _) =
                infix_expression_to_binops(dyvar_definition, finalsym=dyvar_name, gensymbase=generate_dyvar_name(dyvar_name))
            binops_syms = dyvar_exprs_to_symbolic_repr(binops)
            syms = vcat(syms, binops_syms)
        end
    end
    return syms
end

"""
    create_flow_definitions( flows::Vector{Tuple{Symbol,Expr,Symbol}}
                           , dyvar_names :: Vector{Symbol}
                           )

Assemble the flow parameters for a StockAndFlowF data type from definitions of the form
`input_stock => flow_name(flow_equation) => output_stock`.

### Input
- `flows` -- A vector of flow definitions in the form of 3-tuple
             `(input_stock, flow_expr, output_stock)
- `dyvar_names` -- A vector of known dynamic variable names

### Output
A 2-tuple containing in the first cell the flow definitions of the form
`flow_name => dynamic_variable_name`, and in the second cell any dynamic variables
generated in generating the flow definitions.
"""
function create_flow_definitions(flows::Vector{Tuple{Symbol,Expr,Symbol}}, dyvar_names)
    flow_definitions = []
    updated_dyvars = []
    # Start and end objects here can be ignored:
    # In StockAndFlowF, they are factored into the stocks parameter, not here.
    for (_start_object, flow, _end_object) in flows
        (additional_dyvars, flow_definition) = flow_expr_to_symbolic_repr(flow, dyvar_names)
        push!(flow_definitions, flow_definition)
        updated_dyvars = vcat(updated_dyvars, additional_dyvars)
    end
    return (flow_definitions, updated_dyvars)
end

"""
    flow_expr_to_symbolic_repr(flow_expression :: Expr, dyvar_names :: Vector{Symbol})

Converts a flow definition from one of three forms into a format suitable
for input to StockAndFlowF data types:

- `flow_name(name)` -- uses the dynamic variable `name` as the flow equation
- `flow_name(expr)` -- hoists the expr out as a dynamic variable and replaces
                       it with the generated name
- `flow_name(expr, name=name)` -- hoists the expr out as a dynamic variable and replaces
                                  it with the given name

### Input
- `flow_expression` -- a flow definition equation from Stock and Flow syntax
- `dyvar_names` -- a vector of symbols that the flow_expression may be referring to

### Output
The flow_expression's representation StockAndFlowF parameter format:
`flow_name => dynamic_variable`
"""
function flow_expr_to_symbolic_repr(flow_expression, dyvar_names)
    @match flow_expression begin
        :($flow_name($expr)) => begin
            if typeof(expr) <: Symbol
                if expr in dyvar_names
                    return ([], flow_name => expr)
                else
                    error("Unknown dynamic variable referenced " * String(expr))
                end
            else
                (additional_dyvars, var_name) = infix_expression_to_binops(expr, gensymbase=generate_dyvar_name(flow_name))
                dyvs = dyvar_exprs_to_symbolic_repr(additional_dyvars)
                return (dyvs, flow_name => var_name)
            end
        end
        :($flow_name($expr, name=$sym)) => begin
            (additional_dyvars, var_name) = infix_expression_to_binops(expr, finalsym=sym)
            dyvs = dyvar_exprs_to_symbolic_repr(additional_dyvars)
            return (dyvs, flow_name => sym)
        end
        Expr(c, _, _) || Expr(c, _, _, _) => begin
            error("Unhandled expression in flow equation definition " * String(c))
        end
    end
end


"""
  extract_function_name_and_args_expr(flow::Expr)

Given a Julia expression of the form f(a), return the symbol :f and the expression, a,
the function is being called with.

### Input
- `flow` -- a julia expression of the form f(a)

### Output
The name of the function as a symbol, and the expression being passed to it as an argument.

### Examples
```julia-repl
julia> Syntax.extract_flow_name_and_equation(:(infectionRate(a + b + c)))
(:infectionRate, :(a + b + c))
```
"""
function extract_function_name_and_args_expr(flow_equation::Expr)
    @match flow_equation begin
        :($flow_name($expr)) => (flow_name, expr)
        :($flow_name($expr, name=$_)) => (flow_name, expr)
        Expr(en, _, _, _) || Expr(en, _, _) =>
            error("Unhandled expression in flow name definition " * String(en))
    end
end
"""
    default_or_value_or_tuple(arrows :: Vector{Symbol})

Given a vector of arrow names, modify it into suitable input for a StockAndFlowF data type:
default for an empty vector, the value alone for a singleton vector,
and the vector itself if there are multiple arrows given.

### Input
- `arrows` -- A vector of symbols which represents some of the arrows for a stock.
- `default` -- A default symbol if an empty vector is passed in

### Output
given default, a single symbol, or a list of symbols.

### Example
```julia-repl
`
"""
function default_or_value_or_vector(arrows::Vector{Symbol}; default=:F_NONE)
    if isempty(arrows)
        default
    elseif length(arrows) == 1
        arrows[1]
    else
        arrows
    end
end


"""
    fnone_or_tuple(arrows :: Vector{Symbol})

Given a vector of arrow names, modify it into suitable input for a StockAndFlowF data type:
:F_NONE for an empty vector, the value alone for a singleton vector,
and the vector itself if there are multiple arrows given.

### Input
- `arrows` -- A vector of symbols which represents some of the arrows for a stock.

### Output
:F_NONE, a single symbol, or a list of symbols.

### Example
```julia-repl
`
"""
function fnone_value_or_vector(arrows::Vector{Symbol})
    default_or_value_or_vector(arrows)
end

"""
   infix_expression_to_binops( expression :: Expr; gensymbase :: String = ""
                             , finalsym :: Union{Nothing, Symbol} = nothing)

Convert a nested expression of the form (a * b + c - d / e ...) into a series of binary
operations with named values (sym1 = a * b; sym2 = sym1 + c; ...; symn = symn-1 + somevar).
If finalsym is given, symn is replaced with the one given.

### Input
- `expression` -- a single expression containing nested function calls.
- `gensymbase` -- the base name of the generated interim symbols for each binop
- `sym`   -- (optional, default: `nothing`) a name for the final result to replace
             one of the generated symbols. If not given, it will be of the form
             `Symbol("###<NUMBER>")` e.g. `Symbol("###418")`.

### Output
A vector of tuples of variable definitions as symbols and the corresponding
Julia expression that calculates the variable.

### Examples
```julia-repl
julia> infix_expression_to_binops(:(a + b + c))
julia> Syntax.infix_expression_to_binops(:(a + b + c))
( Tuple{Symbol, Expr}[ (Symbol("###495"), :(a + b))
                     , (Symbol("###496"), :(var"###495" + c))
                     ]
, Symbol("###496")
)
```
That is, this converts the expression `a + b + c` into two expressions:
```julia
###495 = a + b
###496 = ###495 + c
```

```
julia> Syntax.infix_expression_to_binops( :(a + b + c)
                                        , gensymbase="generated_variable"
                                        , lastsym=:infectionRate)
( Tuple{Symbol, Expr}[ (Symbol("##generated_variable#1045"), :(a + b))
                     , (:infectionRate, :(var"##generated_variable#1045" + c))
                     ]
, :infectionRate)
```
That is, this converts the expression `a + b + c` into the expressions:
```julia
generated_variable#1045 = a + b
infectionRate = generatedVariable#1045 + c
```
This would be used in the case the original expression was `infectionRate = a + b + c`.
"""
function infix_expression_to_binops(
    expression::Expr;
    gensymbase::String="",
    finalsym::Union{Nothing,Symbol}=nothing
)
    exprs::Vector{Tuple{Symbol,Expr}} = []
    function loop(e)
        @match e begin
            ::Symbol || ::Float32 || ::Float64 || ::Int || ::String => e
            Expr(:call, f, a) => begin
                asym = loop(a)
                varname = gensym(gensymbase)
                push!(exprs, (varname, :($f($asym))))
                varname
            end
            Expr(:call, f, a, b) => begin
                asym = loop(a)
                bsym = loop(b)
                varname = gensym(gensymbase)
                push!(exprs, (varname, :($f($asym, $bsym))))
                varname
            end
            Expr(:call, f, args...) => begin
                if (isempty(args))
                    error("Expression f() cannot be converted into form f(a, b)")
                end
                argsyms = map(loop, args)
                lastsym = gensym(gensymbase)
                a = popfirst!(argsyms)
                b = popfirst!(argsyms)
                symexpr = :($f($a, $b))
                push!(exprs, (lastsym, symexpr))
                for argsym in argsyms
                    currsym = gensym(gensymbase)
                    push!(exprs, (currsym, :($f($lastsym, $argsym))))
                    lastsym = currsym
                end
                lastsym
            end
            Expr(en, _, _, _) || Expr(en, _, _) => begin
                error(
                    "Unhandled expression type " * String(en) * " cannot be converted into form f(a, b)",
                )
            end
        end
    end
    last_generated_sym = loop(expression)
    if finalsym !== nothing
        set_final_binop_varname!(exprs, finalsym)
        return (exprs, finalsym)
    else
        return (exprs, last_generated_sym)
    end
end

"""
    sum_variables(sum_syntax_elements :: Vector{Tuple{Symbol, Expr}})

Extract the names of the sum variables from its syntax elements.

### Input
`sum_syntax_elements` -- A vector of 2-tuples of the sum variable's name and
                         its constituent stocks.

### Output
A vector of sum variable names.
"""
sum_variables(sum_syntax_elements) =
    [sum_name for (sum_name, _sum_definition) in sum_syntax_elements]


"""
    is_binop_or_unary(e :: Expr)

Check if a Julia expression is a call of the form `op(a, b)` or `a op b`

### Input
- `e` -- a Julia expression

### Output
A boolean indicating if the given julia expression is a function call of two parameters.

### Examples
```julia-repl
julia> is_binop_or_unary(:(f()))
false
julia> is_binop_or_unary(:(f(a)))
true
julia> is_binop_or_unary(:(f(a, b)))
true
julia> is_binop_or_unary(:(a * b))
true
julia> is_binop_or_unary(:(f(a, b, c)))
false
```
"""
function is_binop_or_unary(e::Expr)
    @match e begin
        Expr(:call, f::Symbol, a) => true
        Expr(:call, f::Symbol, a, b) => true
        _ => false
    end
end

"""
    set_finaL_binop_varname!(exprs::Vector{Tuple{Symbol, Expr}}, targetsym::Symbol)

Given an expression that is in the form of a series of binary operations
(e.g. from infix_expression_to_binops), replace the final generated symbol with a given one.

### Input
- `exprs` -- A vector of tuples of variable name symbols and their corresponding expression
             definitions.
- `varname` -- The final variable name to set.

### Output
The original collection, with the last element updated to have a new variable name
in its tuple.

### Examples
```julia-repl
julia> (binops, _final_sym_name) = Syntax.infix_expression_to_binops(:(a + b + c))
( Tuple{Symbol, Expr}[ (Symbol("###1415"), :(a + b))
                     , (Symbol("###1416"), :(var"###1415" + c))
                     ]
, Symbol("###1416")
)

julia> binops
2-element Vector{Tuple{Symbol, Expr}}:
 (Symbol("###1415"), :(a + b))
 (Symbol("###1416"), :(var"###1415" + c))

julia> Syntax.set_final_binop_varname!(binops, :infectionRate)

julia> binops
2-element Vector{Tuple{Symbol, Expr}}:
 (Symbol("###1415"), :(a + b))
 (:infectionRate, :(var"###1415" + c))
```
"""
function set_final_binop_varname!(exprs::Vector{Tuple{Symbol,Expr}}, varname::Symbol)
    idx = lastindex(exprs)
    (_oldvarname, expr) = last(exprs)
    exprs[idx] = (varname, expr)
end




"""
    foot(block :: Expr)

Create a foot with S => N syntax, where S is stock, N is sum variable.
```julia
@foot P => Q
@foot S1 => ()
@foot () => N
@foot () => ()
@foot A => N, () => NI
```
"""
macro foot(block::Expr)
    Base.remove_linenums!(block)
    return create_foot(block)
end


"""
    feet(block :: Expr)

Create Vector of feet using same notation for foot macro.
Separated by newlines.
First argument is stock, second is sum variable.

```julia
feetses = @feet begin
    A => B
    () => N
    C => ()
    D => E
    () => ()
    P => NI, R => NI, () => N
end
```
"""
macro feet(block::Expr)
    Base.remove_linenums!(block)
    @match block begin
        quote
          $((block...))
        end => map(create_foot, block) # also matches empty
        Expr(e, _...) => [create_foot(block)] # this also matches the above, so it's necessary this comes second.
    end
end


"""
    create_foot(block :: Expr)

Create foot (StockAndFlow0) using format A => B, where A is a stock and B is a sum variable.  Use () to represent no stock or sum variable.
To have multiple stocks or sum variables, chain together multiple pairs with commas.  Repeated occurences of the same symbol will be treated as the same stock or sum variable.
You cannot create distinct stocks or sum variables with the same name using this format. 

```julia
standard_foot = @foot A => N
emptyfoot = @foot () => ()
all_seir_links = @foot S => N, E => N, I => N, R => N
no_stocks = @foot () => N, () => NI, () => NS
no_sum = @foot A => ()
multiple_links = @foot A => B, A => B # will have two links from A to B.
```

"""
function create_foot(block::Expr)
    @match block.head begin

        :tuple => begin 
            if isempty(block.args) # case for create_foot(:()) 
                error("Cannot create foot with no arguments.")
            end
            foot_s = Vector{Symbol}()
            foot_sv = Vector{Symbol}()
            foot_ssv = Vector{Pair{Symbol, Symbol}}()

            for (s, sv, ssv) ∈ map(match_foot_format, block.args)
                if s != () push!(foot_s, s) end
                if sv != () push!(foot_sv, sv) end
                if ssv != () push!(foot_ssv, ssv) end
            end
            return foot(unique(foot_s), unique(foot_sv), foot_ssv)
        end
        :call => foot(match_foot_format(block)...)
        _ => error("Invalid expression type $(block.head).  Expecting tuple or call.")
    end
end

"""
    match_foot_format(footblock :: Expr)

Takes as argument an expression of the form A => B and returns a tuple in a format acceptable as arguments to create a foot.
Return type is Tuple{Union{Tuple{}, Symbol}, Union{Tuple{}, Symbol}, Union{Tuple{}, Pair{Symbol, Symbol}}}.  The empty tuple represents no stocks, no flows, or no links.

"""
function match_foot_format(footblock::Expr) 
    @match footblock begin
        :(()              => ())              => ((), (), ())
        :($(s :: Symbol)  => ())              => (s, (), ())
        :(()              => $(sv :: Symbol)) => ((), sv, ())
        :($(s :: Symbol)  => $(sv :: Symbol)) => (s, sv, s => sv)
        :($(s :: Symbol)  => sv)              => error("Non-symbolic second argument of foot: $sv")
        :($s              => $(sv :: Symbol)) => error("Non-symbolic first argument of foot: $s")
        :($s              => $sv)             => error("Foot definition requires symbolic names. Received: $s, $sv")
        Expr(:call, name, args...)            => error("Received: $name called with $args. Expected foot definition of form: A => B.")
        _                                     => error("Invalid foot definition.")
    end
end

const TEMP_STRAT_DEFAULT = :_

"""
Take an expression of the form a1, ..., => t <= s1, ..., where every element is a symbol, and return a 2-tuple of dictionaries of form ((a1 => t, a2 => t, ...), (s1 => t, ...))
"""
function interpret_stratification_notation(mapping_pair::Expr)
    @match mapping_pair begin

        :(_ => $t <= _) => return (Dict(TEMP_STRAT_DEFAULT => t), Dict(TEMP_STRAT_DEFAULT => t)) # match literal underscores to act as defaults (temporarily map to TEMP_STRAT_DEFAULT)
        :(_ => $t <= $a) => return (Dict(TEMP_STRAT_DEFAULT => t), Dict(a => t))
        :(_ => $t <= $a, $(atail...)) => return (Dict(TEMP_STRAT_DEFAULT => t), push!(Dict(as => t for as in atail), a => t))


        :($s => $t <= _) => return (Dict(s => t), Dict(TEMP_STRAT_DEFAULT => t))
        :($s => $t <= $a) => return (Dict(s => t), Dict(a => t))
        :($s => $t <= $a, $(atail...)) => return (Dict(s => t), push!(Dict(as => t for as in atail), a => t))


        :($(shead...), $s => $t <= _) => return (push!(Dict(ss => t for ss in shead), s => t), Dict(TEMP_STRAT_DEFAULT => t))
        :($(shead...), $s => $t <= $a) => return (push!(Dict(ss => t for ss in shead), s => t), Dict(a => t))


        # The most annoying case. :($(shead...), $s => $t <= $a, $(atail...))
        # basically just gave up here, if you can figure out a way to match it, godspeed.
        if mapping_pair.head == :tuple end => begin
            
            svec::Vector{Symbol} = []
            sdict::Dict{Symbol, Symbol} = Dict() # I don't think I need to initialize this here, but makes scope less ambiguous.
            adict::Dict{Symbol, Symbol} = Dict()
            found_t = false
            t_val = :NONE # this should never, ever be used.  Just need to initialize with something.

            for val in mapping_pair.args
                @match val begin # TODO: determine if there are more cases to deal with here
                    val::Symbol && if !found_t end => push!(svec, val)
                    val::Symbol && if found_t end => push!(adict, val => t_val)
                    :($s => $t <= $a) => begin
                        push!(svec, s)
                        sdict = Dict(prevs => t for prevs in svec)

                        t_val = t
                        found_t = true

                        push!(adict, a => t)

                    end
                    _ => error("Unknown expression found in stratification notation: $val")
                end
                    # val::Expr => begin # found the s => t <= a
                    #     @match val begin
                            
                    # end

            end
            return (sdict, adict)
                # _ => error("Unknown type found in stratification notation: $val of type $(typeof(val))")
        end
        _ => error("Unknown line format found in stratification notation.") 
    end
end


"""
Take 5 dictionaries:
s: strata dict, symbol => index
t: type dict, symbol => index
a: aggregate dict, symbol => index

ds: new strata, strata symbol => type symbol
da: new aggregate, aggregate symbol => type symbol

convert ds to strata index => type index, da to aggregate index => type index
"""
function substitute_symbols(s, t, a, ds, da) # TODO: add assert that all symbols in s,t and a are in ds and da

    map(x -> println(x), [s, t, a, ds, da])

    new_strata_dict = Dict(s[strata_symbol] => t[type_symbol] for (strata_symbol, type_symbol) in ds)
    new_aggregate_dict = Dict(a[aggregate_symbol] => t[type_symbol] for (aggregate_symbol, type_symbol) in da)
    return new_strata_dict, new_aggregate_dict
end


# TODO: the following
"""
    @stratify (strata, type, aggregate) begin ... end

    Ok, so the general idea here is:
    1. Grab all names from strata, type and aggregate, and create dictionaries which map them to their indices
    2. iterate over each line in the block
        2a. Split each line into a dictionary which maps all strata to that type and all aggregate to that type
        2b. Convert from two Symbol => Symbol dictionaries to two Int => Int dictionaries, using the dictionaries from step 1
        2c. Accumulate respective dictionaries
    3. Initialize an array of 0s for stocks, flows, parameters, dyvars and sums for strata and aggregate
    4. Insert into arrays all (nonwildcard)
    5. Do a once over of the arrays and replace wildcard values, replacing all 0s with index for wildcard, ensure there are no zeroes remaining
    6. Do some magic to infer LS, LSV, etc.
    7. Deal with attributes
    8. Construct strata -> type and aggregate -> type ACSetTransformations (Maybe not ACSet, because we don't care about attributes)
    9. Return pullback

"""
macro stratify(sf, block) # Trying to be very vigilant about catching errors.


    @assert sf.head == :tuple && length(sf.args) == 3
    @assert all(x -> x ∈ names(Main), sf.args) # Maybe want it to be not just in Main?

    strata, type, aggregate = map(x -> getfield(Main, x), sf.args) # Attempting to be clever and get around an eval call
    @assert all(x -> typeof(x) <: AbstractStockAndFlowStructureF, [strata, type, aggregate]) # Unsure if we want to be more or less strict with the type check

    Base.remove_linenums!(block)

    # STEP 1 and 2
       
    strata_snames = Dict(S => i for (i, S) in enumerate(snames(strata)))
    strata_svnames = Dict(S => i for (i, S) in enumerate(svnames(strata)))
    strata_vnames = Dict(S => i for (i, S) in enumerate(vnames(strata)))
    strata_fnames = Dict(S => i for (i, S) in enumerate(fnames(strata)))
    strata_pnames = Dict(S => i for (i, S) in enumerate(pnames(strata)))




    type_snames = Dict(S => i for (i, S) in enumerate(snames(type)))
    type_svnames = Dict(S => i for (i, S) in enumerate(svnames(type)))
    type_vnames = Dict(S => i for (i, S) in enumerate(vnames(type)))
    type_fnames = Dict(S => i for (i, S) in enumerate(fnames(type)))
    type_pnames = Dict(S => i for (i, S) in enumerate(pnames(type)))


    aggregate_snames = Dict(S => i for (i, S) in enumerate(snames(aggregate)))
    aggregate_svnames = Dict(S => i for (i, S) in enumerate(svnames(aggregate)))
    aggregate_vnames = Dict(S => i for (i, S) in enumerate(vnames(aggregate)))
    aggregate_fnames = Dict(S => i for (i, S) in enumerate(fnames(aggregate)))
    aggregate_pnames = Dict(S => i for (i, S) in enumerate(pnames(aggregate)))




    strata_stock_mappings_dict::Dict{Int, Int} = Dict()
    strata_flow_mappings_dict::Dict{Int, Int} = Dict()
    strata_dyvar_mappings_dict::Dict{Int, Int} = Dict()
    strata_param_mappings_dict::Dict{Int, Int} = Dict()
    strata_sum_mappings_dict::Dict{Int, Int} = Dict()
    

    aggregate_stock_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_flow_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_dyvar_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_param_mappings_dict::Dict{Int, Int} = Dict()
    aggregate_sum_mappings_dict::Dict{Int, Int} = Dict()
    
    



    strata_all_names = [strata_snames, strata_svnames, strata_vnames, strata_fnames, strata_pnames]

    # println(strata_all_names)


    @assert all(x -> TEMP_STRAT_DEFAULT ∉ keys(x) && allunique(x), strata_all_names)

    type_all_names = [type_snames, type_svnames, type_vnames, type_fnames, type_pnames]

    @assert all(x -> TEMP_STRAT_DEFAULT ∉ keys(x) && allunique(x), type_all_names)

    aggregate_all_names = [aggregate_snames, aggregate_svnames, aggregate_vnames, aggregate_fnames, aggregate_pnames]
    @assert all(x -> TEMP_STRAT_DEFAULT ∉ keys(x) && allunique(x), aggregate_all_names)

    # Inserting the symbol for wildcard into each dictionary
    map(x -> (push!(x, (TEMP_STRAT_DEFAULT => -1))), strata_all_names)
    map(x -> (push!(x, (TEMP_STRAT_DEFAULT => -1))), type_all_names) # TODO: Make this do something?  Maybe if there's only one stock, flow, etc.
    # ...in the type schema, you can replace it with an underscore?
    map(x -> (push!(x, (TEMP_STRAT_DEFAULT => -1))), aggregate_all_names)

    println(strata_all_names)
    println(strata_snames)


    # use 0 for uninitialized, -1 for map to default
    # indices start at 1, so both these are clearly placeholder values.



    # STEP 3

    # TODO: Extract this stuff out to reuse code instead of copy-pasting functions
    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = p -> begin
                    
                    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(p)
                    current_strata_dict, current_aggregate_dict = substitute_symbols(strata_snames, type_snames, aggregate_snames, current_strata_symbol_dict, current_aggregate_symbol_dict)


                    

                    # for (s,t) in current_strata_dict
                    println(keys(current_aggregate_dict))
                    println(strata_stock_mappings_dict)
                    @assert (all(x -> x ∉ keys(strata_stock_mappings_dict), keys(current_aggregate_dict)))
                    merge!(strata_stock_mappings_dict, current_aggregate_dict)
                    #     if strata_stock_mappings[s] != 0
                    #         error("Strata stock at index $s has already been assigned!")
                    #     else
                    #         strata_stock_mappings[s] = t
                    #     end
                    # end

                    @assert (all(x -> x ∉ keys(aggregate_stock_mappings_dict), keys(current_aggregate_dict)))
                    merge!(aggregate_stock_mappings_dict, current_aggregate_dict)

                    # for (a,t) in current_aggregate_dict
                    #     if aggregate_stock_mappings[a] != 0
                    #         error("Aggregate stock at index $s has already been assigned!")
                    #     else
                    #         aggregate_stock_mappings[a] = t
                    #     end
                    # end

                end
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> begin                
                    
                    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(p)
                    current_strata_dict, current_aggregate_dict = substitute_symbols(strata_pnames, type_pnames, aggregate_pnames, current_strata_symbol_dict, current_aggregate_symbol_dict)

                    @assert (all(x -> x ∉ keys(strata_param_mappings_dict), keys(current_aggregate_dict)))
                    merge!(strata_param_mappings_dict, current_aggregate_dict)

                    @assert (all(x -> x ∉ keys(aggregate_param_mappings_dict), keys(current_aggregate_dict)))
                    merge!(aggregate_param_mappings_dict, current_aggregate_dict)


                    # for (s,t) in current_strata_dict
                    #     if strata_param_mappings[s] != 0
                    #         error("Strata parameter at index $s has already been assigned!")
                    #     else
                    #         strata_param_mappings[s] = t
                    #     end
                    # end

                    # for (a,t) in current_aggregate_dict
                    #     if aggregate_param_mappings[a] != 0
                    #         error("Aggregate parameter at index $s has already been assigned!")
                    #     else
                    #         aggregate_param_mappings[a] = t
                    #     end
                    # end
                end
            end
            QuoteNode(:dynamic_variables) => begin
                    current_phase = p -> begin                
                    
                     
                    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(p)
                    current_strata_dict, current_aggregate_dict = substitute_symbols(strata_vnames, type_vnames, aggregate_vnames, current_strata_symbol_dict, current_aggregate_symbol_dict)

                    @assert (all(x -> x ∉ keys(strata_dyvar_mappings_dict), keys(current_aggregate_dict)))
                    merge!(strata_dyvar_mappings_dict, current_aggregate_dict)

                    @assert (all(x -> x ∉ keys(aggregate_dyvar_mappings_dict), keys(current_aggregate_dict)))
                    merge!(aggregate_dyvar_mappings_dict, current_aggregate_dict)
                    
                    # for (s,t) in current_strata_dict
                    #     if strata_dyvar_mappings[s] != 0
                    #         error("Strata dyvar at index $s has already been assigned!")
                    #     else
                    #         strata_dyvar_mappings[s] = t
                    #     end
                    # end

                    # for (a,t) in current_aggregate_dict
                    #     if aggregate_dyvar_mappings[a] != 0
                    #         error("Aggregate dyvar at index $s has already been assigned!")
                    #     else
                    #         aggregate_dyvar_mappings[a] = t
                    #     end
                    # end
                end
            end            
            QuoteNode(:flows) => begin
                current_phase = p -> begin                
                        
                        
                 
                    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(p)
                    current_strata_dict, current_aggregate_dict = substitute_symbols(strata_fnames, type_fnames, aggregate_fnames, current_strata_symbol_dict, current_aggregate_symbol_dict)

                    @assert (all(x -> x ∉ keys(strata_flow_mappings_dict), keys(current_aggregate_dict)))
                    merge!(strata_flow_mappings_dict, current_aggregate_dict)

                    @assert (all(x -> x ∉ keys(aggregate_flow_mappings_dict), keys(current_aggregate_dict)))
                    merge!(aggregate_flow_mappings_dict, current_aggregate_dict)

                    # for (s,t) in current_strata_dict
                    #     if strata_flow_mappings[s] != 0
                    #         error("Strata flow at index $s has already been assigned!")
                    #     else
                    #         strata_flow_mappings[s] = t
                    #     end
                    # end

                    # for (a,t) in current_aggregate_dict
                    #     if aggregate_flow_mappings[a] != 0
                    #         error("Aggregate flow at index $s has already been assigned!")
                    #     else
                    #         aggregate_flow_mappings[a] = t
                    #     end
                    # end
                end
            end                    
                  
            QuoteNode(:sums) => begin
                current_phase = p -> begin                

                    current_strata_symbol_dict, current_aggregate_symbol_dict = interpret_stratification_notation(p)
                    current_strata_dict, current_aggregate_dict = substitute_symbols(strata_svnames, type_svnames, aggregate_svnames, current_strata_symbol_dict, current_aggregate_symbol_dict)
                    

                    @assert (all(x -> x ∉ keys(strata_sum_mappings_dict), keys(current_aggregate_dict)))
                    merge!(strata_sum_mappings_dict, current_aggregate_dict)

                    @assert (all(x -> x ∉ keys(aggregate_sum_mappings_dict), keys(current_aggregate_dict)))
                    merge!(aggregate_sum_mappings_dict, current_aggregate_dict)
                    
                    # for (s,t) in current_strata_dict
                    #     if strata_sum_mappings[s] != 0
                    #         error("Strata sum at index $s has already been assigned!")
                    #     else
                    #         strata_sum_mappings[s] = t
                    #     end
                    # end
                    
                    # for (a,t) in current_aggregate_dict
                    #     if aggregate_sum_mappings[a] != 0
                    #         error("Aggregate sum at index $s has already been assigned!")
                    #     else
                    #         aggregate_sum_mappings[a] = t
                    #     end
                    # end
                end
            end                    


            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    # STEP 4


    # println(strata_stock_mappings)
    # println("I'm sorry.")
    default_index_strata_stock = -1 ∈ keys(strata_stock_mappings_dict) ? strata_stock_mappings_dict[-1] : 0
    default_index_strata_flow = -1 ∈ keys(strata_flow_mappings_dict) ? strata_flow_mappings_dict[-1] : 0
    default_index_strata_dyvar = -1 ∈ keys(strata_dyvar_mappings_dict) ? strata_dyvar_mappings_dict[-1] : 0
    default_index_strata_param = -1 ∈ keys(strata_param_mappings_dict) ? strata_param_mappings_dict[-1] : 0
    default_index_strata_sum = -1 ∈ keys(strata_sum_mappings_dict) ? strata_sum_mappings_dict[-1] : 0

    default_index_aggregate_stock = -1 ∈ keys(aggregate_stock_mappings_dict) ? aggregate_stock_mappings_dict[-1] : 0
    default_index_aggregate_flow = -1 ∈ keys(aggregate_flow_mappings_dict) ? aggregate_flow_mappings_dict[-1] : 0
    default_index_aggregate_dyvar = -1 ∈ keys(aggregate_dyvar_mappings_dict) ? aggregate_dyvar_mappings_dict[-1] : 0
    default_index_aggregate_param = -1 ∈ keys(aggregate_param_mappings_dict) ? aggregate_param_mappings_dict[-1] : 0
    default_index_aggregate_sum = -1 ∈ keys(aggregate_sum_mappings_dict) ? aggregate_sum_mappings_dict[-1] : 0


    strata_stock_mappings = [get(strata_stock_mappings_dict, i, default_index_strata_stock)  for i in 1:ns(strata)]
    strata_flow_mappings =  [get(strata_flow_mappings_dict, i, default_index_strata_flow)  for i in 1:nf(strata)]
    strata_dyvar_mappings::Vector{Int} = [get(strata_dyvar_mappings_dict, i, default_index_strata_dyvar)  for i in 1:nvb(strata)]
    strata_param_mappings::Vector{Int} = [get(strata_param_mappings_dict, i, default_index_strata_param)  for i in 1:np(strata)]
    strata_sum_mappings::Vector{Int} =  [get(strata_sum_mappings_dict, i, default_index_strata_sum)  for i in 1:nsv(strata)]

    aggregate_stock_mappings = [get(aggregate_stock_mappings_dict, i, default_index_aggregate_stock)  for i in 1:ns(aggregate)]
    aggregate_flow_mappings =  [get(aggregate_flow_mappings_dict, i, default_index_aggregate_flow)  for i in 1:nf(aggregate)]
    aggregate_dyvar_mappings::Vector{Int} = [get(aggregate_dyvar_mappings_dict, i, default_index_aggregate_dyvar)  for i in 1:nvb(aggregate)]
    aggregate_param_mappings::Vector{Int} = [get(aggregate_param_mappings_dict, i, default_index_aggregate_param)  for i in 1:np(aggregate)]
    aggregate_sum_mappings::Vector{Int} =  [get(aggregate_sum_mappings_dict, i, default_index_aggregate_sum)  for i in 1:nsv(aggregate)]
    

    # strata_stock_mappings = replace(strata_stock_mappings, 0 => default_index_strata_stock)
    # aggregate_stock_mappings = replace(aggregate_stock_mappings, 0 => default_index_aggregate_stock)

    # strata_flow_mappings = replace(strata_flow_mappings, 0 => default_index_strata_flow)
    # aggregate_flow_mappings = replace(aggregate_flow_mappings, -1 => default_index_aggregate_flow)

    # strata_dyvar_mappings = replace(strata_dyvar_mappings, -1 => default_index_strata_dyvar)
    # aggregate_dyvar_mappings = replace(aggregate_dyvar_mappings, -1 => default_index_aggregate_dyvar)

    # strata_param_mappings = replace(strata_param_mappings, -1 => default_index_strata_param)
    # aggregate_param_mappings = replace(aggregate_param_mappings, -1 => default_index_aggregate_param)

    # strata_sum_mappings = replace(strata_sum_mappings, -1 => default_index_strata_sum)
    # aggregate_sum_mappings = replace(aggregate_sum_mappings, -1 => default_index_aggregate_sum)



    println(strata_stock_mappings)
    println(aggregate_stock_mappings)

    println(strata_flow_mappings)
    println(aggregate_flow_mappings)

    println(strata_dyvar_mappings)
    println(aggregate_dyvar_mappings)

    println(strata_param_mappings)
    println(aggregate_param_mappings)

    println(strata_sum_mappings)
    println(aggregate_sum_mappings)


    # TODO: Once over at the end to convert placeholders.

    # s = StockAndFlowBlock(stocks, params, dyvars, flows, sums)
    # return s
end


end




  
