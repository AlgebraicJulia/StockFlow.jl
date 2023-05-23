using Test
using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax: is_binop_or_unary, sum_variables, infix_expression_to_binops, fnone_value_or_vector, extract_function_name_and_args_expr

@testset "is_binop_or_unary recognises binops" begin
    @test is_binop_or_unary(:(a + b))
    @test is_binop_or_unary(:(f(a, b)))
end
@testset "is_binop_or_unary recognises non-binops as non-binops" begin
    @test !is_binop_or_unary(:(f()))
    @test !is_binop_or_unary(:(a + b + c))
    @test !is_binop_or_unary(:(f(a, b, c)))
end

@testset "sum_variables" begin
    @test sum_variables([]) == []
    @test sum_variables([(:a, 1)]) == [:a]
    @test sum_variables([(:a, 1), (:b, 2)]) == [:a, :b]
end

@testset "infix_expression_to_binops does nothing to binops and unary exprs" begin
    @test infix_expression_to_binops(:(f(a)))[1][1][2] == :(f(a))
    @test infix_expression_to_binops(:(f(a, b)))[1][1][2] == :(f(a, b))
    @test infix_expression_to_binops(:(a + b))[1][1][2] == :(a + b)
end

@testset "infix_expression_to_binops creates right number of expressions" begin
    @test length(infix_expression_to_binops(:(a + b + c))[1]) == 2
    @test length(infix_expression_to_binops(:(a + b + c + d))[1]) == 3
    @test length(infix_expression_to_binops(:(a + b + c + d + e))[1]) == 4
    @test length(infix_expression_to_binops(:(a + b + c + d + e + f))[1]) == 5
end

@testset "infix_expression_to_binops throws exception when no binops provided" begin
    @test_throws Exception infix_expression_to_binops(:(f()))
end

@testset "infix_expression_to_binops uses final symbol" begin
    @test infix_expression_to_binops(:(f(a, b)); finalsym=:testsym)[2] == :testsym
    @test infix_expression_to_binops(:(a + b); finalsym=:testsym)[2] == :testsym
end

@testset "fnone_value_or_vector" begin
    empty_symbol_vector::Vector{Symbol} = []
    @test fnone_value_or_vector(empty_symbol_vector) == :F_NONE
    @test fnone_value_or_vector([:a]) == :a
    @test fnone_value_or_vector([:a, :b]) == [:a, :b]
end

@testset "extract_function_name_args_expr extracts the flow name and flow definition" begin
    @test extract_function_name_and_args_expr(:(testf(a))) == (:testf, :a)
    @test extract_function_name_and_args_expr(:(testf(a + b + c + d))) == (:testf, :(a + b + c + d))
end

@testset "extract_function_name_args_expr rejects invalid flow expressions" begin
    # Undefined flow equation
    @test_throws Exception extract_function_name_and_args_expr(:(testf()))
    # Multiparameter flows
    @test_throws Exception extract_function_name_and_args_expr(:(testf(a, b)))
end

@testset "model allows uses unary functions like log" begin
    SIR_1_via_macro = @stock_and_flow begin
        :stocks
        S
        I
        R

        :parameters
        c
        beta
        tRec

        :dynamic_variables
        v_prevalence = exp(I)
        v_meanInfectiousContactsPerS = c * v_prevalence
        v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS
        v_newInfections = S * v_perSIncidenceRate
        v_newRecovery = log(I)

        :flows
        S => inf(v_newInfections) => I
        I => rec(v_newRecovery) => R

        :sums
        N = [S, I, R]
    end

    SIR_1_canonical = StockAndFlowF(
        # stocks
        (:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),
        # parameters
        (:c, :beta, :tRec),
        # dynamical variables
        (:v_prevalence => (:I => :exp),
            :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*),
            :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
            :v_newInfections => ((:S, :v_perSIncidenceRate) => :*),
            :v_newRecovery => (:I => :log),
        ),
        # flows
        (:inf => :v_newInfections, :rec => :v_newRecovery),
        # sum dynamical variables
        (:N),
    )
    @test SIR_1_via_macro == SIR_1_canonical
end
@testset "stock_and_flow macro generates the expected StockAndFlowF representations" begin
    SIR_1_via_macro = @stock_and_flow begin
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

    SIR_1_canonical = StockAndFlowF(
        # stocks
        (:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),
        # parameters
        (:c, :beta, :tRec),
        # dynamical variables
        (:v_prevalence => ((:I, :N) => :/),
            :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*),
            :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
            :v_newInfections => ((:S, :v_perSIncidenceRate) => :*),
            :v_newRecovery => ((:I, :tRec) => :/),
        ),
        # flows
        (:inf => :v_newInfections, :rec => :v_newRecovery),
        # sum dynamical variables
        (:N),
    )
    @test SIR_1_via_macro == SIR_1_canonical

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

    # Although the variable names are different
    # due to gensym, the models should structurally be the same.
    for part in SIR_2.parts
        @test SIR_2.parts[part] == SIR_1_canonical.parts[part]
    end
end

@testset "stock_and_flow macro base case" begin
    SIR_1_via_macro = @stock_and_flow begin
    end
end
