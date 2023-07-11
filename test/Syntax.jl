import Pkg
Pkg.activate("MyProject")
Pkg.develop(path="/home/silicon/Documents/Git/StockFlow.jl/")
Pkg.instantiate()

using Base: is_unary_and_binary_operator
using Test
using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax: is_binop_or_unary, sum_variables, infix_expression_to_binops, fnone_value_or_vector, extract_function_name_and_args_expr, is_recursive_dyvar, create_foot

@testset "is_binop_or_unary recognises binops" begin
    @test is_binop_or_unary(:(a + b))
    @test is_binop_or_unary(:(f(a, b)))
    @test is_binop_or_unary(:(1.0 + x))
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

#@testset "non-variable parameters in functions" begin
#  @stock_and_flow begin
#      :stocks
#       A
#       B
#       C
#       :parameters
#       x
#       y
#       z
#       :dynamic_variables
# TODO? f = 1.0 - x
#       :flows
#       A => fname(f) => B
#  end
#end

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

@testset "stock_and_flow macro base cases" begin
    empty_via_macro = @stock_and_flow begin end
    @test empty_via_macro == StockAndFlowF()

    no_sums = @stock_and_flow begin
        :stocks
        A
        B
        C

        :parameters
        p
        q

        :dynamic_variables
        dyvar1 = A + B
        dyvar2 = B * C
        dyvar3 = sqrt(q)
        dyvar4 = exp(p)
        dyvar5 = log(dyvar3, dyvar4)

        :flows
        A => f1(dyvar1) => B
        B => f2(dyvar2) => C
        C => f3(dyvar5) => A
    end
    no_sums_canonical = StockAndFlowF(
        #stocks
        (:A => (:f3, :f1, :SV_NONE), :B => (:f1, :f2, :SV_NONE), :C => (:f2, :f3, :SV_NONE)),
        #params
        (:p, :q),
        # dyvars
        (:dyvar1 => ((:A, :B) => :+),
         :dyvar2 => ((:B, :C) => :*),
         :dyvar3 => (:q => :sqrt),
         :dyvar4 => (:p => :exp),
         :dyvar5 => ((:dyvar3, :dyvar4) => :log)),
        #flows
        (:f1 => :dyvar1,
         :f2 => :dyvar2,
         :f3 => :dyvar5),
        ())
    @test no_sums == no_sums_canonical
end
@testset "is_recursive_dyvar detects recursive dyvars" begin
    @test is_recursive_dyvar(:v, :(v + v))
    @test is_recursive_dyvar(:a, :(b + c + d / (e + f + g / (h + i + j / (a - b * a)))))
end
@testset "is_recursive_dyvar does not flag non-recursive dyvars" begin
    @test !is_recursive_dyvar(:w, :(v + v))
    @test !is_recursive_dyvar(:z, :(b + c + d / (e + f + g / (h + i + j / (a - b * a)))))

end
@testset "recursive definitions should be disallowed" begin
    expr = quote
        :dynamic_variables
        v = v + v
    end
    # When used as a macro -- @stock_and_flow -- this exception is thrown
    # at a point that @test_throws cannot capture it.
    @test_throws Exception stock_and_flow(expr)
end



@testset "foot syntax can create all types of feet" begin
    @test (@foot A => B) == foot(:A, :B, :A => :B)
    @test (@foot P => ()) == foot(:P, (), ())
    @test (@foot () => Q) == foot((), :Q, ())
    @test (@foot () => ()) == foot((),(),())

    @test (@foot =>((), SV)) == foot((),:SV,()) 
    @test (@foot A11 => B22) == foot(:A11, :B22, :A11 => :B22)
end

@testset "foot syntax disallows invalid feet" begin # note, @feet calls create_foot for each line, so this should apply to both @foot and @feet
    @test_throws Exception create_foot(:(A => B => C)) # Invalid syntax for second argument of foot: B => C
    @test_throws Exception create_foot(:(oooo2 + f => C)) # Invalid syntax for first argument of foot: oooo2 + f
    @test_throws Exception create_foot(:(A + B)) # Invalid syntax function for foot: +
    @test_throws Exception create_foot(:(=>)) # no method matching create_foot(::Symbol)
    @test_throws Exception create_foot(:(=>(A, B, C, D)))
    @test_throws Exception create_foot(:())
end

@testset "feet syntax can create feet" begin

    @test (@feet begin
        
        A => B
        C => D
        () => ()

        P => ()
        () => Q

    end) == [foot(:A, :B, :A => :B), foot(:C, :D, :C => :D), foot((),(),()), foot(:P, (),()), foot((),:Q,())]

    @test (@feet P => Q) == [foot(:P, :Q, :P => :Q)]
    @test (@feet begin end) == Vector{StockAndFlow0}()
    @test (@feet begin P => Q; L => R end) == [foot(:P, :Q, :P => :Q), foot(:L, :R, :L => :R)]

    @test (@feet P => P) == [foot(:P, :P, :P => :P)]
end

@testset "feet syntax fails on invalid feet" begin # mostly to check that an exception is thrown even if some of the feet are valid.
    @test_throws Exception @eval @feet A => B => C # eval required so the errors occur at runtime, rather than at compilation
    @test_throws Exception @eval @feet begin A => B; =>(D,E,F) end
    @test_throws Exception @eval @feet begin A => B; 1 => 2; end
end



