using Test
using StockFlow.Syntax: is_binop, sum_variables, infix_expression_to_binops

@testset "is_binop recognises binops" begin
    @test is_binop(:(a + b))
    @test is_binop(:(f(a, b)))
end
@testset "is_binop recognises non-binops as non-binops" begin
    @test !is_binop(:(f()))
    @test !is_binop(:(f(a)))
    @test !is_binop(:(a + b + c))
    @test !is_binop(:(f(a, b, c)))
end

@testset "sum_variables" begin
    @test sum_variables([]) == []
    @test sum_variables([(:a, 1)]) == [:a]
    @test sum_variables([(:a, 1), (:b, 2)]) == [:a, :b]
end

@testset "infix_expression_to_binops does nothing to binops" begin
    @test infix_expression_to_binops(:(f(a, b)))[1][1][2] == :(f(a, b))
    @test infix_expression_to_binops(:(a + b))[1][1][2] == :(a + b)
end

@testset "infix_expression_to_binops throws exception when no binops provided" begin
    @test_throws Exception infix_expression_to_binops(:(f()))
    @test_throws Exception infix_expression_to_binops(:(f(a)))
end

@testset "infix_expression_to_binops uses final symbol" begin
    @test infix_expression_to_binops(:(f(a, b)); finalsym=:testsym)[2] == :testsym
    @test infix_expression_to_binops(:(a + b); finalsym=:testsym)[2] == :testsym
end
