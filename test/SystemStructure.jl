
import Pkg
Pkg.activate("MyProject")
Pkg.develop(path="/home/sodium/Documents/JuliaShare/StockFlow.jl")
Pkg.instantiate()

using Test
using StockFlow
using StockFlow.Syntax
import Base.deepcopy
using Catlab.CategoricalAlgebra.StructuredCospans



empty = @stock_and_flow begin end

p = @stock_and_flow begin 
    :stocks
    S
    I
    
    :dynamic_variables
    v = S + I
    
    :parameters
    p1
    p2
    
    :flows
    S => f1(v) => CLOUD
    
    :sums
    N = [S]
    NI = [I]
end

p_prefixed = @stock_and_flow begin 
    :stocks
    prefS
    prefI
    
    :dynamic_variables
    prefv = prefS + prefI
    
    :parameters
    prefp1
    prefp2
    
    :flows
    prefS => preff1(prefv) => CLOUD
    
    :sums
    prefN = [prefS]
    prefNI = [prefI]
end

p_suffixed = @stock_and_flow begin 
    :stocks
    Ssuf
    Isuf
    
    :dynamic_variables
    vsuf = Ssuf + Isuf
    
    :parameters
    p1suf
    p2suf
    
    :flows
    Ssuf => f1suf(vsuf) => CLOUD
    
    :sums
    Nsuf = [Ssuf]
    NIsuf = [Isuf]
end

empty_foot = foot((),(),())

footA = foot(:S, :N, :S => :N)
footA_pref = foot(:prefS, :prefN, :prefS => :prefN)
footA_suf = foot(:Ssuf, :Nsuf, :Ssuf => :Nsuf)

footB = foot(:prefI, :prefNI, :prefI => :prefNI)

OpenA = Open(deepcopy(p), deepcopy(footA))
OpenA_pref = Open(deepcopy(p_prefixed), deepcopy(footA_pref))
OpenB = Open(deepcopy(p_prefixed), deepcopy(footB))

# function deepcopy(open::StructuredCospan)
#     return Open(deepcopy(apex(open)), deepcopy(feet(open)))
# end


@testset "changing names act the same as if the stock flow/foot/open was created with the changed name" begin
    @test empty == add_suffix!(deepcopy(empty), "AAA") == add_prefix!(deepcopy(empty), "BBB")
    @test add_prefix!(deepcopy(p), "pref") == p_prefixed
    @test add_suffix!(deepcopy(p), "suf") == p_suffixed
    
    @test empty_foot == add_suffix!(deepcopy(empty_foot), "CCC") == add_prefix!(deepcopy(empty_foot), "DDD")
    @test add_prefix!(deepcopy(footA), "pref") == footA_pref
    @test add_suffix!(deepcopy(footA), "suf") == footA_suf

    @test add_prefix!(deepcopy(OpenA), "pref") == OpenA_pref


end

# @testset "is_binop_or_unary recognises binops" begin
#     @test is_binop_or_unary(:(a + b))
#     @test is_binop_or_unary(:(f(a, b)))
#     @test is_binop_or_unary(:(1.0 + x))
# end
# @testset "is_binop_or_unary recognises non-binops as non-binops" begin
#     @test !is_binop_or_unary(:(f()))
#     @test !is_binop_or_unary(:(a + b + c))
#     @test !is_binop_or_unary(:(f(a, b, c)))
# end

# @testset "sum_variables" begin
#     @test sum_variables([]) == []
#     @test sum_variables([(:a, 1)]) == [:a]
#     @test sum_variables([(:a, 1), (:b, 2)]) == [:a, :b]
# end

# @testset "infix_expression_to_binops does nothing to binops and unary exprs" begin
#     @test infix_expression_to_binops(:(f(a)))[1][1][2] == :(f(a))
#     @test infix_expression_to_binops(:(f(a, b)))[1][1][2] == :(f(a, b))
#     @test infix_expression_to_binops(:(a + b))[1][1][2] == :(a + b)
# end

# @testset "infix_expression_to_binops creates right number of expressions" begin
#     @test length(infix_expression_to_binops(:(a + b + c))[1]) == 2
#     @test length(infix_expression_to_binops(:(a + b + c + d))[1]) == 3
#     @test length(infix_expression_to_binops(:(a + b + c + d + e))[1]) == 4
#     @test length(infix_expression_to_binops(:(a + b + c + d + e + f))[1]) == 5
# end

# @testset "infix_expression_to_binops throws exception when no binops provided" begin
#     @test_throws Exception infix_expression_to_binops(:(f()))
# end

# @testset "infix_expression_to_binops uses final symbol" begin
#     @test infix_expression_to_binops(:(f(a, b)); finalsym=:testsym)[2] == :testsym
#     @test infix_expression_to_binops(:(a + b); finalsym=:testsym)[2] == :testsym
# end

# @testset "fnone_value_or_vector" begin
#     empty_symbol_vector::Vector{Symbol} = []
#     @test fnone_value_or_vector(empty_symbol_vector) == :F_NONE
#     @test fnone_value_or_vector([:a]) == :a
#     @test fnone_value_or_vector([:a, :b]) == [:a, :b]
# end

# @testset "extract_function_name_args_expr extracts the flow name and flow definition" begin
#     @test extract_function_name_and_args_expr(:(testf(a))) == (:testf, :a)
#     @test extract_function_name_and_args_expr(:(testf(a + b + c + d))) == (:testf, :(a + b + c + d))
# end

# @testset "extract_function_name_args_expr rejects invalid flow expressions" begin
#     # Undefined flow equation
#     @test_throws Exception extract_function_name_and_args_expr(:(testf()))
#     # Multiparameter flows
#     @test_throws Exception extract_function_name_and_args_expr(:(testf(a, b)))
# end

# #@testset "non-variable parameters in functions" begin
# #  @stock_and_flow begin
# #      :stocks
# #       A
# #       B
# #       C
# #       :parameters
# #       x
# #       y
# #       z
# #       :dynamic_variables
# # TODO? f = 1.0 - x
# #       :flows
# #       A => fname(f) => B
# #  end
# #end

# @testset "model allows uses unary functions like log" begin
#     SIR_1_via_macro = @stock_and_flow begin
#         :stocks
#         S
#         I
#         R

#         :parameters
#         c
#         beta
#         tRec

#         :dynamic_variables
#         v_prevalence = exp(I)
#         v_meanInfectiousContactsPerS = c * v_prevalence
#         v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS
#         v_newInfections = S * v_perSIncidenceRate
#         v_newRecovery = log(I)

#         :flows
#         S => inf(v_newInfections) => I
#         I => rec(v_newRecovery) => R

#         :sums
#         N = [S, I, R]
#     end

#     SIR_1_canonical = StockAndFlowF(
#         # stocks
#         (:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),
#         # parameters
#         (:c, :beta, :tRec),
#         # dynamical variables
#         (:v_prevalence => (:I => :exp),
#          :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*),
#          :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
#          :v_newInfections => ((:S, :v_perSIncidenceRate) => :*),
#          :v_newRecovery => (:I => :log),
#         ),
#         # flows
#         (:inf => :v_newInfections, :rec => :v_newRecovery),
#         # sum dynamical variables
#         (:N),
#     )
#     @test SIR_1_via_macro == SIR_1_canonical
# end
# @testset "stock_and_flow macro generates the expected StockAndFlowF representations" begin
#     SIR_1_via_macro = @stock_and_flow begin
#         :stocks
#         S
#         I
#         R

#         :parameters
#         c
#         beta
#         tRec

#         :dynamic_variables
#         v_prevalence = I / N
#         v_meanInfectiousContactsPerS = c * v_prevalence
#         v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS
#         v_newInfections = S * v_perSIncidenceRate
#         v_newRecovery = I / tRec

#         :flows
#         S => inf(v_newInfections) => I
#         I => rec(v_newRecovery) => R

#         :sums
#         N = [S, I, R]
#     end

#     SIR_1_canonical = StockAndFlowF(
#         # stocks
#         (:S => (:F_NONE, :inf, :N), :I => (:inf, :rec, :N), :R => (:rec, :F_NONE, :N)),
#         # parameters
#         (:c, :beta, :tRec),
#         # dynamical variables
#         (:v_prevalence => ((:I, :N) => :/),
#          :v_meanInfectiousContactsPerS => ((:c, :v_prevalence) => :*),
#          :v_perSIncidenceRate => ((:beta, :v_meanInfectiousContactsPerS) => :*),
#          :v_newInfections => ((:S, :v_perSIncidenceRate) => :*),
#          :v_newRecovery => ((:I, :tRec) => :/),
#         ),
#         # flows
#         (:inf => :v_newInfections, :rec => :v_newRecovery),
#         # sum dynamical variables
#         (:N),
#     )
#     @test SIR_1_via_macro == SIR_1_canonical

#     SIR_2 = @stock_and_flow begin
#         :stocks
#         S
#         I
#         R

#         :parameters
#         c
#         beta
#         tRec

#         # We can leave out dynamic variables and let them be inferred from flows entirely!

#         :flows
#         S => inf(S * beta * (c * (I / N))) => I
#         I => rec(I / tRec) => R

#         :sums
#         N = [S, I, R]
#     end

#     # Although the variable names are different
#     # due to gensym, the models should structurally be the same.
#     for part in SIR_2.parts
#         @test SIR_2.parts[part] == SIR_1_canonical.parts[part]
#     end
# end

# @testset "stock_and_flow macro base cases" begin
#     empty_via_macro = @stock_and_flow begin end
#     @test empty_via_macro == StockAndFlowF()

#     no_sums = @stock_and_flow begin
#         :stocks
#         A
#         B
#         C

#         :parameters
#         p
#         q

#         :dynamic_variables
#         dyvar1 = A + B
#         dyvar2 = B * C
#         dyvar3 = sqrt(q)
#         dyvar4 = exp(p)
#         dyvar5 = log(dyvar3, dyvar4)

#         :flows
#         A => f1(dyvar1) => B
#         B => f2(dyvar2) => C
#         C => f3(dyvar5) => A
#     end
#     no_sums_canonical = StockAndFlowF(
#         #stocks
#         (:A => (:f3, :f1, :SV_NONE), :B => (:f1, :f2, :SV_NONE), :C => (:f2, :f3, :SV_NONE)),
#         #params
#         (:p, :q),
#         # dyvars
#         (:dyvar1 => ((:A, :B) => :+),
#          :dyvar2 => ((:B, :C) => :*),
#          :dyvar3 => (:q => :sqrt),
#          :dyvar4 => (:p => :exp),
#          :dyvar5 => ((:dyvar3, :dyvar4) => :log)),
#         #flows
#         (:f1 => :dyvar1,
#          :f2 => :dyvar2,
#          :f3 => :dyvar5),
#         ())
#     @test no_sums == no_sums_canonical
# end
# @testset "is_recursive_dyvar detects recursive dyvars" begin
#     @test is_recursive_dyvar(:v, :(v + v))
#     @test is_recursive_dyvar(:a, :(b + c + d / (e + f + g / (h + i + j / (a - b * a)))))
# end
# @testset "is_recursive_dyvar does not flag non-recursive dyvars" begin
#     @test !is_recursive_dyvar(:w, :(v + v))
#     @test !is_recursive_dyvar(:z, :(b + c + d / (e + f + g / (h + i + j / (a - b * a)))))

# end
# @testset "recursive definitions should be disallowed" begin
#     expr = quote
#         :dynamic_variables
#         v = v + v
#     end
#     # When used as a macro -- @stock_and_flow -- this exception is thrown
#     # at a point that @test_throws cannot capture it.
#     @test_throws Exception stock_and_flow(expr)
# end
