using Base: is_unary_and_binary_operator
using Test
using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax: is_binop_or_unary, sum_variables, infix_expression_to_binops, fnone_value_or_vector, extract_function_name_and_args_expr, is_recursive_dyvar, create_foot, apply_flags, substitute_symbols

@testset "Stratification DSL" begin
    include("Stratification.jl")
end

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
    for part in keys(SIR_2.parts)
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

    @test (@foot () => B, A => ()) == foot(:A, :B, ())
    @test (@foot A => B, A => C) == foot(:A, (:B, :C), (:A => :B, :A => :C))
    @test (@foot A => B, A => B, A => B) == foot(:A, :B, (:A => :B, :A => :B, :A => :B)) # at present, it deduplicates stocks and sums, but not links.
    @test (@foot P => Q, R => (), () => ()) == foot((:P, :R), (:Q), (:P => :Q))
    @test (@foot () => (), () => ()) == foot((), (), ())
end

@testset "foot syntax disallows invalid feet" begin # note, @feet calls create_foot for each line, so this should apply to both @foot and @feet
    @test_throws Exception create_foot(:(A => B => C)) # Invalid syntax for second argument of foot: B => C
    @test_throws Exception create_foot(:(oooo2 + f => C)) # Invalid syntax for first argument of foot: oooo2 + f
    @test_throws Exception create_foot(:(A + B)) # Invalid syntax function for foot: +
    @test_throws Exception create_foot(:(=>)) # no method matching create_foot(::Symbol)
    @test_throws Exception create_foot(:(=>(A, B, C, D)))
    @test_throws Exception create_foot(:())
    @test_throws Exception create_foot(:(([]) => ()))

    @test_throws Exception create_foot(:(A => B, P => Q, C))
    @test_throws Exception create_foot(:(() => E, () => (K,)))

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

    @test (@feet begin A => NS, A => N; B => NS; C => (), D => () end) == [
        foot((:A), (:NS, :N), (:A => :NS, :A => :N)),
        foot((:B), (:NS), (:B => :NS)),
        foot((:C, :D), (), ())
        ]

    @test (@feet begin
        P => Q, R => ()
        () => ()
        J => K, J => Q
    end) == [
        foot((:P, :R), :Q, :P => :Q),
        foot((), (), ()),
        foot(:J, (:K, :Q), (:J => :K, :J => :Q))
    ]

end

@testset "feet syntax fails on invalid feet" begin # mostly to check that an exception is thrown even if some of the feet are valid.
    @test_throws Exception @eval @feet A => B => C # eval required so the errors occur at runtime, rather than at compilation
    @test_throws Exception @eval @feet begin A => B; =>(D,E,F) end
    @test_throws Exception @eval @feet begin A => B; 1 => 2; end
end

###########################

@testset "infer_links works as expected" begin
    # No prior mappings means no inferred mappings
    @test (infer_links(StockAndFlowF(), StockAndFlowF(), Dict{Symbol, Vector{Int64}}(:S => [], :F => [], :SV => [], :P => [], :V => []))
        == Dict(:LS => [], :LSV => [], :LV => [], :I => [], :O => [], :LPV => [], :LVV => []))

    # S: 1 -> 1 and SV: 1 -> 1 implies LS: 1 -> 1
    @test (infer_links(
        (@stock_and_flow begin; :stocks; A; :sums; NA = [A]; end), 
        (@stock_and_flow begin; :stocks; B; :sums; NB = [B]; end),
        Dict{Symbol, Vector{Int64}}(:S => [1], :F => [], :SV => [1], :P => [], :V => []))
    == Dict(:LS => [1], :LSV => [], :LV => [], :I => [], :O => [], :LPV => [], :LVV => []))

    # annoying exanmple, required me to add code to disambiguate using position
    # that is, vA = A + A, vA -> vB, A -> implies that the As in the vA definition map to the Bs in the vB definition
    # But both As link to the same stock and dynamic variable so just looking at those isn't enough to figure out what it maps to.
    # There will exist cases where it's impossible to tell - eg, when there exist multiple duplicate links, and some positions don't match up.
   
    # It does not currently look at the operator.  You could therefore map vA = A + A -> vB = B * B
    # I can see this being useful, actually, specifically when mapping between + and -, * and /, etc.  Probably logs and powers too.
    # Just need to be aware that it won't say it's invalid.
    @test (infer_links(
        (@stock_and_flow begin; :stocks; A; :dynamic_variables; vA = A + A; end), 
        (@stock_and_flow begin; :stocks; B; :dynamic_variables; vB = B + B; end),
        Dict{Symbol, Vector{Int64}}(:S => [1], :F => [], :SV => [], :P => [], :V => [1]))
    == Dict(:LS => [], :LSV => [], :LV => [1,1], :I => [], :O => [], :LPV => [], :LVV => [])) # If duplicate values, always map to first.

    @test (infer_links(
        (@stock_and_flow begin; :stocks; A; :parameters; pA; :dynamic_variables; vA = A + pA; end), 
        (@stock_and_flow begin; :stocks; B; :parameters; pB; :dynamic_variables; vB = pB + B; end),
        Dict{Symbol, Vector{Int64}}(:S => [1], :F => [], :SV => [], :P => [1], :V => [1]))
    == Dict(:LS => [], :LSV => [], :LV => [1], :I => [], :O => [], :LPV => [1], :LVV => []))

    @test (infer_links(
        (@stock_and_flow begin
            :stocks
            S
            I
            R

            :parameters
            p_inf
            p_rec


            :flows
            S => f_StoI(p_inf * S) => I
            I => f_ItoR(I * p_rec) => R

            :sums
            N = [S,I,R]
            NI = [I]
            NS = [S,I,R]
        end),
        (@stock_and_flow begin
            :stocks
            pop

            :parameters
            p_generic


            :flows
            pop => f_generic(p_generic * pop) => pop

            :sums
            N = [pop]
            NI = [pop]
            NS = [pop]
        end),

        Dict{Symbol, Vector{Int64}}(:S => [1,1,1], :F => [1,1], :SV => [1,2,3], :P => [1,1], :V => [1,1]))
    == Dict(:LS => [1,3,1,2,3,1,3], :LSV => [], :LV => [1,1], :I => [1,1], :O => [1,1], :LPV => [1,1], :LVV => []))

 
end


@testset "Applying flags can correctly find substring matches" begin
    @test apply_flags(:f_, Set([:~]), Vector{Symbol}()) == [] 
    @test apply_flags(:f_, Set([:~]), [:f_death, :f_birth]) == [:f_death, :f_birth] 
    @test apply_flags(:NOMATCH, Set([:~]), [:f_death, :f_birth]) == [] 
    @test apply_flags(:f_birth, Set([:~]), [:f_death, :f_birth]) == [:f_birth] 
    @test apply_flags(:f_birth, Set{Symbol}(), [:f_death, :f_birth]) == [:f_birth]

    # Note, apply_flags is specifically meant to work on vectors without duplicates; the vector which is input are the keys of a dictionary.
    # Regardless, the following will hold:
    @test apply_flags(:f_birth, Set{Symbol}(), [:f_death, :f_birth, :f_birth, :f_birth]) == [:f_birth]
    @test apply_flags(:f_birth, Set{Symbol}([:~]), [:f_death, :f_birth, :f_birth, :f_birth]) == [:f_birth, :f_birth, :f_birth]
end


@testset "substitute_symbols will correctly associate values of the two provided dictionaries based on user defined mappings" begin
    # substitute_symbols(s::Dict{Symbol, Int}, t::Dict{Symbol, Int}, m::Vector{DSLArgument} ; use_flags::Bool=true)::Dict{Int, Int}


    # Note, these dictionaries represent a vector where all the entries are unique, and the values are the original indices.
    # So, both keys and values should be unique.
    # For stratification, first dictionary is strata or aggregate, second is type, and the vector of DSLArgument are the user-defined maps.
    # For homomorphism, first argument is src, second is dest, vector are user-defined maps.
    @test substitute_symbols(Dict{Symbol, Int}(), Dict{Symbol, Int}(), Vector{DSLArgument}()) == Dict{Int, Int}()
    @test substitute_symbols(Dict{Symbol, Int}(), Dict(:B => 2), Vector{DSLArgument}()) == Dict{Int, Int}()

    @test substitute_symbols(Dict(:A => 1), Dict(:B => 1), [DSLArgument(:A, :B, Set{Symbol}())]) == Dict(1 => 1)
    @test substitute_symbols(Dict(:A1 => 1, :A2 => 2), Dict(:B => 1), [DSLArgument(:A1, :B, Set{Symbol}()), DSLArgument(:A2, :B, Set{Symbol}())]) == Dict(1 => 1, 2 => 1)
    @test substitute_symbols(Dict(:A1 => 1), Dict(:B1 => 1, :B2 => 2), [DSLArgument(:A1, :B2, Set{Symbol}())]) == Dict(1 => 2)


    @test substitute_symbols(Dict(:A1 => 1, :A2 => 2), Dict(:B1 => 1, :B2 => 2), [DSLArgument(:A, :B2, Set{Symbol}([:~]))]) == Dict(1 => 2, 2 => 2)

    # 1:100
    # 1:50
    # All multiples x of 14 below 100 go to x % 10 + 1
    @test (substitute_symbols(Dict(Symbol(i) => i for i ∈ 1:100), Dict(Symbol(-i) => i for i ∈ 1:50), [DSLArgument(Symbol(i), Symbol(-((i%10) + 1)), Set{Symbol}()) for i ∈ 1:100 if i % 14 == 0])
    == Dict(14 => 5, 28 => 9, 42 => 3, 56 => 7, 70 => 1, 84 => 5, 98 => 9))

    # Captures everything with a 7 as a digit
    @test (substitute_symbols(Dict(Symbol(i) => i for i ∈ 1:100), Dict(Symbol(-i) => i for i ∈ 1:50), [DSLArgument(Symbol(7), Symbol(-1), Set{Symbol}([:~]))])
    == Dict(7 => 1, 17 => 1, 27 => 1, 37 => 1, 47 => 1, 57 => 1, 67 => 1, 70 => 1, 71 => 1, 72 => 1, 73 => 1, 74 => 1, 75 => 1, 76 => 1, 77 => 1, 78 => 1, 79 => 1, 87 => 1, 97 => 1))

    @test substitute_symbols(Dict(Symbol("~") => 1), Dict(:R => 1), [DSLArgument(Symbol("~"), :R, Set([:~]))], ; use_flags = false) == Dict(1 => 1) # Note, the Set([:~]) is ignored because use_flags is false

end


@testset "non-natural transformations fail infer_links" begin

    # Map both dynamic variables to the same
    # Obviously, this will fail, as the new dynamic variable needs a LVV and one LV, but instead has two LV
    @test_throws KeyError (infer_links(
        (@stock_and_flow begin
        :stocks
        A

        :dynamic_variables
        v1 = A + A
        v2 = v1 + A
    end),
    (@stock_and_flow begin
        :stocks
        A
        
        :dynamic_variables
        v1 = A + A
    end), 
    Dict{Symbol, Vector{Int64}}(:S => [1], :V => [1,1])))

    

    # Mapping it all to I

    # This one fails when trying to figure out the inflow.  Stock maps to 2, and flow maps to 2,
    # But inflows on the target have (1,2) and (2,3)

    # This also wouldn't work if we tried mapping flow to 1 instead.  Outflows expect 1,1 or 2,2,
    # so it fails on (2,1).
     @test_throws  KeyError (infer_links(
        (@stock_and_flow begin
        :stocks
        pop

        :parameters
        p_generic


        :flows
        pop => f_generic(p_generic * pop) => pop

        :sums
        N = [pop]
        NI = [pop]
        NS = [pop]
    end),
    (@stock_and_flow begin
        :stocks
        S
        I
        R

        :parameters
        p_inf
        p_rec


        :flows
        S => f_StoI(p_inf * S) => I
        I => f_ItoR(I * p_rec) => R

        :sums
        N = [S,I,R]
        NI = [I]
        NS = [S,I,R]
    end),
    Dict{Symbol, Vector{Int64}}(:S => [2], :F => [2], :SV => [1,2,3], :P => [2], :V => [2])))

end

@testset "Applying flags throws on invalid inputs" begin
    @test_throws ErrorException apply_flags(:f_, Set([:+]), [:f_death, :f_birth]) # fails because :+ is not a defined operation
    @test_throws ErrorException apply_flags(:f_birth, Set([:~, :+]), [:f_death, :f_birth]) # also fails for same reason

    @test_throws AssertionError apply_flags(:NOMATCH, Set{Symbol}(), Vector{Symbol}()) # fails because it's not looking for substrings, and :NOMATCH isn't in the list of options.
    @test_throws AssertionError apply_flags(:NOMATCH, Set{Symbol}(), [:nomatch]) # same reason

end