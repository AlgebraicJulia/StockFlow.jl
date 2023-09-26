using Test
using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax.Composition
import StockFlow.Syntax.Composition: interpret_composition_notation

@testset "Composition creates expected stock flows" begin
    empty_sf = StockAndFlowF()


    @test (@compose (begin # composing no stock flows returns an empty stock flow.
        ()
    end)) == empty_sf

    @test (@compose empty_sf begin
        (sf,)
    end) == empty_sf

    @test (@compose (@stock_and_flow begin; :stocks; A; end;) (@stock_and_flow begin; :stocks; B; end;) (begin
        (sf1, sf2)
    end)) == (@stock_and_flow begin; :stocks; A; B; end;) # Combining without any composing

    @test (@compose (@stock_and_flow begin; :stocks; A; end;) (@stock_and_flow begin; :stocks; A; end;) (begin
        (sf1, sf2)
    end)) == (@stock_and_flow begin; :stocks; A; A; end;)

    @test (@compose (@stock_and_flow begin; :stocks; A; end;) (@stock_and_flow begin; :stocks; A; end;) (begin
        (sf1, sf2)
        sf1, sf2 ^ A => ()
    end)) == (@stock_and_flow begin; :stocks; A; end;)

    @test ((@compose (@stock_and_flow begin
        :stocks
        A
        B

        :dynamic_variables
        v1 = A + B

        :sums
        N = [A,B]
    end) (@stock_and_flow begin
            :stocks
            B
            C

            :dynamic_variables
            v2 = B + C

            :sums
            N = [B,C]
        end) (begin
            (sfA, sfC)
            sfA, sfC ^ B => N
        end))
        ==
        (@stock_and_flow begin
            :stocks
            A
            B
            C

            :dynamic_variables
            v1 = A + B
            v2 = B + C

            :sums
            N = [A, B, C]
        end))
        


end

@testset "interpret_composition_notation interprets arguments correctly" begin
    # @test interpret_composition_notation(:(() ^ A => N)) == (Vector{Symbol}(), (@foot A => N))
    @test interpret_composition_notation(:(sf ^ A => N)) == ([:sf], (@foot A => N))
    @test interpret_composition_notation(:(sf1, sf2 ^ A => N)) == ([:sf1,:sf2], (@foot A => N))
    @test interpret_composition_notation(:(sf1, sf2 ^ A => N, A => NI)) == ([:sf1,:sf2], (@foot A => N, A => NI))
    @test interpret_composition_notation(:(sf1, sf2, sf3, sf4 ^ () => NI)) == ([:sf1, :sf2, :sf3, :sf4], (@foot () => NI))
    @test interpret_composition_notation(:(sf1, sf2 ^ L => ())) == ([:sf1,:sf2], (@foot L => ()))

    @test interpret_composition_notation(:(sf1, sf2 ^ () => ())) == ([:sf1,:sf2], (@foot () => ()))

end

@testset "invalid composition expressions fail" begin
    @test_throws AssertionError interpret_composition_notation(:(B => C))
    @test_throws AssertionError interpret_composition_notation(:(A, B, C))
    @test_throws AssertionError interpret_composition_notation(:(A ^ B ^ C))
    @test_throws AssertionError interpret_composition_notation(:(A => B => C))
    @test_throws ErrorException interpret_composition_notation(:(A ^ B => C => D)) # caught by create_foot
end

@testset "invalid sfcompose calls fail" begin
    @test_throws AssertionError sfcompose([(@stock_and_flow begin; :stocks; A; end;), (@stock_and_flow begin; :stocks; A; end;)], quote
        (sf1, sf2)
        sf1, sf2 ^ () => ()
    end) # not allowed to map to empty
    @test_throws AssertionError sfcompose([(@stock_and_flow begin; :stocks; A; end;), (@stock_and_flow begin; :stocks; A; end;)], quote
        (sf1, sf2)
        sf1 ^ A => ()
        sf2 ^ A => ()
    end) # not allowed to map to the same foot twice
    @test_throws AssertionError sfcompose([(@stock_and_flow begin; :stocks; A; end;), (@stock_and_flow begin; :stocks; A; end;)], quote
        (sf1,)
        sf1 ^ A => ()
    end) # incorrect number of symbols on the first line in the quote
end
