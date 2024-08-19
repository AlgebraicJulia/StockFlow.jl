using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax.Composition
import StockFlow.Syntax.Composition: interpret_composition_notation
using Catlab.CategoricalAlgebra

function ≅(x,y)
  !isnothing(isomorphism(x,y))
end

@testset "Composition creates expected stock flows" begin
    empty_sf = StockAndFlowF()
    A = @stock_and_flow begin; :stocks; A; end;
    AA = @stock_and_flow begin; :stocks; A; A; end;
    
    B = @stock_and_flow begin; :stocks; B; end;
    AB = @stock_and_flow begin :stocks; A; B; end; 
    BA = @stock_and_flow begin :stocks; B; A; end; 

    @test (@compose empty_sf begin
        (sf,)
    end) == empty_sf

    @test (@compose A begin
        (sf,)
    end) == A

    @test (@compose A B (begin
        (sf1, sf2)
    end)) == AB # Combining without any composing

    @test (@compose A A  (begin
        (sf1, sf2)
    end)) == AA

    @test (@compose A A (begin
        (sf1, sf2)
        (sf1, sf2) ^ A => ()
    end)) == A

    @test (@compose A B (begin
    end)) == AB

    @test (@compose B B A (begin
        (B1, B2)
        (B1, B2) ^ B => ()
        end)) == BA

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
            (sfA, sfC) ^ B => N
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
    @test interpret_composition_notation(:(sf ^ A => N), create_foot) == ([:sf], (@foot A => N))
    @test interpret_composition_notation(:((sf1, sf2) ^ A => N), create_foot) == ([:sf1,:sf2], (@foot A => N))
    @test interpret_composition_notation(:((sf1, sf2) ^ A => N, A => NI), create_foot) == ([:sf1,:sf2], (@foot A => N, A => NI))
    @test interpret_composition_notation(:((sf1, sf2, sf3, sf4) ^ () => NI), create_foot) == ([:sf1, :sf2, :sf3, :sf4], (@foot () => NI))
    @test interpret_composition_notation(:((sf1, sf2) ^ L => ()), create_foot) == ([:sf1,:sf2], (@foot L => ()))

    @test interpret_composition_notation(:((sf1, sf2) ^ () => ()), create_foot) == ([:sf1,:sf2], (@foot () => ()))


    @test interpret_composition_notation(:((sf1, sf2, sf3) ^ () => (), A => N), create_foot) == ([:sf1,:sf2,:sf3], (@foot () => (), A => N))

end

@testset "invalid composition expressions fail" begin
    @test_throws ErrorException interpret_composition_notation(:(B => C), create_foot)
    @test_throws ErrorException interpret_composition_notation(:(A, B, C), create_foot)
    @test_throws ErrorException interpret_composition_notation(:(A ^ B ^ C), create_foot)
    @test_throws ErrorException interpret_composition_notation(:(A => B => C), create_foot)
    @test_throws ErrorException interpret_composition_notation(:(A ^ B => C => D), create_foot)
end

@testset "invalid sfcompose calls fail" begin
    @test_throws AssertionError sfcompose([(@stock_and_flow begin; :stocks; A; end;), (@stock_and_flow begin; :stocks; A; end;)], (quote
        (sf1, sf2)
        (sf1, sf2) ^ () => ()
    end), StockAndFlowF, StockAndFlow0, create_foot) # not allowed to map to empty
    @test_throws AssertionError sfcompose([(@stock_and_flow begin; :stocks; A; end;), (@stock_and_flow begin; :stocks; A; end;)], quote
        (sf1, sf2)
        sf1 ^ A => ()
        sf2 ^ A => ()
    end, StockAndFlowF, StockAndFlow0, create_foot) # not allowed to map to the same foot twice
end

@testset "Causal Loop composition" begin
    CL1 = @cl B => +A, B => +C, C => +A, G;
    CL2 = @cl B => +A, A => +D, D => -E;       
    CL3 = @cl B => +C, C => -F, F => +E, E => + B, G;

    BigCL = @causal_loop begin
        :nodes
        A; B; C; D; E; F; G;

        :edges
        A => +D;
        B => +A; B => +C;
        C => +A; C => -F;
        D => -E;
        E => +B;
        F => +E;
    end

    @test (@compose CL1 CL2 CL3 begin
        (CL1, CL2, CL3)
        (CL1, CL2) ^ B => +A
        (CL1, CL2, CL3) ^ B
        (CL2, CL3) ^ E
        (CL1, CL3) ^ B => +C, G
      end) ≅ BigCL

    AB = @cl A => +B;
    BC = @cl B => -C;
    CA = @cl C => +A;

    ABC = @cl A => +B, B => -C, C => +A;

    @test (@compose AB BC CA begin
        (AB, BC, CA)
        (AB, BC) ^ B
        (BC, CA) ^ C
        (CA, AB) ^ A
    end) ≅ ABC

    # no polarities
    ABC = (@cl A => B, B => C)
    BCD = (@cl B => C, C => D)
    ABCD = (@cl A => B, B => C, C => D)
    
    @test (@compose ABC BCD begin
        (bc, cd)
        (bc, cd) ^ B => C
    end) == ABCD

    @test (@compose ABC ABC begin
        (abc1, abc2)
        (abc1, abc2) ^ A => B, B => C
    end) == ABC
    
      

end