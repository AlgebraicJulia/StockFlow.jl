using Pkg;
Pkg.activate(".")

using Test
using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax.Rewrite

@testset "sfrewrite changes stock flows as expected" begin
    @test sfrewrite(StockAndFlowF(), quote end) == StockAndFlowF()
    @test sfrewrite(StockAndFlowF(), (quote; :stocks; +A; end)) == (@stock_and_flow begin; :stocks; A; end;)
    @test sfrewrite((@stock_and_flow begin; :stocks; A; end;), (quote; :stocks; +A; end)) == (@stock_and_flow begin; :stocks; A; A; end;)
    @test sfrewrite((@stock_and_flow begin; :stocks; A; end;), (quote; :stocks; -A; end)) == StockAndFlowF()
    @test sfrewrite((@stock_and_flow begin; :stocks; A; end;), (quote; :swaps; A => B; :stocks; -A; +B; end)) == (@stock_and_flow begin; :stocks; B; end;)
    @test sfrewrite((@stock_and_flow begin; :stocks; A; end;), (quote; :swaps; A := B; :stocks; -A; +B; end)) == (@stock_and_flow begin; :stocks; B; end;)
    @test sfrewrite((@stock_and_flow begin; :stocks; A; :dynamic_variables; v1 = A + A; end;), (quote; :swaps; v1 => v2; end)) == (@stock_and_flow begin; :stocks; A; :dynamic_variables; v1 = A + A; v2 = A + A;  end;)
    @test sfrewrite((@stock_and_flow begin; :stocks; A; :dynamic_variables; v1 = A + A; end;), (quote; :swaps; v1 := A - A; end)) == (@stock_and_flow begin; :stocks; A; :dynamic_variables; v1 = A - A; end;)

    @test sfrewrite(
    (@stock_and_flow begin
        :stocks
        A

        :dynamic_variables
        v1 = A + A
    end),
    (quote
        :swaps
        v1 => v2

        :dynamic_variables
        +v2 = A - A
        -v1
    end)  
    ) ==
    @stock_and_flow begin
        :stocks
        A

        :dynamic_variables
        v2 = A - A
    end

end