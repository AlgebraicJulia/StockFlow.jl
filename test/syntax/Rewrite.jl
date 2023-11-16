using Test

using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax.Rewrite

@testset "Trivial rewrite examples act as expected." begin
  empty = StockAndFlowF()
  A = @stock_and_flow begin; :stocks; A; end
  p = @stock_and_flow begin; :parameters; p; end

  Ap = @stock_and_flow begin; :stocks; A; :parameters; p; end

  @test (@rewrite empty begin end) == empty
  @test (@rewrite A begin end) == A

  @test (@rewrite empty begin
    :stocks
    +A    
    end) == A

  @test (@rewrite A begin
    :stocks
    -A
    end) == empty

  @test (@rewrite A begin
    :stocks
    -B
    end) == A

  @test (@rewrite p begin
    :parameters
    -p
  end) == empty

  @test (@rewrite Ap begin
    :parameters
    -p
  end) == A

  @test (@rewrite Ap begin
    :stocks
    -A
  end) == p

  @test (@rewrite p begin
    :stocks
    +A
  end) == Ap

  @test (@rewrite A begin
    :parameters
    +p
  end) == Ap

    # @test (@rewrite A begin
    #   :stocks
    #   -B
    #   end) == A

    # delete A, but shouldn't
    # We might want to just add a :removes section and if you put it under
    # a :parameters or :stocks header it's just added.

    # @test (@rewrite A begin 
    #   :parameters
    #   -A
    #   end) == A
    

end


@testset "Basic dynamic variable manipulation" begin
  # You should never swap both variables in a dynamic variable in one rewrite.


  ABv = @stock_and_flow begin
    :stocks
    A
    B
    :dynamic_variables
    v = A + B
  end

  AAv = @stock_and_flow begin
    :stocks
    A
    :dynamic_variables
    v = A + A
  end

  BBv = @stock_and_flow begin
    :stocks
    B
    :dynamic_variables
    v = B + B
  end

  Av = @stock_and_flow begin
    :stocks
    A
    :dynamic_variables
    v = +(A)
  end





  @test (@rewrite ABv begin
    :swaps
    B => A
    :stocks
    -B
  end) == AAv

  @test (@rewrite ABv begin
    :stocks
    -B
    :swaps
    B => A
  end) == AAv

  @test (@rewrite ABv begin
    :redefs
    v := A + A
    :stocks
    -B
  end) == AAv

  @test (@rewrite ABv begin
    :redefs
    v := +(A)
    :stocks
    -B
  end) == Av

  @test (@rewrite ABv begin
    :stocks
    -B
    :redefs
    v := +(A)
  end) == Av

  @test (@rewrite ABv begin
    :stocks
    -B

  end) == Av


  @test (@rewrite Av begin
    :redefs
    v := B + B
    :stocks
    +B
    -A
  end) == BBv


end

