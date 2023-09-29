using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax: NothingFunction
using StockFlow.Syntax.Homomorphism

using Catlab.CategoricalAlgebra

@testset "hom macro creates correct homomorphisms" begin
  empty = @stock_and_flow begin end
  empty_hom = ACSetTransformation(empty, empty; :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction )
  @test (@hom empty empty begin end) == empty_hom

  sfA = @stock_and_flow begin; :stocks; A; end;
  sfB = @stock_and_flow begin; :stocks; B; end;
  @test (@hom sfA sfB begin; :stocks; A => B; end;) == ACSetTransformation(sfA, sfB ; :S => [1], :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)

end