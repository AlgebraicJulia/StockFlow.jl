using Catlab.CategoricalAlgebra

nn = StockFlow.nn
ne = StockFlow.ne

@testset "Empty CausalLoopF" begin
  empty = CausalLoopF()
  @test ( nn(empty) == 0 
    && ne(empty) == 0 
    && isempty(subpart(empty, :epolarity)) )
end

@testset "Simple CausalLoopF" begin
  cl1 = CausalLoopF([:A, :B], [:A => :B], [POL_BALANCING])
  @test nnames(cl1) == [:A, :B] && ne(cl1) == 1 && epol(cl1, 1) == POL_BALANCING
  
  cl2 = CausalLoopF([:A, :B, :C], [:A => :B, :B => :C, :C => :A], [POL_REINFORCING, POL_NOT_WELL_DEFINED, POL_ZERO])
  @test nnames(cl2) == [:A, :B, :C] && ne(cl2) == 3 && epols(cl2) == [POL_REINFORCING, POL_NOT_WELL_DEFINED, POL_ZERO]

  cl3 = CausalLoopF([:A], [:A => :A for _ in 1:5], [POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED])
  @test nnames(cl3) == [:A] && ne(cl3) == 5 && epols(cl3) == [POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED]
end

