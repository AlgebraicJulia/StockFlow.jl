using Catlab.CategoricalAlgebra

nn = StockFlow.nn
ne = StockFlow.ne

@testset "Empty CausalLoopF" begin
  empty = CausalLoopF()
  @test ( nn(empty) == 0 
    && ne(empty) == 0 
    && isempty(subpart(empty, :epolarity)) )
end

@testset "Helper functions" begin
  cl0 = CausalLoopF([:A, :B], [:A => :B, :B => :A], [POL_BALANCING, POL_BALANCING])
  @test epol(cl0, 1) == POL_BALANCING 
  @test epols(cl0)[1] == POL_BALANCING
  @test outgoing_edges(cl0, 1) == [1] && outgoing_edges(cl0, 2) == [2]
  @test incoming_edges(cl0, 1) == [2] && incoming_edges(cl0, 2) == [1]
end

@testset "Simple CausalLoopF" begin
  cl1 = CausalLoopF([:A, :B], [:A => :B], [POL_BALANCING])
  @test nnames(cl1) == [:A, :B] && ne(cl1) == 1 && epol(cl1, 1) == POL_BALANCING
  
  cl2 = CausalLoopF([:A, :B, :C], [:A => :B, :B => :C, :C => :A], [POL_REINFORCING, POL_NOT_WELL_DEFINED, POL_ZERO])
  @test nnames(cl2) == [:A, :B, :C] && ne(cl2) == 3 && epols(cl2) == [POL_REINFORCING, POL_NOT_WELL_DEFINED, POL_ZERO]

  cl3 = CausalLoopF([:A], [:A => :A for _ in 1:5], [POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED])
  @test nnames(cl3) == [:A] && ne(cl3) == 5 && epols(cl3) == [POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED]
end

@testset "Extract Loops" begin
  cl4 = CausalLoopF([:A], [:A => :A], [POL_UNKNOWN])
  @test extract_loops(cl4) == Dict([1] => POL_UNKNOWN)

  # 2 balancing make a reinforcing
  cl5 = CausalLoopF([:A, :B], [:A => :A, :A => :B, :B => :A], [POL_BALANCING, POL_BALANCING, POL_BALANCING])
  @test extract_loops(cl5) == Dict([1] => POL_BALANCING, [2, 3] => POL_REINFORCING)

  # A -> B -> C -> D -> E
  # ^ - - - - |    |    |
  #^ - - - - - - - -    |
  #^ - - - - - - - - - - 

  # B -> C is not well defined, C -> D is unknown, D -> E is zero.
  cl6 = CausalLoopF([:A, :B, :C, :D, :E], 
    [:A => :B, :B => :C, :C => :A, :C => :D, :D => :A, :D => :E, :E => :A], 
    [POL_BALANCING, POL_NOT_WELL_DEFINED, POL_BALANCING, POL_UNKNOWN, POL_BALANCING, POL_ZERO, POL_BALANCING])
  
  # shows that not well defined < unknown < zero
  # run Graph_RB(cl6) to see the cycles
  @test extract_loops(cl6) == Dict([1,2,3] => POL_NOT_WELL_DEFINED, [1,2,4,5] => POL_UNKNOWN, [1,2,4,6,7] => POL_ZERO)
end

@testset "Discard Zero Pol" begin
  cl_disc_zero = CausalLoopF([:A], [:A => :A for _ in 1:10], [POL_ZERO for _ in 1:10])
  cl_disc_zero′ = discard_zero_pol(cl_disc_zero)
  @test nnames(cl_disc_zero′) == [:A] && ne(cl_disc_zero′) == 0 && epols(cl_disc_zero′) == []
end

