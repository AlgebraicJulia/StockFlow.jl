using Pkg;
Pkg.activate(".")

using Test

using StockFlow
using StockFlow.Syntax
using Catlab.CategoricalAlgebra


@testset "Empty CausalLoopPM" begin
    empty = CausalLoopF()
    @test (nvert(empty) == 0
           && np(empty) == 0
           && nm(empty) == 0)

    @test from_clp(to_clp(empty)) == empty

end

@testset "Basic CausalLoop Creation" begin
    cl1 = (@cl A => +B, B => -C, D, C => -A)
    cl2 = (@causal_loop begin
               :nodes #TODO: rename to vertices
               A
               B
               C
               D

               :edges
               A => +B
               B => -C
               C => -A
           end)
    cl3 = CausalLoopF([:A,:B,:C,:D], [:A => :B, :B => :C, :C => :A], [POL_POSITIVE, POL_NEGATIVE, POL_NEGATIVE])
    
    @test cl1 == cl2
    @test cl2 == cl3
    @test cl3 == cl1

    @test from_clp(to_clp(cl1)) == cl1

end


@testset "Cycles" begin
    clc = (@cl A => +B, A => -B, B => -B, B => -A)
    
    @test cl_cycles(CausalLoopF()) == Vector{Vector{Int}}()
    @test cl_cycles(clc) == [[1, 4], [2, 4], [3]]

    cl = (@cl A => +B, C => -D, D => -C, E => +E)
    
    # Unfortunately, this is pretty non-indicative of what they actually are.
    # Edges are ordered first by pos > neg > zero > nwd > unknown, then by order in which they appear.
    # So, the order of edges here is A => +B == 1, E => +E == 2, C => -D == 3,
    # D => -C == 4
    @test cl_cycles(cl) == [[3,4], [2]] 

    @test extract_loops(CausalLoopF()) == Vector{Pair{Vector{Int}, Polarity}}()
    @test extract_loops(clc) == [[1,4] => POL_NEGATIVE, [2, 4] => POL_POSITIVE, [3] => POL_NEGATIVE]

end

@testset "Walk polarity" begin

    @test is_walk(CausalLoopF(), Vector{Int}())
    @test !is_walk((@cl A => +B, B => +C), [2, 1])
    @test !is_walk(CausalLoopF(), [1])

    @test is_circuit((@cl A => +B, B => -A), [1,2])
    @test !is_circuit((@cl A => +B, B => -A), [1])
    @test !is_circuit((@cl A => +B), Vector{Int}())
    @test !is_circuit((@cl A => +B), [2])

    @test walk_polarity(CausalLoopF(), Vector{Int}()) == POL_POSITIVE
    @test walk_polarity((@cl A => -B), [1]) == POL_NEGATIVE 
    @test walk_polarity((@cl A => -B, B => -A), [1,2]) == POL_POSITIVE
    @test walk_polarity((@cl A => -B, B => -A), [1,2,1]) == POL_NEGATIVE
end


@testset "Count of loops a variable is on" begin
    cll = (@cl A => +B, B => +C, C => -D, D => +A, D => -A, E => +E, E => -E, F => +F, G)

    @test num_loops_var_on(cll, :A) == 2
    @test num_indep_loops_var_on(cll, :A) == 1
    @test num_loops_var_on(cll, :G) == num_indep_loops_var_on(cll, :G) == 0
    @test num_loops_var_on(cll, :E) == 2
    @test num_indep_loops_var_on(cll, :E) == 1

    @test_throws ArgumentError num_loops_var_on(cll, :H)



end


@testset "all paths" begin

    # Won't hit same node twice
    @test extract_all_nonduplicate_paths( (@cl )) == Dict([Vector{Int}() => POL_POSITIVE])
    @test extract_all_nonduplicate_paths((@cl A => +B)) == Dict([Vector{Int}() => POL_POSITIVE,
                                                                 [1] => POL_POSITIVE])
    @test extract_all_nonduplicate_paths((@cl A => +B, A => -C, A => +C, B => -C)) == Dict([
                                                                                       Vector{Int}() => POL_POSITIVE, [1] => POL_POSITIVE, [1, 4] => POL_NEGATIVE, [2] => POL_POSITIVE, [3] => POL_NEGATIVE 
                                                                                      ])

end

#ne = StockFlow.ne
#
#@testset "Empty CausalLoopF" begin
  #empty = CausalLoopF()
  #@test ( nvert(empty) == 0 
    #&& ne(empty) == 0 
    #&& isempty(subpart(empty, :epolarity)) )
#end
#
#@testset "Helper functions" begin
  #cl0 = CausalLoopF([:A, :B], [:A => :B, :B => :A], [POL_NEGATIVE, POL_NEGATIVE])
  #@test epol(cl0, 1) == POL_NEGATIVE 
  #@test epols(cl0)[1] == POL_NEGATIVE
  #@test outgoing_edges(cl0, 1) == [1] && outgoing_edges(cl0, 2) == [2]
  #@test incoming_edges(cl0, 1) == [2] && incoming_edges(cl0, 2) == [1]
#end
#
#@testset "Simple CausalLoopF" begin
  #cl1 = CausalLoopF([:A, :B], [:A => :B], [POL_NEGATIVE])
  #@test nnames(cl1) == [:A, :B] && ne(cl1) == 1 && epol(cl1, 1) == POL_NEGATIVE
  #
  #cl2 = CausalLoopF([:A, :B, :C], [:A => :B, :B => :C, :C => :A], [POL_POSITIVE, POL_NOT_WELL_DEFINED, POL_ZERO])
  #@test nnames(cl2) == [:A, :B, :C] && ne(cl2) == 3 && epols(cl2) == [POL_POSITIVE, POL_NOT_WELL_DEFINED, POL_ZERO]
#
  #cl3 = CausalLoopF([:A], [:A => :A for _ in 1:5], [POL_POSITIVE, POL_NEGATIVE, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED])
  #@test nnames(cl3) == [:A] && ne(cl3) == 5 && epols(cl3) == [POL_POSITIVE, POL_NEGATIVE, POL_UNKNOWN, POL_ZERO, POL_NOT_WELL_DEFINED]
#end
#
#@testset "Extract Loops" begin
  #cl4 = CausalLoopF([:A], [:A => :A], [POL_UNKNOWN])
  #@test extract_loops(cl4) == [[1] => POL_UNKNOWN]
#
  ## 2 balancing make a reinforcing
  #cl5 = CausalLoopF([:A, :B], [:A => :A, :A => :B, :B => :A], [POL_NEGATIVE, POL_NEGATIVE, POL_NEGATIVE])
  #@test extract_loops(cl5) == [[1] => POL_NEGATIVE, [2, 3] => POL_POSITIVE]
#
  ## A -> B -> C -> D -> E
  ## ^ - - - - |    |    |
  ##^ - - - - - - - -    |
  ##^ - - - - - - - - - - 
#
  ## B -> C is not well defined, C -> D is unknown, D -> E is zero.
  #cl6 = CausalLoopF([:A, :B, :C, :D, :E], 
    #[:A => :B, :B => :C, :C => :A, :C => :D, :D => :A, :D => :E, :E => :A], 
    #[POL_NEGATIVE, POL_NOT_WELL_DEFINED, POL_NEGATIVE, POL_UNKNOWN, POL_NEGATIVE, POL_ZERO, POL_NEGATIVE])
  #
  ## shows that not well defined < unknown < zero
  ## run Graph_RB(cl6) to see the cycles
  #@test extract_loops(cl6) == [[1,2,3] => POL_NOT_WELL_DEFINED, [1,2,4,5] => POL_UNKNOWN, [1,2,4,6,7] => POL_ZERO]
#end
#
#@testset "Discard Zero Pol" begin
  #cl_disc_zero = CausalLoopF([:A], [:A => :A for _ in 1:10], [POL_ZERO for _ in 1:10])
  #cl_disc_zero′ = discard_zero_pol(cl_disc_zero)
  #@test nnames(cl_disc_zero′) == [:A] && ne(cl_disc_zero′) == 0 && epols(cl_disc_zero′) == []
#end
#
