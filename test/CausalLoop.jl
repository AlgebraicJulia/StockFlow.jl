using StockFlow
using StockFlow.Syntax
using StockFlow.PremadeModels

using Catlab.CategoricalAlgebra    


@testset "Empty CausalLoop" begin
    
    e = CausalLoop() # Graph with names
    @test (nvert(e) == 0
           && nedges(e) == 0)

    e2 = CausalLoopPol()
    @test (nvert(e2) == 0
           && nedges(e2) == 0
           && epols(e2) == Vector{Polarity}())


    empty = CausalLoopPM() # CausalLoopPM
    @test (nvert(empty) == 0
           && np(empty) == 0
           && nm(empty) == 0)

    @test from_clp(to_clp(empty)) == empty

end



@testset "sf to causal loop" begin
    @test (convertToCausalLoop(seir()) 
           == (@causal_loop begin
            :nodes
            S; E; I; R; f_birth; f_incid; f_deathS; f_inf; f_deathE; f_rec; f_deathI; f_deathR; N; NI; NS; v_prevalence; v_meanInfectiousContactsPerS; v_perSIncidenceRate; μ; β; tlatent; trecovery; δ; c;
            :edges
             S => N
             S => NS
             E => N
             E => NS
             I => N
             I => NI
             I => NS
             R => N
             R => NS
             
             NI => v_prevalence
             NS => v_prevalence
             N => f_birth
             
             S => f_incid
             E => f_inf
             I => f_rec
             
             S => f_deathS
             E => f_deathE
             I => f_deathI
             R => f_deathR
             f_birth => S
             f_incid => E
             f_inf => I
             f_rec => R
             f_incid => S 
             f_deathS => S
             f_inf => E
             f_deathE => E
             f_rec => I
             f_deathI => I
             f_deathR => R
 
             c => v_meanInfectiousContactsPerS
             β => v_perSIncidenceRate
             μ => f_birth
             tlatent => f_inf
             trecovery => f_rec
             δ => f_deathS
             δ => f_deathE
             δ => f_deathI
             δ => f_deathR
             
             v_prevalence => v_meanInfectiousContactsPerS
             v_meanInfectiousContactsPerS => v_perSIncidenceRate
             v_perSIncidenceRate => f_incid
 
         
            end)

          )
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
    cl3 = CausalLoopPM([:A,:B,:C,:D], [:A => :B, :B => :C, :C => :A], [POL_POSITIVE, POL_NEGATIVE, POL_NEGATIVE])
    
    @test cl1 == cl2
    @test cl2 == cl3
    @test cl3 == cl1

    @test from_clp(to_clp(cl1)) == cl1

end


@testset "Cycles" begin
    clc = (@cl A => +B, A => -B, B => -B, B => -A)
    
    @test cl_cycles(CausalLoopPM()) == Vector{Vector{Int}}()
    @test cl_cycles(clc) == [[1, 4], [2, 4], [3]]

    cl = (@cl A => +B, C => -D, D => -C, E => +E)
    
    # Unfortunately, this is pretty non-indicative of what they actually are.
    # Edges are ordered first by pos > neg > zero > nwd > unknown, then by order in which they appear.
    # So, the order of edges here is A => +B == 1, E => +E == 2, C => -D == 3,
    # D => -C == 4
    @test cl_cycles(cl) == [[3,4], [2]] 

    @test extract_loops(CausalLoopPM()) == Dict{Vector{Int}, Polarity}()
    @test extract_loops(clc) == Dict([[1,4] => POL_NEGATIVE, [2, 4] => POL_POSITIVE, [3] => POL_NEGATIVE])

end

@testset "Walk polarity" begin
    @test is_walk(CausalLoopPM(), Vector{Int}())
    @test !is_walk((@cl A => +B, B => +C), [2, 1])
    @test !is_walk(CausalLoopPM(), [1])

    @test is_circuit((@cl A => +B, B => -A), [1,2])
    @test !is_circuit((@cl A => +B, B => -A), [1])
    @test !is_circuit((@cl A => +B), Vector{Int}())
    @test !is_circuit((@cl A => +B), [2])

    @test walk_polarity(CausalLoopPM(), Vector{Int}()) == POL_POSITIVE
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
    @test (extract_all_nonduplicate_paths((@cl A => +B, A => -C, A => +C, B => -C))
      == Dict([Vector{Int}() => POL_POSITIVE, [1] => POL_POSITIVE, [1, 4] => POL_NEGATIVE, [2] => POL_POSITIVE, [3] => POL_NEGATIVE, [4] => POL_NEGATIVE ]))

end

@testset "Number of Outputs" begin
    @test num_inputs_outputs(@cl ) == Dict{Symbol, Tuple{Int, Int}}()
    @test num_inputs_outputs(@cl A) == Dict([:A => (0, 0)])
    @test num_inputs_outputs(@cl A => +B, B => +C, C => -D) == Dict([:A => (0, 1), :B => (1, 1), :C => (1, 1), :D => (1, 0)])

    @test (num_inputs_outputs_pols(@cl A => +B, A => -C, A => +D, A => -E, A => -B) == 
    Dict([:A => (0, 2, 0, 3), :B => (1, 0, 1, 0),
     :C => (0, 0, 1, 0), :D => (1, 0, 0, 0), :E => (0, 0, 1, 0)    
    ]))

end

@testset "Shortest path functions" begin

    @testset "Shortest paths" begin
        @test shortest_paths((@cl A), :A, :A) == Vector{Int}[[]] # trivial path
        @test shortest_paths((@cl A => A), :A, :A) == Vector{Int}[[]]
        @test shortest_paths((@cl A => B), :A, :A) == Vector{Int}[[]]
        @test shortest_paths((@cl A => B), :A, :B) == Vector{Int}[[1]]
        @test shortest_paths((@cl A => B), :B, :A) == Vector{Int}[] # no path

        @test shortest_paths((@cl A => B, B => C), :A, :C) == Vector{Int}[[1,2]]
        @test shortest_paths((@cl A => B, B => C, A => C), :A, :C) == Vector{Int}[[3]]

        @test shortest_paths((@cl A => B, B => C, C => D, A => B′, B′ => C′, C′ => D), :A, :D) == Vector{Int}[[1, 2, 3], [4,5,6]]
        @test shortest_paths((@cl A => B, A => B), :A, :B) == Vector{Int}[[1], [2]]

        @test shortest_paths(convertToCausalLoop(seir()), :S, :E) ==  Vector{Int}[[13, 21]] #  S => f_incid => E

    end

    @testset "All shortest paths" begin

        vvi = Vector{Vector{Int}}()
        vve = Vector{Vector{Int}}([[]])

        @test (all_shortest_paths(@cl) == Matrix{Vector{Vector{Int}}}(undef, 0, 0))
        @test (all_shortest_paths(@cl A => A) == fill(vve, (1,1)))
        @test (all_shortest_paths(@cl A => B, B => C)
        == permutedims(stack(
            [[vve, [[1]], [[1,2]]],
             [vvi, vve, [[2]]],
             [vvi,  vvi,    vve]])
        ))

        @test (all_shortest_paths(@cl A => B, B => C, C => A, D) 
        == permutedims(stack(
            [
                [vve, [[1]], [[1,2]], vvi],
                [[[2,3]], vve, [[2]], vvi ],
                [[[3]], [[3,1]], vve, vvi],
                [vvi, vvi, vvi, vve]
            ]

        ))
        )


        @test (all_shortest_paths(@cl A => B, A => C, B => A, B => C, C => A, C => B)
        == permutedims(stack(
            [
                [vve, [[1]], [[2]]],
                [[[3]], vve, [[4]]],
                [[[5]], [[6]], vve]
            ]
 

        )))



        multi = Matrix{Vector{Vector{Int}}}(undef, 2, 2)
        multi[1,1] = vve
        multi[1,2] = [[1],[2]]
        multi[2,1] = [[4]]
        multi[2,2] = vve


        @test (all_shortest_paths(@cl A => B, A => B, B => B, B => A)
        == multi
        )

    end


    @testset "Betweenness" begin
        
        @test betweenness(@cl) == Vector{Float64}()
        @test betweenness(@cl A => B) == [0.0, 0.0]
        @test betweenness(@cl A => B, B => C) == [0.0, 1.0, 0.0]


    end


end