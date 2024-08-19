# # Causal Loop Operations

using StockFlow
using StockFlow.Syntax
using StockFlow.PremadeModels

GraphCL(convertToCausalLoop(seir()))

C2 = CausalLoopPM([:C2], [:C2 => :C2 for _ in 1:5], [POL_NEGATIVE, POL_POSITIVE, POL_NEGATIVE, POL_POSITIVE, POL_NEGATIVE])

C2′ = @causal_loop begin
    :nodes
    C2
    :edges
    C2 => -C2
    C2 => +C2
    C2 => -C2
    C2 => +C2
    C2 => -C2
end;

C2′′ = (@cl C2 => -C2, C2 => +C2, C2 => -C2, C2 => +C2, C2 => -C2);

C2 == C2′ == C2′′

GraphCL(C2)

GraphRB(C2)

CL_ABC = @cl A => +B, B => +C, C => -A, D

GraphCL(CL_ABC)

GraphRB(CL_ABC)

nvert(CL_ABC)

nedges(CL_ABC)

np(CL_ABC)

nm(CL_ABC)

sedge(CL_ABC, 3) # Source node for edge 3, which is C => A, is C, which is index 3

tedge(CL_ABC, 3) # C => A, target is A, with index 1

vnames(CL_ABC)

epol(CL_ABC, 1) # Polarity of edge 1, A => B, is positive

epols(CL_ABC)

outgoing_edges(CL_ABC, 1) # indices of all edges with src 1 (in this case, A)

outgoing_edges(CL_ABC, 4) # D

incoming_edges(CL_ABC, 1) # indices of all edges with src 1

cl_cycles(CL_ABC)

extract_loops(CL_ABC)

is_walk(CL_ABC, [1,2,3,1,2])

is_walk(CL_ABC, [3,2])

is_walk(CL_ABC, Vector{Int}())

is_circuit(CL_ABC, [1,2,3])

is_circuit(CL_ABC, [1,2])

walk_polarity(CL_ABC, [1,2,3,1,2,3])

extract_all_nonduplicate_paths(CL_ABC)

num_loops_var_on(CL_ABC, :D)

num_loops_var_on(CL_ABC, :A)

num_loops_var_on(C2, :C2)

num_indep_loops_var_on(C2, :C2) # Treating each pair of nodes as if there is at most one edge between them

to_graphs_graph(CL_ABC)

betweenness(CL_ABC)

to_graphs_graph(C2) # eliminates duplicate edges!

betweenness(C2)

num_inputs_outputs(C2) # in, out

num_inputs_outputs_pols(C2) # pos in, pos out, neg in, neg out

all_shortest_paths(C2) 

 all_shortest_paths(@cl((A => B, B => C, C => D, A => B′, B′ => C′, C′ => D)))

shortest_paths( (@cl A => B, B => C, C => D, A => B′, B′ => C′, C′ => D), :A, :D )

 betweenness(@cl((A => B, B => C, C => D, A => B′, B′ => C′, C′ => D)))

all_shortest_paths(convertToCausalLoop(seir()))

all_shortest_paths(CL_ABC)

GraphCL(CL_ABC)

shortest_paths(CL_ABC, :A, :C)

cl = @causal_loop begin
        :nodes
        A
        B
        C
        D
        E
        :edges
        A => B
        B => C
        B => C
        B => D
        D => C
    end;

cl2 = @causal_loop begin
        :nodes
        A
        B
        C
        D
        E
        :edges
        A => B
        B => C
        B => D
        D => C
    end;

GraphCL(cl)

betweenness(cl)

betweenness(cl2)

# ### Note, negative polarities always come after positive!

cl_small = @cl A => -B, B => +C 

epols(cl_small)

to_simple_cl(cl_small) == (@cl A, B => C, A => B) # Note the order!

to_simple_cl(cl_small)

using StockFlow.Syntax.Composition

ABC = (@cl A => B, B => C)
BCD = (@cl B => C, C => D)
@compose ABC BCD begin
    (ABC, BCD)
    (ABC, BCD) ^ B => C
end

ABC_pol = (@cl A => +B, B => -C)
BCD_pol = (@cl B => -C, C => +D)
@compose ABC_pol BCD_pol begin
    (ABC, BCD)
    (ABC, BCD) ^ B => -C
end

