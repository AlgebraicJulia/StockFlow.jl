export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop,
nvert, nedges, vname, np, nm,
sedge, tedge, convertToCausalLoop, vnames, epol, epols,
Polarity, POL_POSITIVE, POL_NEGATIVE,
add_node!, add_nodes!, add_edge!, add_edges!,
outgoing_edges, incoming_edges, extract_loops, is_walk, is_circuit, walk_polarity, cl_cycles,
CausalLoopPol, to_clp, from_clp, CausalLoopPM, leg,
extract_all_nonduplicate_paths, num_loops_var_on, num_indep_loops_var_on,
betweenness, to_simple_cl, num_inputs_outputs, num_inputs_outputs_pols,
to_graphs_graph, shortest_paths, all_shortest_paths


using MLStyle

import Graphs
import Graphs: SimpleDiGraph, simplecycles, SimpleEdge, betweenness_centrality, a_star


import Base: *



"""
P - sp - >  
  - tp - >  
            V
M - sm - >
  - tm - >
"""
@present TheoryCausalLoopNameless(FreeSchema) begin
  
  V::Ob
  P::Ob
  M::Ob

  sp::Hom(P, V)
  tp::Hom(P, V)

  sm::Hom(M, V)
  tm::Hom(M, V)

end


"""
P - sp - >    
  - tp - >    
            V - vname - > Name
M - sm - >
  - tm - >
"""
@present TheoryCausalLoopPM <: TheoryCausalLoopNameless begin
  
  Name::AttrType
  vname::Attr(V, Name)

end




@abstract_acset_type AbstractCausalLoop
@acset_type CausalLoopNamelessUntyped(TheoryCausalLoopNameless, index=[:sp,:tp, :sn, :tn]) <: AbstractCausalLoop
@abstract_acset_type AbstractNamedCausalLoop <: AbstractCausalLoop
@acset_type CausalLoopPMUntyped(TheoryCausalLoopPM, index=[:sp,:tp, :sn, :tn]) <: AbstractNamedCausalLoop

const CausalLoopNameless = CausalLoopNamelessUntyped
const CausalLoopPM = CausalLoopPMUntyped{Symbol}


const OpenCausalLoopNamelessOb, OpenCausalLoopNameless = OpenACSetTypes(CausalLoopNamelessUntyped, CausalLoopNamelessUntyped)
const OpenCausalLoopPMOb, OpenCausalLoopPM = OpenACSetTypes(CausalLoopPMUntyped, CausalLoopPMUntyped)




@enum Polarity begin
  POL_POSITIVE
  POL_NEGATIVE
end

"""
  - src -> 
E           V - vname -> Name
  - tgt -> 
"""
@present TheoryCausalLoop <: SchGraph begin
  Name::AttrType
  vname::Attr(V, Name)
end

"""
                      - src -> 
Polarity <- epol - E           V - vname -> Name
                      - tgt -> 
"""
@present TheoryCausalLoopPol <: TheoryCausalLoop begin
  Polarity::AttrType
  epol::Attr(E, Polarity)
end

@abstract_acset_type AbstractSimpleCausalLoop
@acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
@acset_type CausalLoopPolUntyped(TheoryCausalLoopPol, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
const CausalLoop = CausalLoopUntyped{Symbol}
const CausalLoopPol = CausalLoopPolUntyped{Symbol, Polarity}

const OpenCausalLoopPolOb, OpenCausalLoopPol = OpenACSetTypes(CausalLoopPolUntyped, CausalLoopPolUntyped)


vname(c::AbstractSimpleCausalLoop,n) = subpart(c,n,:vname)
vnames(c::AbstractSimpleCausalLoop) = subpart(c, :vname)

vname(c::CausalLoopPM, n) = subpart(c, n, :vname)
vnames(c::CausalLoopPM) = subpart(c, :vname)

ename(c::CausalLoopPol, e) = (vname(c, sedge(c, e)), vname(c, tedge(c, e)), epol(c,e))
enames(c::AbstractSimpleCausalLoop) = [ename(c,e) for e in 1:nedges(c)]

leg(a::CausalLoopPol, x::CausalLoopPol) = OpenACSetLeg(a, E=ntcomponent(enames(a), enames(x)), V=ntcomponent(vnames(a), vnames(x)))
Open(p::CausalLoopPol, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenCausalLoopPol{Symbol,Polarity}(p, legs...)
end

"""
Create a CausalLoop (Graph with named vertices) with a vector of vertices, and
a vector of pairs of vertices.

CausalLoop([:A, :B], [:A => :B, :B => :B]) will create a CausalLoop with 
vertices A and B, an edge A => B and an edge B => B.
"""
function CausalLoop(vs::Vector{Symbol}, es::Vector{Pair{Symbol, Symbol}})
    c = CausalLoop()
    add_parts!(c, :V, length(vs) ; vname = vs)
    vs_idx = state_dict(vs)

    s = map(first, es)
    t = map(last, es)


    add_parts!(c, :E, length(es) ; src=map(x->vs_idx[x], s), tgt=map(x->vs_idx[x], t))

    c
end
"""
Construct a CausalLoop from a StockFlow.
"""
function convertToCausalLoop(p::AbstractStockAndFlowStructure)
    
    sns=snames(p)
    fns=fnames(p)
    svns=svnames(p)
    flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
    vNotf=setdiff(1:nvb(p),flowVariableIndexs)
    vNotfns=[vname(p,v) for v in vNotf]
    
    ns=vcat(sns,fns,svns,vNotfns)

    lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
    lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
    lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
    fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
    foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]

    es=vcat(lses,lsvfes,lfves,fies,foes)

    return CausalLoop(ns,es)
end


"""
Convert a CausalLoopPol to a CausalLoop (forget polarities).
"""
function to_simple_cl(cl::CausalLoopPol)
    c = CausalLoop()
    add_parts!(c, :V, length(vnames(cl)) ; vname = vnames(cl))
    add_parts!(c, :E, nedges(cl) ; src=subpart(cl, :src), tgt=subpart(cl, :tgt))
    c
end

""" Identity function on CausalLoop. """
to_simple_cl(cl::CausalLoop) = cl;

"""
Convert a CausalLoopPM to a CausalLoop (forget polarities).
"""
function to_simple_cl(cl::CausalLoopPM)
    c = CausalLoop()
    add_parts!(c, :V, length(vnames(cl)) ; vname = vnames(cl))
    add_parts!(c, :E, nedges(cl) ; src=vcat(subpart(cl, :sp), subpart(cl, :sm)),
     tgt=vcat(subpart(cl, :tp), subpart(cl, :tm)))
    c
end

"""
Convert StockFlow to CausalLoop.
Nodes: stocks, flows, sum variables, parameters, nonflow dynamic variables
Edges: morphisms in stock flow
"""
function convertToCausalLoop(p::AbstractStockAndFlowStructureF)
    
    sns=snames(p)
    fns=fnames(p)
    svns=svnames(p)
    pns=pnames(p)
    flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
    vNotf=setdiff(1:nvb(p),flowVariableIndexs)
    vNotfns=[vname(p,v) for v in vNotf]
    
    ns=vcat(sns,fns,svns,vNotfns,pns)

    lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
    lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
    lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
    fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
    foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]
    lpvs=[pname(p,subpart(p,lp,:lpvp))=>subpart(p,lp,:lpvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lp,:lpvv),:fv))) : vname(p,subpart(p,lp,:lpvv)) for lp in 1:nlpv(p)]
    lvvs=[subpart(p,lv,:lvsrc) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvsrc),:fv))) : vname(p,subpart(p,lv,:lvsrc))=>subpart(p,lv,:lvtgt) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvtgt),:fv))) : vname(p,subpart(p,lv,:lvtgt)) for lv in 1:nlvv(p)]


    es=vcat(lses,lsvfes,lfves,fies,foes,lpvs,lvvs)

    return CausalLoop(ns,es)
end



"""
Convert a CausalLoopPol to a CausalLoopPM.
"""
function from_clp(cl::CausalLoopPol)
  pols = subpart(cl, :epol)
  names = Dict([i => x for (i, x) in enumerate(subpart(cl, :vname))])
  src = map(x -> names[x], subpart(cl, :src))
  tgt = map(x -> names[x], subpart(cl, :tgt))
  st = Vector{Pair{Symbol, Symbol}}(map(((x,y),) -> x => y, zip(src, tgt)))
  CausalLoopPM(subpart(cl, :vname), st, pols)
end

"""
Create a CausalLoopPol from a vector of node names, and two vectors indicating 
which indices for vertices will act as edges.

to_clp([:A, :B], [1 => 2], Vector{Pair{Int, Int}}()) will create a 
CausalLoopPol with a positive polarity edge from A to B.
"""
function to_clp(nodes::Vector{Symbol}, reinf::Vector{Pair{Int, Int}}, 
  bal::Vector{Pair{Int, Int}})

  ne = length(reinf) + length(bal)
  pols = vcat(
    repeat([POL_POSITIVE], length(reinf)),
    repeat([POL_NEGATIVE], length(bal)),
  )

  edges = vcat(reinf, bal)
  src = map(first, edges)
  tgt = map(last, edges)

  clp = CausalLoopPol()
  add_parts!(clp, :V, length(nodes); vname=nodes)
  add_parts!(clp, :E, ne; src = src, tgt = tgt, epol = pols)

  clp

end

"""
Convert CausalLoopPM to CausalLoopPol.
"""
function to_clp(cl::CausalLoopPM)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    )
end

"""
Polarity multiplication for edge concatenation.

POS * POS == NEG * NEG == POS
POS * NEG == NEG * POS == NEG
"""
function *(p1::Polarity, p2::Polarity)
  if p1 == p2
    return POL_POSITIVE
  else
    return POL_NEGATIVE
  end
end

"""Add vertex to CausalLoopPM."""
add_vertex!(c::AbstractCausalLoop;kw...) = add_part!(c,:V;kw...)
"""Add vertices to CausalLoopPM."""
add_vertices!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:V,n;kw...)

"""Add positive edge to CausalLoopPM."""
add_plus!(c::AbstractCausalLoop;kw) = add_part!(c, :P; kw...)
"""Add positive edges to CausalLoopPM."""
add_pluses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :P, n; kw...)

"""Add negative edge to CausalLoopPM."""
add_minus!(c::AbstractCausalLoop;kw) = add_part!(c, :M; kw...)
"""Add negative edges to CausalLoopPM."""
add_minuses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :M, n; kw...)


"""
Create CausalLoopPM from vector of symbols, vector of pairs of symbols, and
vector of polarities.
"""
CausalLoopPM(ns::Vector{Symbol}, es::Vector{Pair{Symbol, Symbol}}, pols::Vector{Polarity}) = begin
  @assert length(pols) == length(es)

  
  c = CausalLoopPM()
 

  ns = vectorify(ns)
  es = vectorify(es)
  
  ns_idx=state_dict(ns)
  add_vertices!(c, length(ns), vname=ns)

  s=map(first,es)
  t=map(last,es)


  for i in eachindex(pols)
    src = s[i]
    tgt = t[i]

    if pols[i] == POL_POSITIVE
      add_part!(c, :P; sp = ns_idx[src], tp = ns_idx[tgt])
    elseif pols[i] == POL_NEGATIVE
      add_part!(c, :M; sm =  ns_idx[src], tm = ns_idx[tgt])
    end
  end

  c
end




""" Return count of edges of CausalLoop. """
nedges(c::AbstractSimpleCausalLoop) = nparts(c,:E) #edges
""" Return count of vertices of CausalLoop. """
nvert(c::AbstractSimpleCausalLoop) = nparts(c,:V) #vertices

"""Return count of vertices of CausalLoopPM. """
nvert(c::AbstractCausalLoop) = nparts(c, :V)

""" Return count of positive edges of CausalLoopPM. """
np(c::AbstractCausalLoop) = nparts(c, :P)
""" Return count of negative edges of CausalLoopPM. """
nm(c::AbstractCausalLoop) = nparts(c, :M)

""" Return total number of edges of CausalLoopPM (# pos + # neg). """
nedges(c::AbstractCausalLoop) = np(c) + nm(c)





"""
Return source vertex index of an edge of CausalLoopPM by index.
Negative edges come after positive edges.

```julia-repl
julia> using StockFlow; StockFlow.Syntax; 
julia> cl = (@cl A => +B, B => -C, C => +D);
julia> sedge(cl, 3)
2
````
The nodes are ordered A, B, C.
The edges are ordered A => +B, C => +D, B => -C; so, the source index of the 
third edge is B, which has index 2.
"""
sedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :sp) : subpart(c, e - np(c), :sm)
end


"""
Return target vertex index of an edge of CausalLoopPM by index.
Negative edges come after positive edges.
See sedge.
"""
tedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :tp) : subpart(c, e - np(c), :tm)
end

""" CausalLoopPM positive edge source. """
spedge(c::AbstractCausalLoop, e) = subpart(c, e, :sp)
""" CausalLoopPM positive edge target. """
tpedge(c::AbstractCausalLoop, e) = subpart(c, e, :tp)
""" CausalLoopPM positive edge sources. """
spedges(c::AbstractCausalLoop) = subpart(c, :sp)
""" CausalLoopPM positive edge targets. """
tpedges(c::AbstractCausalLoop) = subpart(c, :tp)

""" CausalLoopPM negative edge source. """
smedge(c::AbstractCausalLoop, e) = subpart(c, e, :sm)
""" CausalLoopPM negative edge target. """
tmedge(c::AbstractCausalLoop, e) = subpart(c, e, :tm)
""" CausalLoopPM negative edge sources. """
smedges(c::AbstractCausalLoop) = subpart(c, :sm)
""" CausalLoopPM negative edge targets. """
tmedges(c::AbstractCausalLoop) = subpart(c, :tm)



""" CausalLoopPM, name a positive edge by its source and target.  For composition. """
pname(c::AbstractCausalLoop, e) = (vname(c, spedge(c, e)), vname(c, tpedge(c, e)))
""" CausalLoopPM, pairs of source, target for positive edges.  Used for composition. """
pnames(c::AbstractCausalLoop) = Vector{Tuple{Symbol, Symbol}}([pname(c,e) for e in 1:np(c)])

""" CausalLoopPM, name a negative edge by its source and target.  For composition. """
mname(c::AbstractCausalLoop, e) = (vname(c, smedge(c, e)), vname(c, tmedge(c, e)))
""" CausalLoopPM, pairs of source, target for negative edges.  Used for composition. """
mnames(c::AbstractCausalLoop) = Vector{Tuple{Symbol, Symbol}}([mname(c,e) for e in 1:nm(c)])




""" CausalLoopPol, return edge's name with src number e """
sedge(c::Union{CausalLoopPol,CausalLoop},e) = subpart(c,e,:src)
""" CausalLoopPol, return edge's name with tgt number e """
tedge(c::Union{CausalLoopPol, CausalLoop},e) = subpart(c,e,:tgt)



""" CausalLoopPol, return Polarity of edge e. """
epol(c::CausalLoopPol,e) = subpart(c,e,:epol)

""" CausalLoopPol, return vector of Polarities. """
epols(c::CausalLoopPol) = Vector{Polarity}([epol(c, n) for n in 1:nedges(c)])

""" CausalLoopPol, indices of all edges with src vertex with index n. """
outgoing_edges(c::Union{CausalLoopPol, CausalLoopPM, CausalLoop}, n) = Vector{Int}(collect(filter(i -> sedge(c,i) == n, 1:nedges(c))))
""" CausalLoopPol, indices of all edges with tgt vertex with index n. """
incoming_edges(c::Union{CausalLoopPol,CausalLoopPM, CausalLoop}, n) = Vector{Int}(collect(filter(i -> tedge(c,i) == n, 1:nedges(c))))

""" CausalLoopPM, used for composition. """
leg(a::CausalLoopPM, x::CausalLoopPM) = OpenACSetLeg(a, P=ntcomponent(pnames(a), pnames(x)),  M=ntcomponent(mnames(a), mnames(x)), V=ntcomponent(vnames(a), vnames(x)))

"""Construct an OpenCausalLoopPM with a CausalLoopPM, and any number of additional CausalLoopPM to act as feet. """
Open(p::CausalLoopPM, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenCausalLoopPM{Symbol}(p, legs...)
end



""" Convert a CausalLoopPol to a Graphs package Graph. """
function to_graphs_graph(cl::Union{CausalLoopPol, CausalLoop})
  g = SimpleDiGraph(SimpleEdge.(zip(subpart(cl, :src), subpart(cl, :tgt))))
  Graphs.add_vertices!(g, nvert(cl) - Graphs.nv(g))
  g
end

""" 
  CausalLoopPM, return all cycles of a causal loop as a vector of vectors of int, 
  where positive edges come before negative.

  Each cycle will include each edge at most once.

  Don't count empty vector as a cycle.
"""
function cl_cycles(cl::AbstractCausalLoop)
  cl_cycles(to_clp(cl))
end


""" 
  CausalLoopPol, return all cycles of a causal loop as a vector of vectors of int.

  Each cycle will include each edge at most once.

  Don't count empty vector as a cycle.
"""
function cl_cycles(cl::CausalLoopPol)
  edges = collect(zip(subpart(cl, :src), subpart(cl, :tgt)))
  # Unique are sufficient for making simple graph.
  g = SimpleDiGraph(SimpleEdge.(edges))

  all_cycles = Vector{Vector{Int}}()
  # Edges => Polarity
  for cycle ∈ simplecycles(g)
    cycle_length = length(cycle)
    # Last pair is cycle[end], cycle[1]
    node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
    nonsimple_cycles = Vector{Vector{Int}}()
    for (p1, p2) in node_pairs
      # grab all edges with p1 as start and p2 as end
      push!(nonsimple_cycles, Vector{Int}(intersect(incident(cl, p1, :src), incident(cl, p2, :tgt))))
    end
    # generated_cycles = Vector{Vector{Int}}(Base.Iterators.product(nonsimple_cycles...))
    # For loop instead of comprehension to get around product being multidimensional
    for c in Base.Iterators.product(nonsimple_cycles...)
      push!(all_cycles, collect(c))
    end
  end

  all_cycles

end


"""
Return dict of pairs of edges => polarity.
"""
function extract_loops(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  cycles = cl_cycles(cl)
  Dict{Vector{Int}, Polarity}(map(x -> x => walk_polarity(cl, x), cycles))
end

""" CausalLoopPM, return polarity of edge with index e. """
epol(cl::CausalLoopPM, e) = begin
  @assert e <= nedges(cl)
  e <= np(cl) ? POL_POSITIVE : POL_NEGATIVE
end

"""
CausalLoopPM, polarities of edges.
"""
epols(cl::CausalLoopPM) = vcat(repeat([POL_POSITIVE], np(cl)), repeat([POL_NEGATIVE], nm(cl)))

""" Return polarity of walk.  Empty walk is positive. """
function walk_polarity(cl::K, edges::Vector{Int}) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  @assert is_walk(cl, edges)
  foldl(*, map(x -> epol(cl, x), edges); init = POL_POSITIVE)
end

"""
CausalLoopPol, return all dictionary of path => Polarity, for all paths with nonduplicate edges.
Each path therefore will be at most |E|.
"""
function extract_all_nonduplicate_paths(clp::CausalLoopPol)


  function rec_search!(path, nodes, paths)
    target = tedge(clp, path[end])
    outgoing = Vector{Int}(incident(clp, target, :src)) # all connected edges from target vertex

    for o in outgoing
      outgoing_tgt = tedge(clp, o) # note, clp is defined in outer func
      new_path = [path..., o]
      if outgoing_tgt ∈ nodes
	# If this path is a cycle, we may already have it with a different start.
	# In that case, don't add this one.
	# That is, only add this path if there does not already exist a path with all the same edges.
	if !(any(k -> all(∈(k), new_path), keys(paths)))
	  push!(paths, new_path => paths[path] * epol(clp, o))
        end
        continue
      end

      new_nodes = Set{Int}([nodes..., outgoing_tgt])


      push!(paths, new_path => paths[path] * epol(clp, o))
      rec_search!(new_path, new_nodes, paths)
    end
  end


  paths = Dict{Vector{Int}, Polarity}(Vector{Int}() => POL_POSITIVE)
  for e in 1:nedges(clp)
    nodes = Set{Int}([sedge(clp, e)])
    push!(paths, [e] => epol(clp, e))
    rec_search!([e], nodes, paths)
  end

  paths

end


"""
CausalLoopPM, return all dictionary of path => Polarity, for all paths with nonduplicate edges.
Each path therefore will be at most |E|.
"""
function extract_all_nonduplicate_paths(cl::K) where K <: AbstractCausalLoop
  clp = to_clp(cl)
  extract_all_nonduplicate_paths(clp)
end



"""
Return true if a given list of edges is a walk; that the target for each edge is the
source of the next.  Empty list of edges counts as a walk.

Note, negative edges come after positive edges:
```julia-repl
julia> using StockFlow; using StockFlow.Syntax;
julia> cl = (@cl A => +B, B => -C, C => +D, D => -E);
julia> is_walk(cl, [1,3,2,4])
true
```
"""
function is_walk(cl::Union{CausalLoopPM, CausalLoopPol}, edges::Vector{Int})
    length(edges) == 1 ? only(edges) <= nedges(cl) :
  all(x -> tedge(cl, edges[x]) == sedge(cl, edges[x+1]), eachindex(edges[1:end-1]))
end

"""
Return true if a given list of edges is a walk, and the target of the final edge is
the source of the first.

Empty list of edges does not count as a circuit.
"""
function is_circuit(cl::Union{CausalLoopPM, CausalLoopPol}, edges::Vector{Int})
    length(edges) > 0 && is_walk(cl, edges) && sedge(cl, edges[1]) == tedge(cl, edges[end])
end

"""
Calculate betweenness centrality.
"""
function betweenness(cl::CausalLoop)
  @assert allunique(vnames(cl))
  if nvert(cl) == 0
    # Just deal with this edge case right here.
    Array{Float64}(undef, 0) # tried doing undef, 1, 0, but seemed to turn into 0, 0
  end

  betweenness_cent = fill(0//1, nvert(cl))

  sp = all_shortest_paths(cl)
  # Technically, we should probably also be mapping empty lists to that particular node, but it doesn't affect betweenness
  sp_nodes = map(paths -> (map(path -> (length(path) == 0 ? Vector{Int}() : vcat([sedge(cl, path[1])], (x -> tedge(cl, x)).(path))), paths)), sp) 

  σₛₜ = Matrix{Int}(map(x -> length(x), sp))

  for i in 1:(nvert(cl))
    for j in 1:(nvert(cl))
      if length(sp[i,j]) > 0 && length(sp_nodes[i,j][1]) <= 2 # We don't care about start or end.
        continue
      end
      for path in sp_nodes[i,j]
        for node in path[2:end-1]
          betweenness_cent[node] += (1 // σₛₜ[i,j])
        end
      end
    end
  end

  betweenness_cent
  
end

"""
Calculate betweenness centrality.
"""
function betweenness(cl::Union{CausalLoopPM, CausalLoopPol})
  betweenness(to_simple_cl(cl))
end


"""
Convert CausalLoopPM to a Graphs' library graph.
Note, this removes duplicate edges!
"""
function to_graphs_graph(cl::AbstractCausalLoop)
  to_graphs_graph(to_clp(cl))
end



"""
Count how many loops a variable is on.

A => +B, B => +C, C => +A, C => -A will be treated as 2 loops: [1,2,3] and [1,2,4].
Each variable will be treated as being on 2 loops.

Takes a CausalLoopPM or a CausalLoopPol, and a Symbol.  Throws an error if that Symbol
is the name for more than one variable.
"""
function num_loops_var_on(c::Union{AbstractCausalLoop, CausalLoopPol}, name::Symbol)
  name_index = only(incident(c, name, :vname))
  node_cycles = map(x -> map(y -> tedge(c, y), x), cl_cycles(c)) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), node_cycles)
end


"""
Count how many loops a variable is on, ignoring when two vertices have more than one edge between them.

A => +B, B => +C, C => +A, C => -A will be treated as 1 loop.

Takes a CausalLoopPM or a CausalLoopPol, and a Symbol.  Throws an error if that Symbol
is the name for more than one variable.
"""
function num_indep_loops_var_on(c::K, name::Symbol) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  g = to_graphs_graph(c)
  sc = simplecycles(g)
  name_index = only(incident(c, name, :vname))

  # node_cycles = map(x -> map(y -> tedge(c, y), x),sc) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), sc)
end

"""
Return vector of tuple of (name, num inputs, num outputs)
"""
function num_inputs_outputs(cl::CausalLoopPol)
  @assert allunique(vnames(cl))
  ssvec = Dict{Symbol, Tuple{Int, Int}}()
  for i in 1:nvert(cl)
    push!(ssvec, subpart(cl, i, :vname) => (length(incident(cl, i, :tgt)), length((incident(cl, i, :src))))) # name, num inputs, num outputs
  end
  ssvec
end

"""
Return vector of tuple of (name, num inputs, num outputs)
"""
function num_inputs_outputs(cl::CausalLoopPM)
  num_inputs_outputs(to_clp(cl))
end

"""
Return vector of tuple of (name, num pos inputs, num pos outputs, num neg inputs, num neg outputs)
"""
function num_inputs_outputs_pols(cl::CausalLoopPol)
  @assert allunique(vnames(cl))
  ssvec = Dict{Symbol, Tuple{Int, Int, Int, Int}}() # name, pos in, pos out, neg in, neg out
  for i in 1:nvert(cl)
    push!(ssvec, 
    subpart(cl, i, :vname) => 
    ((count(x -> epol(cl, x) == POL_POSITIVE, incident(cl, i, :tgt))), 
    (count(x -> epol(cl, x) == POL_POSITIVE, incident(cl, i, :src))), 
    (count(x -> epol(cl, x) == POL_NEGATIVE, incident(cl, i, :tgt))), 
    (count(x -> epol(cl, x) == POL_NEGATIVE, incident(cl, i, :src))))
    )
    
  end
  ssvec
end

"""
Return vector of tuple of (name, num pos inputs, num pos outputs, num neg inputs, num neg outputs)
"""
function num_inputs_outputs_pols(cl::CausalLoopPM)
  num_inputs_outputs_pols(to_clp(cl))
end

"""
Return vector of all shortest paths between two nodes.  Takes node names as args.
"""
function shortest_paths(cl::Union{CausalLoopPM, CausalLoopPol, CausalLoop}, s::Symbol, d::Symbol)
  @assert allunique(vnames(cl))
  sindex = only(incident(cl, s, :vname))
  dindex = only(incident(cl, d, :vname))
  shortest_paths(to_simple_cl(cl), sindex, dindex)
end

"""
Return vector of all shortest paths between two nodes.  Takes node indices as args.
"""
function shortest_paths(cl::Union{CausalLoopPM, CausalLoopPol, CausalLoop}, s::Int, d::Int)
  paths = Vector{Vector{Int}}()
  minimum = Inf
 
  function rec_search!(path, nodes)
    if length(path) >= minimum
      return
    end

    outgoing = outgoing_edges(cl, nodes[end])
    for o in outgoing
      if tedge(cl, o) in nodes
        continue
      end
      new_path = [path..., o]
      if tedge(cl, o) == d
        if length(new_path) == minimum
          push!(paths, new_path)
        else # must be shorter
          paths = [new_path]
          minimum = length(new_path)
        end
      else
        rec_search!(new_path, [nodes..., tedge(cl, o)])
      end
    end
  end


  if s == d
    return [Vector{Int}()]
  end

  cl = to_simple_cl(cl)
  rec_search!(Vector{Int}(), [s])

  paths

end


"""
Return Matrix{Vector{Vector{Int}}} of all shortest paths.
First index is src, second is tgt, and it points to a vector of all shortest paths,
represented as a sequence of edge indices.

Shortest path to self is empty.  [Vector{Int}()]

An empty Vector at i, j represents there being no path i -> j.
"""
function all_shortest_paths(cl::Union{CausalLoopPM, CausalLoopPol})
  all_shortest_paths(to_simple_cl(cl))
end

"""
Return Matrix{Vector{Vector{Int}}} of all shortest paths.
First index is src, second is tgt, and it points to a vector of all shortest paths,
represented as a sequence of edge indices.

Shortest path to self is empty.  [Vector{Int}()]

An empty Vector at i, j represents there being no path i -> j.
"""
function all_shortest_paths(cl::CausalLoop)
  all_paths = Matrix{Vector{Vector{Int}}}(undef, (nvert(cl), nvert(cl)))
  for i in 1:(nvert(cl))
    for j in 1:(nvert(cl))
      all_paths[i,j] = Vector{Vector{Int}}()
    end
    push!(all_paths[i,i], Vector{Int}()) # length 0
  end

  for node in 1:nvert(cl) # length 1
    outgoing = outgoing_edges(cl, node)
    for o in outgoing
      target = tedge(cl, o)
      if target == node
        continue
      end
      push!(all_paths[node, target], [o])
    end
  end
  made_change = true
  while made_change
    made_change = false
    for node in 1:nvert(cl)
      for target in 1:nvert(cl)
        if node == target || isempty(all_paths[node,target]) 
          continue
        end
        for tdest in 1:nvert(cl) # O(V^3) but who cares
          if isempty(all_paths[target, tdest])
            continue
          end

          if isempty(all_paths[node, tdest]) || length(all_paths[node, tdest][1]) > length(all_paths[node, target][1]) + length(all_paths[target, tdest][1])
            empty!(all_paths[node, tdest])
            made_change = true
            for p1 in all_paths[node, target]
              for p2 in all_paths[target, tdest]
                push!(all_paths[node, tdest], vcat(p1, p2))
              end
            end
          elseif length(all_paths[node, tdest][1]) == length(all_paths[node, target][1]) + length(all_paths[target, tdest][1])
            test_path = vcat(all_paths[node, target][1], all_paths[target, tdest][1])
            if test_path in all_paths[node, tdest]
              continue
            end
            made_change = true
            for p1 in all_paths[node, target]
              for p2 in all_paths[target, tdest]
                new_path = vcat(p1, p2)
                push!(all_paths[node, tdest], new_path)
              end
            end
          end # the other case is len(A -> B -> C) > len(A -> C), so we don't care.
        end
      end
    end
  end

  all_paths

end

    