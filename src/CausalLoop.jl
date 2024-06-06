export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop, CausalLoopF,
nn, ne, nname,
sedge, tedge, convertToCausalLoop, nnames, CausalLoopF, epol, epols,
Polarity, POL_ZERO, POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_NOT_WELL_DEFINED,
add_node!, add_nodes!, add_edge!, add_edges!, discard_zero_pol,
outgoing_edges, incoming_edges, extract_loops, is_walk, is_circuit, walk_polarity, cl_cycles


using MLStyle

# using DataMigrations

import Graphs: SimpleDiGraph, simplecycles, SimpleEdge

import Catlab.Graphs: nv


import Base: +, *, -, /


@present TheoryCausalLoopNameless(FreeSchema) begin
  
  V::Ob
  P::Ob
  M::Ob

  sp::Hom(P, V)
  tp::Hom(P, V)

  sm::Hom(M, V)
  tm::Hom(M, V)

end


@present TheoryCausalLoopPM <: TheoryCausalLoopNameless begin
  
  Name::AttrType
  vname::Attr(V, Name)

end


@present TheoryCausalLoopZ <: TheoryCausalLoopPM begin
  Z::Ob

  sz::Hom(Z, V)
  tz::Hom(Z, V)
end

@present TheoryCausalLoopFull <: TheoryCausalLoopZ begin
  NWD::Ob
  U::Ob

  snwd::Hom(NWD, V)
  tnwd::Hom(NWD, V)

  su::Hom(U, V)
  tu::Hom(U, V)
end


@abstract_acset_type AbstractCausalLoop
@acset_type CausalLoopNamelessUntyped(TheoryCausalLoopNameless, index=[:sp,:tp, :sn, :tn]) <: AbstractCausalLoop
@acset_type CausalLoopPMUntyped(TheoryCausalLoopPM, index=[:sp,:tp, :sn, :tn]) <: AbstractCausalLoop
@acset_type CausalLoopZUntyped(TheoryCausalLoopZ, index=[:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractCausalLoop
@acset_type CausalLoopFullUntyped(TheoryCausalLoopFull, index = [:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractCausalLoop

const CausalLoopNameless = CausalLoopNamelessUntyped
const CausalLoopPM = CausalLoopPMUntyped{Symbol}
const CausalLoopZ = CausalLoopZUntyped{Symbol}
const CausalLoopFull = CausalLoopFullUntyped{Symbol}








# const CausalLoop = CausalLoopUntyped{Symbol} 
# const CausalLoopF = CausalLoopFUntyped{Symbol, Polarity}




@present TheoryCausalLoop(FreeSchema) begin
  E::Ob
  V::Ob

  s::Hom(E,V)
  t::Hom(E,V)

  # Attributes:
  Name::AttrType
  
  vname::Attr(N, Name)
end

@present TheoryCausalLoopPol <: TheoryCausalLoop begin
  Polarity::AttrType
  epol::Attr(E, Polarity)
end

# TODO: Make subtyping a bit more sensible
@abstract_acset_type AbstractSimpleCausalLoop
@acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:s,:t]) <: AbstractSimpleCausalLoop
@acset_type CausalLoopPolUntyped(TheoryCausalLoopPol, index=[:s,:t]) <: AbstractSimpleCausalLoop
const CausalLoop = CausalLoopUntyped{Symbol}
const CausalLoopPol = CausalLoopPolUntyped{Symbol, Polarity}

function from_clp(cl::CausalLoopPol)
  pols = subpart(cl, :epol)
  CausalLoopF(subpart(cl, :vname), collect(zip(subpart(cl, :s), subpart(cl, :t))), pols)
end

function to_clp(cl::CausalLoopPM)
  p = np(cl)
  m = nm(cl)

  clp = CausalLoopPol()
  add_parts!(clp, :V, nparts(cl, :V))
  # !TODO: Rename reinforcing and balancing; probably to 'positive' and 'negative'
  add_parts!(clp, :E, p; s = subpart(cl, :sp), t = subpart(cl, :tp), epol=repeat([POL_REINFORCING], p))
  add_parts!(clp, :E, m; s = subpart(cl, :sm), t = subpart(cl, :tm), epol=repeat([POL_BALANCING], m))

  clp

end

function to_clp(cl::CausalLoopZero)
  clp = CausalLoopPol()

  p = np(cl)
  m = nm(cl)
  z = nz(cl)

  add_parts!(clp, :V, nparts(cl, :V))
  add_parts!(clp, :E, p; s = subpart(cl, :sp), t = subpart(cl, :tp), epol=repeat([POL_REINFORCING], p))
  add_parts!(clp, :E, m; s = subpart(cl, :sm), t = subpart(cl, :tm), epol=repeat([POL_BALANCING], m))
  add_parts!(clp, :Z, m; s = subpart(cl, :sz), t = subpart(cl, :tz), epol=repeat([POL_ZERO], z))

end

function to_clp(cl::CausalLoopFull)
  clp = CausalLoopPol()

  p = np(cl)
  m = nm(cl)
  z = nz(cl)
  nwd = nnwd(cl)
  u = nu(cl)


  add_parts!(clp, :V, nparts(cl, :V))
  add_parts!(clp, :E, p; s = subpart(cl, :sp), t = subpart(cl, :tp), epol=repeat([POL_REINFORCING], p))
  add_parts!(clp, :E, m; s = subpart(cl, :sm), t = subpart(cl, :tm), epol=repeat([POL_BALANCING], m))
  add_parts!(clp, :E, m; s = subpart(cl, :sz), t = subpart(cl, :tz), epol=repeat([POL_ZERO], z))
  add_parts!(clp, :E, m; s = subpart(cl, :snwd), t = subpart(cl, :tnwd), epol=repeat([POL_NOT_WELL_DEFINED], nwd))
  add_parts!(clp, :E, m; s = subpart(cl, :su), t = subpart(cl, :tu), epol=repeat([POL_UNKNOWN], u))
end


@enum Polarity begin
  POL_ZERO
  POL_REINFORCING
  POL_BALANCING
  POL_UNKNOWN
  POL_NOT_WELL_DEFINED
end
  
function *(p1::Polarity, p2::Polarity)
  if (p1 == POL_ZERO || p2 == POL_ZERO) return POL_ZERO end
  if (p1 == POL_UNKNOWN || p2 == POL_UNKNOWN) return POL_UNKNOWN end
  if (p1 == POL_NOT_WELL_DEFINED || p2 == POL_NOT_WELL_DEFINED) return POL_NOT_WELL_DEFINED end
  if ((p1 == POL_BALANCING && p2 == POL_BALANCING) || (p1 == POL_REINFORCING && p2 == POL_REINFORCING)) return POL_REINFORCING end
  return POL_BALANCING
end

function +(p1::Polarity, p2::Polarity)
  @match (p1, p2) begin
    (POL_ZERO, _) => p2
    (_, POL_ZERO) => p1

    (POL_UNKNOWN, _) || (_, POL_UNKNOWN) => POL_UNKNOWN
    
    (POL_NOT_WELL_DEFINED, _) || (_, POL_NOT_WELL_DEFINED) => POL_NOT_WELL_DEFINED
    (POL_REINFORCING, POL_BALANCING) || (POL_BALANCING, POL_REINFORCING) => POL_NOT_WELL_DEFINED

    (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
    (POL_BALANCING, POL_BALANCING) => POL_BALANCING
  end
end

# -(p::Polarity) = p == POL_ZERO ? POL_ZERO : DomainError(p)

# function -(p1::Polarity, p2::Polarity)
#   @match (p1, p2) begin
#     (_, POL_ZERO) => p1
#     (POL_ZERO, POL_ZERO) => POL_ZERO
#     (POL_ZERO, _) => DomainError((p1, p2))

#     (POL_UNKNOWN, _) => POL_UNKNOWN

#     (POL_NOT_WELL_DEFINED, POL_UNKNOWN) => DomainError((p1, p2))
#     (POL_NOT_WELL_DEFINED, _) => POL_NOT_WELL_DEFINED

#     (POL_BALANCING, POL_BALANCING) => POL_BALANCING
#     (POL_BALANCING, _) => DomainError((p1, p2))

#     (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
#     (POL_REINFORCING, _) => DomainError((p1, p2))

#   end
# end

# """p1 is a path of concatenated pols, p2 is what's being removed from the end """
# function /(p1::Polarity, p2::Polarity)
#   @match (p1, p2) begin
#     (POL_BALANCING, POL_BALANCING) || (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
#     (POL_REINFORCING, POL_BALANCING) || (POL_BALANCING, POL_REINFORCING) => POL_BALANCING
#     (POL_ZERO, POL_ZERO) => POL_UNKNOWN # it's possible the previous value was the only zero
    
#     (_, POL_ZERO) => DivideError()

#     (POL_UNKNOWN, POL_UNKNOWN) => POL_UNKNOWN
#     (_, POL_UNKNOWN) => DivideError()

    
#     (POL_NOT_WELL_DEFINED, _) => POL_NOT_WELL_DEFINED
#     (_, POL_NOT_WELL_DEFINED) => DivideError()
#   end
# end




# @present TheoryCausalLoopF <: TheoryCausalLoop begin
#   Polarity::AttrType
#   epolarity::Attr(E, Polarity)
# end

# @abstract_acset_type AbstractCausalLoop
# @acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:s,:t]) <: AbstractCausalLoop
# @acset_type CausalLoopFUntyped(TheoryCausalLoopF, index=[:s,:t]) <: AbstractCausalLoop
# const CausalLoop = CausalLoopUntyped{Symbol} 
# const CausalLoopF = CausalLoopFUntyped{Symbol, Polarity}

add_vertex!(c::AbstractCausalLoop;kw...) = add_part!(c,:V;kw...) 
add_vertices!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:V,n;kw...)

add_plus!(c::AbstractCausalLoop;kw) = add_part!(c, :P; kw...)
add_pluses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :P, n; kw...)

add_minus!(c::AbstractCausalLoop;kw) = add_part!(c, :M; kw...)
add_minuses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :M, n; kw...)

# add_edge!(c::AbstractCausalLoop,s,t;kw...) = add_part!(c,:E,s=s,t=t;kw...) 
# add_edges!(c::AbstractCausalLoop,n,s,t;kw...) = add_parts!(c,:E,n,s=s,t=t;kw...)

"""
    CausalLoop(ns,es)
Create causal loop diagram from collection of nodes and collection of edges.
"""
# CausalLoop(ns,es) = begin
#     c = CausalLoop()
#     ns = vectorify(ns)
#     es = vectorify(es)
    
#     ns_idx=state_dict(ns)
#     add_nodes!(c, length(ns), nname=ns)

#     s=map(first,es)
#     t=map(last,es)
#     add_edges!(c, length(es), map(x->ns_idx[x], s), map(x->ns_idx[x], t))

#     c
# end



CausalLoopF() = CausalLoopPM()
CausalLoopF(ns, es, pols) = begin
  @assert length(pols) == length(es)

  if (POL_NOT_WELL_DEFINED ∈ pols || POL_UNKNOWN ∈ pols)
    c = CausalLoopFull()
  elseif(POL_ZERO ∈ pols)
    c = CausalLoopZ()
  else
    c = CausalLoopPM()
  end


  # c = CausalLoopF()
  ns = vectorify(ns)
  es = vectorify(es)
  
  ns_idx=state_dict(ns)
  add_vertices!(c, length(ns), nname=ns)

  s=map(first,es)
  t=map(last,es)


  for i in eachindex(pols)
    src = s[i]
    tgt = t[i]
    if pols[i] == POL_REINFORCING
      add_part!(c, :P; sp = ns_idx[src], tp = ns_idx[tgt])
    elseif pols[i] == POL_BALANCING
      add_part!(c, :M; sm =  ns_idx[src], tm = ns_idx[tgt])
    elseif pols[i] == POL_ZERO
      add_part!(c, :Z; sz =  ns_idx[src], tz =ns_idx[tgt])
    elseif pols[i] == POL_NOT_WELL_DEFINED
      add_part!(c, :NWD; snwd =  ns_idx[src], tnwd = ns_idx[tgt])
    elseif pols[i] == POL_UNKNOWN
      add_part!(c, :U; su =  ns_idx[src], tu = ns_idx[tgt])
    end
  end

  c
end



# # type piracy, hooray
# # return the count of each components
# """ return count of nodes of CLD """
# nv(c::AbstractCausalLoop) = 

# """ return count of edges of CLD """
# ne(c::AbstractCausalLoop) = nparts(c,:E) #edges
np(c::AbstractCausalLoop) = nparts(c, :P)
nm(c::AbstractCausalLoop) = nparts(c, :M)
nz(c::CausalLoopZ) = nparts(c, :Z)
nnwd(c::CausalLoopFull) = nparts(c, :NWD)
nu(c::CausalLoopFull) = nparts(c, :U)

nedges(c::AbstractCausalLoop) = np(c) + nm(c)
nedges(c::CausalLoopZ) = np(c) + nm(c) + nz(c)
nedges(c::CausalLoopFull) = np(c) + nm(c) + nz(c) + nnwd(c) + nu(c)

sedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :sp) : subpart(c, e - np(c), :sm)
end



tedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :tp) : subpart(c, e - np(c), :tm)
end



""" return node's name with index n """
nname(c::AbstractCausalLoop,n) = subpart(c,n,:nname) # return the node's name with index of s
# """ return edge's name with target number t """
# sedge(c::AbstractCausalLoop,e) = subpart(c,e,:s)
# """ return edge's name with edge number e """
# tedge(c::AbstractCausalLoop,e) = subpart(c,e,:t)

""" return node names of CLD """
nnames(c::AbstractCausalLoop) = [nname(c, n) for n in 1:nn(c)]

# epol(c::CausalLoopF,e) = subpart(c,e,:epolarity)

# epols(c::CausalLoopF) = [epol(c, n) for n in 1:ne(c)]


# outgoing_edges(c::AbstractCausalLoop, n) = collect(filter(i -> sedge(c,i) == n, 1:ne(c)))
# incoming_edges(c::AbstractCausalLoop, n) = collect(filter(i -> tedge(c,i) == n, 1:ne(c)))



# function convertToCausalLoop(p::AbstractStockAndFlowStructure)
    
#     sns=snames(p)
#     fns=fnames(p)
#     svns=svnames(p)
#     flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
#     vNotf=setdiff(1:nvb(p),flowVariableIndexs)
#     vNotfns=[vname(p,v) for v in vNotf]
    
#     ns=vcat(sns,fns,svns,vNotfns)

#     lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
#     lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
#     lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
#     fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
#     foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]

#     es=vcat(lses,lsvfes,lfves,fies,foes)

#     return CausalLoop(ns,es)
# end

"""
Convert StockFlow to CLD.
Nodes: stocks, flows, sum variables, parameters, nonflow dynamic variables
Edges: morphisms in stock flow
"""
# function convertToCausalLoop(p::AbstractStockAndFlowStructureF)
    
#     sns=snames(p)
#     fns=fnames(p)
#     svns=svnames(p)
#     pns=pnames(p)
#     flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
#     vNotf=setdiff(1:nvb(p),flowVariableIndexs)
#     vNotfns=[vname(p,v) for v in vNotf]
    
#     ns=vcat(sns,fns,svns,vNotfns,pns)

#     lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
#     lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
#     lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
#     fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
#     foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]
#     lpvs=[pname(p,subpart(p,lp,:lpvp))=>subpart(p,lp,:lpvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lp,:lpvv),:fv))) : vname(p,subpart(p,lp,:lpvv)) for lp in 1:nlpv(p)]
#     lvvs=[subpart(p,lv,:lvsrc) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvsrc),:fv))) : vname(p,subpart(p,lv,:lvsrc))=>subpart(p,lv,:lvtgt) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvtgt),:fv))) : vname(p,subpart(p,lv,:lvtgt)) for lv in 1:nlvv(p)]


#     es=vcat(lses,lsvfes,lfves,fies,foes,lpvs,lvvs)

#     return CausalLoop(ns,es)
# end

# function from_catlab_graph(g::Catlab.Graph, p::Vector{Polarity})
#   cl = CausalLoopF()
#   add_parts!(cl, :N, Catlab.nn(g))
#   add_parts!(cl, :E, Catlab.ne(g); s = subpart(g, :src), t = subpart(g, :tgt), epolarity = p)
#   cl
# end

# function from_graphs_graph(g::Graphs.Graph, p::Vector{Polarity})
#   cl = CausalLoopF()
#   add_parts!(cl, :N, Graphs.nn(g))
#   add_parts!(cl, :E, Graphs.ne(g); s = subpart(g, :src), t = subpart(g, :tgt), epolarity = p)
# end


# function simple_cycles(c::)
# end

function to_graphs_graph(cl::CausalLoopPM)
  edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
  g = SimpleDiGraph(SimpleEdge.(edges))
  return (g, np(cl))
end

# function cl_cycles(cl::CausalLoopPM)
#   g, last_pos = to_graphs_graph(cl)
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)

#   end
# end



function cl_cycles(cl::CausalLoopPM) 
  edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
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
      push!(nonsimple_cycles, Vector{Int}(
        union(
          intersect(incident(cl, p1, :sp), incident(cl, p2, :tp)),
          intersect(incident(cl, p1, :sm), incident(cl, p2, :tm)) .+ np(cl)
        )
      )
      )
    end
    # generated_cycles = Vector{Vector{Int}}(Base.Iterators.product(nonsimple_cycles...))
    # For loop instead of comprehension to get around product being multidimensional
    for c in Base.Iterators.product(nonsimple_cycles...)
      push!(all_cycles, collect(c))
    end
  end

  all_cycles

end

function extract_loops(cl::CausalLoopPM)
  cycles = cl_cycles(cl)
  map(x -> x => walk_polarity(cl, x), cycles)
end

epol(cl::CausalLoopPM, e) = begin
  @assert e <= nedges(cl)
  e <= np(cl) ? POL_REINFORCING : POL_BALANCING
end

function walk_polarity(cl::CausalLoopPM, edges::Vector{Int})
  foldl(*, map(x -> epol(cl, x), edges); init = POL_REINFORCING)
end

function extract_all_nonduplicate_paths(cl::CausalLoopPM)
  
end

# function walk_polarity(cl::CausalLoopF, edges::Vector{Int})
#   foldl(*, map(x -> epol(cl, x), edges); init = POL_REINFORCING)
# end

""" 
Cycles are uniquely characterized by sets of edges, not sets of nodes
We construct a simple graph, then for each edge, we check if there exist
multiple edges with the same start and end node
We then product all of them, to get every cycle.

This could be made more efficient, but it should be fine for now.
"""
# function extract_loops(cl::CausalLoopF)

#   cycle_pol = Vector{Pair{Vector{Int}, Polarity}}()


#   for cycle_instance in cl_cycles(cl)
#     collected_cycle = collect(cycle_instance)
#     push!(cycle_pol, collected_cycle => walk_polarity(cl, collected_cycle))
#   end

#   cycle_pol
          
# end

function is_walk(cl::CausalLoopPM, edges::Vector{Int})
  all(x -> tedge(cl, edges[x]) == sedge(cl, edges[x+1]), eachindex(edges[1:end-1]))
end

function is_circuit(cl::CausalLoopPM, edges::Vector{Int})
  is_path(cl, edges) && sedge(cl, edges[1]) == tedge(cl, edges[end])
end

# TODO: How, pray tell, is this a functor
function betweenness(cl::K) where K <: AbstractCausalLoop
  g = to_graphs_graph(cl)
  Graphs.betweenness_centrality(g)
end


# function num_loops(cl::CausalLoopF, name::Symbol)
#   el = cl_cycles(cl)
#   node_num = only(incident(cl, :nname, name))
#   return count(x -> node_num ∈ x, el)
# end






# # Graphs.betweenness_centrality
# # ^ Use this if we don't care about polarities
# function betweenness(cl::CausalLoopF; max_edges = typemax(Int))
#   g = to_graphs_graph(cl)
#   num_graph_vert = vertices(g)[end] # ensures we don't go over if the conversion deleted the end 
#   # deleting the end can only happen if there are no edges which come from or leave the node, so the betweenness centrality value for it is 0


#   # past a certain threshold, we start getting infinity confused with very large.

#   # Though we're using typemax(Int) for this so why does this matter
#   @assert length(edges(g)) < max_edges

#   betweenness_values = zeroes(nn(cl))

#   for node in nn(cl)
#     if node > num_graph_vert
#       break
#     end

#     dij = dijkstra_shortest_paths(g, node ; maxdist=max_edges)
#     preds = dij.predecessors

#     # cache the 




#   end



  
# end


function to_graphs_graph(cl::K) where K <: AbstractCausalLoop
  edges = collect(zip(subpart(cl, :s), subpart(cl, :t)))
  g = SimpleDiGraph(SimpleEdge.(edges)) # Note, this will discard the final nodes if they have no edges
  g 
end

function num_loops_var_on(c::CausalLoopPM, name::Symbol)
  name_index = only(incident(c, name, :nname))
  node_cycles = map(x -> map(y -> tedge(c, y), x), cl_cycles(c)) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), node_cycles)
end

function num_indep_loops_var_on(c::CausalLoopPM, name::Symbol)
  g = to_graphs_graph(c)
  sc = simplecycles(g)
  name_index = only(incident(c, name, :nname))
  node_cycles = map(x -> map(y -> tedge(c, y), x),sc) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), node_cycles)
end

# function to_catlab_graph(cl::K) where K <: AbstractCausalLoop
#   g = Catlab.Graph(nn(cl))
#   add_parts!(g, :E, ne(cl) ; src = subpart(cl, :s), tgt = subpart(cl, :t))
#   g
# end



# ! This works, but we have version conflicts right now
# ! I added DataMigrations to TOML; presumably, just need to wait a few days for updated requirements in used packages

# function to_catlab_graph(cl::CausalLoopNameless)
#   mig = @migration SchGraph TheoryCausalLoopNameless begin
#     E => @cases (p::P; m::M)
#     V => V
#     src => begin
#       p => sp
#       m => sm
#     end
#     tgt => begin
#       p => tp
#       m => tm
#     end 
#   end
#   return migrate(Graph, cl, mig)
# end

# function to_graphs_graph(cl::)
# end


# function discard_zero_pol(c)
#   cl = CausalLoopF()
#   add_vertices!(cl, nn(c) ; nname = nnames(c))
#   for edge in 1:ne(c)
#     pol = epol(c, edge)
#     if pol != POL_ZERO
#       add_edge!(cl, sedge(c, edge), tedge(c, edge) ; epolarity = pol)
#     end
#   end
#   cl
# end



