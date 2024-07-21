export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop, CausalLoopF,
nvert, nedges, vname, np, nm,
sedge, tedge, convertToCausalLoop, nnames, CausalLoopF, epol, epols,
Polarity, POL_POSITIVE, POL_NEGATIVE,
add_node!, add_nodes!, add_edge!, add_edges!, discard_zero_pol,
outgoing_edges, incoming_edges, extract_loops, is_walk, is_circuit, walk_polarity, cl_cycles,
CausalLoopPol, to_clp, from_clp, CausalLoopPM, leg,
extract_all_nonduplicate_paths, num_loops_var_on, num_indep_loops_var_on,
betweenness, to_simple_cl


using MLStyle

using DataMigrations

import Graphs: SimpleDiGraph, simplecycles, SimpleEdge, betweenness_centrality


import Base: *


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


@present TheoryCausalLoop <: SchGraph begin
  Name::AttrType
  vname::Attr(V, Name)
end


@present TheoryCausalLoopPol <: TheoryCausalLoop begin
  Polarity::AttrType
  epol::Attr(E, Polarity)
end

# TODO: Make subtyping a bit more sensible
@abstract_acset_type AbstractSimpleCausalLoop
@acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
@acset_type CausalLoopPolUntyped(TheoryCausalLoopPol, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
const CausalLoop = CausalLoopUntyped{Symbol}
const CausalLoopPol = CausalLoopPolUntyped{Symbol, Polarity}

const OpenCausalLoopPolOb, OpenCausalLoopPol = OpenACSetTypes(CausalLoopPolUntyped, CausalLoopPolUntyped)


vname(c::AbstractSimpleCausalLoop,n) = subpart(c,n,:vname)
vnames(c::AbstractSimpleCausalLoop) = subpart(c, :vname)

vname(c::K, n) where K <: CausalLoopPM = subpart(c, n, :vname)
vnames(c::K) where K <: CausalLoopPM = subpart(c, :vname)

ename(c::CausalLoopPol, e) = (vname(c, sedge(c, e)), vname(c, tedge(c, e)), epol(c,e))
enames(c::AbstractSimpleCausalLoop) = [ename(c,e) for e in 1:nedges(c)]

leg(a::CausalLoopPol, x::CausalLoopPol) = OpenACSetLeg(a, E=ntcomponent(enames(a), enames(x)), V=ntcomponent(vnames(a), vnames(x)))
Open(p::CausalLoopPol, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenCausalLoopPol{Symbol,Polarity}(p, legs...)
end


function CausalLoop(vs, es)
    c = CausalLoop()
    add_parts!(c, :V, length(vs) ; vname = vs)
    vs_idx = state_dict(vs)

    s = map(first, es)
    t = map(last, es)


    add_parts!(c, :E, length(es) ; src=map(x->vs_idx[x], s), tgt=map(x->vs_idx[x], t))

    c
end

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

function to_simple_cl(cl::CausalLoopPol)
    c = CausalLoop()
    add_parts!(c, :V, length(vnames(cl)) ; vname = vnames(cl))
    add_parts!(c, :E, nedges(cl) ; src=subpart(cl, :src), tgt=subpart(cl, :tgt))
    c
end

function to_simple_cl(cl::CausalLoopPM)
    to_simple_cl(to_clp(cl))
end

"""
Convert StockFlow to CLD.
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




function from_clp(cl::CausalLoopPol)
  pols = subpart(cl, :epol)
  names = Dict([i => x for (i, x) in enumerate(subpart(cl, :vname))])
  src = map(x -> names[x], subpart(cl, :src))
  tgt = map(x -> names[x], subpart(cl, :tgt))
  st = Vector{Pair{Symbol, Symbol}}(map(((x,y),) -> x => y, zip(src, tgt)))
  CausalLoopF(subpart(cl, :vname), st, pols)
end


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

function to_clp(cl::CausalLoopPM)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    )
end

# assuming only positive and negative exist
function *(p1::Polarity, p2::Polarity)
  if p1 == p2
    return POL_POSITIVE
  else
    return POL_NEGATIVE
  end
end


add_vertex!(c::AbstractCausalLoop;kw...) = add_part!(c,:V;kw...) 
add_vertices!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:V,n;kw...)

add_plus!(c::AbstractCausalLoop;kw) = add_part!(c, :P; kw...)
add_pluses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :P, n; kw...)

add_minus!(c::AbstractCausalLoop;kw) = add_part!(c, :M; kw...)
add_minuses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :M, n; kw...)



CausalLoopF() = CausalLoopPM()

CausalLoopF(ns::Vector{Symbol}, es::Vector{Pair{Symbol, Symbol}}, pols::Vector{Polarity}) = begin
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




# """ return count of edges of CLD """
# TODO: This should be unnecessary now, since it subtypes graph.  Just use ne and nv like normal.
nedges(c::AbstractSimpleCausalLoop) = nparts(c,:E) #edges
nvert(c::AbstractSimpleCausalLoop) = nparts(c,:V) #vertices


nvert(c::AbstractCausalLoop) = nparts(c, :V)

np(c::AbstractCausalLoop) = nparts(c, :P)
nm(c::AbstractCausalLoop) = nparts(c, :M)


nedges(c::AbstractCausalLoop) = np(c) + nm(c)






sedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :sp) : subpart(c, e - np(c), :sm)
end



tedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :tp) : subpart(c, e - np(c), :tm)
end


spedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :sp)
tpedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :tp)
spedges(c::K) where K <: AbstractCausalLoop = subpart(c, :sp)
tpedges(c::K) where K <: AbstractCausalLoop = subpart(c, :tp)

smedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :sm)
tmedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :tm)
smedges(c::K) where K <: AbstractCausalLoop = subpart(c, :sm)
tmedges(c::K) where K <: AbstractCausalLoop = subpart(c, :tm)




pname(c::K, e) where K <: AbstractCausalLoop = (vname(c, spedge(c, e)), vname(c, tpedge(c, e)))
pnames(c::K) where K <: AbstractCausalLoop = [pname(c,e) for e in 1:np(c)]

mname(c::K, e) where K <: AbstractCausalLoop = (vname(c, smedge(c, e)), vname(c, tmedge(c, e)))
mnames(c::K) where K <: AbstractCausalLoop = [mname(c,e) for e in 1:nm(c)]




""" return node's name with index n """
# nname(c::AbstractCausalLoop,n) = subpart(c,n,:nname) # return the node's name with index of s
""" return edge's name with target number t """
sedge(c::CausalLoopPol,e) = subpart(c,e,:src)
""" return edge's name with edge number e """
tedge(c::CausalLoopPol,e) = subpart(c,e,:tgt)

""" return node names of CLD """
# nnames(c::AbstractCausalLoop) = [nname(c, n) for n in 1:nn(c)]

epol(c::CausalLoopPol,e) = subpart(c,e,:epol)

epols(c::CausalLoopPol) = Vector{Polarity}([epol(c, n) for n in 1:nedges(c)])


outgoing_edges(c::CausalLoopPol, n) = collect(filter(i -> sedge(c,i) == n, 1:ne(c)))
incoming_edges(c::CausalLoopPol, n) = collect(filter(i -> tedge(c,i) == n, 1:ne(c)))


leg(a::CausalLoopPM, x::CausalLoopPM) = OpenACSetLeg(a, P=ntcomponent(pnames(a), pnames(x)),  M=ntcomponent(mnames(a), mnames(x)), V=ntcomponent(vnames(a), vnames(x)))
Open(p::CausalLoopPM, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenCausalLoopPM{Symbol}(p, legs...)
end




function to_graphs_graph(cl::CausalLoopPol)
  g = SimpleDiGraph(SimpleEdge.(zip(subpart(cl, :src), subpart(cl, :tgt))))
  g
end

function cl_cycles(cl::K) where K <: AbstractCausalLoop
  cl_cycles(to_clp(cl))
end

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




function extract_loops(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  cycles = cl_cycles(cl)
  Vector{Pair{Vector{Int}, Polarity}}(map(x -> x => walk_polarity(cl, x), cycles))
end

# TODO: terrible
epol(cl::CausalLoopPM, e) = begin
  @assert e <= nedges(cl)
  e <= np(cl) ? POL_POSITIVE : POL_NEGATIVE
end

function walk_polarity(cl::K, edges::Vector{Int}) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  @assert is_walk(cl, edges)
  foldl(*, map(x -> epol(cl, x), edges); init = POL_POSITIVE)
end


function extract_all_nonduplicate_paths(clp::CausalLoopPol)


  function rec_search!(path, nodes, paths)
    target = tedge(clp, path[end])
    outgoing = Vector{Int}(incident(clp, target, :src)) # all connected edges from target vertex

    for o in outgoing
      outgoing_tgt = tedge(clp, o) # note, clp is defined in outer func
      if outgoing_tgt ∈ nodes
        continue
      end

      new_nodes = Set{Int}([nodes..., outgoing_tgt])

      new_path = [path..., o]
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

function extract_all_nonduplicate_paths(cl::K) where K <: AbstractCausalLoop
  clp = to_clp(cl)
  extract_all_nonduplicate_paths(clp)
end




function is_walk(cl::CausalLoopPM, edges::Vector{Int})
    length(edges) == 1 ? only(edges) <= nedges(cl) :
  all(x -> tedge(cl, edges[x]) == sedge(cl, edges[x+1]), eachindex(edges[1:end-1]))
end

function is_circuit(cl::CausalLoopPM, edges::Vector{Int})
    length(edges) > 0 && is_walk(cl, edges) && sedge(cl, edges[1]) == tedge(cl, edges[end])
end

# TODO: How, pray tell, is this a functor
function betweenness(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  g = to_graphs_graph(cl)
  betweenness_centrality(g)
end


# NOTE: simplecycles returns NODES!!!

#TODO: Deal with nameless AbstractCausalLoop, or create a subtype without nameless

function to_graphs_graph(cl::K) where K <: AbstractCausalLoop
  to_graphs_graph(to_clp(cl))
end

# Acts as if there can be more than one edge between nodes
function num_loops_var_on(c::K, name::Symbol) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  name_index = only(incident(c, name, :vname))
  node_cycles = map(x -> map(y -> tedge(c, y), x), cl_cycles(c)) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), node_cycles)
end

# Acts as if there is at most one edge between nodes
function num_indep_loops_var_on(c::K, name::Symbol) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  g = to_graphs_graph(c)
  sc = simplecycles(g)
  name_index = only(incident(c, name, :vname))

  # node_cycles = map(x -> map(y -> tedge(c, y), x),sc) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), sc)
end


function num_inputs_outputs(cl::CausalLoopPol)
  ssvec = Vector{Tuple{Symbol, Int, Int}}()
  for i in 1:nvert(cl)
    push!(ssvet, (i, incident(cl, i, :src), (incident(cl, i, :tgt))))
  end
  ssvec
end



