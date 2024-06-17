export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop, CausalLoopF,
nn, ne, vname,
sedge, tedge, convertToCausalLoop, nnames, CausalLoopF, epol, epols,
Polarity, POL_ZERO, POL_POSITIVE, POL_NEGATIVE, POL_UNKNOWN, POL_NOT_WELL_DEFINED,
add_node!, add_nodes!, add_edge!, add_edges!, discard_zero_pol,
outgoing_edges, incoming_edges, extract_loops, is_walk, is_circuit, walk_polarity, cl_cycles,
CausalLoopPol, to_clp, from_clp, CausalLoopPM, CausalLoopZ, CausalLoopFull, leg,
nvert


using MLStyle

using DataMigrations

import Graphs: SimpleDiGraph, simplecycles, SimpleEdge, betweenness_centrality

import Catlab.Graphs: nv


import Base: +, *, -, /

# Whenever we map from something with attributes to something without, can think
# of it like mapping the attribute to a type with a single instance, which should
# be equivalent to identity

# Need some more work on the theory to ensure it's functorial, though.

@present TheoryStockAndFlowPol <: TheoryStockAndFlow0 begin
  F::Ob
  I::Ob
  O::Ob

  ifn::Hom(I,F)
  is::Hom(I,S)
  ofn::Hom(O,F)
  os::Hom(O,S)

  P::Ob

  Add::Ob

  Sub::Ob

  Times::Ob

  Pos::Ob
  PosLink::Ob
  Neg::Ob
  NegLink::Ob

  Out::Ob


  la::Hom(Add, Pos)
  ra::Hom(Add, Pos)

  ls::Hom(Sub, Pos)
  rs::Hom(Sub, Neg)

  lt::Hom(Times, Pos)
  rt::Hom(Times, Pos)

  poslinkpos::Hom(PosLink, Pos)
  poslinkout::Hom(PosLink, Out)

  neglinkneg::Hom(NegLink, Neg)
  neglinkout::Hom(NegLink, Out)

  sout::Hom(S, Out)
  svout::Hom(SV, Out)
  fout::Hom(F, Out)
  pout::Hom(P, Out)

  addout::Hom(Add, Out)
  subout::Hom(Sub, Out)
  timesout::Hom(Times, Out)

  fname::Attr(F, Name)
  pname::Attr(P, Name)

  addname::Attr(Add, Name)
  subname::Attr(Sub, Name)
  timesname::Attr(Times, Name)

  # outname::Attr(Out, Name)

end
@abstract_acset_type AbstractStockAndFlowPol <: AbstractStockAndFlow0
@acset_type StockAndFlowUntypedPol(TheoryStockAndFlowPol, index=[:is,:os,:ifn,:ofn, :la, :ra, :ls, :rs, :lt, :rt, :lroutre, :lroutout, :lboutbal, :lboutout, 
  :sout, :svout, :fout, :pout, :aout, :suout, :tout,    :lss,:lssv]) <: AbstractStockAndFlowPol
const StockAndFlowPol = StockAndFlowUntypedPol{Symbol}



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
@abstract_acset_type AbstractNamedCausalLoop <: AbstractCausalLoop
@acset_type CausalLoopPMUntyped(TheoryCausalLoopPM, index=[:sp,:tp, :sn, :tn]) <: AbstractNamedCausalLoop
@acset_type CausalLoopZUntyped(TheoryCausalLoopZ, index=[:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractNamedCausalLoop
@acset_type CausalLoopFullUntyped(TheoryCausalLoopFull, index = [:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractNamedCausalLoop

const CausalLoopNameless = CausalLoopNamelessUntyped
const CausalLoopPM = CausalLoopPMUntyped{Symbol}
const CausalLoopZ = CausalLoopZUntyped{Symbol}
const CausalLoopFull = CausalLoopFullUntyped{Symbol}

const OpenCausalLoopNamelessOb, OpenCausalLoopNameless = OpenACSetTypes(CausalLoopNamelessUntyped, CausalLoopNamelessUntyped)
const OpenCausalLoopPMOb, OpenCausalLoopPM = OpenACSetTypes(CausalLoopPMUntyped, CausalLoopPMUntyped)
const OpenCausalLoopZOb, OpenCausalLoopZ = OpenACSetTypes(CausalLoopZUntyped, CausalLoopZUntyped)
const OpenCausalLoopFullOb, OpenCausalLoopFull = OpenACSetTypes(CausalLoopFullUntyped, CausalLoopFullUntyped)



@enum Polarity begin
  POL_ZERO
  POL_POSITIVE
  POL_NEGATIVE
  POL_UNKNOWN
  POL_NOT_WELL_DEFINED
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






function from_clp(cl::CausalLoopPol)
  pols = subpart(cl, :epol)
  names = Dict([i => x for (i, x) in enumerate(subpart(cl, :vname))])
  src = map(x -> names[x], subpart(cl, :src))
  tgt = map(x -> names[x], subpart(cl, :tgt))
  st = Vector{Pair{Symbol, Symbol}}(map(((x,y),) -> x => y, zip(src, tgt)))
  CausalLoopF(subpart(cl, :vname), st, pols)
end


function to_clp(nodes::Vector{Symbol}, reinf::Vector{Pair{Int, Int}}, 
  bal::Vector{Pair{Int, Int}}, zero::Vector{Pair{Int, Int}}, unknown::Vector{Pair{Int, Int}},
  nwd::Vector{Pair{Int, Int}})

  ne = length(reinf) + length(bal) + length(zero) + length(unknown) + length(nwd)
  pols = vcat(
    repeat([POL_POSITIVE], length(reinf)),
    repeat([POL_NEGATIVE], length(bal)),
    repeat([POL_ZERO], length(zero)),
    repeat([POL_UNKNOWN], length(unknown)),
    repeat([POL_NOT_WELL_DEFINED], length(nwd)),
  )

  edges = vcat(reinf, bal, zero, unknown, nwd)
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
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}()
    )
end

function to_clp(cl::CausalLoopZ)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sz), subpart(cl, :tz)))),
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}()
    )
end

function to_clp(cl::CausalLoopFull)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sz), subpart(cl, :tz)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :su), subpart(cl, :tu)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :snwd), subpart(cl, :tnwd)))),
    )
end
  
function *(p1::Polarity, p2::Polarity)
  if (p1 == POL_ZERO || p2 == POL_ZERO) return POL_ZERO end
  if (p1 == POL_UNKNOWN || p2 == POL_UNKNOWN) return POL_UNKNOWN end
  if (p1 == POL_NOT_WELL_DEFINED || p2 == POL_NOT_WELL_DEFINED) return POL_NOT_WELL_DEFINED end
  if ((p1 == POL_NEGATIVE && p2 == POL_NEGATIVE) || (p1 == POL_POSITIVE && p2 == POL_POSITIVE)) return POL_POSITIVE end
  return POL_NEGATIVE
end

function +(p1::Polarity, p2::Polarity)
  @match (p1, p2) begin
    (POL_ZERO, _) => p2
    (_, POL_ZERO) => p1

    (POL_UNKNOWN, _) || (_, POL_UNKNOWN) => POL_UNKNOWN
    
    (POL_NOT_WELL_DEFINED, _) || (_, POL_NOT_WELL_DEFINED) => POL_NOT_WELL_DEFINED
    (POL_POSITIVE, POL_NEGATIVE) || (POL_NEGATIVE, POL_POSITIVE) => POL_NOT_WELL_DEFINED

    (POL_POSITIVE, POL_POSITIVE) => POL_POSITIVE
    (POL_NEGATIVE, POL_NEGATIVE) => POL_NEGATIVE
  end
end

add_vertex!(c::AbstractCausalLoop;kw...) = add_part!(c,:V;kw...) 
add_vertices!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:V,n;kw...)

add_plus!(c::AbstractCausalLoop;kw) = add_part!(c, :P; kw...)
add_pluses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :P, n; kw...)

add_minus!(c::AbstractCausalLoop;kw) = add_part!(c, :M; kw...)
add_minuses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :M, n; kw...)

# add_edge!(c::AbstractCausalLoop,s,t;kw...) = add_part!(c,:E,s=s,t=t;kw...) 
# add_edges!(c::AbstractCausalLoop,n,s,t;kw...) = add_parts!(c,:E,n,s=s,t=t;kw...)

CausalLoopF()

CausalLoopF() = CausalLoopPM()
# CausalLoopF(g::Graph, pols::Vector{Polarity}) = begin
  
# end
CausalLoopF(ns::Vector{Symbol}, es::Vector{Pair{Symbol, Symbol}}, pols::Vector{Polarity}) = begin
  @assert length(pols) == length(es)

  if (POL_NOT_WELL_DEFINED ∈ pols || POL_UNKNOWN ∈ pols)
    c = CausalLoopFull()
  elseif (POL_ZERO ∈ pols)
    c = CausalLoopZ()
  else
    c = CausalLoopPM()
  end


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




# """ return count of edges of CLD """
# TODO: This should be unnecessary now, since it subtypes graph.  Just use ne and nv like normal.
nedges(c::AbstractSimpleCausalLoop) = nparts(c,:E) #edges
nvert(c::AbstractSimpleCausalLoop) = nparts(c,:V) #vertices




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


spedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :sp)
tpedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :tp)
spedges(c::K) where K <: AbstractCausalLoop = subpart(c, :sp)
tpedges(c::K) where K <: AbstractCausalLoop = subpart(c, :tp)

smedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :sm)
tmedge(c::K, e) where K <: AbstractCausalLoop = subpart(c, e, :tm)
smedges(c::K) where K <: AbstractCausalLoop = subpart(c, :sm)
tmedges(c::K) where K <: AbstractCausalLoop = subpart(c, :tm)

szedge(c::K, e) where K <: Union{CausalLoopZ, CausalLoopFull} =  subpart(c, e, :sz)
tzedge(c::K, e) where K <: Union{CausalLoopZ, CausalLoopFull} =  subpart(c, e, :tz)
szedges(c::K) where K <: Union{CausalLoopZ, CausalLoopFull} =  subpart(c, :sz)
tzedges(c::K) where K <: Union{CausalLoopZ, CausalLoopFull} =  subpart(c, :tz)

snwdedge(c::K, e) where K <: CausalLoopFull =  subpart(c, e, :snwd)
tnwdedge(c::K, e) where K <: CausalLoopFull =  subpart(c, e, :tnwd)
snwdedges(c::K) where K <: CausalLoopFull =  subpart(c, :snwd)
tnwdedges(c::K) where K <: CausalLoopFull =  subpart(c, :tnwd)

suedge(c::K, e) where K <: CausalLoopFull =  subpart(c, e, :su)
tuedge(c::K, e) where K <: CausalLoopFull =  subpart(c, e, :tu)
suedges(c::K) where K <: CausalLoopFull =  subpart(c, :su)
tuedges(c::K) where K <: CausalLoopFull =  subpart(c, :tu)


pname(c::K, e) where K <: AbstractCausalLoop = (vname(c, spedge(c, e)), vname(c, tpedge(c, e)))
pnames(c::K) where K <: AbstractCausalLoop = [pname(c,e) for e in 1:np(c)]

mname(c::K, e) where K <: AbstractCausalLoop = (vname(c, smedge(c, e)), vname(c, tmedge(c, e)))
mnames(c::K) where K <: AbstractCausalLoop = [mname(c,e) for e in 1:nm(c)]

zname(c::K, e) where K <: Union{CausalLoopZ, CausalLoopFull} = (vname(c, szedge(c, e)), vname(c, tzedge(c, e)))
mnames(c::K) where K <: Union{CausalLoopZ, CausalLoopFull} = [zname(c,e) for e in 1:nz(c)]

nwdname(c::K, e) where K <: CausalLoopFull = (vname(c, swdedge(c, e)), vname(c, tnwdedge(c, e)))
nwdnames(c::K) where K <: CausalLoopFull = [nwdname(c,e) for e in 1:nwd(c)]

uname(c::K, e) where K <: CausalLoopFull = (vname(c, suedge(c, e)), vname(c, tuedge(c, e)))
unames(c::K) where K <: CausalLoopFull = [uname(c,e) for e in 1:nu(c)]




""" return node's name with index n """
# nname(c::AbstractCausalLoop,n) = subpart(c,n,:nname) # return the node's name with index of s
""" return edge's name with target number t """
sedge(c::CausalLoopPol,e) = subpart(c,e,:src)
""" return edge's name with edge number e """
tedge(c::CausalLoopPol,e) = subpart(c,e,:tgt)

""" return node names of CLD """
# nnames(c::AbstractCausalLoop) = [nname(c, n) for n in 1:nn(c)]

epol(c::CausalLoopPol,e) = subpart(c,e,:epol)

# epols(c::CausalLoopPol) = [epol(c, n) for n in 1:ne(c)]


outgoing_edges(c::CausalLoopPol, n) = collect(filter(i -> sedge(c,i) == n, 1:ne(c)))
incoming_edges(c::CausalLoopPol, n) = collect(filter(i -> tedge(c,i) == n, 1:ne(c)))


leg(a::CausalLoopPM, x::CausalLoopPM) = OpenACSetLeg(a, P=ntcomponent(pnames(a), pnames(x)),  M=ntcomponent(mnames(a), mnames(x)), V=ntcomponent(vnames(a), vnames(x)))
Open(p::CausalLoopPM, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenCausalLoopPM{Symbol}(p, legs...)
end



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

# function to_graphs_graph(cl::CausalLoopPM)
#   edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
#   g = SimpleDiGraph(SimpleEdge.(edges))
#   g
#   # return (g, np(cl))
# end

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


# function cl_cycles(cl::CausalLoopPol)

#   edges = collect(zip(subpart(cl, :s), subpart(cl, :t)))
#   pair_to_edge = state_dict(edges) # Note, this is unique pairs, not all.
#   # Unique are sufficient for making simple graph.
#   g = SimpleDiGraph(SimpleEdge.(edges))

#   # Edges => Polarity
#   cycle_pol = Vector{Pair{Vector{Int}, Polarity}}()
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     # Last pair is cycle[end], cycle[1]
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
#     nonsimple_cycles = Vector{Vector{Int}}()
#     for (p1, p2) in node_pairs
#       matching_edges = Vector{Int}()
#       # grab all edges with p1 as start and p2 as end
#       # This is the bit that could be made more efficient, it loops over all edges every time
#       for cl_node in 1:ne(cl)
#         if sedge(cl, cl_node) == p1 && tedge(cl, cl_node) == p2
#           push!(matching_edges, cl_node)
#         end
#       end
#       push!(nonsimple_cycles, matching_edges)
#     end

#     for cycle_instance in Base.Iterators.product(nonsimple_cycles...)
#       balancing_count = 0
#       is_unknown = false
#       is_zero = false
#       is_not_well_defined = false
#       for node_number in cycle_instance
#         current_pol = epol(cl, node_number)
#         if current_pol == POL_NEGATIVE
#           balancing_count += 1
#         elseif current_pol == POL_UNKNOWN
#           is_unknown = true
#         elseif current_pol == POL_ZERO
#           is_zero = true
#           break
#         elseif current_pol == POL_NOT_WELL_DEFINED
#           is_not_well_defined = true
#         end
#       end

#       collected_cycle = collect(cycle_instance)

#       if is_zero
#         push!(cycle_pol, collected_cycle => POL_ZERO)
#       elseif is_unknown
#         push!(cycle_pol, collected_cycle => POL_UNKNOWN)
#       elseif is_not_well_defined
#         push!(cycle_pol, collected_cycle => POL_NOT_WELL_DEFINED)
#       elseif iseven(balancing_count)
#         push!(cycle_pol, collected_cycle => POL_POSITIVE)
#       else
#         push!(cycle_pol, collected_cycle => POL_NEGATIVE)
#       end
#     end
#   end

#   cycle_pol
          
# end



# function cl_cycles(cl::CausalLoopPM) 
#   edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
#   # Unique are sufficient for making simple graph.
#   g = SimpleDiGraph(SimpleEdge.(edges))

#   all_cycles = Vector{Vector{Int}}()

#   # Edges => Polarity
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     # Last pair is cycle[end], cycle[1]
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
#     nonsimple_cycles = Vector{Vector{Int}}()
#     for (p1, p2) in node_pairs
#       # grab all edges with p1 as start and p2 as end
#       push!(nonsimple_cycles, Vector{Int}(
#         union(
#           intersect(incident(cl, p1, :sp), incident(cl, p2, :tp)),
#           intersect(incident(cl, p1, :sm), incident(cl, p2, :tm)) .+ np(cl)
#         )
#       )
#       )
#     end
#     # generated_cycles = Vector{Vector{Int}}(Base.Iterators.product(nonsimple_cycles...))
#     # For loop instead of comprehension to get around product being multidimensional
#     for c in Base.Iterators.product(nonsimple_cycles...)
#       push!(all_cycles, collect(c))
#     end
#   end

#   all_cycles

# end



function extract_loops(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  cycles = cl_cycles(cl)
  map(x -> x => walk_polarity(cl, x), cycles)
end

# TODO: terrible
epol(cl::CausalLoopPM, e) = begin
  @assert e <= nedges(cl)
  e <= np(cl) ? POL_POSITIVE : POL_NEGATIVE
end

epol(cl::CausalLoopZ, e) = begin
  @assert e <= nedges(cl)
  if e <= np(cl)
    POL_POSITIVE
  elseif e <= (np(cl) + nm(cl))
    POL_NEGATIVE
  else
    POL_ZERO
  end
end

epol(cl::CausalLoopFull, e) = begin
  @assert e <= nedges(cl)
  if e <= np(cl)
    POL_POSITIVE
  elseif e <= (np(cl) + nm(cl))
    POL_NEGATIVE
  elseif e <= (np(cl) + nm(cl) + nz(cl))
    POL_ZERO
  elseif e <= (np(cl) + nm(cl) + nz(cl) + nu(cl))
    POL_UNKNOWN
  else
    POL_NOT_WELL_DEFINED
  end
end


function walk_polarity(cl::K, edges::Vector{Int}) where K <: Union{AbstractCausalLoop, CausalLoopPol}
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

# function walk_polarity(cl::CausalLoopF, edges::Vector{Int})
#   foldl(*, map(x -> epol(cl, x), edges); init = POL_POSITIVE)
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
function betweenness(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  g = to_graphs_graph(cl)
  betweenness_centrality(g)
end


# NOTE: simplecycles returns NODES!!!

#TODO: Deal with nameless AbstractCausalLoop, or create a subtype without nameless

function to_graphs_graph(cl::K) where K <: AbstractCausalLoop
  to_graphs_graph(to_clp(cl))
  # edges = collect(zip(subpart(cl, :s), subpart(cl, :t)))
  # g = SimpleDiGraph(SimpleEdge.(edges)) # Note, this will discard the final nodes if they have no edges
  # g 
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
  # @show sc
  name_index = only(incident(c, name, :vname))
  # retuer

  # node_cycles = map(x -> map(y -> tedge(c, y), x),sc) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), sc)
end


# function to_catlab_graph(cl::K) where K <: AbstractCausalLoop
#   g = Catlab.Graph(nn(cl))
#   add_parts!(g, :E, ne(cl) ; src = subpart(cl, :s), tgt = subpart(cl, :t))
#   g
# end



# ! This works, but we have version conflicts right now
# ! I added DataMigrations to TOML; presumably, just need to wait a few days for updated requirements in used packages

function to_catlab_graph(cl::CausalLoopNameless)
  mig = @migration SchGraph TheoryCausalLoopNameless begin
    E => @cases (p::P; m::M)
    V => V
    src => begin
      p => sp
      m => sm
    end
    tgt => begin
      p => tp
      m => tm
    end 
  end
  return migrate(Graph, cl, mig)
end

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

F = @finfunctor TheoryStockAndFlowPol TheoryCausalLoopPM begin
  Out => V
  S => V
  F => V
  P => V
  SV => V
  Add => V
  Sub => V
  Times => V
  Pos => V
  Neg => V

  PosLink => P
  I => P
  LS => P

  NegLink => M
  O => M

  lss => sp
  poslinkpos => sp
  is => sp

  lssv => tp
  poslinkout => tp
  ifn => tp

  neglinkneg => sm
  os => sm

  neglinkout => tm
  ofn => tm

  sout => id(V)
  pout => id(V)
  fout => id(V)
  svout => id(V)
  addout => id(V)
  subout => id(V)
  timesout => id(V)
  la => id(V)
  ra => id(V)
  ls => id(V)
  rs => id(V)
  lt => id(V)
  rt => id(V)
  

  Name => Name
  sname => vname
  fname => vname
  pname => vname
  svname => vname
  addname => vname
  subname => vname
  timesname => vname
  
end

ΣF = SigmaMigrationFunctor(F, StockAndFlowPol, CausalLoopPM) 

function migrate_stockflow(sf::StockAndFlowPol)
  
  ΣF(sf)

end


function smallest_required_cl(cl::AbstractCausalLoop) 
  if typeof(cl) == CausalLoopFull
    if length(subpart(cl, :NWD)) == length(subpart(cl, :U)) == 0
      if length(subpart(cl, :Z)) == 0
        clpm = CausalLoopPM()
        add_parts!(clpm, :N, nv(cl))
        add_parts!(clpm, :P, np(cl) ; sp = subpart(cl, :sp), tp = subpart(cl, :tp))
        clpm
      else
        clz = CausalLoopZ()
        add_parts!(clz, :N, nv(cl))
        add_parts!(clz, :P, np(cl) ; sp = subpart(cl, :sp), tp = subpart(cl, :tp))
        add_parts!(clz, :P, np(cl) ; sz = subpart(cl, :sz), tz = subpart(cl, :tz))
        clz
      end
    else
      cl
    end
  end

  if typeof(cl) == CausalLoopZ
    if length(subpart(cl, :Z)) == 0
      clpm = CausalLoopPM()
      add_parts!(clpm, :N, nv(cl))
      add_parts!(clpm, :P, np(cl) ; sp = subpart(cl, :sp), tp = subpart(cl, :tp))
      clpm
    else
      cl
    end
  end

  cl

end


