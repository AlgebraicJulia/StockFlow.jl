export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop, CausalLoopF,
nn, ne, nname,
sedge, tedge, convertToCausalLoop, nnames, CausalLoopF, epol

using MLStyle




@present TheoryCausalLoop(FreeSchema) begin
  E::Ob
  N::Ob

  s::Hom(E,N)
  t::Hom(E,N)

  # Attributes:
  Name::AttrType
  
  nname::Attr(N, Name)
end

@enum Polarity begin
  POL_ZERO
  POL_REINFORCING
  POL_BALANCING
  POL_UNKNOWN
end
  

@present TheoryCausalLoopF <: TheoryCausalLoop begin
  Polarity::AttrType
  epolarity::Attr(E, Polarity)
end

@abstract_acset_type AbstractCausalLoop
@acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:s,:t]) <: AbstractCausalLoop
@acset_type CausalLoopFUntyped(TheoryCausalLoopF, index=[:s,:t]) <: AbstractCausalLoop
const CausalLoop = CausalLoopUntyped{Symbol} 
const CausalLoopF = CausalLoopFUntyped{Symbol, Polarity}

add_node!(c::AbstractCausalLoop;kw...) = add_part!(c,:N;kw...) 
add_nodes!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:N,n;kw...)

add_edge!(c::AbstractCausalLoop,s,t;kw...) = add_part!(c,:E,s=s,t=t;kw...) 
add_edges!(c::AbstractCausalLoop,n,s,t;kw...) = add_parts!(c,:E,n,s=s,t=t;kw...)

"""
    CausalLoop(ns,es)
Create causal loop diagram from collection of nodes and collection of edges.
"""
CausalLoop(ns,es) = begin
    c = CausalLoop()
    ns = vectorify(ns)
    es = vectorify(es)
    
    ns_idx=state_dict(ns)
    add_nodes!(c, length(ns), nname=ns)

    s=map(first,es)
    t=map(last,es)
    add_edges!(c, length(es), map(x->ns_idx[x], s), map(x->ns_idx[x], t))

    c
end


CausalLoopF(ns, es, pols) = begin
  @assert length(pols) == length(es)

  c = CausalLoopF()
  ns = vectorify(ns)
  es = vectorify(es)
  
  ns_idx=state_dict(ns)
  add_nodes!(c, length(ns), nname=ns)

  s=map(first,es)
  t=map(last,es)

  add_edges!(c, length(es), map(x->ns_idx[x], s), map(x->ns_idx[x], t), epolarity=pols)

  c
end

# return the count of each components
""" return count of nodes of CLD """
nn(c::AbstractCausalLoop) = nparts(c,:N) #nodes
""" return count of edges of CLD """
ne(c::AbstractCausalLoop) = nparts(c,:E) #edges

""" return node's name with index n """
nname(c::AbstractCausalLoop,n) = subpart(c,n,:nname) # return the node's name with index of s
""" return edge's name with target number t """
sedge(c::AbstractCausalLoop,e) = subpart(c,e,:s)
""" return edge's name with edge number e """
tedge(c::AbstractCausalLoop,e) = subpart(c,e,:t)

""" return node names of CLD """
nnames(c::AbstractCausalLoop) = [nname(c, n) for n in 1:nn(c)]

epol(c::CausalLoopF,e) = subpart(c,e,:epolarity)


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
