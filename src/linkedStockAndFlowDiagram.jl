export TheoryLinkedBoneStockFlow, LinkedBoneStockFlow, AbstractLinkedBoneStockFlow, nl, add_link!, add_links!, lstock, lflow,
OpenLinkedBoneStockFlowOb_S, OpenLinkedBoneStockFlow_S, OpenLinkedBoneStockFlowOb_F, OpenLinkedBoneStockFlow_F,
TheoryLabelledLinkedBoneStockFlow, LabelledLinkedBoneStockFlow, AbstractLabelledLinkedBoneStockFlow,
OpenLabelledLinkedBoneStockFlowOb_S, OpenLabelledLinkedBoneStockFlow_S, OpenLabelledLinkedBoneStockFlowOb_F, OpenLabelledLinkedBoneStockFlow_F,
TheoryLinkedStockFlow, LinkedStockFlow, AbstractLinkedStockFlow,
OpenLinkedStockFlowOb_S, OpenLinkedStockFlow_S, OpenLinkedStockFlowOb_F, OpenLinkedStockFlow_F,
TheoryLabelledLinkedStockFlow, LabelledLinkedStockFlow, AbstractLabelledLinkedStockFlow,
OpenLabelledLinkedStockFlowOb_S, OpenLabelledLinkedStockFlow_S, OpenLabelledLinkedStockFlowOb_F, OpenLabelledLinkedStockFlow_F,
initialValue, funcFlow, initialValues, funcFlows,
TheoryLabelledLinkedStockFlowC, LabelledLinkedStockFlowC, AbstractLabelledLinkedStockFlowC,
OpenLabelledLinkedStockFlowObC_S, OpenLabelledLinkedStockFlowC_S, OpenLabelledLinkedStockFlowObC_F, OpenLabelledLinkedStockFlowC_F


@present TheoryLinkedBoneStockFlow <: TheoryBoneStockFlow begin
  L::Ob

  ls::Hom(L,S)
  lf::Hom(L,F)
end

const AbstractLinkedBoneStockFlow = AbstractACSetType(TheoryLinkedBoneStockFlow)
const LinkedBoneStockFlow = CSetType(TheoryLinkedBoneStockFlow,index=[:ifn,:is,:ofn,:os,:ls,:lf])

# support compose via stocks where the legs are stocks for the decorated cospan
const OpenLinkedBoneStockFlowOb_S, OpenLinkedBoneStockFlow_S = OpenCSetTypes(LinkedBoneStockFlow,:S) # :S indicates the foot are S, which means composed based on stocks?
Open_S(p::AbstractLinkedBoneStockFlow) = OpenLinkedBoneStockFlow_S(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLinkedBoneStockFlow, legs...) = OpenLinkedBoneStockFlow_S(p, map(l->FinFunction(l, ns(p)), legs)...)
Open_S(n, p::AbstractLinkedBoneStockFlow, m) = Open_S(p, n, m)

# support compose via flows where the legs are flows for the decorated cospan
const OpenLinkedBoneStockFlowOb_F, OpenLinkedBoneStockFlow_F = OpenCSetTypes(LinkedBoneStockFlow,:F)
Open_F(p::AbstractLinkedBoneStockFlow) = OpenLinkedBoneStockFlow_F(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLinkedBoneStockFlow, legs...) = OpenLinkedBoneStockFlow_F(p, map(l->FinFunction(l, nf(p)), legs)...)
Open_F(n, p::AbstractLinkedBoneStockFlow, m) = Open_F(p, n, m)

# no. of flows, a stock: inflows, outflows, linked flows
#LinkedBoneStockFlow(6, (1,(2,4),(1,2)), (2,(3,5),(1,2,3)),(3,6,(1,2)))
# function does not consider empty inlfows or outflows
LinkedBoneStockFlow(n,ts...) = begin
  p = LinkedBoneStockFlow()
  add_flows!(p, n)
  add_stocks!(p, length(ts))
  for (i,(ins,outs,fs)) in enumerate(ts)
    ins = vectorify(ins)
    outs = vectorify(outs)
    fs = vectorify(fs)
    ins = ins[ins .!= FK_FLOW_INDEX]
    outs = outs[outs .!= FK_FLOW_INDEX]
    fs = fs[fs .!= FK_FLOW_INDEX]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), collect(ins))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), collect(outs))
    end
    if length(fs)>0
      add_links!(p, length(fs), repeat([i], length(fs)), collect(fs))
    end
  end
  p
end

ns(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = nparts(p,:S)
nf(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = nparts(p,:F)
ni(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = nparts(p,:I)
no(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = nparts(p,:O)

add_stock!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow};kw...) = add_part!(p,:S;kw...) # before ; are positional arguments, and after ; are keyword arguments
add_stocks!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},n;kw...) = add_parts!(p,:S,n;kw...)

add_flow!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow};kw...) = add_part!(p,:F;kw...)
add_flows!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},n;kw...) = add_parts!(p,:F,n;kw...)

add_inflow!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},s,t;kw...) = add_part!(p,:I;is=s,ifn=t,kw...)
add_inflows!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},n,s,t;kw...) = add_parts!(p,:I,n;is=s,ifn=t,kw...)

add_outflow!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},s,t;kw...) = add_part!(p,:O;os=s,ofn=t,kw...)
add_outflows!(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},n,s,t;kw...) = add_parts!(p,:O,n;ofn=t,os=s,kw...)

sname(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},s) = (1:ns(p))[s]
fname(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},f) = (1:nf(p))[f]

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert -- cite from the package of AlgebraicPetri

#inflows return the index array of in_flows to a specific stock
#outflows return the index array of out_flows to a specific stock
inflows(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},s) = subpart(p,incident(p,s,:is),:ifn) # subtpart return columns of :is, :ifn, :ofn, :os of the cset table, with the row of returned value of incident(p,s,:ofn)
outflows(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow},s) = subpart(p,incident(p,s,:os),:ofn) # incident: return the incident (indexes) of the column :ofn when :ofn = s

# return the total number of links in the model
nl(p::AbstractLinkedBoneStockFlow) = nparts(p,:L)
add_link!(p::AbstractLinkedBoneStockFlow,s,f;kw...) = add_part!(p,:L;ls=s,lf=f,kw...)
add_links!(p::AbstractLinkedBoneStockFlow,n,s,f;kw...) = add_parts!(p,:L,n;ls=s,lf=f,kw...)
lstock(p::AbstractLinkedBoneStockFlow,s::Int) = subpart(p,incident(p,s,:ls),:lf)
lflow(p::AbstractLinkedBoneStockFlow,f::Int) = subpart(p,incident(p,f,:lf),:ls)

struct TransitionMatrices # row represent flows, column represent stocks; and element of 1 of matrix indicates whether there is a connection between the flow and stock; 0 indicates no connection
  inflow::Matrix{Int}
  outflow::Matrix{Int}
  link::Matrix{Int}
  TransitionMatrices(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = begin
    inflow, outflow, link = zeros(Int,(nf(p),ns(p))), zeros(Int,(nf(p),ns(p))), zeros(Int,(nf(p),ns(p)))
    for i in 1:ni(p)
      inflow[subpart(p,i,:ifn),subpart(p,i,:is)] += 1
    end
    for o in 1:no(p)
      outflow[subpart(p,o,:ofn),subpart(p,o,:os)] += 1
    end
    if p isa AbstractLinkedBoneStockFlow
      for l in 1:nl(p)
        link[subpart(p,l,:lf),subpart(p,l,:ls)] += 1
      end
    end
    new(inflow,outflow,link)
  end
end

valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end

vectorfield(pn::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow}) = begin
  tm = TransitionMatrices(pn)
  f(du,u,p,t) = begin
    u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
    p_m = [p[fname(pn, i)] for i in 1:nf(pn)]
    for i in 1:ns(pn)
      du[sname(pn, i)] = 0
      for j in 1:nf(pn)
        if tm.inflow[j,i] == 1
          du[sname(pn, i)] = du[sname(pn, i)] + valueat(p_m[j],u,t)
        end
        if tm.outflow[j,i] == 1
          du[sname(pn, i)] = du[sname(pn, i)] - valueat(p_m[j],u,t)
        end
      end
    end
    return du
  end
  return f
end

@present TheoryLabelledLinkedBoneStockFlow <: TheoryLinkedBoneStockFlow begin
  Name::Data

  fname::Attr(F, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledLinkedBoneStockFlow = AbstractACSetType(TheoryLabelledLinkedBoneStockFlow)
const LabelledLinkedBoneStockFlowUntyped = ACSetType(TheoryLabelledLinkedBoneStockFlow, index=[:ifn,:is,:ofn,:os,:ls,:lf])
const LabelledLinkedBoneStockFlow = LabelledLinkedBoneStockFlowUntyped{Symbol}

# open with stocks
const OpenLabelledLinkedBoneStockFlowObUntyped_S, OpenLabelledLinkedBoneStockFlowUntyped_S = OpenACSetTypes(LabelledLinkedBoneStockFlowUntyped,:S)
const OpenLabelledLinkedBoneStockFlowOb_S, OpenLabelledLinkedBoneStockFlow_S = OpenLabelledLinkedBoneStockFlowObUntyped_S{Symbol}, OpenLabelledLinkedBoneStockFlowUntyped_S{Symbol}
Open_S(p::AbstractLabelledLinkedBoneStockFlow) = OpenLabelledLinkedBoneStockFlow_S(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLabelledLinkedBoneStockFlow, legs...) = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledLinkedBoneStockFlow_S(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)),legs)...)
end
Open_S(n, p::AbstractLabelledLinkedBoneStockFlow, m) = Open_S(p, n, m)

# open with flows
const OpenLabelledLinkedBoneStockFlowObUntyped_F, OpenLabelledLinkedBoneStockFlowUntyped_F = OpenACSetTypes(LabelledLinkedBoneStockFlowUntyped,:F)
const OpenLabelledLinkedBoneStockFlowOb_F, OpenLabelledLinkedBoneStockFlow_F = OpenLabelledLinkedBoneStockFlowObUntyped_F{Symbol}, OpenLabelledLinkedBoneStockFlowUntyped_F{Symbol}
Open_F(p::AbstractLabelledLinkedBoneStockFlow) = OpenLabelledLinkedBoneStockFlow_F(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLabelledLinkedBoneStockFlow, legs...) = begin
  f_idx = Dict(fname(p, f)=>f for f in 1:nf(p))
  OpenLabelledLinkedBoneStockFlow_F(p, map(l->FinFunction(map(i->f_idx[i], l), nf(p)),legs)...)
end
Open_F(n, p::AbstractLabelledLinkedBoneStockFlow, m) = Open_F(p, n, m)

#e.g. LabelledLinkedBoneStockFlow([:birth, :inf, :rec, :deathS, :deathI, :deathR], :S=>(:birth,(:inf,:deathS),(:birth,:inf,:deathS)), :I=>(:inf,(:rec,:deathI),(:birth,:inf,:rec,:deathI)),:R=>(:rec,:deathR),(:birth,:inf,:deathR))
# function consider empty inlfows or outflows
LabelledLinkedBoneStockFlow(n,ts...) = begin
  p = LabelledLinkedBoneStockFlow()
  n = vectorify(n)
  state_idx = state_dict(n)
  add_flows!(p, length(n), fname=n)
  for (name,(ins,outs,links)) in ts
    i = add_stock!(p, sname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    links = vectorify(links)
    # filter out the fake flows named :F_NONE
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    links = links[links .!= FK_FLOW_NAME]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
    end
    if length(links)>0
      add_links!(p, length(links), repeat([i], length(links)), map(x->state_idx[x], collect(links)))
    end
  end
  p
end

@present TheoryLinkedStockFlow <: TheoryLinkedBoneStockFlow begin
  FuncFlow::Data
  InitialValue::Data

  funcFlow::Attr(F, FuncFlow)
  initialValue::Attr(S, InitialValue)
end

const AbstractLinkedStockFlow = AbstractACSetType(TheoryLinkedStockFlow)
const LinkedStockFlow = ACSetType(TheoryLinkedStockFlow, index=[:ifn,:is,:ofn,:os,:ls,:lf])

# open with stocks
const OpenLinkedStockFlowOb_S, OpenLinkedStockFlow_S = OpenACSetTypes(LinkedStockFlow,:S)
Open_S(p::AbstractLinkedStockFlow{FF,I}) where {FF,I} = OpenLinkedStockFlow_S{FF,I}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLinkedStockFlow{FF,I}, legs...) where {FF,I} = OpenLinkedStockFlow_S{FF,I}(p, map(l->FinFunction(l, ns(p)), legs)...)
Open_S(n, p::AbstractLinkedStockFlow, m) = Open_S(p, n, m)

# open with flows
const OpenLinkedStockFlowOb_F, OpenLinkedStockFlow_F = OpenACSetTypes(LinkedStockFlow,:F)
Open_F(p::AbstractLinkedStockFlow{FF,I}) where {FF,I} = OpenLinkedStockFlow_F{FF,I}(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLinkedStockFlow{FF,I}, legs...) where {FF,I} = OpenLinkedStockFlow_F{FF,I}(p, map(l->FinFunction(l, nf(p)), legs)...)
Open_F(n, p::AbstractLinkedStockFlow, m) = Open_F(p, n, m)


# function does not consider empty inlfows or outflows
# "FF" represents the parameteric type for the attributes of functions of the flows
# "I" represents the parameteric type for the attributes of the initial value of the stocks

LinkedStockFlow{FF,I}(n,ts...) where {FF,I} = begin
  p = LinkedStockFlow{FF,I}()
  add_flows!(p, length(n), funcFlow=n)
  for (i, (initialValue,(ins,outs,links))) in enumerate(ts)
    i = add_stock!(p, initialValue=initialValue)
    ins = vectorify(ins)
    outs = vectorify(outs)
    links = vectorify(links)
    ins = ins[ins .!= FK_FLOW_INDEX]
    outs = outs[outs .!= FK_FLOW_INDEX]
    links = links[links .!= FK_FLOW_INDEX]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), collect(ins))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), collect(outs))
    end
    if length(links)>0
      add_links!(p, length(links), repeat([i], length(links)), collect(links))
    end
  end
  p
end

initialValue(p::Union{AbstractStockFlow,AbstractLinkedStockFlow},s) = subpart(p,s,:initialValue)
funcFlow(p::Union{AbstractStockFlow,AbstractLinkedStockFlow},f) = subpart(p,f,:funcFlow)
initialValues(p::Union{AbstractStockFlow,AbstractLinkedStockFlow}) = map(s->initialValue(p, s), 1:ns(p))
funcFlows(p::Union{AbstractStockFlow,AbstractLinkedStockFlow}) = map(f->funcFlow(p, f), 1:nf(p))

# LabelledLinkedStockFlow{Function, Int}((:birth=>f_birth, :inf=>f_inf, :rec=>f_rec, :deathS=>f_deathS, :deathI=>f_deathI, :deathR=>f_deathR), (:S, 990)=>(:birth,(:inf,:deathS),(:birth,:inf,:deathS)), (:I, 10)=>(:inf,(:rec,:deathI),(:birth,:inf,:rec,:deathI)),(:R, 0)=>(:rec,:deathR,(:birth,:inf,:deathR)))
@present TheoryLabelledLinkedStockFlow <: TheoryLinkedStockFlow begin
  Name::Data

  sname::Attr(S, Name)
  fname::Attr(F, Name)
end

const AbstractLabelledLinkedStockFlow = AbstractACSetType(TheoryLabelledLinkedStockFlow)
const LabelledLinkedStockFlowUntyped = ACSetType(TheoryLabelledLinkedStockFlow, index=[:ifn,:is,:ofn,:os,:ls,:lf])
const LabelledLinkedStockFlow{FF,I} = LabelledLinkedStockFlowUntyped{FF,I,Symbol}

# open with stocks
const OpenLabelledLinkedStockFlowObUntyped_S, OpenLabelledLinkedStockFlowUntyped_S = OpenACSetTypes(LabelledLinkedStockFlowUntyped,:S)
const OpenLabelledLinkedStockFlowOb_S{FF,I} = OpenLabelledLinkedStockFlowObUntyped_S{FF,I,Symbol}
const OpenLabelledLinkedStockFlow_S{FF,I} = OpenLabelledLinkedStockFlowUntyped_S{FF,I,Symbol}
Open_S(p::AbstractLabelledLinkedStockFlow{FF,I}) where {FF,I} = OpenLabelledLinkedStockFlow_S{FF,I}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLabelledLinkedStockFlow{FF,I}, legs...) where {FF,I} = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledLinkedStockFlow_S{FF,I}(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)), legs)...)
end
Open_S(n, p::AbstractLabelledLinkedStockFlow, m) = Open_S(p, n, m)


# open with flows
const OpenLabelledLinkedStockFlowObUntyped_F, OpenLabelledLinkedStockFlowUntyped_F = OpenACSetTypes(LabelledLinkedStockFlowUntyped,:F)
const OpenLabelledLinkedStockFlowOb_F{FF,I} = OpenLabelledLinkedStockFlowObUntyped_F{FF,I,Symbol}
const OpenLabelledLinkedStockFlow_F{FF,I} = OpenLabelledLinkedStockFlowUntyped_F{FF,I,Symbol}
Open_F(p::AbstractLabelledLinkedStockFlow{FF,I}) where {FF,I} = OpenLabelledLinkedStockFlow_F{FF,I}(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLabelledLinkedStockFlow{FF,I}, legs...) where {FF,I} = begin
  f_idx = Dict(fname(p, f)=>f for f in 1:nf(p))
  OpenLabelledLinkedStockFlow_F{FF,I}(p, map(l->FinFunction(map(i->f_idx[i], l), nf(p)), legs)...)
end
Open_F(n, p::AbstractLabelledLinkedStockFlow, m) = Open_F(p, n, m)


# function consider empty inlfows or outflows
LabelledLinkedStockFlow{FF,I}(n,ts...) where {FF,I} = begin
  p = LabelledLinkedStockFlow{FF,I}()
  n = vectorify(n)
  flows = map(first, collect(n))
  funcFlows = map(last, collect(n))
  state_idx = state_dict(flows)
  add_flows!(p, length(flows), funcFlow=funcFlows, fname=flows)
  for (i, ((name,initialValue),(ins,outs,links))) in enumerate(ts)
    i = add_stock!(p,initialValue=initialValue, sname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    links = vectorify(links)
    # filter out the fake flows named :F_NONE
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    links = links[links .!= FK_FLOW_NAME]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
    end
    if length(links)>0
      add_links!(p, length(links), repeat([i], length(links)), map(x->state_idx[x], collect(links)))
    end
  end
  p
end

initialValue(p::AbstractLabelledLinkedStockFlow,s) = subpart(p,s,:initialValue)
funcFlow(p::AbstractLabelledLinkedStockFlow,f) = subpart(p,f,:funcFlow)

initialValues(p::AbstractLabelledLinkedStockFlow) = begin
  snames = [sname(p, s) for s in 1:ns(p)]
  LVector(;[(snames[s]=>initialValue(p, s)) for s in 1:ns(p)]...)
end

funcFlows(p::AbstractLabelledLinkedStockFlow) = begin
  fnames = [fname(p, f) for f in 1:nf(p)]
  LVector(;[(fnames[f]=>funcFlow(p, f)) for f in 1:nf(p)]...)
end

#############################################
## some tries to let the model can check the flow function constrains of stocks defined by links
# LabelledLinkedStockFlowC{Function, Int}((:birth=>(f_birth,(:S,:I,:R)), :inf=>(f_inf,(:S,:I,:R)), :rec=>(f_rec,:I), :deathS=>(f_deathS,:S), :deathI=>(f_deathI,:I), :deathR=>(f_deathR,:R)), (:S, 990)=>(:birth,(:inf,:deathS),(:birth,:inf,:deathS)), (:I, 10)=>(:inf,(:rec,:deathI),(:birth,:inf,:rec,:deathI)),(:R, 0)=>(:rec,:deathR,(:birth,:inf,:deathR)))
@present TheoryLabelledLinkedStockFlowC <: TheoryLabelledLinkedStockFlow begin
  Lstocks::Data # the dependant stocks each flow's in the user defined function

  fls::Attr(F, Lstocks)
end

const AbstractLabelledLinkedStockFlowC = AbstractACSetType(TheoryLabelledLinkedStockFlowC)
const LabelledLinkedStockFlowUntypedC = ACSetType(TheoryLabelledLinkedStockFlowC, index=[:ifn,:is,:ofn,:os,:ls,:lf])
const LabelledLinkedStockFlowC{FF,I} = LabelledLinkedStockFlowUntypedC{FF,I,Symbol,Array{Symbol,1}}

# open with stocks
const OpenLabelledLinkedStockFlowObUntypedC_S, OpenLabelledLinkedStockFlowUntypedC_S = OpenACSetTypes(LabelledLinkedStockFlowUntypedC,:S)
const OpenLabelledLinkedStockFlowObC_S{FF,I} = OpenLabelledLinkedStockFlowObUntypedC_S{FF,I,Symbol,Array{Symbol,1}}
const OpenLabelledLinkedStockFlowC_S{FF,I} = OpenLabelledLinkedStockFlowUntypedC_S{FF,I,Symbol,Array{Symbol,1}}
Open_S(p::AbstractLabelledLinkedStockFlowC{FF,I}) where {FF,I} = OpenLabelledLinkedStockFlowC_S{FF,I}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLabelledLinkedStockFlowC{FF,I}, legs...) where {FF,I} = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledLinkedStockFlowC_S{FF,I}(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)), legs)...)
end
Open_S(n, p::AbstractLabelledLinkedStockFlowC, m) = Open_S(p, n, m)


# open with flows
const OpenLabelledLinkedStockFlowObUntypedC_F, OpenLabelledLinkedStockFlowUntypedC_F = OpenACSetTypes(LabelledLinkedStockFlowUntypedC,:F)
const OpenLabelledLinkedStockFlowObC_F{FF,I} = OpenLabelledLinkedStockFlowObUntypedC_F{FF,I,Symbol,Array{Symbol,1}}
const OpenLabelledLinkedStockFlowC_F{FF,I} = OpenLabelledLinkedStockFlowUntypedC_F{FF,I,Symbol,Array{Symbol,1}}
Open_F(p::AbstractLabelledLinkedStockFlowC{FF,I}) where {FF,I} = OpenLabelledLinkedStockFlowC_F{FF,I}(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLabelledLinkedStockFlowC{FF,I}, legs...) where {FF,I} = begin
  f_idx = Dict(fname(p, f)=>f for f in 1:nf(p))
  OpenLabelledLinkedStockFlowC_F{FF,I}(p, map(l->FinFunction(map(i->f_idx[i], l), nf(p)), legs)...)
end
Open_F(n, p::AbstractLabelledLinkedStockFlowC, m) = Open_F(p, n, m)

# return the stock name with index of s
sname(p::Union{AbstractLabelledBoneStockFlow, AbstractLabelledStockFlow, AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},s) = subpart(p,s,:sname)
# return the flow name with index of s
fname(p::Union{AbstractLabelledBoneStockFlow, AbstractLabelledStockFlow, AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},f) = subpart(p,f,:fname)

#given stock index, return the linked flows' name
lstock(p::Union{AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},s::Int) = [fname(p,ls) for ls in subpart(p,incident(p,s,:ls),:lf)]
#given flow index, return the linked stocks' name
lflow(p::Union{AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},f::Int) = [sname(p,lf) for lf in subpart(p,incident(p,f,:lf),:ls)]

#given stock name, return the linked flows' name
lstock(p::Union{AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},s::Symbol) = lstock(p,findfirst(isequal(s), p.tables.S.sname))
#given flow name, return the linked stocks' name
lflow(p::Union{AbstractLabelledLinkedBoneStockFlow, AbstractLabelledLinkedStockFlow, AbstractLabelledLinkedStockFlowC},f::Symbol) = lflow(p,findfirst(isequal(f), p.tables.F.fname))


# function consider empty inlfows or outflows
LabelledLinkedStockFlowC{FF,I}(n,ts...) where {FF,I} = begin
  p = LabelledLinkedStockFlowC{FF,I}()
  n = vectorify(n)
  flows = map(first, collect(n))
  funcs_fls=map(last, collect(n))
  funcFlows = map(first,funcs_fls)
  fls = [vectorify(s) for s in map(last,funcs_fls)]
  state_idx = state_dict(flows)
  add_flows!(p, length(flows), funcFlow=funcFlows, fname=flows, fls=fls)
  for (i, ((name,initialValue),(ins,outs,links))) in enumerate(ts)
    i = add_stock!(p,initialValue=initialValue, sname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    links = vectorify(links)
    # filter out the fake flows named :F_NONE
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    links = links[links .!= FK_FLOW_NAME]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
    end
    if length(links)>0
      add_links!(p, length(links), repeat([i], length(links)), map(x->state_idx[x], collect(links)))
    end
  end
  # assertion: each flow's dependant stocks of user defined in the function and in the Stock and Flow Diagram should be the same
  for f in 1:nf(p)
    @assert isequal(sort(lflow(p,subpart(p,f,:fname))),sort(subpart(p,f,:fls))) "For flow: $(subpart(p,f,:fname)), the user-defined function depend stocks $(sort(subpart(p,f,:fls))) is not equal to the link's depend stocks $(sort(lflow(p,subpart(p,f,:fname))))!"
  end
  p
end
