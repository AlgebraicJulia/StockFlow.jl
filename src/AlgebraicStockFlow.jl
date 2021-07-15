module AlgebraicStockFlow

export TheoryBoneStockFlow, BoneStockFlow, AbstractBoneStockFlow, ns, nf, ni, no,
  add_stock!, add_stocks!, add_flow!, add_flows!,
  add_inflow!, add_inflows!, add_outflow!, add_outflows!, inflows, outflows,
  TransitionMatrices, vectorfield,
  TheoryLabelledBoneStockFlow, LabelledBoneStockFlow, AbstractLabelledBoneStockFlow, sname, fname,
  TheoryStockFlow, StockFlow, AbstractStockFlow, initialValue, initialValues, funcFlow, funcFlows,
  TheoryLabelledStockFlow, LabelledStockFlow, AbstractLabelledStockFlow,
  Open_S, OpenBoneStockFlow_S, OpenLabelledBoneStockFlow_S, OpenStockFlow_S, OpenLabelledStockFlow_S,
  OpenBoneStockFlowOb_S, OpenLabelledBoneStockFlowOb_S, OpenStockFlowOb_S, OpenLabelledStockFlowOb_S, 
  Open_F, OpenBoneStockFlow_F, OpenLabelledBoneStockFlow_F, OpenStockFlow_F, OpenLabelledStockFlow_F,
  OpenBoneStockFlowOb_F, OpenLabelledBoneStockFlowOb_F, OpenStockFlowOb_F, OpenLabelledStockFlowOb_F

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Present
using Catlab.Theories
using LabelledArrays
using LinearAlgebra: mul!
import Base.+,Base.-

# fake inflow or outflow names. Note, we need fake inflow or outflow if a stock do not have any inflow or outflow
const FK_FLOW_NAME=:F_NONE
const FK_FLOW_INDEX=0

+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)


vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# BoneStockFlow: stock and flow diagram only with the stucture (without attributes -> names and functions)
############

@present TheoryBoneStockFlow(FreeSchema) begin
  F::Ob
  S::Ob
  I::Ob
  O::Ob

  ifn::Hom(I,F)
  is::Hom(I,S)
  ofn::Hom(O,F)
  os::Hom(O,S)
end

const AbstractBoneStockFlow = AbstractACSetType(TheoryBoneStockFlow)
const BoneStockFlow = CSetType(TheoryBoneStockFlow,index=[:ifn,:is,:ofn,:os])

# support compose via stocks where the legs are stocks for the decorated cospan
const OpenBoneStockFlowOb_S, OpenBoneStockFlow_S = OpenCSetTypes(BoneStockFlow,:S) # :S indicates the foot are S, which means composed based on stocks?
Open_S(p::AbstractBoneStockFlow) = OpenBoneStockFlow_S(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractBoneStockFlow, legs...) = OpenBoneStockFlow_S(p, map(l->FinFunction(l, ns(p)), legs)...)
Open_S(n, p::AbstractBoneStockFlow, m) = Open_S(p, n, m)

# support compose via flows where the legs are flows for the decorated cospan
const OpenBoneStockFlowOb_F, OpenBoneStockFlow_F = OpenCSetTypes(BoneStockFlow,:F)
Open_F(p::AbstractBoneStockFlow) = OpenBoneStockFlow_F(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractBoneStockFlow, legs...) = OpenBoneStockFlow_F(p, map(l->FinFunction(l, nf(p)), legs)...)
Open_F(n, p::AbstractBoneStockFlow, m) = Open_F(p, n, m)


# function does not consider empty inlfows or outflows
BoneStockFlow(n,ts...) = begin
  p = BoneStockFlow()
  add_flows!(p, n)
  add_stocks!(p, length(ts))
  for (i,(ins,outs)) in enumerate(ts)
    ins = vectorify(ins)
    outs = vectorify(outs)
    ins = ins[ins .!= FK_FLOW_INDEX]
    outs = outs[outs .!= FK_FLOW_INDEX]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), collect(ins))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), collect(outs))
    end
  end
  p
end

ns(p::AbstractBoneStockFlow) = nparts(p,:S)
nf(p::AbstractBoneStockFlow) = nparts(p,:F)
ni(p::AbstractBoneStockFlow) = nparts(p,:I)
no(p::AbstractBoneStockFlow) = nparts(p,:O)

add_stock!(p::AbstractBoneStockFlow;kw...) = add_part!(p,:S;kw...) # before ; are positional arguments, and after ; are keyword arguments
add_stocks!(p::AbstractBoneStockFlow,n;kw...) = add_parts!(p,:S,n;kw...)

add_flow!(p::AbstractBoneStockFlow;kw...) = add_part!(p,:F;kw...)
add_flows!(p::AbstractBoneStockFlow,n;kw...) = add_parts!(p,:F,n;kw...)

add_inflow!(p::AbstractBoneStockFlow,s,t;kw...) = add_part!(p,:I;is=s,ifn=t,kw...)
add_inflows!(p::AbstractBoneStockFlow,n,s,t;kw...) = add_parts!(p,:I,n;is=s,ifn=t,kw...)

add_outflow!(p::AbstractBoneStockFlow,s,t;kw...) = add_part!(p,:O;os=s,ofn=t,kw...)
add_outflows!(p::AbstractBoneStockFlow,n,s,t;kw...) = add_parts!(p,:O,n;ofn=t,os=s,kw...)

sname(p::AbstractBoneStockFlow,s) = (1:ns(p))[s]
fname(p::AbstractBoneStockFlow,f) = (1:nf(p))[f]

# Note: although indexing makes this pretty fast, it is often faster to bulk-convert -- cite from the package of AlgebraicPetri

#inflows return the index array of in_flows to a specific stock
#outflows return the index array of out_flows to a specific stock
inflows(p::AbstractBoneStockFlow,s) = subpart(p,incident(p,s,:is),:ifn) # subtpart return columns of :is, :ifn, :ofn, :os of the cset table, with the row of returned value of incident(p,s,:ofn)
outflows(p::AbstractBoneStockFlow,s) = subpart(p,incident(p,s,:os),:ofn) # incident: return the incident (indexes) of the column :ofn when :ofn = s

struct TransitionMatrices # row represent flows, column represent stocks; and element of 1 of matrix indicates whether there is a connection between the flow and stock; 0 indicates no connection
  inflow::Matrix{Int}
  outflow::Matrix{Int}
  TransitionMatrices(p::AbstractBoneStockFlow) = begin
    inflow, outflow = zeros(Int,(nf(p),ns(p))), zeros(Int,(nf(p),ns(p)))
    for i in 1:ni(p)
      inflow[subpart(p,i,:ifn),subpart(p,i,:is)] += 1
    end
    for o in 1:no(p)
      outflow[subpart(p,o,:ofn),subpart(p,o,:os)] += 1
    end
    new(inflow,outflow)
  end
end


valueat(x::Number, u, t) = x
valueat(f::Function, u, t) = try f(u,t) catch e f(t) end


vectorfield(pn::AbstractBoneStockFlow) = begin
  tm = TransitionMatrices(pn)
  flows = zeros(nf(pn),ns(pn))
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


@present TheoryLabelledBoneStockFlow <: TheoryBoneStockFlow begin
  Name::Data

  fname::Attr(F, Name)
  sname::Attr(S, Name)
end

const AbstractLabelledBoneStockFlow = AbstractACSetType(TheoryLabelledBoneStockFlow)
const LabelledBoneStockFlowUntyped = ACSetType(TheoryLabelledBoneStockFlow, index=[:ifn,:is,:ofn,:os])
const LabelledBoneStockFlow = LabelledBoneStockFlowUntyped{Symbol}

# open with stocks
const OpenLabelledBoneStockFlowObUntyped_S, OpenLabelledBoneStockFlowUntyped_S = OpenACSetTypes(LabelledBoneStockFlowUntyped,:S)
const OpenLabelledBoneStockFlowOb_S, OpenLabelledBoneStockFlow_S = OpenLabelledBoneStockFlowObUntyped_S{Symbol}, OpenLabelledBoneStockFlowUntyped_S{Symbol}
Open_S(p::AbstractLabelledBoneStockFlow) = OpenLabelledBoneStockFlow_S(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLabelledBoneStockFlow, legs...) = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledBoneStockFlow_S(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)),legs)...)
end
Open_S(n, p::AbstractLabelledBoneStockFlow, m) = Open_S(p, n, m)

# open with flows
const OpenLabelledBoneStockFlowObUntyped_F, OpenLabelledBoneStockFlowUntyped_F = OpenACSetTypes(LabelledBoneStockFlowUntyped,:F)
const OpenLabelledBoneStockFlowOb_F, OpenLabelledBoneStockFlow_F = OpenLabelledBoneStockFlowObUntyped_F{Symbol}, OpenLabelledBoneStockFlowUntyped_F{Symbol}
Open_F(p::AbstractLabelledBoneStockFlow) = OpenLabelledBoneStockFlow_F(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLabelledBoneStockFlow, legs...) = begin
  f_idx = Dict(fname(p, f)=>f for f in 1:nf(p))
  OpenLabelledBoneStockFlow_F(p, map(l->FinFunction(map(i->f_idx[i], l), nf(p)),legs)...)
end
Open_F(n, p::AbstractLabelledBoneStockFlow, m) = Open_F(p, n, m)

# function consider empty inlfows or outflows
LabelledBoneStockFlow(n,ts...) = begin
  p = LabelledBoneStockFlow()
  n = vectorify(n)
  state_idx = state_dict(n)
  add_flows!(p, length(n), fname=n)
  for (name,(ins,outs)) in ts
    i = add_stock!(p, sname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    # filter out the fake flows named :F_NONE
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
    end
  end
  p
end


# Stock and flow diagram with attributes - including stocks and flows names, and related initial values (stocks), and functions for flows

@present TheoryStockFlow <: TheoryBoneStockFlow begin
  FuncFlow::Data
  InitialValue::Data

  funcFlow::Attr(F, FuncFlow)
  initialValue::Attr(S, InitialValue)
end

const AbstractStockFlow = AbstractACSetType(TheoryStockFlow)
const StockFlow = ACSetType(TheoryStockFlow, index=[:ifn,:is,:ofn,:os])

# open with stocks
const OpenStockFlowOb_S, OpenStockFlow_S = OpenACSetTypes(StockFlow,:S)
Open_S(p::AbstractStockFlow{FF,I}) where {FF,I} = OpenStockFlow_S{FF,I}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractStockFlow{FF,I}, legs...) where {FF,I} = OpenStockFlow_S{FF,I}(p, map(l->FinFunction(l, ns(p)), legs)...)
Open_S(n, p::AbstractStockFlow, m) = Open_S(p, n, m)

# open with flows
const OpenStockFlowOb_F, OpenStockFlow_F = OpenACSetTypes(StockFlow,:F)
Open_F(p::AbstractStockFlow{FF,I}) where {FF,I} = OpenStockFlow_F{FF,I}(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractStockFlow{FF,I}, legs...) where {FF,I} = OpenStockFlow_F{FF,I}(p, map(l->FinFunction(l, nf(p)), legs)...)
Open_F(n, p::AbstractStockFlow, m) = Open_F(p, n, m)


# function does not consider empty inlfows or outflows
StockFlow{FF,I}(n,ts...) where {FF,I} = begin
  p = StockFlow{FF,I}()
  add_flows!(p, length(n), funcFlow=n)
  for (i, (initialValue,(ins,outs))) in enumerate(ts)
    i = add_stock!(p, initialValue=initialValue)
    ins = vectorify(ins)
    outs = vectorify(outs)
    ins = ins[ins .!= FK_FLOW_INDEX]
    outs = outs[outs .!= FK_FLOW_INDEX]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), collect(ins))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), collect(outs))
    end
  end
  p
end


initialValue(p::AbstractStockFlow,s) = subpart(p,s,:initialValue)
funcFlow(p::AbstractStockFlow,f) = subpart(p,f,:funcFlow)

initialValues(p::AbstractStockFlow) = map(s->initialValue(p, s), 1:ns(p))
funcFlows(p::AbstractStockFlow) = map(f->funcFlow(p, f), 1:nf(p))

@present TheoryLabelledStockFlow <: TheoryStockFlow begin
  Name::Data

  sname::Attr(S, Name)
  fname::Attr(F, Name)
end

const AbstractLabelledStockFlow = AbstractACSetType(TheoryLabelledStockFlow)
const LabelledStockFlowUntyped = ACSetType(TheoryLabelledStockFlow, index=[:ifn,:is,:ofn,:os])
const LabelledStockFlow{FF,I} = LabelledStockFlowUntyped{FF,I,Symbol}

# open with stocks
const OpenLabelledStockFlowObUntyped_S, OpenLabelledStockFlowUntyped_S = OpenACSetTypes(LabelledStockFlowUntyped,:S)
const OpenLabelledStockFlowOb_S{FF,I} = OpenLabelledStockFlowObUntyped_S{FF,I,Symbol}
const OpenLabelledStockFlow_S{FF,I} = OpenLabelledStockFlowUntyped_S{FF,I,Symbol}
Open_S(p::AbstractLabelledStockFlow{FF,I}) where {FF,I} = OpenLabelledStockFlow_S{FF,I}(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open_S(p::AbstractLabelledStockFlow{FF,I}, legs...) where {FF,I} = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenLabelledStockFlow_S{FF,I}(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)), legs)...)
end
Open_S(n, p::AbstractLabelledStockFlow, m) = Open_S(p, n, m)


# open with flows
const OpenLabelledStockFlowObUntyped_F, OpenLabelledStockFlowUntyped_F = OpenACSetTypes(LabelledStockFlowUntyped,:F)
const OpenLabelledStockFlowOb_F{FF,I} = OpenLabelledStockFlowObUntyped_F{FF,I,Symbol}
const OpenLabelledStockFlow_F{FF,I} = OpenLabelledStockFlowUntyped_F{FF,I,Symbol}
Open_F(p::AbstractLabelledStockFlow{FF,I}) where {FF,I} = OpenLabelledStockFlow_F{FF,I}(p, map(x->FinFunction([x], nf(p)), 1:nf(p))...)
Open_F(p::AbstractLabelledStockFlow{FF,I}, legs...) where {FF,I} = begin
  f_idx = Dict(fname(p, f)=>f for f in 1:nf(p))
  OpenLabelledStockFlow_F{FF,I}(p, map(l->FinFunction(map(i->f_idx[i], l), nf(p)), legs)...)
end
Open_F(n, p::AbstractLabelledStockFlow, m) = Open_F(p, n, m)


# function consider empty inlfows or outflows
LabelledStockFlow{FF,I}(n,ts...) where {FF,I} = begin
  p = LabelledStockFlow{FF,I}()
  n = vectorify(n)
  flows = map(first, collect(n))
  funcFlows = map(last, collect(n))
  state_idx = state_dict(flows)
  add_flows!(p, length(flows), funcFlow=funcFlows, fname=flows)
  for (i, ((name,initialValue),(ins,outs))) in enumerate(ts)
    i = add_stock!(p,initialValue=initialValue, sname=name)
    ins = vectorify(ins)
    outs = vectorify(outs)
    # filter out the fake flows named :F_NONE
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    if length(ins)>0
      add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->state_idx[x], collect(ins)))
    end
    if length(outs)>0
      add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->state_idx[x], collect(outs)))
    end
  end
  p
end


sname(p::Union{AbstractLabelledBoneStockFlow, AbstractLabelledStockFlow},s) = subpart(p,s,:sname)
fname(p::Union{AbstractLabelledBoneStockFlow, AbstractLabelledStockFlow},f) = subpart(p,f,:fname)


initialValue(p::AbstractLabelledStockFlow,s) = subpart(p,s,:initialValue)
funcFlow(p::AbstractLabelledStockFlow,f) = subpart(p,f,:funcFlow)

initialValues(p::AbstractLabelledStockFlow) = begin
  snames = [sname(p, s) for s in 1:ns(p)]
  LVector(;[(snames[s]=>initialValue(p, s)) for s in 1:ns(p)]...)
end

funcFlows(p::AbstractLabelledStockFlow) = begin
  fnames = [fname(p, f) for f in 1:nf(p)]
  LVector(;[(fnames[f]=>funcFlow(p, f)) for f in 1:nf(p)]...)
end

include("visualization.jl")

end
