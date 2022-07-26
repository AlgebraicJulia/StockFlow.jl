module AlgebraicStockFlow

export TheoryStockAndFlow0, TheoryStockAndFlow, AbstractStockAndFlow0, AbstractStockAndFlow, StockAndFlow0, StockAndFlow, 
add_flow!, add_flows!, add_stock!, add_stocks!, add_variable!, add_variables!, add_svariable!, add_svariables!, 
add_inflow!, add_inflows!, add_outflow!, add_outflows!, add_Vlink!, add_Vlinks!, add_Slink!, add_Slinks!, add_SVlink!, 
add_SVlinks!, ns, nf, ni, no, nvb, nsv, nls, nlv, nlsv, sname, fname, svname, vname, inflows, outflows, svStocks, 
initialValue, initialValues, funcDynam, flowVariableIndex, funcFlow, funcFlows, funcSV, funcSVs, TransitionMatrices, 
vectorfield, funcFlowsRaw, funcFlowRaw, inflowsAll, outflowsAll,instock,outstock, stockssv, stocksv, svsv, svsstock,
vsstock, vssv, svsstockAllF, vsstockAllF, vssvAllF, StockAndFlowUntyped, StockAndFlowUntyped0, Open, snames, fnames, svnames, vnames,
object_shift_right, foot, leg, lsnames

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
using Catlab.Present
using Catlab.Theories
using LabelledArrays
using LinearAlgebra: mul!
import Base.+,Base.-

# fake names for missed part for declaring the components of stock and flow diagrams
const FK_FLOW_NAME=:F_NONE #fake name of inflows or outflows. e.g., if a stock does not have any inflow or any outflow
const FK_VARIABLE_NAME=:V_NONE #fake name of the auxiliary variables. e.g., if a stock does not link to any (not sum) auxiliary variable
const FK_SVARIABLE_NAME=:SV_NONE #fake name of the sum auxiliary variables. e.g., of a stock does not link to any sum auxiliary variable
const FK_SVVARIABLE_NAME=:SVV_NONE #fake name of the auxiliary variables that links to sum auxiliary variables. e.g., a sum auxiliary variable does not link to any auxiliary variable

+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)
-(f::Function, g::Function) = (x...) -> f(x...) - g(x...)


vectorify(n::Vector) = collect(n)
vectorify(n::Tuple) = length(n) == 1 ? [n] : collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

# define the sub-schema of c0, which includes the three objects: stocks(S), sum-auxiliary-variables(SV), and the linkages between them (LS) to be composed
@present TheoryStockAndFlow0(FreeSchema) begin
  S::Ob
  SV::Ob
  LS::Ob

  lss::Hom(LS,S)
  lssv::Hom(LS,SV)

# Attributes:
  Name::AttrType
  InitialValue::AttrType
  
  initialValue::Attr(S, InitialValue)
  sname::Attr(S, Name)
  svname::Attr(SV, Name)
end

@abstract_acset_type AbstractStockAndFlow0
@acset_type StockAndFlowUntyped0(TheoryStockAndFlow0, index=[:lss,:lssv]) <: AbstractStockAndFlow0
# constrains the attributes data type to be: 
# 1. InitialValue: Real
# 2. Name: Symbol
const StockAndFlow0 = StockAndFlowUntyped0{Symbol,Number} 


# only have stocks
StockAndFlow0(s, initialvalues) = begin

  p0 = StockAndFlow0()
  s = vectorify(s)
  initialvalues = vectorify(initialvalues)
  add_stocks!(p0, length(s), initialValue=initialvalues, sname=s)
  p0

end

# have all the components of stocks, sum auxiliary variables and linkages between them
# Note: it is not possible that discrete sum-auxiliary-variables exists in the C0
StockAndFlow0(s, ssv, initialvalues) = begin

  p0 = StockAndFlow0(s, initialvalues)

  s = vectorify(s)
  ssv = vectorify(ssv)

  sv=unique(map(last, ssv))
  add_svariables!(p0, length(sv), svname=sv)

  sv_idx = state_dict(sv)
  s_idx = state_dict(s)
  add_Slinks!(p0, length(ssv), map(x->s_idx[x], map(first,ssv)), map(x->sv_idx[x], sv))
  p0

end


# define the schema of a general stock and flow diagram
@present TheoryStockAndFlow <: TheoryStockAndFlow0 begin
  F::Ob
  I::Ob
  O::Ob
  V::Ob
  LV::Ob
  LSV::Ob


# TODO: should give constrains that make the ifn and ofn to be monomorphisms!!
  ifn::Hom(I,F)
  is::Hom(I,S)
  ofn::Hom(O,F)
  os::Hom(O,S)
  fv::Hom(F,V)
  lvs::Hom(LV,S)
  lvv::Hom(LV,V)
  lsvsv::Hom(LSV,SV)
  lsvv::Hom(LSV,V)


# Attributes:
  FuncDynam::AttrType

  funcDynam::Attr(V, FuncDynam)
  fname::Attr(F, Name)
  vname::Attr(V, Name)
end

@abstract_acset_type AbstractStockAndFlow <: AbstractStockAndFlow0
@acset_type StockAndFlowUntyped(TheoryStockAndFlow, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv]) <: AbstractStockAndFlow
# constrains the attributes data type to be: 
# 1. InitialValue: Real
# 2. Name: Symbol
# 3. FuncDynam: Function
# Note: those three (or any subgroups) attributes' datatype can be defined by the users. See the example of the PetriNet which allows the Reactionrate and Concentration defined by the users
const StockAndFlow = StockAndFlowUntyped{Symbol,Number,Function} 

# functions of adding components of the model schema
add_flow!(p::AbstractStockAndFlow,v;kw...) = add_part!(p,:F;fv=v,kw...)
add_flows!(p::AbstractStockAndFlow,v,n;kw...) = add_parts!(p,:F,n;fv=v,kw...)

add_stock!(p::AbstractStockAndFlow0;kw...) = add_part!(p,:S;kw...) 
add_stocks!(p::AbstractStockAndFlow0,n;kw...) = add_parts!(p,:S,n;kw...)

add_variable!(p::AbstractStockAndFlow;kw...) = add_part!(p,:V;kw...) 
add_variables!(p::AbstractStockAndFlow,n;kw...) = add_parts!(p,:V,n;kw...)

add_svariable!(p::AbstractStockAndFlow0;kw...) = add_part!(p,:SV;kw...) 
add_svariables!(p::AbstractStockAndFlow0,n;kw...) = add_parts!(p,:SV,n;kw...)

add_inflow!(p::AbstractStockAndFlow,s,t;kw...) = add_part!(p,:I;is=s,ifn=t,kw...)
add_inflows!(p::AbstractStockAndFlow,n,s,t;kw...) = add_parts!(p,:I,n;is=s,ifn=t,kw...)

add_outflow!(p::AbstractStockAndFlow,s,t;kw...) = add_part!(p,:O;os=s,ofn=t,kw...)
add_outflows!(p::AbstractStockAndFlow,n,s,t;kw...) = add_parts!(p,:O,n;ofn=t,os=s,kw...)

# links from Stock to dynamic variable
add_Vlink!(p::AbstractStockAndFlow,s,v;kw...) = add_part!(p,:LV;lvs=s,lvv=v,kw...)
add_Vlinks!(p::AbstractStockAndFlow,n,s,v;kw...) = add_parts!(p,:LV,n;lvs=s,lvv=v,kw...)

# links from Stock to sum dynamic variable
add_Slink!(p::AbstractStockAndFlow0,s,sv;kw...) = add_part!(p,:LS;lss=s,lssv=sv,kw...)
add_Slinks!(p::AbstractStockAndFlow0,n,s,sv;kw...) = add_parts!(p,:LS,n;lss=s,lssv=sv,kw...)

# links from sum dynamic variable to dynamic varibale
add_SVlink!(p::AbstractStockAndFlow,sv,v;kw...) = add_part!(p,:LSV;lsvsv=sv,lsvv=v,kw...)
add_SVlinks!(p::AbstractStockAndFlow,n,sv,v;kw...) = add_parts!(p,:LSV,n;lsvsv=sv,lsvv=v,kw...)

# return the count of each components
ns(p::AbstractStockAndFlow0) = nparts(p,:S) #stocks
nf(p::AbstractStockAndFlow) = nparts(p,:F) #flows
ni(p::AbstractStockAndFlow) = nparts(p,:I) #inflows
no(p::AbstractStockAndFlow) = nparts(p,:O) #outflows
nvb(p::AbstractStockAndFlow) = nparts(p,:V) #auxiliary variables
nsv(p::AbstractStockAndFlow0) = nparts(p,:SV) #sum auxiliary variables
nls(p::AbstractStockAndFlow0) = nparts(p,:LS) #links from Stock to sum dynamic variable
nlv(p::AbstractStockAndFlow) = nparts(p,:LV) #links from Stock to dynamic variable
nlsv(p::AbstractStockAndFlow) = nparts(p,:LSV) #links from sum dynamic variable to dynamic varibale


#EXAMPLE:
#sir_StockAndFlow=StockAndFlow(((:S, 990)=>(:birth,(:inf,:deathS),(:v_inf,:v_deathS),:N), (:I, 10)=>(:inf,(:rec,:deathI),(:v_rec,:v_deathI,:v_fractionNonS),:N),(:R, 0)=>(:rec,:deathR,(:v_deathR,:v_fractionNonS),:N)),
#  (:birth=>:v_birth,:inf=>:v_inf,:rec=>:v_rec,:deathS=>:v_deathS,:deathI=>:v_deathI,:deathR=>:v_deathR),
#  (:v_birth=>f_birth,:v_inf=>f_inf,:v_rec=>f_rec,:v_deathS=>f_deathS,:v_deathI=>f_deathI,:v_deathR=>f_deathR),
#  (:N=>(:v_birth,:v_inf))))

StockAndFlow(s,f,v,sv) = begin

    p = StockAndFlow()

    s = vectorify(s)
    f = vectorify(f)
    v = vectorify(v)
    sv = vectorify(sv)

    sname = map(first,map(first,s))
    fname = map(first, f)
    vname = map(first, v)
    svname = map(first, sv)
    s_idx = state_dict(sname)
    f_idx = state_dict(fname)
    v_idx = state_dict(vname)
    sv_idx = state_dict(svname)

    # adding the objects that do not have out-morphisms firstly
    add_variables!(p, length(vname), vname=vname, funcDynam=map(last, v))    # add objects :V (auxiliary variables)
    add_svariables!(p, length(svname), svname=svname)    # add objects :SV (sum auxiliary variables)
    add_flows!(p,map(x->v_idx[x], map(last,f)),length(fname),fname=fname)    # add objects :F (flows)

    # Parse the elements included in "s" -- stocks
    for (i, ((name, initialValue),(ins,outs,vs,svs))) in enumerate(s)
      i = add_stock!(p,initialValue=initialValue, sname=name) # add objects :S (stocks)
      ins=vectorify(ins) # inflows of each stock
      outs=vectorify(outs) # outflows of each stock
      vs=vectorify(vs) # auxiliary variables depends on the stock
      svs=vectorify(svs) # sum auxiliary variables depends on the stock
      # filter out the fake (empty) elements
      ins = ins[ins .!= FK_FLOW_NAME]
      outs = outs[outs .!= FK_FLOW_NAME]
      vs = vs[vs .!= FK_VARIABLE_NAME]
      svs = svs[svs .!= FK_SVARIABLE_NAME]
      if length(ins)>0
        add_inflows!(p, length(ins), repeat([i], length(ins)), map(x->f_idx[x], ins)) # add objects :I (inflows)
      end
      if length(outs)>0
        add_outflows!(p, length(outs), repeat([i], length(outs)), map(x->f_idx[x], outs)) # add objects :O (outflows)
      end
      if length(vs)>0
        add_Vlinks!(p, length(vs), repeat([i], length(vs)), map(x->v_idx[x], vs)) # add objects :LV (links from Stock to dynamic variable)
      end
      if length(svs)>0
        add_Slinks!(p, length(svs), repeat([i], length(svs)), map(x->sv_idx[x], svs)) # add objects :LS (links from Stock to sum dynamic variable)
      end
    end

    # Parse the elements included in "sv" -- sum auxiliary vairables
    for (i, (svname,vs)) in enumerate(sv)
      vs=vectorify(vs)
      vs = vs[vs .!= FK_SVVARIABLE_NAME]
      if length(vs)>0
        add_SVlinks!(p, length(vs), repeat(collect(sv_idx[svname]), length(vs)), map(x->v_idx[x], collect(vs)))
      end
    end
    p
end


sname(p::AbstractStockAndFlow0,s) = subpart(p,s,:sname) # return the stocks name with index of s
fname(p::AbstractStockAndFlow,f) = subpart(p,f,:fname) # return the flows name with index of f
svname(p::AbstractStockAndFlow0,sv) = subpart(p,sv,:svname) # return the sum auxiliary variables name with index of sv
vname(p::AbstractStockAndFlow,v) = subpart(p,v,:vname) # return the auxiliary variables name with index of v

snames(p::AbstractStockAndFlow0) = [sname(p, s) for s in 1:ns(p)]
fnames(p::AbstractStockAndFlow) = [fname(p, f) for f in 1:nf(p)]
svnames(p::AbstractStockAndFlow0) = [svname(p, sv) for sv in 1:nsv(p)]
vnames(p::AbstractStockAndFlow) = [vname(p, v) for v in 1:nv(p)]

# return the pair of names of (stock, sum-auxiliary-variable) for all linkages between them
lsnames(p::AbstractStockAndFlow0) = begin
    s = map(x->subpart(p,x,:lss),collect(1:nls(p)))
    sv = map(x->subpart(p,x,:lssv),collect(1:nls(p)))
    sn = map(x->sname(p,x),s)
    svn = map(x->svname(p,x),sv)
    pssv = collect(zip(sn, svn))
end

# return inflows of stock index s
inflows(p::AbstractStockAndFlow,s) = subpart(p,incident(p,s,:is),:ifn) 
# return outflows of stock index s
outflows(p::AbstractStockAndFlow,s) = subpart(p,incident(p,s,:os),:ofn) 
# return stocks of flow index f flow in
# TODO: add assertion that the length(instock)=1
instock(p::AbstractStockAndFlow,f) = subpart(p,incident(p,f,:ifn),:is) 
# return stocks of flow index f flow out
# TODO: add assertion that the length(outstock)=1
outstock(p::AbstractStockAndFlow,f) = subpart(p,incident(p,f,:ofn),:os) 
# return stocks of sum variable index sv link to
stockssv(p::AbstractStockAndFlow0,sv) = subpart(p,incident(p,sv,:lssv),:lss) 
# return stocks of auxiliary variable index v link to
stocksv(p::AbstractStockAndFlow,v) = subpart(p,incident(p,v,:lvv),:lvs) 
# return sum variables of auxiliary variable index v link to
svsv(p::AbstractStockAndFlow,v) = subpart(p,incident(p,v,:lsvv),:lsvsv) 
# return sum auxiliary variables a stock s link 
svsstock(p::AbstractStockAndFlow,s) = subpart(p,incident(p,s,:lss),:lssv)
# return auxiliary variables a stock s link 
vsstock(p::AbstractStockAndFlow,s) = subpart(p,incident(p,s,:lvs),:lvv)
# return auxiliary variables a sum auxiliary variable link 
vssv(p::AbstractStockAndFlow,sv) = subpart(p,incident(p,sv,:lsvsv),:lsvv)


# return sum auxiliary variables all stocks link (frequency)
svsstockAllF(p::AbstractStockAndFlow) = [((svsstock(p, s) for s in 1:ns(p))...)...]
# return auxiliary variables all stocks link (frequency)
vsstockAllF(p::AbstractStockAndFlow) = [((vsstock(p, s) for s in 1:ns(p))...)...]
# return auxiliary variables all sum auxiliary variables link (frequency)
vssvAllF(p::AbstractStockAndFlow) = [((vssv(p, sv) for sv in 1:nsv(p))...)...]

# return all inflows
inflowsAll(p::AbstractStockAndFlow) = [((inflows(p, s) for s in 1:ns(p))...)...]
# return all outflows
outflowsAll(p::AbstractStockAndFlow) = [((outflows(p, s) for s in 1:ns(p))...)...]

# return the initial value given an index s
initialValue(p::AbstractStockAndFlow0,s) = subpart(p,s,:initialValue)
# return the LVector of pairs: sname => initialValues
initialValues(p::AbstractStockAndFlow0) = begin
  sn = snames(p)
  LVector(;[(sn[s]=>initialValue(p, s)) for s in 1:ns(p)]...)
end
# return the functions of variables give index v
funcDynam(p::AbstractStockAndFlow,v) = subpart(p,v,:funcDynam)
# return the auxiliary variable's index that related to the flow with index of f
flowVariableIndex(p::AbstractStockAndFlow,f) = subpart(p,f,:fv)
# return the functions (not substitutes the function of sum variables yet) of flow index f
funcFlowRaw(p::AbstractStockAndFlow,f)=funcDynam(p,flowVariableIndex(p,f))
# return the LVector of pairs: fname => function (raw: not substitutes the function of sum variables yet)
funcFlowsRaw(p::AbstractStockAndFlow) = begin
  fnames = [fname(p, f) for f in 1:nf(p)]
  LVector(;[(fnames[f]=>funcFlowRaw(p, f)) for f in 1:nf(p)]...)
end
# generate the function substituting sum variables in with flow index fn
funcFlow(p::AbstractStockAndFlow, fn) = begin
    func=funcFlowRaw(p,fn)
    f(u,t) = begin
        uN=funcSVs(p)
        return functioneat(func,u,uN,t)
    end
end
# return the LVector of pairs: fname => function (with function of sum variables substitue in)
funcFlows(p::AbstractStockAndFlow)=begin
    fnames = [fname(p, f) for f in 1:nf(p)]
    LVector(;[(fnames[f]=>funcFlow(p, f)) for f in 1:nf(p)]...)
end


# generate the function of a sum auxiliary variable (index sv) with the sum of all stocks links to it
funcSV(p::AbstractStockAndFlow0,sv) = begin
    uN(u,t) = begin
        sumS = 0
        for i in stockssv(p,sv)
            sumS=sumS+u[sname(p,i)]
        end
        return sumS        
    end
    return uN       
end
# return the LVector of pairs: svname => function 
funcSVs(p::AbstractStockAndFlow0) = begin
  svnames = [svname(p, sv) for sv in 1:nsv(p)]
  LVector(;[(svnames[sv]=>funcSV(p, sv)) for sv in 1:nsv(p)]...)
end


# given a stock and flow diagram in schema "StockAndFlow", return a stock and flow diagram in schema "StockAndFlow0"
object_shift_right(p::StockAndFlow) = begin
    s = snames(p)
    initialvalues = map(x->initialValues(p)[x], s)
    ssv = map(y->map(x->(sname(p,y),x),svname(p,svsstock(p,y))),collect(1:ns(p)))
    ssv = vcat(ssv...)
    StockAndFlow0(s,ssv,initialvalues)
end

# create open acset, as the structured cospan
const OpenStockAndFlowOb, OpenStockAndFlow = OpenACSetTypes(StockAndFlow,StockAndFlow0)

foot(x::StockAndFlow, s) = begin
    s = vectorify(s)
    initialvalues = map(i->initialValues(x)[i], s)
    StockAndFlow0(s, initialvalues)
end

foot(x::StockAndFlow, s, ssv) = begin
    s = vectorify(s)
    ssv = vectorify(ssv)
    initialvalues = map(i->initialValues(x)[i], s)
    StockAndFlow0(s, ssv, initialvalues)
end

ntcomponent(a, x0) = map(x->state_dict(x0)[x], a)

leg(a::A, x0::A) where {A <: Union{StockAndFlow0, StockAndFlow}} = begin
    ϕs = ntcomponent(snames(a), snames(x0))
    if nsv(a) > 0  # if have sum-auxiliary-variable and links between stocks and sum-auxiliary-variables
      ϕsv = ntcomponent(svnames(a), svnames(x0))
      ϕls = ntcomponent(lsnames(a), lsnames(x0))
    else
      ϕsv = Int[]
      ϕls = Int[]
    end
    ACSetTransformation((S=ϕs, LS=ϕls, SV=ϕsv), a, x0)
end

#TODO: temporarily ignore the OpenACSetTypes waiting for the Catlab update for supporting open for C0
#const OpenStockAndFlowOb, OpenStockAndFlow = OpenACSetTypes(StockAndFlow,StockAndFlow0)

Open(p::StockAndFlow, feet...) = begin
  p0 = object_shift_right(p)
  legs = map(x->leg(x, p0), feet)
  OpenStockAndFlow(p, legs...)
end 

struct TransitionMatrices 
# row represent flows, column represent stocks; and element of 1 of matrix indicates whether there is a connection between the flow and stock; 0 indicates no connection
  inflow::Matrix{Int}
  outflow::Matrix{Int}
  TransitionMatrices(p::AbstractStockAndFlow) = begin
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
functioneat(f::Number, u, uN, t) = x
valueat(f::Function, u, t) = f(u,t)
functioneat(f::Function, u, uN, t)=f(u,uN,t)


# generate the function of ODEs with the correct format used in package "OrdinaryDiffEq"
vectorfield(pn::AbstractStockAndFlow) = begin
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

include("visualization.jl")
# The implementations in this file is specific for the Primitive schema of stock and flow diagram in the ACT paper
include("PrimitiveStockFlowInPaper.jl")

end










































