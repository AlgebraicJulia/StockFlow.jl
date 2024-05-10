module StockFlow

export TheoryStockAndFlow0, TheoryStockAndFlow, TheoryStockAndFlowStructure, TheoryStockAndFlowStructureF, TheoryStockAndFlowF, AbstractStockAndFlow0, AbstractStockAndFlow, AbstractStockAndFlowStructure,
AbstractStockAndFlowStructureF, AbstractStockAndFlowF, StockAndFlow0, StockAndFlow, StockAndFlowStructure, StockAndFlowStructureF, StockAndFlowF, add_flow!, add_flows!, add_stock!, add_stocks!,
add_variable!, add_variables!, add_svariable!, add_svariables!, add_parameter!, add_parameters!, add_VVlink!, add_VVlinks!,
add_inflow!, add_inflows!, add_outflow!, add_outflows!, add_Vlink!, add_Vlinks!, add_Slink!, add_Slinks!, add_SVlink!, add_Plink!, add_Plinks!,
add_SVlinks!, ns, nf, ni, no, nvb, nsv, nls, nlv, nlsv, np, nlpv, nlvv, sname, fname, svname, svnames, vname, inflows, outflows,
funcDynam, flowVariableIndex, funcFlow, funcFlows, funcSV, funcSVs, TransitionMatrices,
vectorfield, funcFlowsRaw, funcFlowRaw, inflowsAll, outflowsAll,instock,outstock, stockssv, stocksv, svsv, svsstock,
vsstock, vssv, svsstockAllF, vsstockAllF, vssvAllF, StockAndFlowUntyped, StockAndFlowFUntyped, StockAndFlowStructureUntyped, StockAndFlowStructureFUntyped, StockAndFlowUntyped0, Open, snames, fnames, svnames, vnames,
object_shift_right, foot, leg, lsnames, OpenStockAndFlow, OpenStockAndFlowOb, fv, fvs, nlvv, nlpv, vtgt, vsrc, vpsrc, vptgt, pname, pnames, make_v_expr,
vop, lvvposition, lvtgtposition, lsvvposition, lpvvposition, recreate_stratified, set_snames!, set_fnames!, set_svnames!, set_vnames!, set_pnames!, set_sname!, set_fname!, set_svname!, set_vname!, set_pname!,
get_lss, get_lssv, get_lsvsv, get_lsvv, get_lvs, get_lvv, get_is, get_ifn, get_os, get_ofn, get_lpvp, get_lpvv, get_lvsrc, get_lvtgt, get_links,
make_v_expr_nonrecursive, get_lpvpposition, get_lvsrcposition, get_lsvsvposition, get_lvsposition


using Catlab
using Catlab.CategoricalAlgebra
using Catlab.CategoricalAlgebra.FinSets
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
vectorify(n::SubArray) = collect(n)
vectorify(n) = [n]

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))

#= operators definition
# Operators:

# binary operators:

x + y binary plus performs addition
x - y binary minus  performs subtraction
x * y times performs multiplication
x / y divide  performs division
x ÷ y integer divide  x / y, truncated to an integer
x ^ y power raises x to the yth power
x % y remainder equivalent to rem(x,y)
log(b,x)  base b logarithm of x

# unary:

sqrt(x) square root of x
exp(x)  natural exponential function at x
log(x)  natural logarithm of x

can refer to:
https://docs.julialang.org/en/v1/manual/mathematical-operations/
=#

Operators = Dict(2 => [:+, :-, :*, :/, :÷, :^, :%, :log, Symbol("")],
                 1 => [:log, :exp, :sqrt, Symbol("")]) #:NN is for NONE, which is used in create the special diagram using graph rewriting


math_expr(op, op1, op2) = Expr(:call, op, op1, op2)
math_expr(op, op1) = Expr(:call, op, op1)
#math_expr(op) = Expr(:call, op)


""" define the sub-schema of c0, which includes the three objects: stocks(S), sum-auxiliary-variables(SV), and the linkages between them (LS) to be composed """
@present TheoryStockAndFlow0(FreeSchema) begin
  S::Ob
  SV::Ob
  LS::Ob

  lss::Hom(LS,S)
  lssv::Hom(LS,SV)

# Attributes:
  Name::AttrType

  sname::Attr(S, Name)
  svname::Attr(SV, Name)
end

@abstract_acset_type AbstractStockAndFlow0
@acset_type StockAndFlowUntyped0(TheoryStockAndFlow0, index=[:lss,:lssv]) <: AbstractStockAndFlow0
# constrains the attributes data type to be:
# 2. Name: Symbol
const StockAndFlow0 = StockAndFlowUntyped0{Symbol}

"""
    StockAndFlow0(s,sv,ssv)
    
for an instance of the sub-schema, the program supports only have stocks, or only have sum auxiliary variables, or both stocks
 and sum auxiliary variables, or have both """
StockAndFlow0(s,sv,ssv) = begin
  p0 = StockAndFlow0()
  s = vectorify(s)
  sv = vectorify(sv)
  ssv = vectorify(ssv)

  if length(s)>0
    s_idx=state_dict(s)
    add_stocks!(p0, length(s), sname=s)
  end

  if length(sv)>0
    sv_idx=state_dict(sv)
    add_svariables!(p0, length(sv), svname=sv)
  end

  if length(ssv)>0
    svl=map(last, ssv)
    sl=map(first, ssv)
    add_Slinks!(p0, length(ssv), map(x->s_idx[x], sl), map(x->sv_idx[x], svl))
  end
  p0
end

""" define the schema of a general stock and flow diagram, not includes the functions of auxiliary variables """
@present TheoryStockAndFlowStructure <: TheoryStockAndFlow0 begin
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
  fname::Attr(F, Name)
  vname::Attr(V, Name)
end

@abstract_acset_type AbstractStockAndFlowStructure <: AbstractStockAndFlow0
@acset_type StockAndFlowStructureUntyped(TheoryStockAndFlowStructure, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv]) <: AbstractStockAndFlowStructure
# constrains the attributes data type to be:
# 1. InitialValue: Real
# 2. Name: Symbol
# Note: those three (or any subgroups) attributes' datatype can be defined by the users. See the example of the PetriNet which allows the Reactionrate and Concentration defined by the users
const StockAndFlowStructure = StockAndFlowStructureUntyped{Symbol}


###### TODO #### delete??
# define the schema of a general stock and flow diagram
@present TheoryStockAndFlow <: TheoryStockAndFlowStructure begin
# Attributes:
  FuncDynam::AttrType
  funcDynam::Attr(V, FuncDynam)
end

@abstract_acset_type AbstractStockAndFlow <: AbstractStockAndFlowStructure
@acset_type StockAndFlowUntyped(TheoryStockAndFlow, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv]) <: AbstractStockAndFlow
const StockAndFlow = StockAndFlowUntyped{Symbol,Function}


""" define the schema of a general stock and flow diagram """
@present TheoryStockAndFlowStructureF <: TheoryStockAndFlowStructure begin
#Objects
  P::Ob
  LVV::Ob
  LPV::Ob
#Morphisms
  lvsrc::Hom(LVV,V)
  lvtgt::Hom(LVV,V)
  lpvp::Hom(LPV,P)
  lpvv::Hom(LPV,V)

# Attributes:
  pname::Attr(P, Name)
end

@abstract_acset_type AbstractStockAndFlowStructureF <: AbstractStockAndFlowStructure
@acset_type StockAndFlowStructureFUntyped(TheoryStockAndFlowStructureF, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp,:lpvv]) <: AbstractStockAndFlowStructureF
const StockAndFlowStructureF = StockAndFlowStructureFUntyped{Symbol}


""" define the schema of a general stock and flow diagram """
@present TheoryStockAndFlowF <: TheoryStockAndFlowStructureF begin

# Attributes:
  Op::AttrType # arithmetic operators
  Position::AttrType # argument's position of expressions

  vop::Attr(V, Op)

  lvsposition::Attr(LV, Position)
  lsvsvposition::Attr(LSV, Position)
  lvsrcposition::Attr(LVV, Position)
  lpvpposition::Attr(LPV, Position)
end

@abstract_acset_type AbstractStockAndFlowF <: AbstractStockAndFlowStructureF
@acset_type StockAndFlowFUntyped(TheoryStockAndFlowF, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp,:lpvv]) <: AbstractStockAndFlowF
const StockAndFlowF = StockAndFlowFUntyped{Symbol,Symbol,Int8}



@present TheoryStockAndFlowStructureC <: TheoryStockAndFlowStructureF begin
    C::Ob
    LC::Ob

    lcc::Hom(LC, C)
    lcv::Hom(LC, V)
    cname::Attr(C, Name)
end

@abstract_acset_type AbstractStockAndFlowStructureC <: AbstractStockAndFlowStructureF
@acset_type StockAndFlowStructureCUntyped(TheoryStockAndFlowStructureC, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp,:lpvv, :lvv, :lcv]) <: AbstractStockAndFlowStructureC
const StockAndFlowStructureC = StockAndFlowStructureCUntyped{Symbol}

@present TheoryStockAndFlowC <: TheoryStockAndFlowStructureC begin

# Attributes:
  Op::AttrType # arithmetic operators
  Position::AttrType # argument's position of expressions
  Constant::AttrType

  vop::Attr(V, Op)

  lvsposition::Attr(LV, Position)
  lsvsvposition::Attr(LSV, Position)
  lvsrcposition::Attr(LVV, Position)
  lpvpposition::Attr(LPV, Position)
  con::Attr(C, Constant)
end

@abstract_acset_type AbstractStockAndFlowC <: AbstractStockAndFlow
@acset_type StockAndFlowCUntyped(TheoryStockAndFlowC, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp,:lpvv, :lvv, :lcv]) <: AbstractStockAndFlowC
const StockAndFlowC = StockAndFlowCUntyped{Symbol, Symbol, Int8, Float64}


# functions of adding components of the model schema
add_flow!(p::AbstractStockAndFlowStructure,v;kw...) = add_part!(p,:F;fv=v,kw...)
add_flows!(p::AbstractStockAndFlowStructure,v,n;kw...) = add_parts!(p,:F,n;fv=v,kw...)

add_stock!(p::AbstractStockAndFlow0;kw...) = add_part!(p,:S;kw...)
add_stocks!(p::AbstractStockAndFlow0,n;kw...) = add_parts!(p,:S,n;kw...)

add_variable!(p::AbstractStockAndFlowStructure;kw...) = add_part!(p,:V;kw...)
add_variables!(p::AbstractStockAndFlowStructure,n;kw...) = add_parts!(p,:V,n;kw...)

add_svariable!(p::AbstractStockAndFlow0;kw...) = add_part!(p,:SV;kw...)
add_svariables!(p::AbstractStockAndFlow0,n;kw...) = add_parts!(p,:SV,n;kw...)

add_inflow!(p::AbstractStockAndFlowStructure,s,t;kw...) = add_part!(p,:I;is=s,ifn=t,kw...)
add_inflows!(p::AbstractStockAndFlowStructure,n,s,t;kw...) = add_parts!(p,:I,n;is=s,ifn=t,kw...)

add_outflow!(p::AbstractStockAndFlowStructure,s,t;kw...) = add_part!(p,:O;os=s,ofn=t,kw...)
add_outflows!(p::AbstractStockAndFlowStructure,n,s,t;kw...) = add_parts!(p,:O,n;ofn=t,os=s,kw...)

# links from Stock to dynamic variable
add_Vlink!(p::AbstractStockAndFlowStructure,s,v;kw...) = add_part!(p,:LV;lvs=s,lvv=v,kw...)
add_Vlinks!(p::AbstractStockAndFlowStructure,n,s,v;kw...) = add_parts!(p,:LV,n;lvs=s,lvv=v,kw...)

# links from Stock to sum dynamic variable
add_Slink!(p::AbstractStockAndFlow0,s,sv;kw...) = add_part!(p,:LS;lss=s,lssv=sv,kw...)
add_Slinks!(p::AbstractStockAndFlow0,n,s,sv;kw...) = add_parts!(p,:LS,n;lss=s,lssv=sv,kw...)

# links from sum dynamic variable to dynamic varibale
add_SVlink!(p::AbstractStockAndFlowStructure,sv,v;kw...) = add_part!(p,:LSV;lsvsv=sv,lsvv=v,kw...)
add_SVlinks!(p::AbstractStockAndFlowStructure,n,sv,v;kw...) = add_parts!(p,:LSV,n;lsvsv=sv,lsvv=v,kw...)

# return the count of each components
""" stock count """
ns(p::AbstractStockAndFlow0) = nparts(p,:S) #stocks
""" flow count """
nf(p::AbstractStockAndFlowStructure) = nparts(p,:F) #flows
""" inflow count """
ni(p::AbstractStockAndFlowStructure) = nparts(p,:I) #inflows
""" outflow count """
no(p::AbstractStockAndFlowStructure) = nparts(p,:O) #outflows
""" dynamic variable count """
nvb(p::AbstractStockAndFlowStructure) = nparts(p,:V) #auxiliary variables
""" sum variable count """
nsv(p::AbstractStockAndFlow0) = nparts(p,:SV) #sum auxiliary variables
""" link stock sum count """
nls(p::AbstractStockAndFlow0) = nparts(p,:LS) #links from Stock to sum dynamic variable
""" link stock dynamic variable count """
nlv(p::AbstractStockAndFlowStructure) = nparts(p,:LV) #links from Stock to dynamic variable
""" link sum variable dynamic variable count """
nlsv(p::AbstractStockAndFlowStructure) = nparts(p,:LSV) #links from sum dynamic variable to dynamic variable
""" link dynamic variable dynamic variable count """
nlvv(p::AbstractStockAndFlowStructureF) = nparts(p,:LVV) #links from dynamic variable to dynamic variable
""" link dynamic variable parameter count """
nlpv(p::AbstractStockAndFlowStructureF) = nparts(p,:LPV) #links from dynamic variable to parameters
""" parameter count """
np(p::AbstractStockAndFlowStructureF) = nparts(p,:P) #parameters


get_lss(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lss].m))
get_lssv(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lssv].m))
get_lsvsv(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lsvsv].m))
get_lsvv(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lsvv].m))
get_lvs(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lvs].m))
get_lvv(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lvv].m))
get_is(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:is].m))
get_ifn(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:ifn].m))
get_os(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:os].m))
get_ofn(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:ofn].m))
get_lpvp(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lpvp].m))
get_lpvv(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lpvv].m))
get_lvsrc(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lvsrc].m))
get_lvtgt(sf::AbstractStockAndFlowF) = collect(values(sf.subparts[:lvtgt].m))

get_vop(sf::AbstractStockAndFlowF) = subpart(sf, :vop)

get_lpvpposition(sf::AbstractStockAndFlowF) = subpart(sf, :lpvpposition)
get_lvsrcposition(sf::AbstractStockAndFlowF) = subpart(sf, :lvsrcposition)
get_lsvsvposition(sf::AbstractStockAndFlowF) = subpart(sf, :lsvsvposition)
get_lvsposition(sf::AbstractStockAndFlowF) = subpart(sf, :lvsposition)

#EXAMPLE:
#sir_StockAndFlow=StockAndFlow(((:S, 990)=>(:birth,(:inf,:deathS),(:v_inf,:v_deathS),:N), (:I, 10)=>(:inf,(:rec,:deathI),(:v_rec,:v_deathI,:v_fractionNonS),:N),(:R, 0)=>(:rec,:deathR,(:v_deathR,:v_fractionNonS),:N)),
#  (:birth=>:v_birth,:inf=>:v_inf,:rec=>:v_rec,:deathS=>:v_deathS,:deathI=>:v_deathI,:deathR=>:v_deathR),
#  (:N=>(:v_birth,:v_inf))))

StockAndFlowStructure(s,f,sv) = begin

    p = StockAndFlowStructure()

    s = vectorify(s)
    f = vectorify(f)
    sv = vectorify(sv)

    sname = map(first,s)
    fname = map(first, f)
    vname = map(last, f)
    svname = map(first, sv)
    s_idx = state_dict(sname)
    f_idx = state_dict(fname)
    v_idx = state_dict(vname)
    sv_idx = state_dict(svname)

    # adding the objects that do not have out-morphisms firstly
    add_variables!(p, length(vname), vname=vname)    # add objects :V (auxiliary variables)
    add_svariables!(p, length(svname), svname=svname)    # add objects :SV (sum auxiliary variables)
    add_flows!(p,map(x->v_idx[x], map(last,f)),length(fname),fname=fname)    # add objects :F (flows)

    # Parse the elements included in "s" -- stocks
    for (i, (name,(ins,outs,vs,svs))) in enumerate(s)
      i = add_stock!(p,sname=name) # add objects :S (stocks)
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

###### TODO #### delete??

#EXAMPLE: Add a tuple with (auxiliary variable => function)
#sir_StockAndFlow=StockAndFlow(((:S, 990)=>(:birth,(:inf,:deathS),(:v_inf,:v_deathS),:N), (:I, 10)=>(:inf,(:rec,:deathI),(:v_rec,:v_deathI,:v_fractionNonS),:N),(:R, 0)=>(:rec,:deathR,(:v_deathR,:v_fractionNonS),:N)),
#  (:birth=>:v_birth,:inf=>:v_inf,:rec=>:v_rec,:deathS=>:v_deathS,:deathI=>:v_deathI,:deathR=>:v_deathR),
#  (:v_birth=>f_birth,:v_inf=>f_inf,:v_rec=>f_rec,:v_deathS=>f_deathS,:v_deathI=>f_deathI,:v_deathR=>f_deathR),
#--  (:N=>(:v_birth,:v_inf))))
#  (:N)))

StockAndFlow(s,f,v,sv) = begin

    p = StockAndFlow()

    s = vectorify(s)
    f = vectorify(f)
    v = vectorify(v)
    sv = vectorify(sv)

    sname = map(first,s)
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
    for (i, (name,(ins,outs,vs,svs))) in enumerate(s)
      i = add_stock!(p,sname=name) # add objects :S (stocks)
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

add_parameter!(p::AbstractStockAndFlowStructureF;kw...) = add_part!(p,:P;kw...)
add_parameters!(p::AbstractStockAndFlowStructureF,n;kw...) = add_parts!(p,:P,n;kw...)

add_variable!(p::AbstractStockAndFlowStructureF;kw...) = add_part!(p,:V;kw...)
add_variables!(p::AbstractStockAndFlowStructureF,n;kw...) = add_parts!(p,:V,n;kw...)

add_VVlink!(p::AbstractStockAndFlowStructureF,vs,vt;kw...) = add_part!(p,:LVV;lvsrc=vs,lvtgt=vt,kw...)
add_VVlinks!(p::AbstractStockAndFlowStructureF,n,vs,vt;kw...) = add_parts!(p,:LVV,n;lvsrc=vs,lvtgt=vt,kw...)

add_Plink!(sf::AbstractStockAndFlowStructureF,p,v;kw...) = add_part!(sf,:LPV;lpvp=p,lpvv=v,kw...)
add_Plinks!(sf::AbstractStockAndFlowStructureF,n,p,v;kw...) = add_parts!(sf,:LPV,n;lpvp=p,lpvv=v,kw...)


# define the function to construct system structure diagram with function transparency: StockAndFlowStructureF
StockAndFlowStructureF(s,p,v,f,sv) = begin
  sf = StockAndFlowStructureF()

  s = vectorify(s)
  f = vectorify(f)
  v = vectorify(v)
  p = vectorify(p)
  sv = vectorify(sv)

  sname = map(first,s)
  fname = map(first,f)
  vname = map(first,v)

  s_idx = state_dict(sname)
  f_idx = state_dict(fname)
  v_idx = state_dict(vname)
  p_idx = state_dict(p)
  sv_idx = state_dict(sv)

  add_parameters!(sf,length(p),pname=p)
  add_svariables!(sf, length(sv), svname=sv)
  add_variables!(sf, length(vname), vname=vname)
  add_flows!(sf,map(x->v_idx[x], map(last,f)),length(fname),fname=fname)


  # Parse the elements included in "s" -- stocks
  for (i, (name,(ins,outs,svs))) in enumerate(s)
    i = add_stock!(sf,sname=name) # add objects :S (stocks)
    ins=vectorify(ins) # inflows of each stock
    outs=vectorify(outs) # outflows of each stock
    svs=vectorify(svs) # sum auxiliary variables depends on the stock
    # filter out the fake (empty) elements
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    svs = svs[svs .!= FK_SVARIABLE_NAME]

    if length(ins)>0
      add_inflows!(sf, length(ins), repeat([i], length(ins)), map(x->f_idx[x], ins)) # add objects :I (inflows)
    end
    if length(outs)>0
      add_outflows!(sf, length(outs), repeat([i], length(outs)), map(x->f_idx[x], outs)) # add objects :O (outflows)
    end
    if length(svs)>0
      add_Slinks!(sf, length(svs), repeat([i], length(svs)), map(x->sv_idx[x], svs)) # add objects :LS (links from Stock to sum dynamic variable)
    end
  end

  # Parse the elements included in "v" -- auxiliary vairables
  for (vn,args) in v

    args=vectorify(args)
#    @assert op in Operators[length(args)]

    for (i,arg) in enumerate(args)
      if arg in sname
        add_Vlink!(sf, s_idx[arg], v_idx[vn])
      elseif arg in vname
        add_VVlink!(sf, v_idx[arg], v_idx[vn])
      elseif arg in p
        add_Plink!(sf, p_idx[arg], v_idx[vn])
      elseif arg in sv
        add_SVlink!(sf, sv_idx[arg], v_idx[vn])
      else
        error("arg: " * string(arg) * " does not belong to any stock, auxilliary variable, constant parameter or sum auxilliary variable!")
      end
    end
  end
  sf
end

"""
    StockAndFlowF(s,p,v,f,sv)

EXAMPLE: when define the dynamical variables, need to define with the correct order
```julia
sir_StockAndFlow=StockAndFlowF((:S=>(:F_NONE,:inf,:N), :I=>(:inf,:F_NONE,:N)),
 (:c, :beta),
 (:v_prevalence=>(:I,:N,:/),:v_meanInfectiousContactsPerS=>(:c,:v_prevalence,:*),:v_perSIncidenceRate=>(:beta,:v_meanInfectiousContactsPerS,:*),:v_newInfetions=>(:S,:v_perSIncidenceRate,:*)),
 (:inf=>:v_newInfetions),
 (:N))
```
"""
StockAndFlowF(s,p,v,f,sv) = begin
  sf = StockAndFlowF()

  s = vectorify(s)
  f = vectorify(f)
  v = vectorify(v)
  p = vectorify(p)
  sv = vectorify(sv)

  sname = map(first,s)
  fname = map(first,f)
  vname = map(first,v)

  op=map(last,map(last,v))

  s_idx = state_dict(sname)
  f_idx = state_dict(fname)
  v_idx = state_dict(vname)
  p_idx = state_dict(p)
  sv_idx = state_dict(sv)

  add_parameters!(sf,length(p),pname=p)
  add_svariables!(sf, length(sv), svname=sv)
  add_variables!(sf, length(vname), vname=vname, vop=op)
  add_flows!(sf,map(x->v_idx[x], map(last,f)),length(fname),fname=fname)


  # Parse the elements included in "s" -- stocks
  for (i, (name,(ins,outs,svs))) in enumerate(s)
    i = add_stock!(sf,sname=name) # add objects :S (stocks)
    ins=vectorify(ins) # inflows of each stock
    outs=vectorify(outs) # outflows of each stock
    svs=vectorify(svs) # sum auxiliary variables depends on the stock
    # filter out the fake (empty) elements
    ins = ins[ins .!= FK_FLOW_NAME]
    outs = outs[outs .!= FK_FLOW_NAME]
    svs = svs[svs .!= FK_SVARIABLE_NAME]

    if length(ins)>0
      add_inflows!(sf, length(ins), repeat([i], length(ins)), map(x->f_idx[x], ins)) # add objects :I (inflows)
    end
    if length(outs)>0
      add_outflows!(sf, length(outs), repeat([i], length(outs)), map(x->f_idx[x], outs)) # add objects :O (outflows)
    end
    if length(svs)>0
      add_Slinks!(sf, length(svs), repeat([i], length(svs)), map(x->sv_idx[x], svs)) # add objects :LS (links from Stock to sum dynamic variable)
    end
  end

  # Parse the elements included in "v" -- auxiliary vairables
  for (vn,(args,op)) in v

    args=vectorify(args)
#    @assert op in Operators[length(args)]

    for (i,arg) in enumerate(args)
      if arg in sname
        add_Vlink!(sf, s_idx[arg], v_idx[vn], lvsposition=i)
      elseif arg in vname
        add_VVlink!(sf, v_idx[arg], v_idx[vn], lvsrcposition=i)
      elseif arg in p
        add_Plink!(sf, p_idx[arg], v_idx[vn], lpvpposition=i)
      elseif arg in sv
        add_SVlink!(sf, sv_idx[arg], v_idx[vn], lsvsvposition=i)
      else
        error("arg: " * string(arg) * " does not belong to any stock, auxilliary variable, constant parameter or sum auxilliary variable!")
      end
    end
  end
  sf
end

""" return the stocks name with index of s """
sname(p::AbstractStockAndFlow0,s) = subpart(p,s,:sname) 
""" return the flows name with index of f """
fname(p::AbstractStockAndFlowStructure,f) = subpart(p,f,:fname)
""" return the sum auxiliary variables name with index of sv """
svname(p::AbstractStockAndFlow0,sv) = subpart(p,sv,:svname) 
""" return the auxiliary variables name with index of v """
vname(p::AbstractStockAndFlowStructure,v) = subpart(p,v,:vname) 
""" return the auxiliary variables name with index of v """
pname(sf::AbstractStockAndFlowStructureF,p) = subpart(sf,p,:pname)

""" return stock names """
snames(p::AbstractStockAndFlow0) = [sname(p, s) for s in 1:ns(p)]
""" return flow names """
fnames(p::AbstractStockAndFlowStructure) = [fname(p, f) for f in 1:nf(p)]
""" return sum variable names """
svnames(p::AbstractStockAndFlow0) = [svname(p, sv) for sv in 1:nsv(p)]
""" return variable names """ 
vnames(p::AbstractStockAndFlowStructure) = [vname(p, v) for v in 1:nvb(p)]
""" return parameter names """
pnames(sf::AbstractStockAndFlowStructureF) = [pname(sf,p) for p in 1:np(sf)]

""" set stock names to vector of symbols"""
set_snames!(p::AbstractStockAndFlow0, names) = set_subpart!(p, :sname, names)
""" set flow names to vector of symbols"""
set_fnames!(p::AbstractStockAndFlowStructure, names) = set_subpart!(p, :fname, names)
""" set sum variable names to vector of symbols"""
set_svnames!(p::AbstractStockAndFlow0, names) = set_subpart!(p, :svname, names)
""" set dynamic variable names to vector of symbols"""
set_vnames!(p::AbstractStockAndFlowStructure, names) = set_subpart!(p, :vname, names)
""" set parameter names to vector of symbols"""
set_pnames!(p::AbstractStockAndFlowStructure, names) = set_subpart!(p, :pname, names)


# Could pull out the right side into another function
set_sname!(p::AbstractStockAndFlow0, index, newname) = set_subpart!(p, :sname, [i == index ? newname : prevnames for (i, prevnames) in enumerate(snames(p))])
set_fname!(p::AbstractStockAndFlowStructure, index, newname) = set_subpart!(p, :fname, [i == index ? newname : prevnames for (i, prevnames) in enumerate(fnames(p))])
set_svname!(p::AbstractStockAndFlow0, index, newname) = set_subpart!(p, :svname, [i == index ? newname : prevnames for (i, prevnames) in enumerate(svnames(p))])
set_vname!(p::AbstractStockAndFlowStructure, index, newname) = set_subpart!(p, :vname, [i == index ? newname : prevnames for (i, prevnames) in enumerate(vnames(p))])
set_pname!(p::AbstractStockAndFlowStructure, index, newname) = set_subpart!(p, :pname, [i == index ? newname : prevnames for (i, prevnames) in enumerate(pnames(p))])

set_vops!(p::AbstractStockAndFlowF, names) = set_subparts!(p, :vop, names)
set_lvvpositions!(p::AbstractStockAndFlowF, names) = set_subparts!(p, :lvvposition, names)
set_lvtgtpositions!(p::AbstractStockAndFlowF, names) = set_subparts!(p, :lvtgtposition, names)
set_lsvvpositions!(p::AbstractStockAndFlowF, names) = set_subparts!(p, :lsvvposition, names)
set_lpvvpositions!(p::AbstractStockAndFlowF, names) = set_subparts!(p, :lpvvposition, names)


set_vop!(p::AbstractStockAndFlowF, index, newname) = set_subpart!(p, :vop, [i == index ? newname : prevnames for (i, prevnames) in enumerate(vop(p))])
set_lvvposition!(p::AbstractStockAndFlowF, index, newname) = set_subpart!(p, :lvvposition, [i == index ? newname : prevnames for (i, prevnames) in enumerate(lvvposition(p))])
set_lvtgtposition!(p::AbstractStockAndFlowF, index, newname) = set_subpart!(p, :lvtgtposition, [i == index ? newname : prevnames for (i, prevnames) in enumerate(lvtgtposition(p))])
set_lsvvposition!(p::AbstractStockAndFlowF, index, newname) = set_subpart!(p, :lsvvposition, [i == index ? newname : prevnames for (i, prevnames) in enumerate(lsvvposition(p))])
set_lpvvposition!(p::AbstractStockAndFlowF, index, newname) = set_subpart!(p, :lpvvposition, [i == index ? newname : prevnames for (i, prevnames) in enumerate(lpvvposition(p))])





""" get flow variable of flow """
fv(p::AbstractStockAndFlowStructure,f) = subpart(p,f,:fv)
""" get ordered list of indices of dynamic variable for each flow """
fvs(p::AbstractStockAndFlowStructure)=[fv(p,f) for f in 1:nf(p)]

""" return the pair of names of (stock, sum-auxiliary-variable) for all linkages between them """
lsnames(p::AbstractStockAndFlow0) = begin
    s = map(x->subpart(p,x,:lss),collect(1:nls(p)))
    sv = map(x->subpart(p,x,:lssv),collect(1:nls(p)))
    sn = map(x->sname(p,x),s)
    svn = map(x->svname(p,x),sv)
    pssv = collect(zip(sn, svn))
end

""" return inflows of stock index s """
inflows(p::AbstractStockAndFlowStructure,s) = subpart(p,incident(p,s,:is),:ifn)
""" return outflows of stock index s """
outflows(p::AbstractStockAndFlowStructure,s) = subpart(p,incident(p,s,:os),:ofn)

# TODO: add assertion that the length(instock)=1
""" return stocks of flow index f flow in """
instock(p::AbstractStockAndFlowStructure,f) = subpart(p,incident(p,f,:ifn),:is)

# TODO: add assertion that the length(outstock)=1
""" return stocks that flow index f flow out from """
outstock(p::AbstractStockAndFlowStructure,f) = subpart(p,incident(p,f,:ofn),:os)
""" return stocks that sum variable index sv link to """
stockssv(p::AbstractStockAndFlow0,sv) = subpart(p,incident(p,sv,:lssv),:lss)
""" return stocks that auxiliary variable index v link to """
stocksv(p::AbstractStockAndFlowStructure,v) = subpart(p,incident(p,v,:lvv),:lvs)
""" return sum variables that auxiliary variable index v link to """
svsv(p::AbstractStockAndFlowStructure,v) = subpart(p,incident(p,v,:lsvv),:lsvsv)
""" return sum auxiliary variables that stock s links to """
svsstock(p::AbstractStockAndFlowStructure,s) = subpart(p,incident(p,s,:lss),:lssv)
""" return auxiliary variables that stock s links to """
vsstock(p::AbstractStockAndFlowStructure,s) = subpart(p,incident(p,s,:lvs),:lvv)
""" return auxiliary variables that sum auxiliary variable sv links to """
vssv(p::AbstractStockAndFlowStructure,sv) = subpart(p,incident(p,sv,:lsvsv),:lsvv)


""" return auxiliary variable's source auxiliary variables that it links to """
vtgt(p::AbstractStockAndFlowStructureF,v) = subpart(p,incident(p,v,:lvsrc),:lvtgt)
""" return auxiliary variable's target auxiliary variables it links to """
vsrc(p::AbstractStockAndFlowStructureF,v) = subpart(p,incident(p,v,:lvtgt),:lvsrc)
""" return auxiliary variable's source constant parameters """
vpsrc(p::AbstractStockAndFlowStructureF,v) = subpart(p,incident(p,v,:lpvv),:lpvp)
""" return constant parameters's link auxiliary variables """
vptgt(p::AbstractStockAndFlowStructureF,pp) = subpart(p,incident(p,pp,:lpvp),:lpvv)


""" return operator of an auxiliary variable """
vop(p::AbstractStockAndFlowF,v) = subpart(p,v,:vop)
""" return argument position of source stock variable of an auxiliary variable """
lvvposition(p::AbstractStockAndFlowF,v) = subpart(p,incident(p,v,:lvv),:lvsposition)
""" return argument position of source auxiliary variable of an auxiliary variable """
lvtgtposition(p::AbstractStockAndFlowF,v) = subpart(p,incident(p,v,:lvtgt),:lvsrcposition)
""" return argument position of source sum dynamical variable of an auxiliary variable """
lsvvposition(p::AbstractStockAndFlowF,v) = subpart(p,incident(p,v,:lsvv),:lsvsvposition)
""" return argument position of source constant parameter of an auxiliary variable """
lpvvposition(p::AbstractStockAndFlowF,v) = subpart(p,incident(p,v,:lpvv),:lpvpposition)




# create a dictionary
make_dict(ks, vs) = begin
    @assert length(ks)==length(vs)
    dic=()
    for (k,v) in zip(ks,vs)
        dic=(dic...,(k,v))
    end
    return Dict(dic)
end
""" create expresision of an auxiliary variable v """
function make_v_expr(p::AbstractStockAndFlowF,v)

    op = vop(p,v)
    srcsv=map(i->sname(p,i),stocksv(p,v))
    srcsvv=map(i->svname(p,i),svsv(p,v))
    srcpv=map(i->pname(p,i),vpsrc(p,v))
    srcvv=map(i->vname(p,i),vsrc(p,v))

    lvvp=lvvposition(p,v)
    lvtgtp=lvtgtposition(p,v)
    lsvvp=lsvvposition(p,v)
    lpvvp=lpvvposition(p,v)

    if length(srcvv)>0
        srcvv=map(x->make_v_expr(p,vsrc(p,v)[x]),1:length(vsrc(p,v)))
    end

    # create dictionary of (key=position, value=symbole of source argument)
    position_src=merge(make_dict(lvvp,srcsv),make_dict(lsvvp,srcsvv),make_dict(lpvvp,srcpv),make_dict(lvtgtp,srcvv))
    ordered_position_src=sort(collect(position_src), by = x->x[1])
    srcs=map(x->last(x),ordered_position_src)

    return math_expr(op,srcs...)
end



"""
Generates the expression of an auxiliary variable, only going one layer deep
If other dynamic variables make up its definition, they will not be substituted.
"""
function make_v_expr_nonrecursive(p::AbstractStockAndFlowF,v)

  op = vop(p,v)
  srcsv=map(i->sname(p,i),stocksv(p,v))
  srcsvv=map(i->svname(p,i),svsv(p,v))
  srcpv=map(i->pname(p,i),vpsrc(p,v))
  srcvv=map(i->vname(p,i),vsrc(p,v))

  lvvp=lvvposition(p,v)
  lvtgtp=lvtgtposition(p,v)
  lsvvp=lsvvposition(p,v)
  lpvvp=lpvvposition(p,v)

  # create dictionary of (key=position, value=symbole of source argument)
  position_src=merge(make_dict(lvvp,srcsv),make_dict(lsvvp,srcsvv),make_dict(lpvvp,srcpv),make_dict(lvtgtp,srcvv))
  ordered_position_src=sort(collect(position_src), by = x->x[1])
  srcs=map(x->last(x),ordered_position_src)

  return math_expr(op,srcs...)
end





# generate an array of all arguments of an expression
generate_expr_args(expr) = begin
    args=expr.args
    argsarray=[]
    for arg in args
        if arg isa Expr
            argsarray=vcat(argsarray,generate_expr_args(arg))
        else
            argsarray=vcat(argsarray,arg)
        end
    end
    ops=vcat(collect(values(Operators))...)
    return setdiff(unique(argsarray),ops)
end
# evaluate an expression to a function
eval_function(expr,s,sv,p,us,uNs,ps) = begin
    f=eval(Expr(:->, Expr(:tuple, s..., sv..., p...), Expr(:block,:(()),expr)))
    return @eval $f($(us...), $(uNs...), $(ps...))
end


""" return sum auxiliary variables all stocks link (frequency) """
svsstockAllF(p::AbstractStockAndFlowStructure) = [((svsstock(p, s) for s in 1:ns(p))...)...]
""" return auxiliary variables all stocks link (frequency) """
vsstockAllF(p::AbstractStockAndFlowStructure) = [((vsstock(p, s) for s in 1:ns(p))...)...]
""" return auxiliary variables all sum auxiliary variables link (frequency) """
vssvAllF(p::AbstractStockAndFlowStructure) = [((vssv(p, sv) for sv in 1:nsv(p))...)...]

""" return all inflows """
inflowsAll(p::AbstractStockAndFlowStructure) = [((inflows(p, s) for s in 1:ns(p))...)...]
""" return all outflows """
outflowsAll(p::AbstractStockAndFlowStructure) = [((outflows(p, s) for s in 1:ns(p))...)...]



""" 
    funcDynam(p::AbstractStockAndFlow,v)
return the functions of variables give index v """
funcDynam(p::AbstractStockAndFlow,v) = subpart(p,v,:funcDynam)

""" 
    funcDynam(sf::AbstractStockAndFlowF,v)
return the functions of variables give index v """
funcDynam(sf::AbstractStockAndFlowF,v) = begin
    expr=make_v_expr(sf,v)
    args=generate_expr_args(expr)

    args_s=args[findall(in(snames(sf)),args)]
    args_sv=args[findall(in(svnames(sf)),args)]
    args_p=args[findall(in(pnames(sf)),args)]

    f(u,uN,p,t)=begin
        us=map(i->u[i],args_s)
        uNs=map(i->uN[i](u,t),args_sv)
        ps=map(i->p[i],args_p)
        return eval_function(expr,args_s,args_sv,args_p,us,uNs,ps)
    end
    return f
end

""" return the auxiliary variable's index that related to the flow with index of f """
flowVariableIndex(p::AbstractStockAndFlowStructure,f) = subpart(p,f,:fv)
""" return the functions (not substitutes the function of sum variables yet) of flow index f """
funcFlowRaw(p::Union{AbstractStockAndFlow,AbstractStockAndFlowF},f)=funcDynam(p,flowVariableIndex(p,f))
""" return the LVector of pairs: fname => function (raw: not substitutes the function of sum variables yet) """
funcFlowsRaw(p::Union{AbstractStockAndFlow,AbstractStockAndFlowF}) = begin
  fnames = [fname(p, f) for f in 1:nf(p)]
  LVector(;[(fnames[f]=>funcFlowRaw(p, f)) for f in 1:nf(p)]...)
end
""" generate the function substituting sum variables in with flow index fn """
funcFlow(pn::Union{AbstractStockAndFlow,AbstractStockAndFlowF}, fn) = begin
    func=funcFlowRaw(pn,fn)
    f(u,p,t) = begin
        uN=funcSVs(pn)
        return valueat(func,u,uN,p,t)
    end
end
""" return the LVector of pairs: fname => function (with function of sum variables substitue in) """
funcFlows(p::Union{AbstractStockAndFlow,AbstractStockAndFlowF})=begin
    fnames = [fname(p, f) for f in 1:nf(p)]
    LVector(;[(fnames[f]=>funcFlow(p, f)) for f in 1:nf(p)]...)
end


""" generate the function of a sum auxiliary variable (index sv) with the sum of all stocks links to it """
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
""" return the LVector of pairs: svname => function """
funcSVs(p::AbstractStockAndFlow0) = begin
  svnames = [svname(p, sv) for sv in 1:nsv(p)]
  LVector(;[(svnames[sv]=>funcSV(p, sv)) for sv in 1:nsv(p)]...)
end


####################### functions of composition of stock and flow diagram####################
# Note: the compositional functions currently are only implemented for the StockFlow diagram (including both structure and functions)

# given a stock and flow diagram in schema "StockAndFlow", return a stock and flow diagram in schema "StockAndFlow0"
object_shift_right(p::StockAndFlowStructure) = begin
    s = snames(p)
    sv = svnames(p)
    ssv = map(y->map(x->(sname(p,y),x),svname(p,svsstock(p,y))),collect(1:ns(p)))
    ssv = vcat(ssv...)
    StockAndFlow0(s,sv,ssv)
end

# create open acset, as the structured cospan
const OpenStockAndFlowStructureOb, OpenStockAndFlowStructure = OpenACSetTypes(StockAndFlowStructureUntyped,StockAndFlowUntyped0)
const OpenStockAndFlowOb, OpenStockAndFlow = OpenACSetTypes(StockAndFlowUntyped,StockAndFlowUntyped0)
const OpenStockAndFlowStructureFOb, OpenStockAndFlowStructureF = OpenACSetTypes(StockAndFlowStructureFUntyped,StockAndFlowUntyped0)
const OpenStockAndFlowFOb, OpenStockAndFlowF = OpenACSetTypes(StockAndFlowFUntyped,StockAndFlowUntyped0)


foot(s, sv, ssv) = StockAndFlow0(s, sv, ssv)

ntcomponent(a, x0) = map(x->state_dict(x0)[x], a)

leg(a::StockAndFlow0, x::Union{StockAndFlowStructure,StockAndFlow,StockAndFlowStructureF,StockAndFlowF}) = begin
    if ns(a)>0 # if have stocks
      ϕs = ntcomponent(snames(a), snames(x))
    else
      ϕs = Int[]
    end

    if nsv(a) > 0  # if have sum-auxiliary-variable
      ϕsv = ntcomponent(svnames(a), svnames(x))
    else
      ϕsv = Int[]
    end

    if nls(a)>0 # if have links between stocks and sum-auxiliary-variables
      ϕls = ntcomponent(lsnames(a), lsnames(x))
    else
      ϕls = Int[]
    end

    result = OpenACSetLeg(a, S=ϕs, LS=ϕls, SV=ϕsv)

    result
end

Open(p::StockAndFlow, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenStockAndFlow{Symbol,Function}(p, legs...)
end

Open(p::StockAndFlowF, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenStockAndFlowF{Symbol,Symbol,Int8}(p, legs...)
end

Open(p::StockAndFlowStructure, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenStockAndFlowStructure{Symbol}(p, legs...)
end


############# functions of generating ODEs ###############

struct TransitionMatrices
# row represent flows, column represent stocks; and element of 1 of matrix indicates whether there is a connection between the flow and stock; 0 indicates no connection
  inflow::Matrix{Int}
  outflow::Matrix{Int}
  TransitionMatrices(p::AbstractStockAndFlowStructure) = begin
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

valueat(x::Number, u, p, t) = x
valueat(f::Number, u, uN, p, t) = x
valueat(f::Function, u, p, t) = f(u,p,t)
valueat(f::Function, u, uN, p, t)=f(u,uN,p,t)


# test argumenterror -- stocks in function of flow "fn" are not linked!
# TODO: find method to generate the exact wrong stocks' names and output in error message
ftest(f::Function, u, p, fn) =
try
  f(u,p,0)
catch e
  if isa(e, ArgumentError)
    println(string("Stocks used in the function of flow ", fn, " are not linked but used!"))
    rethrow(e)
  end
end

# if the function f runs to the end, then throw an ErrorException error!
ferror(f::Function, u, p, fn, umissed) = begin
  f(u,p,0)
  throw(ErrorException(string("stocks ", umissed, " in the function of flow", fn, " are linked but not used!")))
end

# test stocks in function of flow "fn" are missed!
fmisstest(f::Function, u, p, fn, umissed) =
try
  ferror(f,u,p,fn, umissed)
catch e
  if isa(e, ErrorException)
    rethrow(e)
  end
end

# the parameters are the labelled vector of the functions of flows
#vectorfield(pn::AbstractStockAndFlowStructure) = begin
#  tm = TransitionMatrices(pn)
#  f(du,u,p,t) = begin
#    u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
#    ϕ_m = [p[fname(pn, i)] for i in 1:nf(pn)]
#    for i in 1:ns(pn)
#      du[sname(pn, i)] = 0
#      for j in 1:nf(pn)
#        if tm.inflow[j,i] == 1
#          du[sname(pn, i)] = du[sname(pn, i)] + valueat(ϕ_m[j],u,p,t)
#        end
#        if tm.outflow[j,i] == 1
#          du[sname(pn, i)] = du[sname(pn, i)] - valueat(ϕ_m[j],u,p,t)
#        end
#      end
#    end
#    return du
#  end
#  return f
#end


vectorfield(pn::Union{AbstractStockAndFlow,AbstractStockAndFlowF}) = begin
  ϕ=funcFlows(pn)
  tm = TransitionMatrices(pn)
  f(du,u,p,t) = begin
    u_m = [u[sname(pn, i)] for i in 1:ns(pn)]
    ϕ_m = [ϕ[fname(pn, i)] for i in 1:nf(pn)]
    for i in 1:ns(pn)
      du[sname(pn, i)] = 0
      for j in 1:nf(pn)
        if tm.inflow[j,i] == 1
          du[sname(pn, i)] = du[sname(pn, i)] + valueat(ϕ_m[j],u,p,t)
        end
        if tm.outflow[j,i] == 1
          du[sname(pn, i)] = du[sname(pn, i)] - valueat(ϕ_m[j],u,p,t)
        end
      end
    end
    return du
  end
  return f
end

include("Syntax.jl")

include("CausalLoop.jl")
include("SystemStructure.jl")

include("visualization.jl")
# The implementations in this file is specific for the Primitive schema of stock and flow diagram in the ACT paper
include("PrimitiveStockFlowInPaper.jl")

include("PremadeModels.jl")

end
