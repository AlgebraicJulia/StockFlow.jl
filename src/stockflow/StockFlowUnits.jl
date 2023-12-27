module StockFlowUnits

using ..StockFlow
using Catlab.GATs.Presentations, Catlab.CategoricalAlgebra
using MLStyle

import ...StockFlow: StockAndFlowF, state_dict, ns, np, vectorify, ntcomponent

export add_unit!, add_units!, add_cunit!, add_cunits!, add_UCUlink!, add_UCUlinks!,
uname, unames, set_unames!, set_exps!, convert, StockAndFlowU, set_scu!, set_pcu!,
FtoU, StockAndFlow0U, footU, Open


@present TheoryStockAndFlow0U <: TheoryStockAndFlow0 begin
    Exponent::AttrType

    U::Ob
    CU::Ob
    LU::Ob


    uname::Attr(U, Name)
    exp::Attr(LU, Exponent)
    cuname::Attr(CU, Name)


    luu::Hom(LU, U)
    lucu::Hom(LU, CU)

    scu::Hom(S, CU)

end
@acset_type StockAndFlowUntyped0U(TheoryStockAndFlow0U, index=[:lss,:lssv, :luu, :lucu, :scu]) <: AbstractStockAndFlow0
const StockAndFlow0U = StockAndFlowUntyped0U{Symbol, Float64}



nu(p::StockAndFlow0U) = nparts(p, :U)
ncu(p::StockAndFlow0U) = nparts(p, :CU)
nlu(p::Union{StockAndFlowU, StockAndFlow0U}) = nparts(p, :LU)


add_unit!(p::StockAndFlow0U;kw...) = add_part!(p,:U;kw...)
add_units!(p::StockAndFlow0U,n;kw...) = add_parts!(p,:U,n;kw...)

add_cunit!(p::StockAndFlow0U;kw...) = add_part!(p,:CU;kw...)
add_cunits!(p::StockAndFlow0U,n;kw...) = add_parts!(p,:CU,n;kw...)

add_UCUlink!(p::StockAndFlow0U,u,cu;kw...) = add_part!(p,:LU;luu=u,lucu=cu,kw...)
add_UCUlinks!(p::StockAndFlow0U,n,u,cu;kw...) = add_parts!(p,:LU,n;luu=u,lucu=cu, kw...)

uname(p::StockAndFlow0U,u) = subpart(p,u,:uname) 
unames(p::StockAndFlow0U) = subpart(p,:uname) 

set_unames!(p::StockAndFlow0U, names) = set_subpart!(p, :uname, names)
set_cunames!(p::StockAndFlow0U, names) = set_subpart!(p, :cuname, names)
set_exps!(p::StockAndFlow0U, exps) = set_subpart!(p, :exp, exps)


set_scu!(p::StockAndFlow0U, cu) = set_subpart!(p, :scu, cu)


function extract_exponents(exp)
    @match exp begin
        ::Symbol => 
            Dict(exp => 1)
        :(1/$B) => Dict(k => -v for (k, v) in extract_exponents(B))
        Expr(:call, :*, A, B...) =>
            merge!(+, extract_exponents(A), [extract_exponents(b) for b ∈ B]...)
        Expr(:call, :/, A, B...) =>
            merge!(+, extract_exponents(A), [Dict(k => -v for (k, v) in extract_exponents(b)) for b ∈ B]...)
        :($A^$n) =>
            Dict(A => n)
    end
end

StockAndFlow0U(s,sv,ssv,cu,u) = begin
    p0 = StockAndFlow0U()
    s = vectorify(s)
    sv = vectorify(sv)
    ssv = vectorify(ssv)
    cu = vectorify(cu)
    u = vectorify(u)
  


    cu_vec = collect(Set(map(last, cu)))
    cu_idx = state_dict(cu_vec)


    exps = Dict(cu => extract_exponents(cu) for cu in filter(x -> x != :NONE, keys(cu_idx)))

    u_vec = unique([u ; [unit for inner_dict in values(exps) for unit in keys(inner_dict)]])
    u_idx = state_dict(u_vec)

    stock_cu = Dict(cu)

    add_cunits!(p0, length(cu_vec), cuname=cu_vec)
    add_units!(p0, length(u_vec), uname=u_vec)
      
    for cunit in keys(exps)
        for (unit, power) in exps[cunit]
            
            add_UCUlink!(p0, u_idx[unit], cu_idx[cunit], exp=power)
        end
    end
  
    if length(s)>0
      s_idx=state_dict(s)
      add_stocks!(p0, length(s), sname=s, scu=map(x -> cu_idx[stock_cu[x]] , s))
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

footU(s,sv,ssv,cu,u) = StockAndFlow0U(s,sv,ssv,cu,u)



@present TheoryStockAndFlowU <: TheoryStockAndFlowF begin

  Exponent::AttrType

  U::Ob
  CU::Ob
  LU::Ob


  uname::Attr(U, Name)
  exp::Attr(LU, Exponent)
  cuname::Attr(CU, Name)
  

  luu::Hom(LU, U)
  lucu::Hom(LU, CU)

  scu::Hom(S, CU)
  pcu::Hom(P, CU)


end

@acset_type StockAndFlowUUntyped(TheoryStockAndFlowU, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp, :lpvv, :luu, :lucu, :scu, :pcu]) <: AbstractStockAndFlowF
const StockAndFlowU = StockAndFlowUUntyped{Symbol,Symbol,Int8,Float64}

add_unit!(p::StockAndFlowU;kw...) = add_part!(p,:U;kw...)
add_units!(p::StockAndFlowU,n;kw...) = add_parts!(p,:U,n;kw...)

add_cunit!(p::StockAndFlowU;kw...) = add_part!(p,:CU;kw...)
add_cunits!(p::StockAndFlowU,n;kw...) = add_parts!(p,:CU,n;kw...)

add_UCUlink!(p::StockAndFlowU,u,cu;kw...) = add_part!(p,:LU;luu=u,lucu=cu,kw...)
add_UCUlinks!(p::StockAndFlowU,n,u,cu;kw...) = add_parts!(p,:LU,n;luu=u,lucu=cu, kw...)

uname(p::StockAndFlowU,u) = subpart(p,u,:uname) 
unames(p::StockAndFlowU) = subpart(p,:uname) 

cuname(p::Union{StockAndFlowU, StockAndFlow0U}, u) = subpart(p,u, :cuname) 
cunames(p::Union{StockAndFlowU, StockAndFlow0U}) = subpart(p,:cuname) 


set_unames!(p::StockAndFlowU, names) = set_subpart!(p, :uname, names)
set_cunames!(p::StockAndFlowU, names) = set_subpart!(p, :cuname, names)
set_exps!(p::StockAndFlowU, exps) = set_subpart!(p, :exp, exps)


set_scu!(p::StockAndFlowU, cu) = set_subpart!(p, :scu, cu)
set_pcu!(p::StockAndFlowU, cu) = set_subpart!(p, :pcu, cu)




lunames(p::Union{StockAndFlowU, StockAndFlow0U}) = begin
    luu = map(x->subpart(p,x,:luu),collect(1:nlu(p)))
    lucu = map(x->subpart(p,x,:lucu),collect(1:nlu(p)))
    unit_names = map(x->uname(p,x),luu)
    cunit_names = map(x->cuname(p,x),lucu)
    plu = collect(zip(unit_names, cunit_names))
end






function convert(::Type{StockAndFlowU}, sff::K) where {K <: AbstractStockAndFlowF}
    sfg = StockAndFlowU()
    
    add_cunit!(sfg, cuname=:NONE)
    
    add_stocks!(sfg, ns(sff), sname=snames(sff), scu=ones(Int, ns(sff)))
    add_svariables!(sfg, nsv(sff), svname=svnames(sff))
    add_variables!(sfg, nvb(sff), vname=vnames(sff))
    add_flows!(sfg, (fv(sff, i) for i in 1:nf(sff)), nf(sff), fname=fnames(sff))
    add_parameters!(sfg, np(sff), pname=pnames(sff), pcu=ones(Int, np(sff)))
    
    add_inflows!(sfg, ni(sff), subpart(sff, :is), subpart(sff, :ifn))
    add_outflows!(sfg, no(sff), subpart(sff, :os), subpart(sff, :ofn))
    add_Vlinks!(sfg, nlv(sff), subpart(sff, :lvs), subpart(sff, :lvv))
    add_Slinks!(sfg, nls(sff), subpart(sff, :lss), subpart(sff, :lssv))
    add_SVlinks!(sfg, nsv(sff), subpart(sff, :lsvsv), subpart(sff, :lsvv))
    add_VVlinks!(sfg, nlvv(sff), subpart(sff, :lvsrc), subpart(sff, :lvtgt))
    add_Plinks!(sfg, nlpv(sff), subpart(sff, :lpvp), subpart(sff, :lpvv))
    

    vop = subpart(sff, :vop)
    set_subpart!(sfg, :vop, vop)

    lvv = subpart(sff, :lvsposition)
    set_subpart!(sfg, :lvsposition, lvv)


    lvtgt = subpart(sff, :lvsrcposition)
    set_subpart!(sfg, :lvsrcposition, lvtgt)


    lsvv = subpart(sff, :lsvsvposition)
    set_subpart!(sfg, :lsvsvposition, lsvv)


    lpvv = subpart(sff, :lpvpposition)
    set_subpart!(sfg, :lpvpposition, lpvv)

    
    return sfg

end

StockAndFlowU(sff::K) where {K <: AbstractStockAndFlowF} = convert(StockAndFlowU, sff)


function FtoU(sff::K, sunits, punits) where {K <: AbstractStockAndFlowF}
    sfg = StockAndFlowU()
    
    

    cunits = collect(setdiff(Set(sunits) ∪ Set(punits)))
    cunit_dict = state_dict(cunits)

      

    cunit_exp_dict = Dict{Union{Symbol, Expr}, Dict{Symbol, Float64}}()
    unit_set = Set{Symbol}()
    
    for cu in cunits
        if cu == :NONE
            continue
        end
        cunit_exp_dict[cu] = extract_exponents(cu)
        union!(unit_set, keys(cunit_exp_dict[cu]))
    end

    
    units = collect(unit_set)
    
    unit_dict = state_dict(units)
    

    
    add_units!(sfg, length(units), uname=([Symbol(x) for x in units]))
    add_cunits!(sfg, length(cunits), cuname=([Symbol(x) for x in cunits]))
    
    for cu in cunits
        if cu == :NONE
            continue
        end
        for (u, exp) in pairs(cunit_exp_dict[cu])
            add_UCUlink!(sfg, unit_dict[u], cunit_dict[cu], exp=exp)
        end
    end
    


    scu = collect(map(x -> cunit_dict[x], sunits))
    pcu = collect(map(x -> cunit_dict[x], punits))
   
    
    add_stocks!(sfg, ns(sff), sname=snames(sff), scu=scu)
    add_svariables!(sfg, nsv(sff), svname=svnames(sff))
    add_variables!(sfg, nvb(sff), vname=vnames(sff))
    add_flows!(sfg, (fv(sff, i) for i in 1:nf(sff)), nf(sff), fname=fnames(sff))
    add_parameters!(sfg, np(sff), pname=pnames(sff), pcu=pcu)
    
    add_inflows!(sfg, ni(sff), subpart(sff, :is), subpart(sff, :ifn))
    add_outflows!(sfg, no(sff), subpart(sff, :os), subpart(sff, :ofn))
    add_Vlinks!(sfg, nlv(sff), subpart(sff, :lvs), subpart(sff, :lvv))
    add_Slinks!(sfg, nls(sff), subpart(sff, :lss), subpart(sff, :lssv))
    add_SVlinks!(sfg, nsv(sff), subpart(sff, :lsvsv), subpart(sff, :lsvv))
    add_VVlinks!(sfg, nlvv(sff), subpart(sff, :lvsrc), subpart(sff, :lvtgt))
    add_Plinks!(sfg, nlpv(sff), subpart(sff, :lpvp), subpart(sff, :lpvv))
    

    vop = subpart(sff, :vop)
    set_subpart!(sfg, :vop, vop)

    lvv = subpart(sff, :lvsposition)
    set_subpart!(sfg, :lvsposition, lvv)


    lvtgt = subpart(sff, :lvsrcposition)
    set_subpart!(sfg, :lvsrcposition, lvtgt)


    lsvv = subpart(sff, :lsvsvposition)
    set_subpart!(sfg, :lsvsvposition, lsvv)


    lpvv = subpart(sff, :lpvpposition)
    set_subpart!(sfg, :lpvpposition, lpvv)

    
    return sfg
end






# const OpenStockAndFlowStructureUOb, OpenStockAndFlowStructureU = OpenACSetTypes(StockAndFlowUUntyped,StockAndFlowUntyped0U)
const OpenStockAndFlowUOb, OpenStockAndFlowU = OpenACSetTypes(StockAndFlowUUntyped,StockAndFlowUntyped0U)



leg(a::StockAndFlow0U, x::StockAndFlowU) = begin
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

    if nu(a)>0 # if have links between stocks and sum-auxiliary-variables
        ϕu = ntcomponent(unames(a), unames(x))
      else
        ϕu = Int[]
      end

      if ncu(a)>0 # if have links between stocks and sum-auxiliary-variables
        ϕcu = ntcomponent(cunames(a), cunames(x))
      else
        ϕcu = Int[]
      end

      if nlu(a)>0 # if have links between stocks and sum-auxiliary-variables        
        ϕlu = ntcomponent(lunames(a), lunames(x))
      else
        ϕlu = Int[]
      end


    result = OpenACSetLeg(a, S=ϕs, LS=ϕls, SV=ϕsv, U=ϕu, CU=ϕcu, LU=ϕlu)

    result
end


Open(p::StockAndFlowU, feet...) = begin
    legs = map(x->leg(x, p), feet)
    OpenStockAndFlowU{Symbol,Symbol,Int8,Float64}(p, legs...)
  end





end
