module StockFlowUnits

using ..StockFlow
using Catlab.GATs, Catlab.CategoricalAlgebra
using MLStyle

import ..StockFlow: StockAndFlowF, state_dict, ns, np, vectorify, ntcomponent,
FK_FLOW_NAME, FK_VARIABLE_NAME, FK_SVARIABLE_NAME, FK_SVVARIABLE_NAME

export add_unit!, add_units!, add_dunit!, add_dunits!, add_UDUlink!, add_UDUlinks!,
uname, unames, duname, dunames, set_unames!, set_exps!, convert, StockAndFlowU, set_sdu!, set_pdu!,
FtoU, StockAndFlow0U, footU, Open, get_sdu, get_pdu, extract_exponents


@present TheoryStockAndFlow0U <: TheoryStockAndFlow0 begin
    Exponent::AttrType

    U::Ob
    DU::Ob
    LU::Ob


    uname::Attr(U, Name)
    exp::Attr(LU, Exponent)
    duname::Attr(DU, Name)


    luu::Hom(LU, U)
    ludu::Hom(LU, DU)

    sdu::Hom(S, DU)

end
@acset_type StockAndFlowUntyped0U(TheoryStockAndFlow0U, index=[:lss,:lssv, :luu, :ludu, :sdu]) <: AbstractStockAndFlow0
const StockAndFlow0U = StockAndFlowUntyped0U{Symbol, Float64}



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



StockAndFlow0U(s,sv,ssv,du,u) = begin
    p0 = StockAndFlow0U()
    s = vectorify(s)
    sv = vectorify(sv)
    ssv = vectorify(ssv)
    du = vectorify(du)
    u = vectorify(u)
  


    du_vec = collect(Set(map(last, du)))
    du_idx = state_dict(du_vec)


    exps = Dict(du => extract_exponents(du) for du in filter(x -> x != :NONE, keys(du_idx)))

    u_vec = unique([u ; [unit for inner_dict in values(exps) for unit in keys(inner_dict)]])
    u_idx = state_dict(u_vec)

    stock_du = Dict(du)

    add_dunits!(p0, length(du_vec), duname=du_vec)
    add_units!(p0, length(u_vec), uname=u_vec)
      
    for dunit in keys(exps)
        for (unit, power) in exps[dunit]
            
            add_UDUlink!(p0, u_idx[unit], du_idx[dunit], exp=power)
        end
    end
  
    if length(s)>0
      s_idx=state_dict(s)
      add_stocks!(p0, length(s), sname=s, sdu=map(x -> du_idx[stock_du[x]] , s))
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

footU(s,sv,ssv,du,u) = StockAndFlow0U(s,sv,ssv,du,u)



@present TheoryStockAndFlowU <: TheoryStockAndFlowF begin

  Exponent::AttrType

  U::Ob
  DU::Ob
  LU::Ob


  uname::Attr(U, Name)
  exp::Attr(LU, Exponent)
  duname::Attr(DU, Name)
  

  luu::Hom(LU, U)
  ludu::Hom(LU, DU)

  sdu::Hom(S, DU)
  pdu::Hom(P, DU)


end

@acset_type StockAndFlowUUntyped(TheoryStockAndFlowU, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp, :lpvv, :luu, :ludu, :sdu, :pdu]) <: AbstractStockAndFlowF
const StockAndFlowU = StockAndFlowUUntyped{Symbol,Symbol,Int8,Float64}




nu(p::StockAndFlow0U) = nparts(p, :U)
ndu(p::StockAndFlow0U) = nparts(p, :DU)
nlu(p::Union{StockAndFlowU, StockAndFlow0U}) = nparts(p, :LU)


add_unit!(p::StockAndFlow0U;kw...) = add_part!(p,:U;kw...)
add_units!(p::StockAndFlow0U,n;kw...) = add_parts!(p,:U,n;kw...)

add_dunit!(p::StockAndFlow0U;kw...) = add_part!(p,:DU;kw...)
add_dunits!(p::StockAndFlow0U,n;kw...) = add_parts!(p,:DU,n;kw...)

add_UDUlink!(p::StockAndFlow0U,u,du;kw...) = add_part!(p,:LU;luu=u,ludu=du,kw...)
add_UDUlinks!(p::StockAndFlow0U,n,u,du;kw...) = add_parts!(p,:LU,n;luu=u,ludu=du, kw...)

uname(p::StockAndFlow0U,u) = subpart(p,u,:uname) 
unames(p::StockAndFlow0U) = subpart(p,:uname) 

set_unames!(p::StockAndFlow0U, names) = set_subpart!(p, :uname, names)
set_dunames!(p::StockAndFlow0U, names) = set_subpart!(p, :duname, names)
set_exps!(p::StockAndFlow0U, exps) = set_subpart!(p, :exp, exps)


set_sdu!(p::StockAndFlow0U, du) = set_subpart!(p, :sdu, du)



add_unit!(p::StockAndFlowU;kw...) = add_part!(p,:U;kw...)
add_units!(p::StockAndFlowU,n;kw...) = add_parts!(p,:U,n;kw...)

add_dunit!(p::StockAndFlowU;kw...) = add_part!(p,:DU;kw...)
add_dunits!(p::StockAndFlowU,n;kw...) = add_parts!(p,:DU,n;kw...)

add_UDUlink!(p::StockAndFlowU,u,du;kw...) = add_part!(p,:LU;luu=u,ludu=du,kw...)
add_UDUlinks!(p::StockAndFlowU,n,u,du;kw...) = add_parts!(p,:LU,n;luu=u,ludu=du, kw...)

uname(p::StockAndFlowU,u) = subpart(p,u,:uname) 
unames(p::StockAndFlowU) = subpart(p,:uname) 

duname(p::Union{StockAndFlowU, StockAndFlow0U}, u) = subpart(p,u, :duname) 
dunames(p::Union{StockAndFlowU, StockAndFlow0U}) = subpart(p,:duname) 


set_unames!(p::StockAndFlowU, names) = set_subpart!(p, :uname, names)
set_dunames!(p::StockAndFlowU, names) = set_subpart!(p, :duname, names)
set_exps!(p::StockAndFlowU, exps) = set_subpart!(p, :exp, exps)


set_sdu!(p::StockAndFlowU, du) = set_subpart!(p, :sdu, du)
set_pdu!(p::StockAndFlowU, du) = set_subpart!(p, :pdu, du)

get_sdu(p::Union{StockAndFlowU, StockAndFlow0U}, s) = subpart(p, :sdu)[s]
get_pdu(p::StockAndFlowU, param) = subpart(p, :pdu)[param]

# get_sdu(p::StockAndFlowU, du) = subpart(p, :pdu, )


lunames(p::Union{StockAndFlowU, StockAndFlow0U}) = begin
    luu = map(x->subpart(p,x,:luu),collect(1:nlu(p)))
    ludu = map(x->subpart(p,x,:ludu),collect(1:nlu(p)))
    unit_names = map(x->uname(p,x),luu)
    dunit_names = map(x->duname(p,x),ludu)
    plu = collect(zip(unit_names, dunit_names))
end






function convert(::Type{StockAndFlowU}, sff::K) where {K <: AbstractStockAndFlowF}
    sfg = StockAndFlowU()
    
    add_dunit!(sfg, duname=:NONE)
    
    add_stocks!(sfg, ns(sff), sname=snames(sff), sdu=ones(Int, ns(sff)))
    add_svariables!(sfg, nsv(sff), svname=svnames(sff))
    add_variables!(sfg, nvb(sff), vname=vnames(sff))
    add_flows!(sfg, (fv(sff, i) for i in 1:nf(sff)), nf(sff), fname=fnames(sff))
    add_parameters!(sfg, np(sff), pname=pnames(sff), pdu=ones(Int, np(sff)))
    
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
    
    

    dunits = collect(setdiff(Set(sunits) ∪ Set(punits)))
    dunit_dict = state_dict(dunits)

      

    dunit_exp_dict = Dict{Union{Symbol, Expr}, Dict{Symbol, Float64}}()
    unit_set = Set{Symbol}()
    
    for du in dunits
        if du == :NONE
            continue
        end
        dunit_exp_dict[du] = extract_exponents(du)
        union!(unit_set, keys(dunit_exp_dict[du]))
    end

    
    units = collect(unit_set)
    
    unit_dict = state_dict(units)
    

    
    add_units!(sfg, length(units), uname=([Symbol(x) for x in units]))
    add_dunits!(sfg, length(dunits), duname=([Symbol(x) for x in dunits]))
    
    for du in dunits
        if du == :NONE
            continue
        end
        for (u, exp) in pairs(dunit_exp_dict[du])
            add_UDUlink!(sfg, unit_dict[u], dunit_dict[du], exp=exp)
        end
    end
    


    sdu = collect(map(x -> dunit_dict[x], sunits))
    pdu = collect(map(x -> dunit_dict[x], punits))
   
    
    add_stocks!(sfg, ns(sff), sname=snames(sff), sdu=sdu)
    add_svariables!(sfg, nsv(sff), svname=svnames(sff))
    add_variables!(sfg, nvb(sff), vname=vnames(sff))
    add_flows!(sfg, (fv(sff, i) for i in 1:nf(sff)), nf(sff), fname=fnames(sff))
    add_parameters!(sfg, np(sff), pname=pnames(sff), pdu=pdu)
    
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

      if ndu(a)>0 # if have links between stocks and sum-auxiliary-variables
        ϕdu = ntcomponent(dunames(a), dunames(x))
      else
        ϕdu = Int[]
      end

      if nlu(a)>0 # if have links between stocks and sum-auxiliary-variables        
        ϕlu = ntcomponent(lunames(a), lunames(x))
      else
        ϕlu = Int[]
      end


    result = OpenACSetLeg(a, S=ϕs, LS=ϕls, SV=ϕsv, U=ϕu, DU=ϕdu, LU=ϕlu)

    result
end


Open(p::StockAndFlowU, feet...) = begin
  legs = map(x->leg(x, p), feet)
  OpenStockAndFlowU{Symbol,Symbol,Int8,Float64}(p, legs...)
end











StockAndFlowU(s,p,v,f,sv,du,u) = begin  
  sf = StockAndFlowU()

  s = vectorify(s)
  f = vectorify(f)
  v = vectorify(v)
  p = vectorify(p)
  sv = vectorify(sv)
  du = vectorify(du)
  u = vectorify(u)

  sname = map(x -> first(first(x)), s)
  fname = map(first,f)
  vname = map(first,v)
  pname = map(first, p)
  duname = map(first, du)

  op=map(last,map(last,v))

  s_idx = state_dict(sname)
  f_idx = state_dict(fname)
  v_idx = state_dict(vname)
  p_idx = state_dict(p)
  sv_idx = state_dict(sv)
  du_idx = state_dict(duname)
  u_idx = state_dict(u)

  add_units!(sf, length(u), uname=u)
  add_dunits!(sf, length(du), duname=duname)

  luu = Vector{Int}()
  ludu = Vector{Int}()
  exp = Vector{Float64}()

  for (derived_unit, unit_exponent_pairs) in du
    append!(ludu, fill(du_idx[derived_unit], length(unit_exponent_pairs)))
    append!(luu, [u_idx[u] for u in first.(unit_exponent_pairs)])
    append!(exp, [power for power in last.(unit_exponent_pairs)])
  end


  add_UDUlinks!(sf, length(ludu), luu, ludu ; exp=exp)

  add_parameters!(sf,length(p),pname=pname, pdu=[du_idx[derived_unit] for derived_unit in map(last, p)])
  add_svariables!(sf, length(sv), svname=sv)
  add_variables!(sf, length(vname), vname=vname, vop=op)
  add_flows!(sf,map(x->v_idx[x], map(last,f)),length(fname),fname=fname)


  # Parse the elements included in "s" -- stocks
  for (i, ((name, dunit),(ins,outs,svs))) in enumerate(s)
    i = add_stock!(sf,sname=name, sdu=du_idx[dunit]) # add objects :S (stocks)
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



end