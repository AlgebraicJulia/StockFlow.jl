module StockFlowUnits

using ..StockFlow
using Catlab.GATs.Presentations, Catlab.CategoricalAlgebra
using MLStyle

import ...StockFlow: StockAndFlowF, state_dict, ns, np

export add_unit!, add_units!, add_cunit!, add_cunits!, add_UCUlink!, add_UCUlinks!,
uname, unames, set_unames!, set_exps!, convert, StockAndFlowU, set_scu!, set_pcu!,
FtoU

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
unames(p::StockAndFlowU,u) = subpart(p,:uname) 

set_unames!(p::StockAndFlowU, names) = set_subpart!(p, :uname, names)
set_cunames!(p::StockAndFlowU, names) = set_subpart!(p, :cuname, names)
set_exps!(p::StockAndFlowU, exps) = set_subpart!(p, :exp, exps)


set_scu!(p::StockAndFlowU, cu) = set_subpart!(p, :scu, cu)
set_pcu!(p::StockAndFlowU, cu) = set_subpart!(p, :pcu, cu)








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
    
    

    cunits = collect(Set(sunits) ∪ Set(punits))
    cunit_dict = state_dict(cunits)

      

    cunit_exp_dict = Dict{Union{Symbol, Expr}, Dict{Symbol, Float64}}()
    unit_set = Set{Symbol}()
    
    for cu in cunits
        cunit_exp_dict[cu] = extract_exponents(cu)
        union!(unit_set, keys(cunit_exp_dict[cu]))
    end

    
    units = collect(unit_set)
    
    unit_dict = state_dict(units)
    

    
    add_units!(sfg, length(units), uname=([Symbol(x) for x in units]))
    add_cunits!(sfg, length(cunits), cuname=([Symbol(x) for x in cunits]))
    
    for cu in cunits
        for (u, exp) in pairs(cunit_exp_dict[cu])
            add_UCUlink!(sfg, unit_dict[u], cunit_dict[cu], exp=exp)
        end
    end
    


    scu = collect(map(x -> cunit_dict[x], sunits))
    pcu = collect(map(x -> cunit_dict[x], punits))
    
    #set_subpart!(sf, :scu, scu)
    #set_subpart!(sf, :pcu, pcu)

    #for i in 1:ns(sf)
     #   sf.subparts.scu[i] = scu[i]
    #end
    #for i in 1:np(sf)
    #    sf.subparts.pcu[i] = pcu[i]
    #end
    #set_scu!(sf, scu)
    #set_pcu!(sf, pcu)
    
    
    
    
    #add_cunit!(sfg, cuname=:NONE)
    
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



end
