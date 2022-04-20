export TheoryStockAndFlowp, AbstractStockAndFlowp, StockAndFlowp, OpenStockAndFlowpOb, OpenStockAndFlowp


# define the preliminary schema (including attributes)
@present TheoryStockAndFlowp(FreeSchema) begin

# Objects:
  Flow::Ob
  Stock::Ob
  Link::Ob

# Morphisms:
  u::Hom(Flow, Stock)
  d::Hom(Flow, Stock)
  s::Hom(Link, Stock)
  t::Hom(Link, Flow)

# Attributes:
  Name::AttrType
  FuncFlow::AttrType
  
  sname::Attr(Stock, Name)
  fname::Attr(Flow, Name)
  ϕf::Attr(Flow, FuncFlow)
end

@abstract_acset_type AbstractStockAndFlowp
@acset_type StockAndFlowpUntyped(TheoryStockAndFlowp, index=[:u,:d,:s,:t]) <: AbstractStockAndFlowp

const StockAndFlowp = StockAndFlowpUntyped{Symbol,Function} 
# open stock and flow diagram
const OpenStockAndFlowpObUntyped, OpenStockAndFlowpUntyped = OpenACSetTypes(StockAndFlowpUntyped,:Stock)
const OpenStockAndFlowpOb = OpenStockAndFlowpObUntyped{Symbol,Function}
const OpenStockAndFlowp = OpenStockAndFlowpUntyped{Symbol,Function}

Open(p::AbstractStockAndFlowp) = OpenStockAndFlowp(p, map(x->FinFunction([x], ns(p)), 1:ns(p))...)
Open(p::AbstractStockAndFlowp, legs...) = begin
  s_idx = Dict(sname(p, s)=>s for s in 1:ns(p))
  OpenStockAndFlowp(p, map(l->FinFunction(map(i->s_idx[i], l), ns(p)), legs)...)
end
Open(n, p::AbstractStockAndFlowp, m) = Open(p, n, m)


StockAndFlowp(s,f) = begin
	d = StockAndFlowp()

    s = vectorify(s)
    add_stocks!(d,length(s),sname=s)

    s_idx = state_dict(s)
    for (i, ((fattr,uds),ls)) in enumerate(f)
      fn = first(fattr)
      ff = last(fattr)
    	sui = s_idx[first(uds)]
    	sdi = s_idx[last(uds)]
    	ls = vectorify(ls)
    	add_flow!(d,sui,sdi,fname=fn,ϕf=ff)
    	add_links!(d,map(x->s_idx[x],ls),repeat([i], length(ls)), length(ls))
    end
    d
end

add_flow!(p::AbstractStockAndFlowp,su,sd;kw...) = add_part!(p,:Flow;u=su,d=sd,kw...)
add_flows!(p::AbstractStockAndFlowp,su,sd,n;kw...) = add_parts!(p,:Flow,n;u=su,d=sd,kw...)

add_stock!(p::AbstractStockAndFlowp;kw...) = add_part!(p,:Stock;kw...) 
add_stocks!(p::AbstractStockAndFlowp,n;kw...) = add_parts!(p,:Stock,n;kw...)

add_link!(p::AbstractStockAndFlowp,ss,ft;kw...) = add_part!(p,:Link;s=ss,t=ft,kw...)
add_links!(p::AbstractStockAndFlowp,ss,ft,n;kw...) = add_parts!(p,:Link,n;s=ss,t=ft,kw...)  

ns(p::AbstractStockAndFlowp) = nparts(p,:Stock) #number of stocks
nf(p::AbstractStockAndFlowp) = nparts(p,:Flow) #number of flows
nl(p::AbstractStockAndFlowp) = nparts(p,:Link) #number of links

sname(p::AbstractStockAndFlowp,s) = subpart(p,s,:sname) # return the stocks name with index of s
fname(p::AbstractStockAndFlowp,f) = subpart(p,f,:fname) # return the flows name with index of f

funcFlow(p::AbstractStockAndFlowp,f) = subpart(p,f,:ϕf)
funcFlows(p::AbstractStockAndFlowp)=begin
    fnames = [fname(p, f) for f in 1:nf(p)]
    LVector(;[(fnames[f]=>funcFlow(p, f)) for f in 1:nf(p)]...)
end

# return inflows of stock index s
inflows(p::AbstractStockAndFlowp,s) = incident(p,s,:d)
# return outflows of stock index s
outflows(p::AbstractStockAndFlowp,s) = incident(p,s,:u)

valueat(x::Number, u, p, t) = x
valueat(f::Function, u, p, t) = f(u,p,t)

# take (F: stock and flow diagram with attributes) and return ODEs of the dynamical system
# Note: ϕ: functions of flows are simply an attribute in F
vectorfield(ps::AbstractStockAndFlowp) = begin
  ϕ=funcFlows(ps)
  f(du,u,p,t) = begin
    for i in 1:ns(ps)
      du[sname(ps, i)] = 0
      for m in 1:length(inflows(ps,i))
        du[sname(ps, i)] = du[sname(ps, i)] + valueat(ϕ[fname(ps,inflows(ps,i)[m])],u,p,t)
      end
      for n in 1:length(outflows(ps,i))
        du[sname(ps, i)] = du[sname(ps, i)] - valueat(ϕ[fname(ps,outflows(ps,i)[n])],u,p,t)
      end      
    end
    return du
  end
  return f
end




































