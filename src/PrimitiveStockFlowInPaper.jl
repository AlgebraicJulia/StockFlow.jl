export TheoryStockAndFlowp, AbstractStockAndFlowp, StockAndFlowp, OpenStockAndFlowpOb, OpenStockAndFlowp, checkfls,
upstock, downstock, linkstock, linkflow

using Combinatorics
# define the primitive schema (including attributes)
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

# return stock of flow f out from (upstream)
upstock(p::AbstractStockAndFlowp,f) = subpart(p,f,:u)
# return stock of flow f stream in (downstream)
downstock(p::AbstractStockAndFlowp,f) = subpart(p,f,:d)

# return the source of a link l, which is a stock
linkstock(p::AbstractStockAndFlowp,l) = subpart(p,l,:s)
# return the target of a link l, which is a flow
linkflow(p::AbstractStockAndFlowp,l) = subpart(p,l,:t)


# return stocks of flow f link to
flinks(p::AbstractStockAndFlowp,f) = subpart(p,incident(p,f,:t),:s)

valueat(x::Number, u, p, t) = x
valueat(f::Function, u, p, t) = f(u,p,t)

# test argumenterror -- stocks in function of flow "fn" are not linked!
# TODO: find method to generate the exact wrong stocks' names and output in error message
ftest(f::Function, u, p, fn) = 
try
  f(u,p,0)
catch e
  if isa(e, ArgumentError)
    println(string("stocks used in the function of flow ", fn, " are not linked but used!"))
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

# check the dependency of links of flow index f
checkfl(ps::AbstractStockAndFlowp, u, p, f) = begin
    s = flinks(ps,f)
    uslabel = map(x->sname(ps,x),s)
    usvalue = map(x->u[x],uslabel)
    # generate the labelled array only includes the stocks the flow linked to
    uf = @LArray usvalue tuple(uslabel...)

    functionf = funcFlow(ps,f)
    namef = fname(ps,f)
    # check used but not linked stocks
    ftest(functionf, uf, p, namef)
    # check linked but not used stocks

    #TODO: what about the situation: the flow function simply return a value, but the flow has stocks linked to??
    # 1. check if a flow and stocks linked, but no stocks used in the function
#    if length(uslabel)>0
#      fmisstest(functionf, [], p, namef, uslabel)
#    end
    # 2. check at least 1 stock used in the flow function
    sub_uslabels=setdiff(collect(combinations(uslabel)),[uslabel]) # generate the sub-arrays of uslabel except itself
    for sub_uslabel in sub_uslabels
      sub_usvalue = map(x->u[x],sub_uslabel)
      sub_uf = @LArray sub_usvalue tuple(sub_uslabel...)
      fmisstest(functionf, sub_uf, p, namef, setdiff(uslabel, sub_uslabel))
    end
end

# check the dependency of links of all flows
checkfls(ps::AbstractStockAndFlowp, u, p) = begin
    for f in 1:nf(ps)
      checkfl(ps, u, p, f)
    end
    "Great! Links dependency check passed!"
end

# take (F: stock and flow diagram with attributes) and return ODEs of the dynamical system
# Note: ϕ: functions of flows are simply an attribute in F
vectorfield(ps::AbstractStockAndFlowp) = begin
  # return ODEs
  ϕ=funcFlows(ps)
  f(du,u,p,t) = begin
    for i in 1:ns(ps)
#      println(i)
      du[sname(ps, i)] = 0
#      println("in")
      for m in 1:length(inflows(ps,i))
#        println(fname(ps,inflows(ps,i)[m]))
        du[sname(ps, i)] = du[sname(ps, i)] + valueat(ϕ[fname(ps,inflows(ps,i)[m])],u,p,t)
      end
#      println("out")
      for n in 1:length(outflows(ps,i))
#        println(fname(ps,outflows(ps,i)[n]))
        du[sname(ps, i)] = du[sname(ps, i)] - valueat(ϕ[fname(ps,outflows(ps,i)[n])],u,p,t)
      end      
    end
    return du
  end
  return f
end

Graph(p::AbstractStockAndFlowp, rd::String="LR") = begin
  stockNodes = [Node(string("$(sname(p, s))"), Attributes(:shape=>"square", :color=>"black", :style=>"filled", :fillcolor=>"#9ACEEB")) for s in 1:ns(p)]
  flowNodes = [Node(string("fn_$(fname(p, f))"), Attributes(:shape=>"invtriangle", :color=>"#9ACEEB", :style=>"filled", :label=>"", :width=>"0.1", :height=>"0.2")) for f in 1:nf(p)]

  # edges of flow
  edges_flows = map(1:nf(p)) do k
    flow_name=fname(p,k)
    stock_index_outfrom = upstock(p,k)
    stock_index_into = downstock(p,k)
    [Edge(["$(sname(p,stock_index_outfrom))", "fn_$flow_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
     Edge(["fn_$flow_name", "$(sname(p,stock_index_into))"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
  end |> flatten |> collect

  # edges of links
  edges_links=map(1:nl(p)) do k
    flow_index=linkflow(p,k)
    stock_index=linkstock(p,k)
    [Edge(["$(sname(p, stock_index))", "fn_$(fname(p,flow_index))"])]
  end |> flatten |> collect

  stmts_nodes = vcat(stockNodes, flowNodes)
  stmts_edges = vcat(edges_flows,edges_links)
  stmts = vcat(stmts_nodes, stmts_edges)

  graph_attrs = Attributes(:rankdir=>rd)
  edge_attrs  = Attributes(:splines=>"ortho")

  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
  return g
end




































