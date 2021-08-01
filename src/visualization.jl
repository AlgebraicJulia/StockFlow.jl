using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Edge
import Base.Iterators: flatten
using StatsBase

export Graph, LGraph

# get the stock index of :is or :os with the flow index of findex from the
# stock inflow and outflow tables: p.tables.I or p.tables.O
# it is notable that a flow only relates to one inflow and/or one outflow.
get_stock_index(sa::Array, fa::Array, findex::Int) = begin
   for i in 1:length(fa)
        if fa[i] == findex
            return sa[i]
        end
    end
end

# only plot stocks and flows, not plot links
function Graph(p::Union{AbstractBoneStockFlow,AbstractLinkedBoneStockFlow})
  innerFlows=intersect(p.tables.I.ifn,p.tables.O.ofn)
  edgeInFlows=symdiff(p.tables.I.ifn,innerFlows)
  edgeOutFlows=symdiff(p.tables.O.ofn,innerFlows)

  stockNodes = [Node(string("$(sname(p, s))"), Attributes(:shape=>"square", :color=>"black")) for s in 1:ns(p)]
  boundInFlowNodes = [Node(string("fs_$(edgeInFlows[s])"),Attributes(:shape=>"point", :color=>"white")) for s in 1:length(edgeInFlows)]
  boundOutFlowNodes = [Node(string("fs_$(edgeOutFlows[s])"),Attributes(:shape=>"point", :color=>"white")) for s in 1:length(edgeOutFlows)]

  stmts_nodes = vcat(stockNodes, boundInFlowNodes, boundOutFlowNodes)

  edges_inner=map(1:length(innerFlows)) do k
    flow_index=innerFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_outfrom = get_stock_index(p.tables.O.os, p.tables.O.ofn, flow_index)
    stock_index_into = get_stock_index(p.tables.I.is, p.tables.I.ifn, flow_index)
    [Edge(["$(sname(p,stock_index_outfrom))", "$(sname(p,stock_index_into))"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
  end |> flatten |> collect

  edges_inflow=map(1:length(edgeInFlows)) do k
    flow_index=edgeInFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_into = get_stock_index(p.tables.I.is, p.tables.I.ifn, flow_index)
    [Edge(["fs_$flow_index", "$(sname(p,stock_index_into))"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
  end |> flatten |> collect

  edges_outflow=map(1:length(edgeOutFlows)) do k
    flow_index=edgeOutFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_outfrom = get_stock_index(p.tables.O.os, p.tables.O.ofn, flow_index)
    [Edge(["$(sname(p,stock_index_outfrom))", "fs_$flow_index"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
  end |> flatten |> collect

  stmts_edges = vcat(edges_inner,edges_inflow,edges_outflow)
  stmts = vcat(stmts_nodes, stmts_edges)

  graph_attrs = Attributes(:rankdir=>"LR")
#  node_attrs  = Attributes(:shape=>"plain", :style=>"filled", :color=>"white")
  edge_attrs  = Attributes(:splines=>"splines")

  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
  return g
end

# plot all properties of stock and flow diagrams, including links
function LGraph(p::AbstractLinkedBoneStockFlow)
  innerFlows=intersect(p.tables.I.ifn,p.tables.O.ofn)
  edgeInFlows=symdiff(p.tables.I.ifn,innerFlows)
  edgeOutFlows=symdiff(p.tables.O.ofn,innerFlows)

  stockNodes = [Node(string("$(sname(p, s))"), Attributes(:shape=>"square", :color=>"black", :style=>"filled", :fillcolor=>"#9ACEEB")) for s in 1:ns(p)]
  flowNodes = [Node(string("fn_$(fname(p, f))"), Attributes(:shape=>"invtriangle", :color=>"#9ACEEB", :style=>"filled", :label=>"", :width=>"0.1", :height=>"0.2")) for f in 1:nf(p)]
  boundInFlowNodes = [Node(string("fs_$(edgeInFlows[s])"),Attributes(:shape=>"point", :color=>"white")) for s in 1:length(edgeInFlows)]
  boundOutFlowNodes = [Node(string("fs_$(edgeOutFlows[s])"),Attributes(:shape=>"point", :color=>"white")) for s in 1:length(edgeOutFlows)]

  stmts_nodes = vcat(stockNodes, boundInFlowNodes, boundOutFlowNodes,flowNodes)

  edges_inner=map(1:length(innerFlows)) do k
    flow_index=innerFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_outfrom = get_stock_index(p.tables.O.os, p.tables.O.ofn, flow_index)
    stock_index_into = get_stock_index(p.tables.I.is, p.tables.I.ifn, flow_index)
    [Edge(["$(sname(p,stock_index_outfrom))", "fn_$flow_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
     Edge(["fn_$flow_name", "$(sname(p,stock_index_into))"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
  end |> flatten |> collect

  edges_inflow=map(1:length(edgeInFlows)) do k
    flow_index=edgeInFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_into = get_stock_index(p.tables.I.is, p.tables.I.ifn, flow_index)
    [Edge(["fs_$flow_index", "fn_$flow_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
     Edge(["fn_$flow_name", "$(sname(p,stock_index_into))"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
  end |> flatten |> collect

  edges_outflow=map(1:length(edgeOutFlows)) do k
    flow_index=edgeOutFlows[k]
    flow_name=fname(p,flow_index)
    stock_index_outfrom = get_stock_index(p.tables.O.os, p.tables.O.ofn, flow_index)
    [Edge(["$(sname(p,stock_index_outfrom))", "fn_$flow_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
     Edge(["fn_$flow_name", "fs_$flow_index"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
  end |> flatten |> collect

  Ltable=p.tables.L
  edges_link=map(1:length(Ltable)) do k
    flow_index=Ltable.lf[k]
    stock_index=Ltable.ls[k]
    [Edge(["$(sname(p, stock_index))", "fn_$(fname(p,flow_index))"])]
  end |> flatten |> collect


  stmts_edges = vcat(edges_inner,edges_inflow,edges_outflow)
  stmts = vcat(stmts_nodes, stmts_edges)

  stmts = vcat(stmts, edges_link)

  graph_attrs = Attributes(:rankdir=>"LR")

  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs)
  return g
end

function Graph(op::Union{OpenBoneStockFlow_S, OpenLabelledBoneStockFlowUntyped_S, OpenStockFlow_S, OpenLabelledStockFlowUntyped_S, OpenBoneStockFlow_F, OpenLabelledBoneStockFlowUntyped_F, OpenStockFlow_F, OpenLabelledStockFlowUntyped_F})
    Graph(apex(op))
end
