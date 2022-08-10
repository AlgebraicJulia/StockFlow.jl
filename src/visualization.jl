using Catlab.CategoricalAlgebra
using Catlab.Graphics.Graphviz
import Catlab.Graphics.Graphviz: Graph, Subgraph
import Base.Iterators: flatten
using StatsBase
using Catlab.Graphics

export Graph, display_uwd

# schema: "c":  the full schema
#         "c0": the simple schema
# the type parameter is only used for schema C, which means if schema=C0, all the component will be plot out anyway, since C0 is quite simple
# type: "SFVL": include all component
#       "SF": only include stocks and flows
#       "SFV": only include stocks, flows, variables (include both auxiliary variables and sum auxiliary variables)

function Graph(p::AbstractStockAndFlow0; schema::String="C", type::String="SFVL", rd::String="LR")

# only full schema C has Flows
  if schema == "C" begin
      inflows=inflowsAll(p)
      outflows=outflowsAll(p)
      innerFlows=intersect(inflows,outflows)
      edgeInFlows=symdiff(inflows,innerFlows)
      edgeOutFlows=symdiff(outflows,innerFlows)
    end
  end

#  stockNodes = [Node(string("$(sname(p, s))"), Attributes(:shape=>"square", :color=>"black", :style=>"filled", :fillcolor=>"#9ACEEB")) for s in 1:ns(p)]
  stockNodes = [Node("s$s", Attributes(:label=>"$(sname(p, s))",:shape=>"square", :color=>"black", :style=>"filled", :fillcolor=>"#9ACEEB")) for s in 1:ns(p)]

  if occursin("V", type)
    if schema == "C"
      vNodes = [Node("v$v", Attributes(:label=>"$(vname(p, v))",:shape=>"plaintext", :color=>"black")) for v in 1:nvb(p)]
    end
    svNodes = [Node("sv$sv", Attributes(:label=>"$(svname(p, sv))",:shape=>"circle", :color=>"black",:fillcolor=>"cornflowerblue", :style=>"filled")) for sv in 1:nsv(p)]
  end
  
  if schema == "C" begin
      boundInFlowNodes = [Node(string("fs_$(edgeInFlows[s])"),Attributes(:label=>"",:shape=>"point", :color=>"white")) for s in 1:length(edgeInFlows)]
      boundOutFlowNodes = [Node(string("fs_$(edgeOutFlows[s])"),Attributes(:label=>"",:shape=>"point", :color=>"white")) for s in 1:length(edgeOutFlows)]
    end
  end
  stmts_nodes = if schema == "C"
                  if occursin("V", type) # SFVL or SFV
                    vcat(stockNodes, boundInFlowNodes, boundOutFlowNodes, vNodes, svNodes)
                  else # SF
                    vcat(stockNodes, boundInFlowNodes, boundOutFlowNodes)
                  end
                else # schema == C0
                  vcat(stockNodes, svNodes)
                end
                  
  # define edges of flows
  if schema == "C" begin 
    edges_inner=map(1:length(innerFlows)) do k
      flow_index=innerFlows[k]
      flow_name=fname(p,flow_index)
      fv_name=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        [Edge(["s$stock_index_outfrom", "v$fv_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
         Edge(["v$fv_name", "s$stock_index_into"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
      else
        [Edge(["s$stock_index_outfrom", "s$stock_index_into"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
      end
    end |> flatten |> collect

    edges_inflow=map(1:length(edgeInFlows)) do k
      flow_index=edgeInFlows[k]
      flow_name=fname(p,flow_index)
      fv_name=flowVariableIndex(p,flow_index)
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        [Edge(["fs_$flow_index", "v$fv_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
         Edge(["v$fv_name", "s$stock_index_into"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
      else
        [Edge(["fs_$flow_index", "s$stock_index_into"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
      end
    end |> flatten |> collect

    edges_outflow=map(1:length(edgeOutFlows)) do k
      flow_index=edgeOutFlows[k]
      flow_name=fname(p,flow_index)
      fv_name=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      if occursin("V", type)
        [Edge(["s$stock_index_outfrom", "v$fv_name"],Attributes(:label=>"", :labelfontsize=>"6", :color=>"black:invis:black", :arrowhead=>"none", :splines=>"ortho")),
         Edge(["v$fv_name", "fs_$flow_index"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black", :splines=>"ortho"))]
      else
        [Edge(["s$stock_index_outfrom", "fs_$flow_index"],Attributes(:label=>"$flow_name", :labelfontsize=>"6", :color=>"black:invis:black"))]
      end
    end |> flatten |> collect

    stmts_edges_flows = if schema == "C"
                           vcat(edges_inner,edges_inflow,edges_outflow)
                         end
    end 
  end

  # linkages from S to sv
  edges_lsv()=begin
    subEdges = Vector{Edge}[]  
    for k in 1:nsv(p)
        subEdges = vcat(map(stockssv(p,k)) do m
                [Edge(["s$m", "sv$k"])]
                end, subEdges)
    end
    return subEdges |> flatten |> collect
  end

  if schema == "C" && occursin("L", type) begin
    # plot the linkages
    # linkages from S to v
    edges_lv()=begin
      subEdges = Vector{Edge}[]  
      for k in 1:nvb(p)
        subEdges = vcat(map(stocksv(p,k)) do m
                [Edge(["s$m", "v$k"])]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end

    # linkages from sv to v
    edges_lsvv()=begin
      subEdges = Vector{Edge}[]  
      for k in 1:nvb(p)
        subEdges = vcat(map(svsv(p,k)) do m
                [Edge(["sv$m", "v$k"])]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end
   end
  end 

  stmts_edges_links = if schema == "C0"
                         edges_lsv()
                      else
                        if occursin("L", type)
                          vcat(edges_lv(), edges_lsv(), edges_lsvv())
                        end
                      end

  stmts = if schema == "C0"
             vcat(stmts_nodes, stmts_edges_links)
          else
             if occursin("L", type)
              vcat(stmts_nodes, stmts_edges_flows, stmts_edges_links)
             else
              vcat(stmts_nodes, stmts_edges_flows)
             end
          end

  graph_attrs = Attributes(:rankdir=>rd)
  edge_attrs  = Attributes(:splines=>"splines")

  g = Graphviz.Digraph("G", stmts; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
  return g

end


display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"))
