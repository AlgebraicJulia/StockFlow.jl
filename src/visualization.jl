using Catlab.CategoricalAlgebra
using Catlab.Graphics.Graphviz: Attributes, Node, Edge, Digraph
using Base.Iterators: flatten
using StatsBase

export display_uwd, GraphF, GraphRB, GraphCL, Graph

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"))

def_stock(p, s) = ("s$s", Attributes(:label=>"$(sname(p, s))",
                                     :shape=>"square",
                                     :color=>"black",
                                     :style=>"filled",
                                     :fillcolor=>"#9ACEEB"))

def_parameter(p, pp) = ("p$pp", Attributes(:label=>"$(pname(p, pp))",
                                     :shape=>"circle",
                                     :color=>"black"))

def_auxiliaryV(p, v) = ("v$v", Attributes(:label=>"$(vname(p, v))",
                                          :shape=>"plaintext",
                                          :fontcolor=>"black"))

#def_auxiliaryVF(p, v) = ("v$v", Attributes(:label=>"$(vname(p, v))"*"="*"$(make_v_expr(p,v))",
#                                          :shape=>"plaintext",
#                                          :fontcolor=>"black"))

def_auxiliaryVF(p, v) = ("v$v", Attributes(:label=>p isa AbstractStockAndFlowF ? "$(make_v_expr(p,v))" : "$(vname(p, v))",
                                          :shape=>"plaintext",
                                          :fontcolor=>"black"))

def_sumV(p, sv) = ("sv$sv", Attributes(:label=>"$(svname(p, sv))",
                                       :shape=>"circle",
                                       :color=>"black",
                                       :fillcolor=>"cornflowerblue",
                                       :style=>"filled"))

# currently, we use invisible fake stock nodes as outer clouds
# type: string of "u" or "d"
#       "u" indicates this cloud node is an upstream cloud to a half-edge inflow
#       "d" indicates this cloud node is a downstream cloud to a half-edge outflow
def_cloud(c, type) = ("fs_$c"*type, Attributes(:label=>"",
                                               :shape=>"point",
                                               :color=>"white"))

# function of define the flow with auxiliary variables are also output, e.g., type = "SFVL" or "SFV"
# us: string of up stream stock name (including clouds) of the flow
# ds: string of down stream stock nme (including clouds) of the flow
# v: the index of the auxiliary variable of the flow
# fname: the index of the flow
def_flow_V(p, us, ds, v, f) = begin
    labelfontsize = "6"
    color = "black:invis:black"
    arrowhead = "none"
    splines = "ortho"
    return ([us, "v$v"],Attributes(:label=>"", :labelfontsize=>labelfontsize, :color=>color, :arrowhead=>arrowhead, :splines=>splines)),
           (["v$v", ds],Attributes(:label=>"$(fname(p,f))", :labelfontsize=>labelfontsize, :color=>color, :splines=>splines))
end

# function of define the flow without auxiliary variables, e.g., type = "SF"
def_flow_noneV(p, us, ds, f) = ([us, ds],Attributes(:label=>"$(fname(p,f))", :labelfontsize=>"6", :color=>"black:invis:black"))

# s: string of the source name
# t: string of the target name
def_link(s,t) = ([s, t])

"""
Graph a Causal Loop diagram (no polarities).
"""
function GraphCL(c::CausalLoop)

  NNodes = [Node("n$n", Attributes(:label=>"$(vname(c, n))",:shape=>"plaintext")) for n in 1:nvert(c)]

  Edges=map(1:nedges(c)) do k
    [Edge(["n$(sedge(c,k))", "n$(tedge(c,k))"],Attributes(:color=>"blue"))]
  end |> flatten |> collect

  stmts=vcat(NNodes,Edges)

  g = Digraph("G", stmts;graph_attrs=Attributes(:rankdir=>"LR"))
  return g

end

"""
Graph a CausalLoopPM.
"""
function GraphCL(c::CausalLoopPM; schema="BASE")
  GraphCL(to_clp(c) ; schema=schema)
end

"""
Graph a CausalLoopPol.
"""
function GraphCL(c::CausalLoopPol; schema="BASE")

  NNodes = [Node("n$n", Attributes(:label=>"$(vname(c, n))",:shape=>"plaintext")) for n in 1:nvert(c)]

  if occursin("BASE", schema)
    Edges=map(1:nedges(c)) do k
      pol_int = epol(c,k)
      if pol_int == POL_POSITIVE
        pol = :+
      elseif pol_int == POL_NEGATIVE
        pol = :-
      else
        error("Unknown Polarity $pol_int.")
      end
      [Edge(["n$(sedge(c,k))", "n$(tedge(c,k))"],Attributes(:color=>"blue",:label=>"$(pol)"))]
    end |> flatten |> collect
  elseif schema == "C0"
    Edges=map(1:nedges(c)) do k
      [Edge(["n$(sedge(c,k))", "n$(tedge(c,k))"],Attributes(:color=>"blue"))]
    end |> flatten |> collect
  end

  stmts=vcat(NNodes,Edges)

  g = Digraph("G", stmts;graph_attrs=Attributes(:rankdir=>"LR"))
  return g

end

# schema: "c":  the full schema
#         "c0": the simple schema
# the type parameter is only used for schema C, which means if schema=C0, all the component will be plot out anyway, since C0 is quite simple
# type: "SFVL": include all component
#       "SF": only include stocks and flows
#       "SFV": only include stocks, flows, variables (include both auxiliary variables and sum auxiliary variables)

function Graph(p::AbstractStockAndFlow0; make_stock::Function=def_stock, make_auxiliaryV::Function=def_auxiliaryV,
                                         make_sumV::Function=def_sumV, make_cloud::Function=def_cloud,
                                         make_flow_V::Function=def_flow_V, make_flow_noneV::Function=def_flow_noneV,
                                         make_link::Function=def_link,
                                         schema::String="C", type::String="SFVL", rd::String="LR")

# only full schema C has Flows
  if schema == "C" begin
      inflows=inflowsAll(p)
      outflows=outflowsAll(p)
      innerFlows=intersect(inflows,outflows)
      edgeInFlows=symdiff(inflows,innerFlows)
      edgeOutFlows=symdiff(outflows,innerFlows)
    end
  end

  stockNodes = [Node(make_stock(p,s)...) for s in 1:ns(p)]


  if occursin("V", type)
    if schema == "C"
      vNodes = [Node(make_auxiliaryV(p,v)...) for v in 1:nvb(p)]

    end
    svNodes = [Node(make_sumV(p,sv)...) for sv in 1:nsv(p)]

  end

  if schema == "C" begin
      boundInFlowNodes = [Node(make_cloud(edgeInFlows[s],"u")...) for s in 1:length(edgeInFlows)] #id is e.g., "fs_1u"
      boundOutFlowNodes = [Node(make_cloud(edgeOutFlows[s],"d")...) for s in 1:length(edgeOutFlows)] #id is e.g., "fs_1d"
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
      fv=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"s$stock_index_outfrom", "s$stock_index_into", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"s$stock_index_outfrom", "s$stock_index_into", flow_index)...)]
      end
    end |> flatten |> collect

    edges_inflow=map(1:length(edgeInFlows)) do k
      flow_index=edgeInFlows[k]
      fv=flowVariableIndex(p,flow_index)
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"fs_$flow_index"*"u", "s$stock_index_into", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"fs_$flow_index"*"u", "s$stock_index_into", flow_index)...)]
      end
    end |> flatten |> collect

    edges_outflow=map(1:length(edgeOutFlows)) do k
      flow_index=edgeOutFlows[k]
      fv=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"s$stock_index_outfrom", "fs_$flow_index"*"d", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"s$stock_index_outfrom", "fs_$flow_index"*"d", flow_index)...)]
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
                [Edge(make_link("s$m", "sv$k")...)]
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
                [Edge(make_link("s$m", "v$k")...)]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end

    # linkages from sv to v
    edges_lsvv()=begin
      subEdges = Vector{Edge}[]
      for k in 1:nvb(p)
        subEdges = vcat(map(svsv(p,k)) do m
                 [Edge(make_link("sv$m", "v$k")...)]
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

  g = Digraph("G", stmts; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
  return g

end


# schema: "c":  the full schema
#         "c0": the simple schema
# the type parameter is only used for schema C, which means if schema=C0, all the component will be plot out anyway, since C0 is quite simple
# type: "SFVL": include all component
#       "SF": only include stocks and flows
#       "SFV": only include stocks, flows, variables (include both auxiliary variables and sum auxiliary variables)
function GraphF(p::AbstractStockAndFlow0; make_stock::Function=def_stock, make_auxiliaryV::Function=def_auxiliaryVF,
                                         make_sumV::Function=def_sumV, make_cloud::Function=def_cloud,
                                         make_flow_V::Function=def_flow_V, make_flow_noneV::Function=def_flow_noneV,
                                         make_link::Function=def_link, make_parameter::Function=def_parameter,
                                         schema::String="C", type::String="SFVL", rd::String="LR")

# only full schema C has Flows
  if schema == "C" begin
      inflows=inflowsAll(p)
      outflows=outflowsAll(p)
      innerFlows=intersect(inflows,outflows)
      edgeInFlows=symdiff(inflows,innerFlows)
      edgeOutFlows=symdiff(outflows,innerFlows)

      parameterNodes = [Node(make_parameter(p,pp)...) for pp in 1:np(p)]
    end
  end

  stockNodes = [Node(make_stock(p,s)...) for s in 1:ns(p)]


  if occursin("V", type)
    if schema == "C"
      vNodes = [Node(make_auxiliaryV(p,v)...) for v in 1:nvb(p)]

    end
    svNodes = [Node(make_sumV(p,sv)...) for sv in 1:nsv(p)]

  end

  if schema == "C" begin
      boundInFlowNodes = [Node(make_cloud(edgeInFlows[s],"u")...) for s in 1:length(edgeInFlows)] #id is e.g., "fs_1u"
      boundOutFlowNodes = [Node(make_cloud(edgeOutFlows[s],"d")...) for s in 1:length(edgeOutFlows)] #id is e.g., "fs_1d"
    end
  end
  stmts_nodes = if schema == "C"
                  if occursin("V", type) # SFVL or SFV
                    vcat(stockNodes, parameterNodes, boundInFlowNodes, boundOutFlowNodes, vNodes, svNodes)
                  else # SF
                    vcat(stockNodes, parameterNodes, boundInFlowNodes, boundOutFlowNodes)
                  end
                else # schema == C0
                  vcat(stockNodes, svNodes)
                end

  # define edges of flows
  if schema == "C" begin
    edges_inner=map(1:length(innerFlows)) do k
      flow_index=innerFlows[k]
      fv=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"s$stock_index_outfrom", "s$stock_index_into", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"s$stock_index_outfrom", "s$stock_index_into", flow_index)...)]
      end
    end |> flatten |> collect

    edges_inflow=map(1:length(edgeInFlows)) do k
      flow_index=edgeInFlows[k]
      fv=flowVariableIndex(p,flow_index)
      stock_index_into = first(instock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"fs_$flow_index"*"u", "s$stock_index_into", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"fs_$flow_index"*"u", "s$stock_index_into", flow_index)...)]
      end
    end |> flatten |> collect

    edges_outflow=map(1:length(edgeOutFlows)) do k
      flow_index=edgeOutFlows[k]
      fv=flowVariableIndex(p,flow_index)
      stock_index_outfrom = first(outstock(p,flow_index))
      if occursin("V", type)
        flow=make_flow_V(p,"s$stock_index_outfrom", "fs_$flow_index"*"d", fv, flow_index)
        [Edge(flow[1]...),Edge(flow[2]...)]
      else
        [Edge(make_flow_noneV(p,"s$stock_index_outfrom", "fs_$flow_index"*"d", flow_index)...)]
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
                [Edge(make_link("s$m", "sv$k")...)]
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
                [Edge(make_link("s$m", "v$k")...)]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end

    # linkages from sv to v
    edges_lsvv()=begin
      subEdges = Vector{Edge}[]
      for k in 1:nvb(p)
        subEdges = vcat(map(svsv(p,k)) do m
                 [Edge(make_link("sv$m", "v$k")...)]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end

    # linkages from v to v
    edges_lvv()=begin
      subEdges = Vector{Edge}[]
      for k in 1:nvb(p)
        subEdges = vcat(map(vsrc(p,k)) do m
                 [Edge(make_link("v$m", "v$k")...)]
                end, subEdges)
      end
      return subEdges |> flatten |> collect
    end

    # linkages from p to v ###################################
    edges_lpv()=begin
      subEdges = Vector{Edge}[]
      for k in 1:nvb(p)
        subEdges = vcat(map(vpsrc(p,k)) do m
                 [Edge(make_link("p$m", "v$k")...)]
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
                          vcat(edges_lv(), edges_lsv(), edges_lsvv(), edges_lvv(), edges_lpv())
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

  g = Digraph("G", stmts; graph_attrs=graph_attrs, edge_attrs=edge_attrs)
  return g

end




"""
Graph reinforcing and balancing loops of a causal loop diagram.
"""
function GraphRB(c::Union{CausalLoopPM, CausalLoopPol} ; cycle_color=:yellow, edge_label_color=:lightblue)
    NNodes = [Node("n$n", Attributes(:label=>"$(vname(c, n))",:shape=>"square")) for n in 1:nvert(c)]
    Edges = Vector{Edge}()
    edge_to_intermediate_node = Dict{Int, Int}()

    for k in 1:nedges(c)
      pol_int = epol(c,k)
      if pol_int == POL_POSITIVE
        pol = :+
      elseif pol_int == POL_NEGATIVE
        pol = :-
      else
        error("Unknown Polarity $pol_int.")
      end
      new_node_index = length(NNodes) + 1

      # Graphviz doesn't allow us to have an edge from an edge, so we put a node
      # in the middle of each edge and use that instead.

      # Intermediate node in edge, with label
      push!(NNodes, Node("n$new_node_index", Attributes(:label => "$pol", :penwidth => "0", :style => "filled", :shape => "cds", :fillcolor => "$edge_label_color", )))
      # edge to intermediate
      push!(Edges, Edge(["n$(sedge(c,k))", "n$new_node_index"],Attributes(:color=>"blue", :arrowhead=>"none")))
      # edge from intermediate
      push!(Edges, Edge(["n$new_node_index", "n$(tedge(c,k))"],Attributes(:color=>"blue")))
      # mapping of edge to its intermediate node
      push!(edge_to_intermediate_node, k => new_node_index)
    end

    for (edges, polarity) ∈ extract_loops(c)

      if polarity == POL_POSITIVE
        label = "R" # reinforcing
      elseif polarity == POL_NEGATIVE
        label = "B" # balancing
      else
        error("Unknown cycle polarity $polarity")
      end
      new_node_index = length(NNodes) + 1
      # node to represent cycle polarity
      push!(NNodes, Node("n$new_node_index", Attributes(:label => "$label", :shape => "circle", :fillcolor => "$cycle_color", :style => "filled")))
      for edge ∈ edges
          edge_node = edge_to_intermediate_node[edge]
          push!(Edges, Edge(["n$edge_node", "n$new_node_index"], Attributes(:color=>"$cycle_color", :style=>"dashed")))
      end

    end


    stmts=vcat(NNodes,Edges)

    g = Digraph("G", stmts;graph_attrs=Attributes(:rankdir=>"LR"))
    return g
end






# function collapse_on_RB(c)
#   all_loops = extract_loops(c)
#   all_loop_sequences = [edges for (edges, _) in all_loops]
#   for (edges, polarity) in all_loops
#     if
#     # if an edge is only in one loop, we can collapse it entirely.
#     # If an edge is in multiple loops, we need to keep this subgraph.

#   end
# end
