module Rewrite

export sfrewrite, @rewrite

# regex to search for uncommented print:  [^#]\s*@show

using ...StockFlow

using ..Syntax
import ..Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, 
parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments,
StockArgT, ParamArgT, create_flow_definitions, Binop, Ref, get


using Catlab.ACSets.ACSetInterface
using Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite


using MLStyle
using MLStyle.Modules.AST

"""
Convert a stockflow block to a stockflow
"""
function block_to_sf(block)
  args = stock_and_flow_syntax_to_arguments(block)
  sf = StockAndFlowF(args.stocks, args.params,  
    map(kv -> kv.first => StockFlow.Syntax.get(kv.second), args.dyvars), 
    args.flows, args.sums)
  return sf
end


struct SFPointer
  type::Symbol
  index::Int
end


get_dyvar_args(dv) = (@capture ($v = $f($(args...))) dv)

function add_link_if_not_already!(connect_dict, link, link_object)
  if !(link in connect_dict[link_object])
    push!(connect_dict[link_object], link)
  end
end

function add_part_if_not_already!(check_set, sf, operand, operand_type, operand_index, sf_block)
  if !(operand in check_set)
    push!(check_set, operand)
    add_object_of_type!(sf, operand, operand_type, operand_index, sf_block)
  end
end






function add_operand_link!(sf, operand, operand_type, dyvar_index, operand_index)
  if operand_type == :S
    stock_index = findfirst(==(operand), snames(sf))
    add_part!(sf, :LV ; lvs = stock_index, lvv = dyvar_index, lvsposition = operand_index)
  elseif operand_type == :P
    param_index = findfirst(==(operand), pnames(sf))
    add_part!(sf, :LPV ; lpvp = param_index, lpvv = dyvar_index, lpvpposition = operand_index)
  elseif operand_type == :V
    dyvar2_index = findfirst(==(operand), vnames(sf))
    add_part!(sf, :LVV ; lvsrc = dyvar2_index, lvtgt = dyvar_index, lvsrcposition = operand_index)
  elseif operand_type == :SV
    sum_index = findfirst(==(operand), svnames(sf))
    add_part!(sf, :LSV ; lsvsv = sum_index, lsvv = dyvar_index, lsvsvposition = operand_index)
  end
end

function add_object_of_type!(sf, object, object_type, object_index, sf_block)
  if object_type == :S
    add_stock!(sf ; sname = object)
  elseif object_type == :P
    add_parameter!(sf ; pname = object)
  elseif object_type == :V
    dyvar_op = sf_block.dyvars[object_index][2].args[1]
    add_variable!(sf ; vname = object, vop = dyvar_op)
  elseif object_type == :SV
    add_svariable!(sf ; svname = object)
  elseif object_type == :F
    flow_dyvar = sf_block.flows[object_index][2].args[2]
    flow_dyvar_index = findfirst(==(flow_dyvar), vnames(sf))
    @assert !(isnothing(flow_dyvar_index)) "Tried adding a flow before its dyvar!"
    add_flow!(sf, flow_dyvar_index ; fname = object)
  end
end


function remove_from_sums!(stock_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for sum ∈ sf_block.sums
    sum_name = sum[1]
    sum_stocks = sum[2]
    # for loop because it's technically possible for a stock
    # to be linked to a sum twice.
    # TODO: lmao, replace
    for sum_stock in sum_stocks
      if sum_stock == stock_name
        stock_index = name_dict[stock_name].index
        add_part_if_not_already!(L_set, L, sum_name, :SV, stock_index, sf_block)
        sum_index = findfirst(==(sum_name), svnames(L))
        # ! WE'RE ASSUMING THERE IS AT MAX ONE LINK BETWEEN A SUM AND A STOCK
        link = (stock_name, sum_name)
        add_link_if_not_already!(L_connect_dict, link, :LS)
        add_link_if_not_already!(remove_connect_dict, link, :LS)
      end
    end
  end
end


function dyvar_link_name_from_object_type(object_type)
  if object_type == :S
    return :LV
  elseif object_type == :SV
    return :LSV
  elseif object_type == :P
    return :LPV
  elseif object_type == :V
    return :LVV
  end
end

function add_links_from_dict!(sf, connect_dict)
  for ((lvs, lvv, lvsposition)) in connect_dict[:LV]
    lvs = findfirst(==(lvs), snames(sf))
    lvv = findfirst(==(lvv), vnames(sf))

    add_part!(sf, :LV; lvs = lvs, lvv = lvv, lvsposition = lvsposition)
  end
  for ((lpvp, lpvv, lpvpposition)) in connect_dict[:LPV]
    lpvp = findfirst(==(lpvp), pnames(sf))
    lpvv = findfirst(==(lpvv), vnames(sf))

    add_part!(sf, :LPV; lpvp = lpvp, lpvv = lpvv, lpvpposition = lpvpposition)
  end
  for ((lsvsv, lsvv, lsvsvposition)) in connect_dict[:LSV]
    lsvsv = findfirst(==(lsvsv), svnames(sf))
    lsvv = findfirst(==(lsvv), vnames(sf))

    add_part!(sf, :LSV; lsvsv = lsvsv, lsvv = lsvv, lsvsvposition = lsvsvposition)
  end
  for ((lvsrc, lvtgt, lvsrcposition)) in connect_dict[:LVV]

    lvsrc = findfirst(==(lvsrc), vnames(sf))
    lvtgt = findfirst(==(lvtgt), vnames(sf))

    add_part!(sf, :LVV; lvsrc = lvsrc, lvtgt = lvtgt, lvsrcposition = lvsrcposition)
  end
  for ((lss, lssv)) in connect_dict[:LS]
    lss = findfirst(==(lss), snames(sf))
    lssv = findfirst(==(lssv), svnames(sf))

    add_part!(sf, :LS; lss = lss, lssv = lssv)
  end
  for ((is, ifn)) in connect_dict[:I]
    is = findfirst(==(is), snames(sf))
    ifn = findfirst(==(ifn), fnames(sf))

    add_part!(sf, :I; is = is, ifn = ifn)
  end
  for ((os, ofn)) in connect_dict[:O]
    os = findfirst(==(os), snames(sf))
    ofn = findfirst(==(ofn), fnames(sf))
    add_part!(sf, :O; os = os, ofn = ofn)
  end
end

function remove_from_dyvars!(object_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for dyvar ∈ sf_block.dyvars
    dyvar_name = dyvar[1]
    dyvar_expression = dyvar[2].args
    dyvar_op = dyvar_expression[1]
    for (operand_index, operand) in enumerate(dyvar_expression[2:end])
      if operand == object_name

        object_type = name_dict[object_name].type
      
        add_part_if_not_already!(L_set, L, dyvar_name, :V, name_dict[dyvar_name].index, sf_block)

        dyvar_index = findfirst(==(dyvar_name), vnames(L))

        link_name = dyvar_link_name_from_object_type(object_type)
        link = (object_name, dyvar_name, operand_index)
        add_link_if_not_already!(L_connect_dict, link, link_name)
        add_link_if_not_already!(remove_connect_dict, link, link_name)
      end # if
    end # for oper
  end # for dyvar
end #func



function add_redefintions!(L, L_redef_queue, R_dyvar_queue, R_sum_queue, R_flow_queue, name_dict, L_set, sf_block, L_connect_dict, remove_connect_dict, removed_set)
  for object in L_redef_queue
    if length(object.args) == 2
      object_name = object.args[1]
    elseif length(object.args) == 3 # flow
      object_name = object.args[3].args[2].args[1]
    end
    object_pointer = name_dict[object_name]
    object_type = object_pointer.type
    object_index = object_pointer.index

    if object_type == :V

      dyvar_operands = object.args[2].args
      dyvar_op = dyvar_operands[1]
      push!(R_dyvar_queue, object)

      add_part_if_not_already!(L_set, L, object_name, :V, object_index, sf_block)

      dyvar_index = findfirst(==(object_name), vnames(L))

      for (operand_index, operand) in enumerate(sf_block.dyvars[object_index][2].args[2:end])

        if operand in keys(name_dict)
          # need to add all links to L
          operand_pointer = name_dict[operand]
          operand_type = operand_pointer.type
          operand_original_index = operand_pointer.index

          add_part_if_not_already!(L_set, L, operand, operand_type, operand_original_index, sf_block)



          original_def = sf_block.dyvars[name_dict[object_name].index]

          
          dyvar_link_type = dyvar_link_name_from_object_type(operand_type)

          link = (operand, object_name, operand_index)

          add_link_if_not_already!(L_connect_dict, link, dyvar_link_type)
          # we want to remove all links such that they aren't in the redefinition
          # operand is in old, we want to check if it's in the new

          if (length(object.args[2].args[2:end]) < operand_index || object.args[2].args[operand_index + 1] != operand )
            add_link_if_not_already!(remove_connect_dict, link, dyvar_link_type)
          end

        end # if operand in namedict
      end # for each operand (except operator)

    elseif object_type == :SV
      add_part_if_not_already!(L_set, L, object_name, :SV, name_dict[object_name].index, sf_block)

      push!(R_sum_queue, object)

      sum_index = findfirst(==(object_name), svnames(L))
      for operand in object.args[2].args
        if operand in keys(name_dict)
          add_part_if_not_already!(L_set, L, operand, :S, name_dict[operand].index, sf_block)
          
          stock_index = findfirst(==(operand), snames(L))
          link = (operand, object_name)

          add_link_if_not_already!(L_connect_dict, link, :LS)

          # Don't remove link if it's retained
          original_links = sf_block.sums[name_dict[object_name].index][2]

          if !(operand in original_links)
            add_link_if_not_already!(remove_connect_dict, link, :LS)
          end          
        end
      end
    elseif object_type == :F
      # If the dyvar has changed, need to have it in L, nothing in I, re-add new to R
      original_flow_var = sf_block.flows[name_dict[object_name].index][2].args[2]    
      original_flow_index = name_dict[original_flow_var].index
      add_part_if_not_already!(L_set, L, original_flow_var, :V, original_flow_index, sf_block)
      add_part_if_not_already!(L_set, L, object_name, :F, name_dict[object_name].index, sf_block)


      flow_index = findfirst(==(object_name), fnames(L))

      original_inflow = sf_block.flows[flow_index][3]
      original_outflow = sf_block.flows[flow_index][1]


    

      if (original_inflow != :F_NONE)
        add_part_if_not_already!(L_set, L, original_inflow, :S, name_dict[original_inflow].index, sf_block)
        link = (original_inflow, object_name)
        add_link_if_not_already!(L_connect_dict, link, :I)
        add_link_if_not_already!(remove_connect_dict, link, :I)
      end 

      if (original_outflow != :F_NONE)
        add_part_if_not_already!(L_set, L, original_outflow, :S, name_dict[original_outflow].index, sf_block)
        link = (original_outflow, object_name)
        add_link_if_not_already!(L_connect_dict, link, :O)
        add_link_if_not_already!(remove_connect_dict, link, :O)
      end 

      # different dynamic variables
      if !(original_flow_var == object.args[3].args[2].args[2])
      
        push!(removed_set, object_name)

        if (original_inflow != :F_NONE)
          link = (original_inflow, object_name)
          add_link_if_not_already!(L_connect_dict, link, :I)
          add_link_if_not_already!(remove_connect_dict, link, :I)
        end
        if (original_outflow != :F_NONE)
          link = (original_outflow, object_name)
          add_link_if_not_already!(L_connect_dict, link, :O)
          add_link_if_not_already!(remove_connect_dict, link, :O)
        end
      end

      # needs to have dyvar added...
      push!(R_flow_queue, object)
      # end


    end 


  end
end




function remove_from_flows!(object_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for flow in sf_block.flows
    inflow_stock = flow[3]
    outflow_stock = flow[1]
    flow_name = flow[2].args[1]
    flow_dyvar = flow[2].args[2]

    link = (object_name, flow_name)

    if (inflow_stock == object_name) || (outflow_stock == object_name)
      add_part_if_not_already!(L_set, L, flow_dyvar, :V, name_dict[flow_dyvar].index, sf_block)
      add_part_if_not_already!(L_set, L, flow_name, :F, name_dict[flow_name].index, sf_block)
    end

    if (inflow_stock == object_name)
      add_link_if_not_already!(L_connect_dict, link, :I)
      add_link_if_not_already!(remove_connect_dict, link, :I)
    end

    if (outflow_stock == object_name)
      add_link_if_not_already!(L_connect_dict, link, :O)
      add_link_if_not_already!(remove_connect_dict, link, :O)
    end


  end

end


function remove_block!(removed, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)

  if removed isa Expr && removed.args[1] == :~
    for k in keys(name_dict)
      if  occursin(String(removed.args[2]), String(k))
        remove_block!(k, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
      end
    end
    return
  end
  
  object_pointer = name_dict[removed]
  object_type = object_pointer.type
  push!(removed_set, removed)

  add_part_if_not_already!(L_set, L, removed, object_type, object_pointer.index, sf_block)



  #TODO: factor out

  if object_type == :S
    object_definition = sf_block.stocks[object_pointer.index]
    remove_from_sums!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    remove_from_dyvars!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    remove_from_flows!(removed, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    

  elseif object_type == :P
    object_definition = sf_block.params[object_pointer.index]
    remove_from_dyvars!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  elseif object_type == :V
    object_definition = sf_block.dyvars[object_pointer.index]
    vname = object_definition[1]

    for (operand_index, operand) in enumerate(object_definition[2].args[2:end])

      add_part_if_not_already!(L_set, L, operand, name_dict[operand].type, name_dict[operand].index, sf_block)

      operand_link_type = dyvar_link_name_from_object_type(name_dict[operand].type)
      link = (operand, removed, operand_index)
      add_link_if_not_already!(L_connect_dict, link, operand_link_type)
      add_link_if_not_already!(remove_connect_dict, link, operand_link_type)
    end

    remove_from_dyvars!(vname, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  elseif object_type == :SV
    object_definition = sf_block.sums[object_pointer.index]
    svname = object_definition[1]

    for sum_stock in object_definition[2]

      add_part_if_not_already!(L_set, L, sum_stock, :S, name_dict[sum_stock].index, sf_block)

      link = (sum_stock, removed)
      add_link_if_not_already!(L_connect_dict, link, :LS)
      add_link_if_not_already!(remove_connect_dict, link, :LS)
    end

    remove_from_dyvars!(svname, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  elseif object_type == :F
    flow = sf_block.flows[object_pointer.index]
    inflow_stock = flow[3]
    outflow_stock = flow[1]
    flow_name = flow[2].args[1]
    flow_dyvar = flow[2].args[2]

    inflow_link = (inflow_stock, flow_name)
    outflow_link = (outflow_stock, flow_name)

    add_part_if_not_already!(L_set, L, flow_dyvar, :V, name_dict[flow_dyvar].index, sf_block)

   

    if !(inflow_stock == :F_NONE)
      add_part_if_not_already!(L_set, L, inflow_stock, :S, name_dict[inflow_stock].index, sf_block)
      add_link_if_not_already!(L_connect_dict, inflow_link, :I)
      add_link_if_not_already!(remove_connect_dict, inflow_link, :I)
    end


    if !(inflow_stock == :F_NONE)
      add_part_if_not_already!(L_set, L, outflow_stock, :S, name_dict[outflow_stock].index, sf_block)
      add_link_if_not_already!(L_connect_dict, outflow_link, :O)
      add_link_if_not_already!(remove_connect_dict, outflow_link, :O)
    end

  end
end



function remove_links_block!(l, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
  @match l begin
    
    :($(src::Symbol) => $(tgt::Symbol)) => begin
      tgt_type = name_dict[tgt].type

      add_part_if_not_already!(L_set, L, tgt, tgt_type, name_dict[tgt].index, sf_block)
      add_part_if_not_already!(L_set, L, src, name_dict[src].type, name_dict[src].index, sf_block)
      

      if tgt_type == :SV
        link = (src, tgt)
        add_link_if_not_already!(L_connect_dict, link, :LS)
        add_link_if_not_already!(remove_connect_dict,  link, :LS)
      elseif tgt_type == :V # tgt must be dyvar
        positions = findall(==(src), sf_block.dyvars[name_dict[tgt].index][2].args[2:end])
        link_type = dyvar_link_name_from_object_type(name_dict[src].type)
        for pos in positions
          link = (src, tgt, pos)
          add_link_if_not_already!(L_connect_dict, link, link_type)
          # TODO: note, in previous version, remove_connect_dict didn't have a check that it wasn't already in
          # I don't think it matters that much either way, really.
          add_link_if_not_already!(remove_connect_dict, link, link_type)
        end
      elseif tgt_type == :F
        tgt_pointer = name_dict[tgt]
        flow_definintion = sf_block.flows[tgt_pointer.index]
        flow_dyvar = flow_definintion[2].args[2]


        add_part_if_not_already!(L_set, L, flow_dyvar, :V, name_dict[flow_dyvar].index, sf_block)

        if flow_definintion[3] == src
          link = (src, tgt)
          add_link_if_not_already!(L_connect_dict, link, :I)
          # also didn't have the check for this one...
          add_link_if_not_already!(remove_connect_dict, link, :I)
        end
        if flow_definintion[1] == src
          add_link_if_not_already!(L_connect_dict, link, :O)
          # ...or this one
          add_link_if_not_already!(remove_connect_dict, link, :O)
        end
      end # if tgt_type == 
    end # :($(src:: ... => begin

    :($(src::Symbol) => $(tgt::Symbol), $(position::Int)) => begin
      tgt_type = name_dict[tgt].type


      add_part_if_not_already!(L_set, L, tgt, tgt_type, name_dict[tgt].index, sf_block)
      add_part_if_not_already!(L_set, L, src, name_dict[src].type, name_dict[src].index, sf_block)

      if tgt_type == :SV
        link = (src, tgt)
        add_link_if_not_already!(L_connect_dict, link, :LS)
        add_link_if_not_already!(remove_connect_dict, link, :LS)
      
      elseif tgt_type == :V
        link_type = dyvar_link_name_from_object_type(name_dict[src].type)
        link = (src, tgt, position)
        add_link_if_not_already!(L_connect_dict, link, link_type)
        add_link_if_not_already!(remove_connect_dict, link, link_type)
      elseif tgt_type == :F
        link = (src, tgt)
        if position == 2
          add_link_if_not_already!(L_connect_dict, link, :I)
          add_link_if_not_already!(remove_connect_dict, link, :I)
        elseif position == 1
          add_link_if_not_already!(L_connect_dict, link, :O)
          add_link_if_not_already!(remove_connect_dict, link, :O)
        end
      end #if tgt_type              

      # TODO: error case

    end # ... $(position::Int)) => begin
  end # match
end # current_phase



function add_links_block!(l, L_set, name_dict, L, sf_block, R_link_vector)
  src = nothing
  tgt = nothing
  position = nothing
  @match l begin 
    :($(srca::Symbol) => $(tgta::Symbol), $(positiona::Int)) => begin
      src = srca
      tgt = tgta
      position = positiona
    end
    :($(srca::Symbol) => $(tgta::Symbol)) => begin
      src = srca
      tgt = tgta
      position = 0
    end
  end
  if src in keys(name_dict)
    add_part_if_not_already!(L_set, L, src, name_dict[src].type, name_dict[src].index, sf_block)
  end
  if tgt in keys(name_dict)
    add_part_if_not_already!(L_set, L, tgt, name_dict[tgt].type, name_dict[tgt].index, sf_block)
  end
  push!(R_link_vector, (src => tgt, position))      
end 


function dyvar_swaps_block!(dw, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict, R_link_vector)
  capture_dict = (@capture $old > $new dw)
  new = capture_dict[:new]
  old = capture_dict[:old]


  add_part_if_not_already!(L_set, L, old, name_dict[old].type, name_dict[old].index, sf_block)
  

  for ((dyvar_name, dyvar_expr)) in sf_block.dyvars
    dyvar_op = dyvar_expr.args[1]
    dyvar_operands = dyvar_expr.args[2:end]
    matching_indices = findall(==(old), dyvar_operands)
    if !(isempty(matching_indices))
      add_part_if_not_already!(L_set, L, dyvar_name, :V, name_dict[dyvar_name].index, sf_block)

      for index in matching_indices
        old_link = (old, dyvar_name, index)

        old_link_type = dyvar_link_name_from_object_type(name_dict[old].type)
        add_link_if_not_already!(L_connect_dict, old_link, old_link_type)
        add_link_if_not_already!(remove_connect_dict, old_link, old_link_type)


        
        if new in keys(name_dict)
          add_part_if_not_already!(L_set, L, new, name_dict[new].type, name_dict[new].index, sf_block)
        end

        push!(R_link_vector, (new => dyvar_name, index))

      end
    end

  end
end





function sfrewrite(sf::K, block::Expr) where {K <: AbstractStockAndFlowF}

  Base.remove_linenums!(block)
  name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]
  @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"

  # Unfortunately, this makes it un-extendable to general schemas...it'll only work for StockAndFlowF...
  sf_block::StockAndFlowBlock = sf_to_block(sf)

  L = StockAndFlowF()

  name_dict = Dict{Symbol, SFPointer}([
    [s => SFPointer(:S, i) for (i, s) ∈ enumerate(snames(sf))] ;
    [sv => SFPointer(:SV, i) for (i, sv) ∈ enumerate(svnames(sf))] ;
    [v => SFPointer(:V, i) for (i, v) ∈ enumerate(vnames(sf))] ;
    [f => SFPointer(:F, i) for (i, f) ∈ enumerate(fnames(sf))] ;
    [p => SFPointer(:P, i) for (i, p) ∈ enumerate(pnames(sf))]
  ])

  # symbols added to L
  L_set = Set{Symbol}()

  L_connect_dict = Dict{Symbol, Vector{Tuple}}(
    :LV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LPV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LVV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LSV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LS => Vector{Tuple{Symbol, Symbol}}(),
    :I => Vector{Tuple{Symbol, Symbol}}(),
    :O => Vector{Tuple{Symbol, Symbol}}(),
  )

  remove_connect_dict = Dict{Symbol, Vector{Tuple}}(
    :LV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LPV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LVV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LSV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LS => Vector{Tuple{Symbol, Symbol}}(),
    :I => Vector{Tuple{Symbol, Symbol}}(),
    :O => Vector{Tuple{Symbol, Symbol}}(),
  )

  L_redef_queue = Vector{Expr}()

  # used to construct I.  Take L_set - removed_set
  removed_set = Set{Symbol}()

  R_stock_queue = Vector{Symbol}()
  R_param_queue = Vector{Symbol}()
  R_sum_queue = Vector{Expr}()
  R_dyvar_queue = Vector{Expr}()
  R_flow_queue = Vector{Expr}()

  R_link_vector = Vector{Tuple}()


  current_phase = (_, _) -> ()
  for statement in block.args
    @match statement begin
      QuoteNode(:removes) => begin
        current_phase = removed -> remove_block!(removed, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
      end

      # TODO: figure out if I remove links at this stage.  If I do, then remove_links needs to come after
      # all the parts which are defined in it
      QuoteNode(:remove_links) => begin
       current_phase = l -> remove_links_block!(l, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
      end # quotenode

      QuoteNode(:add_links) => begin
        current_phase = l -> add_links_block!(l, L_set, name_dict, L, sf_block, R_link_vector)
      end
      QuoteNode(:dyvar_swaps) => begin
        current_phase = dw -> dyvar_swaps_block!(dw, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict, R_link_vector)
      end



      QuoteNode(:redefs) => begin
        current_phase = r -> push!(L_redef_queue, r)
      end
      
      QuoteNode(:stocks) => begin
        current_phase = s -> push!(R_stock_queue, s)
      end
      QuoteNode(:parameters) => begin
        current_phase = p -> push!(R_param_queue, p)
      end
      QuoteNode(:dynamic_variables) => begin
        current_phase = d -> push!(R_dyvar_queue, d)
      end
      QuoteNode(:flows) => begin
        current_phase = f -> push!(R_flow_queue, f)
      end
      QuoteNode(:sums) => begin
        current_phase = s -> push!(R_sum_queue, s)
      end
      QuoteNode(kw) =>
        error("Unknown block type for rewrite syntax: " * String(kw))
      _ => current_phase(statement)

    end # match
  end # for


  # extract links from redefinitions

  add_redefintions!(L, L_redef_queue, R_dyvar_queue, R_sum_queue, R_flow_queue, name_dict, L_set, sf_block, L_connect_dict, remove_connect_dict, removed_set)

  add_links_from_dict!(L, L_connect_dict)
  



  parsed_flows = Vector{Tuple{Symbol, Expr, Symbol}}(parse_flow.(R_flow_queue))
  updated_flows, new_dyvars = create_flow_definitions(parsed_flows, vcat(map(x -> first(x.args), R_dyvar_queue), vnames(L)))
  new_dyvar_names, new_dyvar_parsedform = zip(new_dyvars...)
  new_dyvar_unwrapped = collect(zip(new_dyvar_names, map(Syntax.get, new_dyvar_parsedform)))
    # ! TODO: Make this not suck!!!

  new_dyvar_definitions = map(kv -> kv[2][1] isa Tuple ? Expr(:(=), kv[1], Expr(:call, kv[2][2], kv[2][1][1], (kv[2][1][2]))) : :Expr(:(=), kv[1], Expr(:call, kv[2][2], (kv[2][1]))),  new_dyvar_unwrapped)


  for (i, flow) in enumerate(R_flow_queue)
    flow.args[3].args[2].args = [updated_flows[i][1], updated_flows[i][2]]
  end



 

  append!(R_dyvar_queue, new_dyvar_definitions)
  

  flow_stocks_out, _, flow_stocks_in = zip(parsed_flows...)


  for dyvar in R_dyvar_queue
    for operand in dyvar.args[2].args[2:end]
      if operand in keys(name_dict)
        add_part_if_not_already!(L_set, L, operand, name_dict[operand].type, name_dict[operand].index, sf_block)
      end
    end
  end


  for sum in R_sum_queue
    for stock_sum in sum.args[2].args
      if stock_sum in keys(name_dict)
        add_part_if_not_already!(L_set, L, stock_sum, name_dict[stock_sum].type, name_dict[stock_sum].index, sf_block)
      end
    end
  end
  for flow in R_flow_queue
    outflow = flow.args[2]
    inflow = flow.args[3].args[3]
    dyvar = flow.args[3].args[2].args[2]
   
    
    for value in [outflow, inflow, dyvar]
      if value in keys(name_dict)
        add_part_if_not_already!(L_set, L, value, name_dict[value].type, name_dict[value].index, sf_block)
      end
    end
  end

  
  # L should now be complete
  

  I_set = setdiff(L_set, removed_set)
  I_connect_dict = Dict(key => [value for value in filter(v -> !(v in remove_connect_dict[key]), L_connect_dict[key])] for key in keys(L_connect_dict))
  I = StockAndFlowF()
  
  I_set_flows = filter(x -> name_dict[x].type == :F, I_set)


  # TODO: This isn't preserving any links.  Verify if I need to...
  # probably need to preserve inflows and outflows, at least...
  for object in I_set
    if object in I_set_flows
      continue
    end
    object_pointer = name_dict[object]
    object_type = object_pointer.type
    object_index = object_pointer.index

    add_object_of_type!(I, object, object_type, object_index, sf_block)

  end

  for object in I_set_flows
    add_object_of_type!(I, object, :F, name_dict[object].index, sf_block)
  end


  add_links_from_dict!(I, I_connect_dict)

  R = deepcopy(I)
  R_name_dict = Dict{Symbol, SFPointer}([
    [s => SFPointer(:S, i) for (i, s) ∈ enumerate(snames(R))] ;
    [sv => SFPointer(:SV, i) for (i, sv) ∈ enumerate(svnames(R))] ;
    [v => SFPointer(:V, i) for (i, v) ∈ enumerate(vnames(R))] ;
    [f => SFPointer(:F, i) for (i, f) ∈ enumerate(fnames(R))] ;
    [p => SFPointer(:P, i) for (i, p) ∈ enumerate(pnames(R))]
  ])

  not_yet_added_stocks = collect(filter(x -> !(x in snames(R)), R_stock_queue))
  not_yet_added_params = collect(filter(x -> !(x in pnames(R)), R_param_queue))

  

  type_index_counter = ns(R) + 1
  (x -> (push!(R_name_dict, x => SFPointer(:S, type_index_counter)); type_index_counter+=1)).(not_yet_added_stocks)

  type_index_counter = np(R) + 1
  (x -> (push!(R_name_dict, x => SFPointer(:P, type_index_counter)); type_index_counter+=1)).(not_yet_added_params)



  add_stocks!(R, length(not_yet_added_stocks) ; sname = not_yet_added_stocks)
  add_parameters!(R, length(not_yet_added_params) ; pname = not_yet_added_params)



  sum_names, sum_lists = zip(parse_sum.( R_sum_queue)...)
  
  dyvar_names, dyvar_definitions = zip(parse_dyvar.(R_dyvar_queue)...)
  dyvar_definitions = collect(dyvar_definitions)
  dyvar_names = collect(dyvar_names)
  dyvar_ops = collect((x -> x.args[1]).(dyvar_definitions))












  # add dyvars and sums to R_name_dict, so we know their index when we need to
  # link them to a dyvar
  type_index_counter = nvb(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:V, type_index_counter))); type_index_counter+=1)).(dyvar_names)

  type_index_counter = nsv(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:SV, type_index_counter))); type_index_counter+=1)).(sum_names)


  filtered_dyvars = collect(filter(x -> !(x.args[1] in vnames(R)), R_dyvar_queue))
  
  filtered_dyvar_names = map(x -> first(x.args), filtered_dyvars)

  filtered_dyvar_ops = map(x -> x.args[2].args[1], filtered_dyvars)


  filtered_sum_names = map(x -> x.args[1], collect(filter(x -> !(x.args[1] in svnames(R)), R_sum_queue  )    ))

  add_svariables!(R, length(filtered_sum_names), svname = filtered_sum_names)
  add_variables!(R, length(filtered_dyvar_names), vname = filtered_dyvar_names, vop = filtered_dyvar_ops)


  filtered_flow_names = filter(x -> !(x in fnames(R)), first.(updated_flows))
  flow_variables_indices = map(x -> findfirst(==(x), vnames(R)), map(x -> x[2], updated_flows))

  add_flows!(R, flow_variables_indices, length(filtered_flow_names), fname = filtered_flow_names)

 

  for sum_index in eachindex(R_sum_queue)
    real_sum_index = findfirst(==(R_sum_queue[sum_index].args[1]), svnames(R))
    for stock in sum_lists[sum_index]
      stock_index = findfirst(==(stock), snames(R))
      link = (stock, R_sum_queue[sum_index].args[1])
      if !(link in I_connect_dict[:LS])
        add_part!(R, :LS ; lss = stock_index, lssv = real_sum_index)
      end
    end
  end


  for dyvar_index in eachindex(R_dyvar_queue)
    real_dyvar_index = findfirst(==(R_dyvar_queue[dyvar_index].args[1]), vnames(R))
    dyvar_definition = R_dyvar_queue[dyvar_index].args[2]
    for (operand_index, operand) in enumerate(dyvar_definition.args[2:end])
      link = (operand, R_dyvar_queue[dyvar_index].args[1], operand_index)
      if !(link in I_connect_dict[dyvar_link_name_from_object_type(R_name_dict[operand].type)])
        add_operand_link!(R, operand, R_name_dict[operand].type, real_dyvar_index, operand_index)
      end
    end
  end

  no_flow_names = Set([:CLOUD, :☁])
  for flow_index in eachindex(R_flow_queue)
    stock_in = flow_stocks_in[flow_index]
    if !(stock_in in no_flow_names)
      stock_in_index = R_name_dict[stock_in].index
      link = (stock_in, R_flow_queue[flow_index].args[3].args[2].args[1])
      if !(link in I_connect_dict[:I])
        add_part!(R, :I ; is = stock_in_index, ifn = flow_index)
      end
    end

    stock_out = flow_stocks_out[flow_index]
    if !(stock_out in no_flow_names)
      stock_out_index = R_name_dict[stock_out].index
      link = (stock_out, R_flow_queue[flow_index].args[3].args[2].args[1])
      if !(link in I_connect_dict[:O])
        add_part!(R, :O ; os = stock_out_index, ofn = flow_index)
      end
    end 
  end


  for link in R_link_vector
    src = link[1][1]
    tgt = link[1][2]
    src_pointer = R_name_dict[src]
    tgt_pointer = R_name_dict[tgt]

    src_index = src_pointer.index
    tgt_index = tgt_pointer.index

    tgt_type = tgt_pointer.type

    if tgt_type == :SV
      add_part!(R, :LS ; lss = src_index, lssv = tgt_index)
    elseif tgt_type == :F
      position = link[2]
      if position == 1
        add_part!(R, :O ; is = src_index, ifn = tgt_index)
      elseif position == 2
        add_part!(R, :I ; os = src_index, ofn = tgt_index)
      end
    elseif tgt_type == :V
      add_operand_link!(R, src, src_pointer.type, tgt_index, link[2])
    end


  end

 

  hom1 = homomorphism(I,L)
  hom2 = homomorphism(I,R)


  rule = Rule(hom1, hom2)


  sf_rewritten = rewrite(rule, sf)


  if isnothing(sf_rewritten)
    @show sf
    @show homomorphism(L, sf)
    @show hom1
    @show hom2

    @show L
    @show I
    @show R
    return error("Failed to apply rule!")
  end
  return sf_rewritten

end





"""
Given a stockflow and a syntax block, create a new stockflow by using a rewrite
rule to apply modified homomorphisms to it.

Define three new blocks L, I and R, such that sf ⊇ L ⊇ I and R ⊇ I.  Define
homomorphisms I -> L and I -> R.  Apply this change to the original stockflow.

Use :stocks, :flows, :sums, :dynamic_variables and :parameters to add those
particular objects.  Same format as @stock_and_flow.

Use :dyvar_swaps to replace all instances of an object in a dynamic variable
with another object.  Does not delete the original object.

Use :redefs to provide a new definition for an existing sum or dyvar.

Use :removes to delete objects.

```julia
@rewrite aged_sir begin

  :redefs
  v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC_post
  v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA_post


  :parameters
  fcc
  fca
  fac
  faa

  :dynamic_variables

  v_CCContacts = fcc * v_prevalencev_INC
  v_CAContacts = fca * v_prevalencev_INA
  
  v_ACContacts = fac * v_prevalencev_INC
  v_AAContacts = faa * v_prevalencev_INA
  
  v_prevalencev_INC_post = v_CCContacts + v_CAContacts
  v_prevalencev_INA_post = v_ACContacts + v_AAContacts

end
```
"""
macro rewrite(sf, block)
  escaped_block = Expr(:quote, block)
  sf = esc(sf)
  quote
    sfrewrite($sf, $escaped_block)
  end
end
  


end 