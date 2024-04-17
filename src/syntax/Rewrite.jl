module Rewrite

export sfrewrite, @rewrite



using ...StockFlow

using ..Syntax
import ..Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, 
parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments,
StockArgT, ParamArgT


using Catlab.ACSets.ACSetInterface
using Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite


using MLStyle

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


"""
Takes two sfs, the latter a subset of the former, and makes the positions of the latter match those of the former.
Assumes that the pair of morphisms which identify each stock, param, etc as linking to a dyvar are unique.
If they aren't, then only the first is reset.
"""
function reset_positions!(sfold, sfnew)
  # This only resets the first match.
  # Unless we're dealing with equations with more than two variables, this should be
  # all that's necessary.


  # 1: select all links to dv which have a position:
  oldlv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lvs(sfold), get_lvv(sfold), get_lvsposition(sfold)))
  oldlsv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lsvsv(sfold), get_lsvv(sfold), get_lsvsvposition(sfold)))
  oldlvv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lvsrc(sfold), get_lvtgt(sfold), get_lvsrcposition(sfold)))
  oldlpv = Dict((i, j) => k for (i, j, k) ∈ zip(get_lpvp(sfold), get_lpvv(sfold), get_lpvpposition(sfold)))



  
  newlv = zip(get_lvs(sfnew), get_lvv(sfnew))
  newlsv = zip(get_lsvsv(sfnew), get_lsvv(sfnew))
  newlvv = zip(get_lvsrc(sfnew), get_lvtgt(sfnew))
  newlpv = zip(get_lpvp(sfnew), get_lpvv(sfnew))


  newsnames = snames(sfnew)
  newpnames = pnames(sfnew)
  newvnames = vnames(sfnew)
  newsvnames = svnames(sfnew)

  oldsnames = snames(sfold)
  oldpnames = pnames(sfold)
  oldvnames = vnames(sfold)
  oldsvnames = svnames(sfold)

  used_lvs = Set{Tuple{Int, Int}}()
  used_lsvsv = Set{Tuple{Int, Int}}()
  used_lvsrc = Set{Tuple{Int, Int}}()
  used_lpvp = Set{Tuple{Int, Int}}()
  
  
  #2: iterate over all new values, find what they link to in new, use that to find the corresponding symbol in old


  # could do zero arrays, but doing it this way ensures positions for new objects aren't overwritten
  restored_lvsposition = subpart(sfnew, :lvsposition)
  restored_lsvsvposition = subpart(sfnew, :lsvsvposition)
  restored_lvsrcposition = subpart(sfnew, :lvsrcposition)
  restored_lpvpposition = subpart(sfnew, :lpvpposition)

  
  
  for (lv, (lvs, lvv)) ∈ enumerate(newlv)
    n1 = newsnames[lvs]
    n2 = newvnames[lvv]

    if n1 ∉ oldsnames || n2 ∉ oldvnames
      continue
    end
    old = (findfirst(==(n1), oldsnames), findfirst(==(n2), oldvnames))
    if old ∉ keys(oldlv) || old ∈ used_lvs
      continue
    end
    push!(used_lvs, old) # only first match will be set
    old_pos = oldlv[old]
    restored_lvsposition[lv] = old_pos
  end

    
  for (lsv, (lsvsv, lsvv)) ∈ enumerate(newlsv)
    n1 = newsvnames[lsvsv]
    n2 = newvnames[lsvv]

    if n1 ∉ oldsvnames || n2 ∉ oldvnames
      continue
    end
    old = (findfirst(==(n1), oldsvnames), findfirst(==(n2), oldvnames))
    if old ∉ keys(oldlsv) || old ∈ used_lsvsv
      continue
    end
    push!(used_lsvsv, old)
    old_pos = oldlsv[old]
    restored_lsvsvposition[lsv] = old_pos
  end

  
  for (lvv, (lvsrc, lvtgt)) ∈ enumerate(newlvv)
       n1 = newvnames[lvsrc]
    n2 = newvnames[lvtgt]

    if n1 ∉ oldvnames || n2 ∉ oldvnames
      continue
    end

    old = (findfirst(==(n1), oldvnames), findfirst(==(n2), oldvnames))

    if old ∉ keys(oldlvv) || old ∈ used_lvsrc
      continue
    end
    push!(used_lvsrc, old)
    old_pos = oldlvv[old]
    restored_lvsrcposition[lvv] = old_pos
  end

  for (lpv, (lpvp, lpvv)) ∈ enumerate(newlpv)
      n1 = newpnames[lpvp]
    n2 = newvnames[lpvv]

    if n1 ∉ oldpnames || n2 ∉ oldvnames
      continue
    end
    old = (findfirst(==(n1), oldpnames), findfirst(==(n2), oldvnames))
    if old ∉ keys(oldlpv) || old ∈ used_lpvp
      continue
    end
    push!(used_lpvp, old)
    old_pos = oldlpv[old]
    restored_lpvpposition[lpv] = old_pos
  end

  set_subpart!(sfnew, :lvsposition, restored_lvsposition)
  set_subpart!(sfnew, :lsvsvposition, restored_lsvsvposition)
  set_subpart!(sfnew, :lvsrcposition, restored_lvsrcposition)
  set_subpart!(sfnew, :lpvpposition, restored_lpvpposition)

  return sfnew
end

"""
Given a stockflow block, an index, and the type of object one is retrieving,
return the corresponding object of that type at the index.
"""
function get_from_correct_vector(block, index, object_type)
  if object_type == :S
    return block.stocks[index]
  elseif object_type == :P
    return block.params[index]
  elseif object_type == :F
    return block.flows[index]
  elseif object_type == :V
    return block.dyvars[index]
  elseif object_type == :SV
    return block.sums[index]
  end
end


"""
Given a stockflow block, an object definition and the type of object,
add that object definition to the vector corresponding to that object type
"""
function add_to_correct_vector!(block, object_definition, object_type)
  if object_type == :S
    push!(block.stocks, object_definition)
  elseif object_type == :P
    push!(block.params, object_definition)
  elseif object_type == :F
    push!(block.flows, object_definition)
  elseif object_type == :V
    push!(block.dyvars, object_definition)
  elseif object_type == :SV
    push!(block.sums, object_definition)
  end
end


"""
Adds sum and all stocks to L if not already there
"""
function remove_from_sums!(stock_name, sf_block, L_block, L_set, name_dict)
  for sum ∈ sf_block.sums
    sum_name = sum[1]
    sum_stocks = sum[2]
    if stock_name ∈ sum_stocks && sum_name ∉ L_set
        push!(L_block.sums, sum)
        push!(L_set, sum_name)
    end

  end
end

"""
Given a dyvar, recursively add all dyvars contained in it to L.
"""
function recursively_add_dyvars_L!(current_dyvar, sf_block, new_block, new_set, name_dict)
  dyvar_name = current_dyvar[1]
  dyvar_copy = deepcopy(current_dyvar)
  dyvar_definition = dyvar_copy[2]
    
  if dyvar_name ∉ new_set
    push!(new_block.dyvars, dyvar_copy)
    push!(new_set, dyvar_name => dyvar_copy)
  end
    
  # For everything not yet in the dict, add to dict
  for object ∈ filter(∉(new_set), dyvar_definition.args[2:end])
    object_type = name_dict[object][1]
    object_index = name_dict[object][2]
    object_definition = get_from_correct_vector(sf_block, object_index, object_type)
    push!(new_set, object)
    add_to_correct_vector!(new_block, object_definition, object_type)


    if object_type == :V
      recursively_add_dyvars_L!(object_definition, sf_block, new_block, new_set, name_dict)
    end
    
  end
  
end
        
        

"""
For every dyvar in the original stockflow which contains a link to object_name:
- If it hasn't been already, add the original dyvar definition to L's dyvars
- Add all objects which are a part of the dyvar to L. 
  (Deal with recursively adding dyvars at end)
"""
function remove_from_dyvars!(object_name, sf_block, L_block, L_set, name_dict)
  for dyvar ∈ sf_block.dyvars
    dyvar_name = dyvar[1]
    dyvar_expression = dyvar[2]
    if object_name ∉ dyvar_expression.args || dyvar_name ∈ L_set
      continue
    end
  
    push!(L_block.dyvars, dyvar)
    push!(L_set, dyvar_name)
    for particular_object ∈ dyvar_expression.args[2:end]
      if particular_object ∉ L_set
        

        object_type = name_dict[particular_object][1]
        object_index = name_dict[particular_object][2]
        object_definition = get_from_correct_vector(sf_block, object_index, object_type)

          push!(L_set, particular_object)
          add_to_correct_vector!(L_block, object_definition, object_type)
      end
    end
  end
end

"""
For every dyvar in the original stockflow which contains a link to src:
- If it hasn't been already, add the original dyvar definition to L's dyvars
- Add all objects which are a part of the dyvar to L. 
  (Deal with recursively adding dyvars at end)
- Create a new dyvar with src replaced with tgt and add to R
"""
function swap_from_dyvars!(src, tgt, sf_block, L_block, L_set, R_block, R_dict, name_dict)
  # check every dyvar to see if it contains src 
  for dyvar ∈ sf_block.dyvars
    dyvar_name = dyvar[1]
    dyvar_expression = dyvar[2]
    if src ∈ dyvar_expression.args
      if dyvar_name ∉ L_set
        push!(L_block.dyvars, dyvar)
        push!(L_set, dyvar_name)
        for particular_object ∈ dyvar_expression.args[2:end]
          if particular_object ∉ L_set
            
            
            object_type = name_dict[particular_object][1]
            object_index = name_dict[particular_object][2]
            object_definition = get_from_correct_vector(sf_block, object_index, object_type)
            push!(L_set, particular_object)
            add_to_correct_vector!(L_block, object_definition, object_type)
          end
        end
      end

      if tgt ∈ keys(name_dict) && tgt ∉ L_set
        tgt_type = name_dict[tgt][1]
        tgt_index = name_dict[tgt][2]
        tgt_definition = get_from_correct_vector(sf_block, tgt_index, tgt_type)
        add_to_correct_vector!(L_block, tgt_definition, tgt_type)
        push!(L_set, tgt)
      end


      if dyvar_name ∉ keys(R_dict)
        new_dyvar = deepcopy(dyvar)
        push!(R_block.dyvars, new_dyvar)
        push!(R_dict, dyvar_name => new_dyvar)
      else
        new_dyvar = R_dict[dyvar_name]
      end
      replace!(new_dyvar[2].args, src => tgt)
    end
  end
end

"""
Given a dyvar, recursively add all dyvars contained in it to R.
"""
function recursively_add_dyvars_R!(current_dyvar, sf_block, new_block, new_dict, name_dict)
  dyvar_name = current_dyvar[1]
  dyvar_copy = deepcopy(current_dyvar)
  dyvar_definition = dyvar_copy[2]
  if dyvar_name ∉ keys(new_dict)
    push!(new_block.dyvars, dyvar_copy)
    push!(new_dict, dyvar_name => dyvar_copy)
  end

  # For everything not yet in the dict, add to dict
  for object ∈ filter(∉(keys(new_dict)), dyvar_definition.args[2:end])
    object_type = name_dict[object][1]
    object_index = name_dict[object][2]
    object_definition = get_from_correct_vector(sf_block, object_index, object_type)
    if object_definition isa Symbol  
      push!(new_dict, object => (object_definition,))
    else
      push!(new_dict, object => object_definition)
    end

    add_to_correct_vector!(new_block, object_definition, object_type)


    if object_type == :V
      recursively_add_dyvars_R!(object_definition, sf_block, new_block, new_dict, name_dict)
    end
    
  end
  
end

struct SFPointer
  type::Symbol
  index::Int
end



function remove_from_sums2!(stock_name, sf_block, L, L_set, name_dict)
  for sum ∈ sf_block.sums
    sum_name = sum[1]
    sum_stocks = sum[2]
    # for loop because it's technically possible for a stock
    # to be linked to a sum twice.
    for sum_stock in sum_stocks
      if sum_stock == stock_name
        stock_index = name_dict[stock_name].index
        if sum_name ∉ L_set
          add_part!(L, :SV ; svname = sum_name)
          # add_part!(L, :LS ; lvs = stock_index, lvv = nsv(L))
          sum_index = nsv(L)
        else
          sum_index = findfirst(==(sum_name), svnames(L))
        end
        add_part!(L, :LS ; lvs = stock_index, lvv = sum_index)
      end
    end
  end
end

function remove_from_dyvars2!(object_name, sf_block, L, L_set, name_dict)
  for dyvar ∈ sf_block.dyvars
    dyvar_name = dyvar[1]
    dyvar_expression = dyvar[2].args
    dyvar_op = dyvar_expression[1]
    for (operand_index, operand) in enumerate(dyvar_expression[2:end])
      if operand == object_name

        object_type = name_dict[object_name].type
      

        if dyvar_name ∉ L_set
          add_part!(L, :V ; vname = dyvar_name, vop = dyvar_op)
          push!(L_set, dyvar_name)
          dyvar_index = nvb(L)
        else
          dyvar_index = findfirst(==(dyvar_name), vnames(L))
        end

        if object_type == :S
          stock_index = findfirst(==(object_name), snames(L))
          add_part!(L, :LV ; lvv = dyvar_index, lvs = stock_index, lvsposition = operand_index)

        elseif object_type == :P
          param_index = findfirst(==(object_name), pnames(L))
          add_part!(L, :LPV ; lpvv = dyvar_index, lpvp = param_index, lpvpposition = operand_index)

        elseif object_type == :V
          dyvar2_index = findfirst(==(object_name), vnames(L))
          add_part!(L, :LVV ; lvtgt = dyvar_index, lvsrc = dyvar2_index, lvsrcposition = operand_index)
        
        elseif object_type == :SV
          sum_index = findfirst(==(object_name), svnames(L))
          add_part!(L, :LSV ; lsvv = dyvar_index, lsvsv = sum_index, lsvsvposition = operand_index)
        
        # could have a check if it's a F as well...

        end #if/elseif
      end # if
    end # for oper
  end # for dyvar
end #func

# TODO
function remove_from_flow!(object_name, sf_block, L, L_set, name_dict) end

function sfrewrite2(sf::K, block::Expr) where {K <: AbstractStockAndFlowF}

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

  L_redef_queue = Vector{Expr}()

  # used to construct I.  Take L_set - removed_set
  removed_set = Set{Symbol}()

  R_stock_queue = Vector{Symbol}()
  R_param_queue = Vector{Symbol}()
  R_sum_queue = Vector{Expr}()
  R_dyvar_queue = Vector{Expr}()
  R_flow_queue = Vector{Expr}()

  # R_name_dict = Dict{Symbol, SFPointer}()

  current_phase = (_, _) -> ()
  for statement in block.args
    @match statement begin
      QuoteNode(:removes) => begin
        current_phase = removed -> begin
          push!(removed_set, removed)
          if (removed in L_set)
            return
          end

          push!(L_set, removed)

          object_pointer = name_dict[removed]
          object_type = object_pointer.type
          if object_type == :S
            object_definition = sf_block.stocks[object_pointer.index]
            add_part!(L, :S ; sname = object_definition.val)
            remove_from_sums2!(object_definition.val, sf_block, L, L_set, name_dict)
            remove_from_dyvars2!(object_definition.val, sf_block, L, L_set, name_dict)
          elseif object_type == :P
            object_definition = sf_block.params[object_pointer.index]
            add_part!(L, :P ; pname = object_definition.val)
            remove_from_dyvars2!(object_definition.val, sf_block, L, L_set, name_dict)
          elseif object_type == :V
            object_definition = sf_block.dyvars[object_pointer.index]
            dyvar_op = object_definition[2].args[1]
            vname = object_definition[1]
            add_part!(L, :V ; vname = vname, vop = dyvar_op)
            remove_from_dyvars2!(vname, sf_block, L, L_set, name_dict)
          elseif object_type == :SV
            object_definition = sf_block.sums[object_pointer.index]
            svname = object_definition[1]
            add_part!(L, :SV ; svname = svname)
            remove_from_dyvars2!(svname, sf_block, L, L_set, name_dict)
          elseif object_type == :F
            object_definition = sf_block.flows[object_pointer.index]
            #TODO: need to add dyvar before adding flow.
            # fname = object_definition.args[3].args[2].args[1]
            # remove_from_flow!(object_name, sf_block, L, L_set, name_dict)
          end
        end
      end

      QuoteNode(:remove_links) => begin
        current_phase = l -> begin
          src = nothing
          tgt = nothing
          pos = 0
          @match l begin
            
            :($(srca::Symbol) => $(tgta::Symbol)) => begin
              src = srca
              tgt = tgta
            end

            :($(srca::Symbol) => $(tgta::Symbol), $(position::Int)) => begin
              src = srca
              tgt = tgta
              pos = position
            end

            # TODO: error case

          end


          src_pointer = name_dict[src]
          tgt_pointer = name_dict[tgt]

          src_type = src_pointer.type
          tgt_type = tgt_pointer.type
          # @show L_set
          # @show L
          if tgt_type == :V
            # @show tgt
            if !(tgt in L_set)
              push!(L_set, tgt)
              dyvar_op = sf_block.dyvars[tgt_pointer.index][2].args[1]
              add_variable!(L ; vname = tgt, vop = dyvar_op)

            end

            tgt_index = findfirst(==(tgt), vnames(L))
            original_tgt_index = tgt_pointer.index


            if pos == 0
              pos_matches = findall(==(src), sf_block.dyvars[original_tgt_index][2].args[2:end])
            else
              pos_matches = [pos]
            end


            if src_type == :S
              if !(src in L_set)
                push!(L_set, src)
                add_stock!(L ; sname = src)

              end
              
              src_L_index = findfirst(==(src), snames(L))

              for p in pos_matches
                add_part!(L, :LV ; lvs = src_L_index, lvv = tgt_index, lvsposition = p)
              end


            elseif src_type == :P

              if !(src in L_set)
                push!(L_set, src)
                add_parameter!(L ; pname = src)
              end

              src_L_index = findfirst(==(src), pnames(L))

              for p in pos_matches
                add_part!(L, :LPV ; lpvp = src_L_index, lpvv = tgt_index, lpvpposition = p)
              end
            elseif src_type == :SV
              if !(src in L_set)
                push!(L_set, src)
                add_svariable!(L ; svname = src)
              end

              src_L_index = findfirst(==(src), svnames(L))

              for p in pos_matches
                add_part!(L, :LSV ; lsvsv = src_L_index, lsvv = tgt_index, lsvsvposition = p)
              end
            elseif src_type == :V
              if !(src in L_set)
                push!(L_set, src)
                add_variable!(L ; vname = src, vop = sf_block.dyvars[src_pointer.index][2].args[1])
              end
              src_L_index = findfirst(==(src), vnames(L))

              for p in pos_matches
                add_part!(L, :LVV ; lvsrc = src_L_index, lvtgt = tgt_index, lvsrcposition = p)
              end
            #TODO: check for flows?
            end # if src_type

          elseif tgt_type == :SV
            # src must be stock
            tgt_index = findfirst(==(tgt), svnames(L))
            src_L_index = findfirst(==(src), snames(L))


            if !(tgt in L_set)
              push!(L_set, tgt)
              add_variable!(L ; svname = tgt)
            end

            # number of links
            if pos == 0
              link_num = 1
            else
              link_num = pos
            end

            if !(src in L_set)
              push!(L_set, src)
              add_stock!(L ; sname = src)
            end

            # for most sane implementations, link_num should be 1
            add_parts!(L, :LS, link_num ; lss = repeat([src_L_index], link_num), lssv = repeat([tgt_index], link_num))
          
          elseif tgt_type == :F
            tgt_index = findfirst(==(tgt), fnames(L))
            src_L_index = findfirst(==(src), snames(L))

            if !(tgt in L_set)
              push!(L_set, tgt)
              flow_dyvar = sf_block.flows[tgt_pointer.index][2].args[2]
              if !(flow_dyvar in L_set)
                flow_dyvar_op = sf_block.dyvars[name_dict[flow_dyvar].index][2].args[1]
                push!(L_set, flow_dyvar)
                add_variable!(L ; vname = flow_dyvar, vop = flow_dyvar_op)
              end
              flow_dyvar_index = findfirst(==(flow_dyvar), vnames(L))

              add_flow!(L ; fname = tgt, vop = flow_dyvar_index)
            end

            if !(src in L_set)
              push!(L_set, src)
              add_stock!(L ; sname = src)
            end

            if pos == 1
              add_part!(L, :I ; is = src_L_index, ifn = tgt_index)
            elseif pos == 2
              add_part!(L, :O ; os = src_L_index, ofn = tgt_index)
            end

          end # if tgt_type
          
        end # current_phase
      end # quote node

      QuoteNode(:redefs) => begin # only allow vars and sums
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



  for object in L_redef_queue
    object_name = object.args[1]
    object_pointer = name_dict[object_name]
    object_type = object_pointer.type
    object_index = object_pointer.index

    if object_type == :V
      @show L_set

      dyvar_operands = object.args[2].args
      dyvar_op = dyvar_operands[1]
      push!(R_dyvar_queue, object)
      if !(object_name in L_set) 
        add_variable!(L ; vname = object_name, vop = dyvar_op)
        # dyvar_index = nvb(L)
      end
      dyvar_index = findfirst(==(object_name), vnames(L))


      for (operand_index, operand) in enumerate(dyvar_operands[2:end])

        if operand in keys(name_dict)
          # need to add all links to L
          operand_pointer = name_dict[operand]
          operand_type = operand_pointer.type
          operand_original_index = operand_pointer.index


          if !(operand in L_set)
            # note, IT IS NOT REMOVED!
            # It could be removed separately, but as is, this should
            # place it in I

            # TODO: User's responsibility to add the values in the redef under
            # the correct header.  If we have a redef N = [S, I],
            # user needs to add S and I under a :stocks
            # Is this reasonable?

            push!(L_set, operand)
            if operand_type == :S
              operand_definition = sf_block.stocks[operand_original_index]
              add_part!(L, :S ; sname = operand_definition.val)
            elseif operand_type == :P
              operand_definition = sf_block.params[operand_original_index]
              add_part!(L, :P ; pname = operand_definition.val)
            
            elseif operand_type == :V
              operand_definition = sf_block.dyvars[operand_original_index]
              dyvar_op = operand_definition[2].args[1]
              vname = object_definition[1]
              add_part!(L, :V ; vname = vname, vop = dyvar_op)
            elseif operand_type == :SV
              operand_definition = sf_block.sums[operand_original_index]
              svname = operand_definition[1]
              add_part!(L, :SV ; svname = svname)
            end
          end # if !operand

          # only want to add link if there was one in that position originally,
          # hence the continue.

          original_def = sf_block.dyvars[findfirst(x -> x[1] == object_name, sf_block.dyvars)]
          if original_def[2].args[operand_index + 1] != operand
            continue
          end


          if operand_type == :S
            stock_index = findfirst(==(operand), snames(L))
            add_part!(L, :LV ; lvs = stock_index, lvv = dyvar_index, lvsposition = operand_index)
          elseif operand_type == :P
            param_index = findfirst(==(operand), pnames(L))

            add_part!(L, :LPV ; lpvp = param_index, lvv = dyvar_index, lpvpposition = operand_index)
          elseif operand_type == :V
            dyvar2_index = findfirst(==(operand), vnames(L))
            add_part!(L, :LVV ; lvsrc = dyvar2_index, lvtgt = dyvar_index, lvsrcposition = operand_index)
          elseif operand_type == :SV
            sum_index = findfirst(==(operand), svnames(L))
            add_part!(L, :LSV ; lsvsv = sum_index, lsvv = dyvar_index, lsvsvposition = operand_index)
          end
        end # if operand in namedict
      end # for each operand (except operator)

    elseif object_type == :SV
      if !(object_name in L_set)
        push!(R_sum_queue, object)
        add_svariable!(L ; svname = object_name)
      end
      sum_index = findfirst(==(object_name), svnames(L))
      for operand in object.args[2].args
        if operand in keys(name_dict)
          if !(operand in L_set)
            push!(L_set, operand)
            add_stock!(L ; sname = operand)
          end
          
          stock_index = findfirst(==(operand), snames(L))
          add_part!(L, :LS ; lss = stock_index, lssv = sum_index)
        end
    end


    end 


  end


  I_set = setdiff(L_set, removed_set)
  I = StockAndFlowF()
  
  # TODO: This isn't preserving any links.  Verify if I need to...
  # probably need to preserve inflows and outflows, at least...
  for object in I_set
    object_pointer = name_dict[object]
    object_type = object_pointer.type
    object_index = object_pointer.index

    if object_type == :S
      add_stock!(I ; sname = object)
    elseif object_type == :P
      add_parameter!(I ; pname = object)
    elseif object_type == :V
      dyvar_op = sf_block.dyvars[object_index][2].args[1]
      add_variable!(I ; vname = object, vop = dyvar_op)
    elseif object_type == :SV
      add_svariable!(I ; svname = object)
    #TODO: deal with flows.  They need their variable to be added by this point.
    # elseif object_type == :F
    #   flow_var = sf_block.flows[object_index][2].args[1]
      
    #   add_flow!(I, flow_var ; )
    end

  end

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



  if (length(R_sum_queue) > 0)
    sum_names, sum_lists = parse_sum.(collect(filter(x -> !(x.args[1] in svnames(R)), R_sum_queue)))
  else
    sum_names = Vector{Symbol}()
    sum_lists = Vector{Expr}()
  end

  if (length(R_dyvar_queue) > 0)
    dyvars = parse_dyvar.(collect(filter(x -> !(x.args[1] in vnames(R)), R_dyvar_queue)))
    dyvar_names = map(x -> getindex(x, 1), dyvars)
    dyvar_definitions = map(x -> getindex(x, 2), dyvars)
    dyvar_ops = (x -> x.args[1]).(dyvar_definitions)
  else
    dyvar_names = Vector{Symbol}()
    dyvar_definitions = Vector{Expr}()
    dyvar_ops = Vector{Symbol}()
  end

  if (length(R_sum_queue) > 0)
    flow_stocks_in, flows, flow_stocks_out = parse_flow.(collect(filter(x -> !(x.args[2].args[1] in fnames), R_flow_queue)))
    flow_names, flow_variables = (x -> (x.args[1], x.args[2])).(flows)
  else
    flow_stocks_in = Vector{Symbol}()
    flow_stocks_out = Vector{Symbol}()
    flows = Vector{Expr}()
    flow_names = Vector{Symbol}()
    flow_variables = Vector{Symbol}()
  end


  # add dyvars and sums to R_name_dict, so we know their index when we need to
  # link them to a dyvar
  type_index_counter = nvb(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:V, type_index_counter))); type_index_counter+=1)).(dyvar_names)

  type_index_counter = nsv(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:SV, type_index_counter))); type_index_counter+=1)).(sum_names)
  # (x -> push!(R_name_dict, x -> SFPointer(:F, nf(R)))).flow_names



  add_svariables!(R, length(sum_names), svname = sum_names)
  add_variables!(R, length(dyvar_names), vname = dyvar_names, vop = dyvar_ops)
  add_flows!(R, flow_variables, length(flow_names), fnames = flow_names)





  for sum_index in eachindex(R_sum_queue)
    for stock in sum_list.args[1]
      stock_index = findfirst(==(stock), snames(R))
      add_part!(R, :LS ; lss = stock_index, lssv = sum_index)
    end
  end

  # @show R_dyvar_queue

  for dyvar_index in eachindex(R_dyvar_queue)
    dyvar_definition = R_dyvar_queue[dyvar_index].args[2]
    for (operand_index, operand) in enumerate(dyvar_definition.args[2:end])
      object_pointer = R_name_dict[operand]
      object_type = object_pointer.type
      object_index = object_pointer.index
      if object_type == :S
        add_part!(R, :LV ; lvs = object_index, lvv = dyvar_index, lvsposition = operand_index)
      elseif object_type == :P
        add_part!(R, :LPV ; lpvp = object_index, lpvv = dyvar_index, lpvpposition = operand_index)
      elseif object_type == :V
        add_part!(R, :LVV ; lvsrc = object_index, lvtgt = dyvar_index, lvsrcposition = operand_index)
      elseif object_type == :SV
        add_part!(R, :LSV ; lsvsv = object_index, lsvv = dyvar_index, lsvsvposition = operand_index)
      end # could check for flow as well :P
    end
  end

  for flow_index in eachindex(R_flow_queue)
    stock_in = flow_stocks_in[flow_index]
    if stock_in != :F_NONE
      stock_in_index = R_name_dict[stock_in].index
      add_part!(R, :I ; is = stock_in_index, ifn = flow_index)
    end

    stock_out = flow_stocks_out[flow_index]
    if stock_out != :F_NONE
      stock_in_index = R_name_dict[stock_out].index
      add_part!(R, :O ; os = stock_out_index, ofn = flow_index)
    end 
  end


  @show L
  @show I
  @show R

  hom1 = homomorphism(I,L)
  hom2 = homomorphism(I,R)


  rule = Rule(hom1, hom2)


  sf_rewritten = rewrite(rule, sf)


  if isnothing(sf_rewritten)
    println(L)
    println(I)
    println(R)
    return error("Failed to apply rule!")
  end
  return sf_rewritten


  



end



        
"""
Function call to create a new stockflow,
using an existing stockflow and a block describing the modifications.
"""
function sfrewrite(sf::K, block::Expr) where {K <: AbstractStockAndFlowF}

  Base.remove_linenums!(block)
  name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]
  @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"

  sf_block::StockAndFlowBlock = sf_to_block(sf)


  L_stocks::Vector{StockArgT} = []
  L_params::Vector{ParamArgT} = []
  L_dyvars::Vector{Tuple{Symbol,Expr}} = []
  L_flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
  L_sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []

  L_block = StockAndFlowBlock(L_stocks, L_params, L_dyvars, L_flows, L_sums)
  
  R_stocks::Vector{StockArgT} = []
  R_params::Vector{ParamArgT} = []
  R_dyvars::Vector{Tuple{Symbol,Expr}} = []
  R_flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
  R_sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []

  R_block = StockAndFlowBlock(R_stocks, R_params, R_dyvars, R_flows, R_sums)


  # name -> type, index
  name_dict = Dict{Any, Tuple{Symbol, Int}}([[s => (:S, i) for (i, s) ∈ enumerate(snames(sf))] ;
  [sv => (:SV, i) for (i, sv) ∈ enumerate(svnames(sf))] ;
  [v => (:V, i) for (i, v) ∈ enumerate(vnames(sf))] ;
  [f => (:F, i) for (i, f) ∈ enumerate(fnames(sf))] ;
  [p => (:P, i) for (i, p) ∈ enumerate(pnames(sf))]])



  removed_set = Set{Any}()
  modified_dict = Dict{Any, Any}()


  # Ensures we don't add the same object multiple times
  L_set = Set{Any}()


  # Allows me to update in R repeatedly.
  R_dict = Dict{Any, Any}()


  
  current_phase = (_, _) -> ()
  for statement in block.args
    @match statement begin
      QuoteNode(:removes) => begin
        current_phase = removed -> begin
          push!(removed_set, removed)
          if removed ∉ L_set && removed ∈ keys(name_dict)
            push!(L_set, removed)
            object_definitition = get_from_correct_vector(sf_block, name_dict[removed][2], name_dict[removed][1])
            add_to_correct_vector!(L_block, object_definitition, name_dict[removed][1])
            remove_from_sums!(removed, sf_block, L_block, L_set, name_dict)
            remove_from_dyvars!(removed, sf_block, L_block, L_set, name_dict)
            if name_dict[removed][1] == :F
              index = name_dict[removed][2] 
              definition = deepcopy(sf_block.flows[index])
              stock1 = definition[1]
              if stock1 != :F_NONE && stock1 ∉ L_set && stock1 ∈ keys(name_dict)
                push!(L_block.stocks, StockArgUnitSymbol(stock1, :NONE))
                push!(L_set, stock1)
              end
              stock2 = definition[3]
              if stock2 != :F_NONE && stock2 ∉ L_set && stock2 ∈ keys(name_dict)
                push!(L_block.stocks, StockArgUnitSymbol(stock2, :NONE))
                push!(L_set, stock2)
              end
              flow_var = definition[2].args[2]
              flow_var_index = name_dict[flow_var][2]
              flow_var_definition = deepcopy(sf_block.dyvars[flow_var_index])
              if flow_var ∉ L_set && flow_var ∈ keys(name_dict)
                push!(L_block.dyvars, flow_var_definition)
                push!(L_set, flow_var)
              end
            end
          end
        end
      end
      QuoteNode(:redefs) => begin
        current_phase = redef -> begin
          @match redef begin
            Expr(:(=), src, tgt) => begin
              

              src_type = name_dict[src][1]
              src_index = name_dict[src][2]
              if src_type == :V
                object_definition = parse_dyvar(Expr(:(=), src, tgt))
                original_definition = deepcopy(sf_block.dyvars[src_index])

                intersection_equation_symbols = (object_definition[2].args[2:end] ∩ original_definition[2].args[2:end])
                equation_operator = original_definition[2].args[1]


                @assert equation_operator == object_definition[2].args[1] "Different operator in new equation!\
                 new: $object_definition old: $original_definition"


                push!(modified_dict, src => parse_dyvar(Expr(:(=), src, Expr(:call, equation_operator, intersection_equation_symbols...)))) # wow.
                

                if src ∉ keys(R_dict)
                  push!(R_dict, src => object_definition)
                  push!(R_block.dyvars, object_definition)
                end
                if src ∉ L_set
                  push!(L_set, src)
                  push!(L_block.dyvars, original_definition)

                end
              
              elseif src_type == :SV
                object_definition = parse_sum(Expr(:(=), src, tgt))
                original_definition = deepcopy(sf_block.sums[src_index])
                I_definition = parse_sum(Expr(:(=), src, Expr(:vect, (object_definition[2] ∩ original_definition[2])...)))


                push!(modified_dict, src => I_definition)
                if src ∉ keys(R_dict)
                  push!(R_dict, src => object_definition)
                  push!(R_block.sums, object_definition)
                end
                if src ∉ L_set
                  push!(L_set, src)
                  push!(L_block.sums, original_definition)
                end

              elseif src_type == :F
                object_definition = parse_flow(tgt)
                original_definition = deepcopy(sf_block.flows[src_index])

                # push!(modified_dict, src => object_definition)
                if src ∉ keys(R_dict)
                  push!(R_dict, src => object_definition)
                  push!(R_block.flows, object_definition)
                end
                if src ∉ L_set
                  push!(L_set, src)
                  push!(L_block.flows, original_definition)
                end
              end
            end
          end
        end
      end



      QuoteNode(:dyvar_swaps) => begin 
        current_phase = swap -> begin
          @match swap begin
            Expr(:call, :(=>), src, tgt) => begin
              swap_from_dyvars!(src, tgt, sf_block, L_block, L_set, R_block, R_dict, name_dict)
            end
          end
        end
      end
      QuoteNode(:stocks) => begin 
        current_phase = s -> begin 
          @match s begin
            val::Symbol => begin # stocks, params
              object_definition = parse_stock(val)
              # You're adding it, so it's gonna be added.  No checking if it's already there.
              push!(R_dict, object_definition.val => object_definition)
              push!(R_block.stocks, object_definition)
            end
          end
        end
      end
      
      QuoteNode(:parameters) => begin 
        current_phase = p -> begin 
          @match p begin
            val::Symbol => begin # stocks, params
              object_definition = parse_param(val)
              # You're adding it, so it's gonna be added.  No checking if it's already there.
              push!(R_dict, object_definition.val => object_definition)
              push!(R_block.params, object_definition)
            end
          end
        end
      end
      
      QuoteNode(:dynamic_variables) => begin
        current_phase = v -> begin 
          @match v begin
            :($tgt = $def) => begin
              object_definition = parse_dyvar(v)
              # You're adding it, so it's gonna be added.  No checking if it's already there.
              push!(R_dict, object_definition[1] => object_definition)
              push!(R_block.dyvars, object_definition)
            end
          end
        end
      end
      
      QuoteNode(:flows) => begin
        current_phase = f -> begin 
          @match f begin
             Expr(:call, :(=>), S1::Symbol, rest) => begin # flows
              object_definition = parse_flow(Expr(:call, :(=>), S1, rest))
              object_name = object_definition[2].args[1]
              # You're adding it, so it's gonna be added.  No checking if it's already there.
              push!(R_dict, object_name => object_definition)
              push!(R_block.flows, object_definition)
            end
          end
        end
      end
      
      QuoteNode(:sums) => begin
        current_phase = sv -> begin 
          @match sv begin
            Expr(:(=), tgt::Symbol, def) => begin
              object_definition = parse_sum(Expr(:(=), tgt, def))
              # You're adding it, so it's gonna be added.  No checking if it's already there.
              push!(R_dict, object_definition[1] => object_definition)
              push!(R_block.sums, object_definition)
            end
          end
        end
      end
      
      QuoteNode(kw) =>
        error("Unknown block type for rewrite syntax: " * String(kw))
      _ => current_phase(statement)
      end
    end




  # add all that is in dyvar's definition to L
  for current_dyvar ∈ L_block.dyvars

    # dyvar_name = current_dyvar[1]
    # dyvar_copy = deepcopy(current_dyvar)
    # dyvar_definition = dyvar_copy[2]
      
    # if dyvar_name ∉ new_set
    #   push!(new_block.dyvars, dyvar_copy)
    #   push!(L_set, dyvar_name => dyvar_copy)
    # end

    # for object ∈ filter(∉(L_set), dyvar_definition.args[2:end])
    #   object_type = name_dict[object][1]
    #   object_index = name_dict[object][2]
    #   object_definition = get_from_correct_vector(sf_block, object_index, object_type)
    #   if object_type == :SV
    #     sum_name = object_definition.args[1]
    #     sum_definition = :(sum_name = [])
    #     push(L_set, sum_name => object_definition)
    #     push!(L_block.sums, )
    #   elseif object_type == :V
    #     v_name = object_definition.args[1]
    #     v_op = object_definition.args[2].args[1]

    #   end

    #   # push!(new_set, object)
    #   add_to_correct_vector!(new_block, object_definition, object_type)
  


    recursively_add_dyvars_L!(current_dyvar, sf_block, L_block, L_set, name_dict)
  end
  
  for sum ∈ L_block.sums
    for stock ∈ sum[2]
      if stock ∉ L_set
        push!(L_set, stock)
        push!(L_block.stocks, StockArgUnitSymbol(stock, :NONE))
      end
    end
  end

  for flow ∈ sf_block.flows
    if flow[2].args[1] ∈ L_set
      continue
    end
    stock1 = flow[1]
    stock2 = flow[3]
    if stock1 ∈ removed_set || stock2 ∈ removed_set
      push!(L_block.flows, deepcopy(flow))
      push!(L_set, flow[2].args[1])
      if stock1 ∉ L_set
        push!(L_set, stock1)
        push!(L_block.stocks, StockArgUnitSymbol(stock1, :NONE))
      end
      if stock2 ∉ L_set
        push!(L_set, stock2)
        push!(L_block.stocks, StockArgUnitSymbol(stock2, :NONE))
      end
      if flow[2].args[2] ∉ L_set
        push!(L_set, flow[2].args[2])
        flow_dyvar_index = name_dict[flow[2].args[2]][2]
        flow_dyvar_definition = deepcopy(sf_block.dyvars[flow_dyvar_index])
        push!(L_block.dyvars, flow_dyvar_definition)
      end
    end
  end
        

  
  I_block = deepcopy(L_block)
  


  filter!(s -> s.val ∉ removed_set, I_block.stocks)
  filter!(p -> p.val ∉ removed_set, I_block.params)
  filter!(x -> x[1] ∉ removed_set, I_block.sums)
  filter!(x -> x[1] ∉ removed_set, I_block.dyvars)
  filter!(x -> x[2].args[1] ∉ removed_set, I_block.flows)





  for (i, dyvar) ∈ enumerate(I_block.dyvars)
    if dyvar[1] in keys(modified_dict)
      I_block.dyvars[i] = modified_dict[dyvar[1]]
    else
      filter!(∉(removed_set), dyvar[2].args)
    end

  end

  
  
  for (i, sum) ∈ enumerate(I_block.sums)
    if sum[1] in keys(modified_dict)
      I_block.sums[i] = modified_dict[sum[1]]
    else
      filter!(∉(removed_set), sum[2])
    end
  end

  for (i, flow) ∈ enumerate(I_block.flows)
    if flow[2].args[1] in keys(modified_dict)
      I_block.flows[i] = modified_dict[flow[2].args[1]]
    else
      if flow[1] ∈ removed_set
        I_block.flows[i] = tuple(replace(collect(flow), flow[1] => :F_NONE)...)
      end
      if flow[2].args[2] ∈ removed_set
        deleteat!(flow[2].args[2], 2) # Now has no flow variable, good job :P
      end
      if flow[3] ∈ removed_set
        I_block.flows[i] = tuple(replace(collect(flow), flow[3] => :F_NONE)...)
      end
    end
  end

  for s in I_block.stocks
    if s.val ∉ keys(R_dict)
      push!(R_dict, s.val => s)
      push!(R_block.stocks, s)
    end
  end
  for p in I_block.params
    if p.val ∉ keys(R_dict)
      push!(R_dict, p.val => (p, ))
      push!(R_block.params, p)
    end
  end
  for sv in I_block.sums
    if sv[1] ∉ keys(R_dict)
      sum_copy = deepcopy(sv)
      push!(R_dict, sum_copy[1] => sum_copy)
      push!(R_block.sums, sum_copy)
    end
  end
  for v in I_block.dyvars
    if v[1] ∉ keys(R_dict)
      dyvar_copy = deepcopy(v)
      push!(R_dict, dyvar_copy[1] => dyvar_copy)
      push!(R_block.dyvars, dyvar_copy)
    end
  end
  for f in I_block.flows
    if f[2].args[1] ∉ keys(R_dict)
      flow_copy = deepcopy(f)
      push!(R_dict, flow_copy[2].args[1] => flow_copy)
      push!(R_block.flows, flow_copy)
    end
  end


  

  for current_dyvar ∈ R_block.dyvars
    recursively_add_dyvars_R!(current_dyvar, sf_block, R_block, R_dict, name_dict)
  end



  for sum ∈ R_block.sums
    for stock ∈ sum[2]
      if stock ∉ keys(R_dict)
        push!(R_dict, stock => (stock,))
        push!(R_block.stocks, StockArgUnitSymbol(stock, :NONE))
      end
    end
  end

  # maintain order
  # order by indices in original, or don't care if not in original
  # terrible way to do this though.

  max_val = max(length(R_block.stocks), length(R_block.params), length(R_block.dyvars), length(R_block.sums), length(R_block.flows)) + 1
  mv_array = [max_val,max_val]


  sort!(R_block.stocks, by = x-> get(name_dict, x.val, mv_array)[2])
  sort!(R_block.params, by = x-> get(name_dict, x.val, mv_array)[2])
  sort!(R_block.dyvars, by = x-> get(name_dict, x[1], mv_array)[2])
  sort!(R_block.sums, by = x-> get(name_dict, x[1], mv_array)[2])
  sort!(R_block.flows, by = x-> get(name_dict, x[2].args[1], mv_array)[2])
  
  L_args = stock_and_flow_syntax_to_arguments(L_block)

  
  L = StockAndFlowF(L_args.stocks, L_args.params, map(kv -> kv.first => StockFlow.Syntax.get(kv.second), L_args.dyvars), L_args.flows, L_args.sums)

  I_args = stock_and_flow_syntax_to_arguments(I_block)
  

  I = StockAndFlowF(I_args.stocks, I_args.params, map(kv -> kv.first => StockFlow.Syntax.get(kv.second), I_args.dyvars), I_args.flows, I_args.sums)

  R_args = stock_and_flow_syntax_to_arguments(R_block)


  R = StockAndFlowF(R_args.stocks, R_args.params,  map(kv -> kv.first => StockFlow.Syntax.get(kv.second), R_args.dyvars), R_args.flows, R_args.sums)

  


  # Now we line the positions up with the original stockflow



  reset_positions!(sf, L)
  reset_positions!(sf, I)
  reset_positions!(sf, R)
  
  hom = homomorphism
  hom1 = hom(I,L)
  hom2 = hom(I,R)


  rule = Rule(hom1, hom2)


  sf_rewritten = rewrite(rule, sf)


  if isnothing(sf_rewritten)
    println(L)
    println(I)
    println(R)
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
    sfrewrite2($sf, $escaped_block)
  end
end
  


end