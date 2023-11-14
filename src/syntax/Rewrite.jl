module Rewrite

export sfrewrite, @rewrite



using ...StockFlow

using ..Syntax
import ..Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, 
parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments


using Catlab.ACSets.ACSetInterface
using Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite


using MLStyle



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
"""
function reset_positions!(sfold, sfnew)
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
        if old ∉ keys(oldlv)
            continue
        end
        # if (lvs, lvv) ∉ keys(oldlv)
        #     continue
        # end
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
                if old ∉ keys(oldlsv)
            continue
        end
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

        if old ∉ keys(oldlvv)
            continue
        end
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
                if old ∉ keys(oldlpv)
            continue
        end
      old_pos = oldlpv[old]
      restored_lpvpposition[lpv] = old_pos
  end

  set_subpart!(sfnew, :lvsposition, restored_lvsposition)
  set_subpart!(sfnew, :lsvsvposition, restored_lsvsvposition)
  set_subpart!(sfnew, :lvsrcposition, restored_lvsrcposition)
  set_subpart!(sfnew, :lpvpposition, restored_lpvpposition)

  return sfnew
end


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
        
        if stock_name ∈ sum_stocks
            if sum_name ∉ L_set
                push!(L_block.sums, sum)
              
            end
        end
    end
end

function recursively_add_dyvars_L!(current_dyvar, sf_block, new_block, new_set, name_dict)
    dyvar_name = current_dyvar[1]
     dyvar_copy = deepcopy(current_dyvar)
        dyvar_definition = dyvar_copy[2]
        
    if dyvar_name ∉ new_set
       
        push!(new_block.dyvars, dyvar_copy)
        push!(new_set, dyvar_name => dyvar_copy)
    end
        
    # For everything not yet in the dict, add t dict
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
                
                

function remove_from_dyvars!(object_name, sf_block, L_block, L_set, name_dict)
    for dyvar ∈ sf_block.dyvars
        dyvar_name = dyvar[1]
        dyvar_expression = dyvar[2]
        
        if object_name ∈ dyvar_expression.args

            
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
                        
                        if object_type == :SV
                            for stock ∈ object_definition[2]
                                if stock ∉ L_set
                                    push!(L_set, stock)
                                    push!(L_block.stocks, stock)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


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
                        
                        if object_type == :SV
                            for stock ∈ object_definition[2]
                                if stock ∉ L_set
                                    push!(L_set, stock)
                                    push!(L_block.stocks, stock)
                                end
                            end
                        end
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
                
                

        
function sfrewrite(sf::K, block::Expr) where {K <: AbstractStockAndFlowF}
    Base.remove_linenums!(block)
    name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]
    @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"

    sf_block::StockAndFlowBlock = sf_to_block(sf)
    
	L_stocks::Vector{Symbol} = []
    L_params::Vector{Symbol} = []
	L_dyvars::Vector{Tuple{Symbol,Expr}} = []
	L_flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
	L_sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []

    L_block = StockAndFlowBlock(L_stocks, L_params, L_dyvars, L_flows, L_sums)
    
	R_stocks::Vector{Symbol} = []
    R_params::Vector{Symbol} = []
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
    


    # Ensures we don't add the same object multiple times
    L_set = Set{Any}()


    # Allows me to update in R repeatedly.
    R_dict = Dict{Any, Tuple}()


    
    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:swaps) => begin 
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
                        Expr(:call, :+, val) => begin # stocks, params
                            object_definition = parse_stock(val)
                            # You're adding it, so it's gonna be added.  No checking if it's already there.
                            push!(R_dict, object_definition => (object_definition,))
                            push!(R_block.stocks, object_definition)
                        end
                        Expr(:call, :-, val) => begin
                            object_definition = parse_stock(val)
                            push!(removed_set, object_definition)

                            remove_from_sums!(object_definition, sf_block, L_block, L_set, name_dict)
                            remove_from_dyvars!(object_definition, sf_block, L_block, L_set, name_dict)
                        end
                    end
                end
            end
          
          QuoteNode(:parameters) => begin 
                current_phase = p -> begin 
                    @match p begin
                        Expr(:call, :+, val) => begin # stocks, params
                            object_definition = parse_param(val)
                            # You're adding it, so it's gonna be added.  No checking if it's already there.
                            push!(R_dict, object_definition => (object_definition,))
                            push!(R_block.params, object_definition)
                        end
                        Expr(:call, :-, val) => begin
                            object_definition = parse_param(val)

                            push!(removed_set, object_definition)
                            remove_from_dyvars!(object_definition, sf_block, L_block, L_set, name_dict)
                        end
                    end
                end
            end
            
          QuoteNode(:dynamic_variables) => begin
                current_phase = v -> begin 
                    @match v begin
                        Expr(:(=), Expr(:call, :+, tgt), Expr(:block, def)) => begin
                            object_definition = parse_dyvar(Expr(:(=), tgt, def))
                            # You're adding it, so it's gonna be added.  No checking if it's already there.
                            push!(R_dict, object_definition[1] => object_definition)
                            push!(R_block.dyvars, object_definition)
                        end
                        Expr(:(=), Expr(:call, :-, tgt), Expr(:block, def)) => begin
                            object_definition = parse_dyvar(Expr(:(=), tgt, def))

                            push!(removed_set, object_definition[1])
                            remove_from_dyvars!(object_definition[1], sf_block, L_block, L_set, name_dict)
                        end
                    end
                end
            end
            
          QuoteNode(:flows) => begin
                current_phase = sv -> begin 
                    @match sv begin
                       Expr(:call, :(=>), Expr(:call, :+, S1), rest) => begin # flows
                            object_definition = parse_flow(Expr(:call, :(=>), S1, rest))
                            object_name = object_definition[2].args[1]
                            # You're adding it, so it's gonna be added.  No checking if it's already there.
                            push!(R_dict, object_name => object_definition)
                            push!(R_block.flows, object_definition)
                        end
                       Expr(:call, :(=>), Expr(:call, :-, S1), rest) => begin # flows
                            object_definition = parse_flow(Expr(:call, :(=>), S1, rest))
                            object_name = object_definition[2].args[1]
                            push!(removed_set, object_name)
                        end
                    end
                end
            end
            
          QuoteNode(:sums) => begin
                current_phase = sv -> begin 
                    @match sv begin
                        Expr(:(=), Expr(:call, :+, tgt), Expr(:block, def)) => begin
                            object_definition = parse_sum(Expr(:(=), tgt, def))
                            # You're adding it, so it's gonna be added.  No checking if it's already there.
                            push!(R_dict, object_definition[1] => object_definition)
                            push!(R_block.sums, object_definition)
                        end
                        Expr(:(=), Expr(:call, :-, tgt), Expr(:block, def)) => begin
                            object_definition = parse_sum(Expr(:(=), tgt, def))

                            push!(removed_set, object_definition[1])
                            remove_from_dyvars!(object_definition[1], sf_block, L_block, L_set, name_dict)
                        end
                    end
                end
            end
            
          QuoteNode(kw) =>
              error("Unknown block type for rewrite syntax: " * String(kw))
          _ => current_phase(statement)
          end
      end





    for current_dyvar ∈ L_block.dyvars
        recursively_add_dyvars_L!(current_dyvar, sf_block, L_block, L_set, name_dict)
    end
    
    for sum ∈ L_block.sums
        for stock ∈ sum[2]
            if stock ∉ L_set
                push!(L_set, stock)
                push!(L_block.stocks, stock)
            end
        end
    end


    
    I_block = deepcopy(L_block)


    filter!(∉(removed_set), I_block.stocks)
    filter!(∉(removed_set), I_block.params)
    filter!(x -> x[1] ∉ removed_set, I_block.sums)
    filter!(x -> x[1] ∉ removed_set, I_block.dyvars)
    filter!(x -> x[2].args[1] ∉ removed_set, I_block.flows)





    for dyvar ∈ I_block.dyvars
        filter!(∉(removed_set), dyvar[2].args)
    end

    
    
    for sum ∈ I_block.sums
        filter!(∉(removed_set), sum[2])
    end
    
    for flow ∈ I_block.flows
        if flow.args[1] ∈ removed_set
            flow.args[1] = :CLOUD
        end
        if flow.args[2].args[2] ∈ removed_set
            deleteat!(flow.args[2].args[2], 2) # Now has no flow variable, good job :P
        end
        if flow.args[3] ∈ removed_set
            flow.args[3] = :CLOUD
        end
    end




    # I'm not 100% sure that this is actually correct; here's hoping

  

    for current_dyvar ∈ R_block.dyvars
        recursively_add_dyvars_R!(current_dyvar, sf_block, R_block, R_dict, name_dict)
    end



    for sum ∈ R_block.sums
        for stock ∈ sum[2]
            if stock ∉ keys(R_dict)
                push!(R_dict, stock => (stock,))
                push!(R_block.stocks, stock)
            end
        end
    end


    
    L_args = stock_and_flow_syntax_to_arguments(L_block)

    
    L = StockAndFlowF(L_args.stocks, L_args.params, map(kv -> kv.first => StockFlow.Syntax.get(kv.second), L_args.dyvars), L_args.flows, L_args.sums)

    I_args = stock_and_flow_syntax_to_arguments(I_block)
    I = StockAndFlowF(I_args.stocks, I_args.params, map(kv -> kv.first => StockFlow.Syntax.get(kv.second), I_args.dyvars), I_args.flows, I_args.sums)

    R_args = stock_and_flow_syntax_to_arguments(R_block)
    R = StockAndFlowF(R_args.stocks, R_args.params,  map(kv -> kv.first => StockFlow.Syntax.get(kv.second), R_args.dyvars), R_args.flows, R_args.sums)

    
    # # Now we line the positions up with the original stockflow


    reset_positions!(sf, L)
    reset_positions!(sf, I)
    reset_positions!(sf, R)



    
    hom = homomorphism
    rule = Rule(hom(I,L), hom(I,R))

    sf_rewritten = rewrite(rule, sf)

    return sf_rewritten
        
end


macro rewrite(sf, block)
    escaped_block = Expr(:quote, block)
    sf = esc(sf)
    quote
      sfrewrite($sf, $escaped_block)
    end
  end
  


end