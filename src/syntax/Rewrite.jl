module Rewrite


export sfrewrite


using ...StockFlow
using ..StockFlow.Syntax
import ..StockFlow.Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments
using MLStyle


function replace_symbol(A::Symbol, old_new::Dict)
    return get(old_new, A, A)
end
function replace_symbol(A::Expr, old_new::Dict)
    expr_args = [x for x ∈ A.args]
    replace!(expr_args, old_new...)
    return Expr(A.head, expr_args...)
    # return get(old_new, A, A)
end



function add_swap!(sw, sfblock, all_names, modified_stocks, modified_flows, modified_dyvars, modified_params, modified_sums)
    @match sw begin
        :($A => $B) => begin # stock => stock, flow => flow, dv => dv, etc
            typeA, otherA... = all_names[A]
            typeB, otherB... = all_names[B]
            @assert (typeA == typeB) "Cannot swap a $typeA for a $typeB !"
            # TODO: Add the ability to swap S, P and SV, implying you can move links from one to the other?
            # DV will be more complicated.
            @match typeA begin
                :S => push!(modified_stocks, A => B) #TODO: Cascade.  All instances of A need to be swapped with B.
                :P => push!(modified_params, A => B)
                :SV => begin
                    original_index = all_names[A][2]
                    original_definition = sfblock.sums[original_index]
                    new_definition = (B, original_definition[2])
                    push!(modified_sums, original_definition => new_definition)
                end
                :V => begin
                    original_index = all_names[A][2]
                    original_definition = sfblock.dyvars[original_index]
                    new_definition = (B, original_definition[2])
                    push!(modified_dyvars, original_definition => new_definition)
                end
                :F => begin
                    original_index = all_names[A][2]
                    original_definition = sfblock.flows[original_index]

                    original_flow = original_definition[2] # out, flow(dv), in
                    original_fv = original_flow.args[2] # We converted a sf in such a way that flows will always take the form :(f(dv)), or equivalently, Expr(:call, :f, :dv)

                    new_flow = Expr(:call, B, original_fv)
                    new_definition = (original_definition[1], new_flow, original_definition[2])
                    push!(modified_flows, original_definition => new_definition)
                end
            end


        end
        Expr(:(:=), tgt, new) => begin # sum := new, dv := new
            @match all_names[tgt][1] begin # tells us what type it is
            # Three cases: :V, :SV, :F
            # S and P are just the same, since they're atomic

                :S => push!(modified_stocks, tgt => new)
                :P => push!(modified_params, tgt => new)

                :V => begin
                    original_index = all_names[tgt][2]
                    original_definition = sfblock.dyvars[original_index]
                    new_definition = parse_dyvar(Expr(:(=), tgt, new))
                    push!(modified_dyvars, original_definition => new_definition)
                end

                :SV => begin
                    original_index = all_names[tgt][2]
                    original_definition = sfblock.sums[original_index]

                    new_definition = parse_sum(Expr(:(=), tgt, new))
                    push!(modified_sums, original_definition => new_definition)
                end

                :F => begin
                    original_index = all_names[tgt][2]
                    original_definition = sfblock.flows[original_index]

                    new_definition = parse_flow(Expr(:(=), tgt, new))
                    push!(modified_flows, original_definition => new_definition)
                end
            end
        end
        _ => error("Unknown expression in swaps section.")
    end
end


"""
ok, so some explanation
when adding, you need to specify the whole definition, like +v1 = A + B
when removing, you just need to specify the name. -v2. We already have the defintion in packed_sf, which we can
 get from finding the index in all_names
"""

function modify_rewrite_dict!(line, added_objects, removed_objects, parse_definition_function!)
    Base.remove_linenums!(line)
    @match line begin 
        Expr(:call, :+, val) => begin # stocks, params
            object_definition = parse_definition_function!(val)
            push!(added_objects, object_definition)
        end
        Expr(:call, :-, val) => push!(removed_objects, val) 



        Expr(:(=), Expr(:call, :+, tgt), Expr(:block, rest)) => begin # dyvars, sums,
            object_definition = parse_definition_function!(Expr(:(=), tgt, rest))
            push!(added_objects, object_definition)
        end
        Expr(:call, :(=>), Expr(:call, :+, S1), rest) => begin # flows


        object_definition = parse_definition_function!(Expr(:call, :(=>), S1, rest))
        push!(added_objects, object_definition)
    end
        
        _ => error("Unknown expression $line ")
    end
end


"""
After a day and a half of work on this, I realized it made more sense to operate on StockAndFlowBlocks rather than StockAndFlows
Links will be automatically dealt with when creating new StockFlows
"""
function sfrewrite(sf::AbstractStockAndFlowF, block)

    Base.remove_linenums!(block)
    sfblock = sf_to_block(sf)


    added_params::Vector{Symbol} = []
    added_stocks::Vector{Symbol} = []
    added_dyvars::Vector{Tuple{Symbol,Expr}} = []
    added_flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
    added_sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []

    # We just need the names
    removed_params::Vector{Symbol} = [] 
    removed_stocks::Vector{Symbol} = []
    removed_dyvars::Vector{Symbol} = []
    removed_flows::Vector{Symbol} = []
    removed_sums::Vector{Symbol} = []

    modified_params::Dict{Symbol, Symbol} = Dict()
    modified_stocks::Dict{Symbol, Symbol} = Dict()
    modified_dyvars::Dict{Tuple{Symbol,Expr}, Tuple{Symbol,Expr}} = Dict()
    modified_flows::Dict{Tuple{Symbol,Expr,Symbol}, Tuple{Symbol,Expr,Symbol}} = Dict()
    modified_sums::Dict{Tuple{Symbol,Vector{Symbol}}, Tuple{Symbol,Vector{Symbol}}} = Dict()




    sf_snames = Dict(S => i for (i, S) in enumerate(snames(sf)))
    sf_svnames = Dict(S => i for (i, S) in enumerate(svnames(sf)))
    sf_vnames = Dict(S => i for (i, S) in enumerate(vnames(sf)))
    sf_fnames = Dict(S => i for (i, S) in enumerate(fnames(sf)))
    sf_pnames = Dict(S => i for (i, S) in enumerate(pnames(sf)))

    
    

    name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]

    # println(name_vector)
    @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"
    # This may be a bit annoying when after a composition, you can have things that have the same name, but this makes rewriting a lot simpler
    # Unless you'd rather specify indices than names, I guess
    # Can always run this through set_snames! before calling.
    

    flow_variable_indices = Set(fvs(sf))


    all_names = Dict() # note, because we've confirmed all the names are unique, there will be no overwriting
    # Doubtless there's a more graceful way to do this than this, though
    # Basically just creating a dictionary such that, when we encounter a symbol, we know both what type it represents, and what its index is
    merge!(all_names, Dict(s => (:S, i) for (s, i) ∈ sf_snames))
    merge!(all_names, Dict(sv => (:SV, i) for (sv, i) ∈ sf_svnames))
    merge!(all_names, Dict(v => (:V, i, vop(sf, i), i ∈ flow_variable_indices) for (v, i) ∈ sf_vnames)) # Symbol => (:V, index, operator, is flow variable?)
    merge!(all_names, Dict(f => (:F, i) for (f, i) ∈ sf_fnames))
    merge!(all_names, Dict(p => (:P, i) for (p, i) ∈ sf_pnames))


    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:swaps) => begin
                current_phase = sw -> add_swap!(sw, sfblock, all_names, modified_stocks, modified_flows, modified_dyvars, modified_params, modified_sums) # this one will be different than the others
            end
            QuoteNode(:stocks) => begin
                current_phase = s -> modify_rewrite_dict!(s, added_stocks, removed_stocks, x -> parse_stock(x))
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> modify_rewrite_dict!(p, added_params, removed_params, x -> parse_param(x))
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = v -> modify_rewrite_dict!(v, added_dyvars, removed_dyvars, x -> parse_dyvar(x))
            end
            QuoteNode(:flows) => begin
                current_phase = f -> modify_rewrite_dict!(f, added_flows, removed_flows, x -> parse_flow(x))
            end
            QuoteNode(:sums) => begin
                current_phase = sv -> modify_rewrite_dict!(sv, added_sums, removed_sums, x -> parse_sum(x))
            end
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end


    # 😇 (good): Create I, L and R, define and apply rewrite rule
    # 😈 (evil): Create a new stock flow directly from information gathered.

    # Going evil route for now

    new_stocks = add_remove_modify(sfblock.stocks, added_stocks, removed_stocks, modified_stocks)
    new_flows = add_remove_modify(sfblock.flows, added_flows, removed_flows, modified_flows)
    new_dyvars = add_remove_modify(sfblock.dyvars, added_dyvars, removed_dyvars, modified_dyvars)
    new_params = add_remove_modify(sfblock.params, added_params, removed_params, modified_params)
    new_sums = add_remove_modify(sfblock.sums, added_sums, removed_sums, modified_sums)


    # new_stocks_dict = Dict(stock => i for (i, stock) ∈ enumerate(new_stocks))
    # new_params_dict = Dict(param => i for (i, param) ∈ enumerate(new_params))
    # new_flows_dict = Dict(flow[1] => i for (i, flow) ∈ enumerate(new_flows))
    # new_dyvars_dict = Dict(dyvar[1] => i for (i, dyvar) ∈ enumerate(new_dyvars))
    # new_sums_dict = Dict(sum => i for (i, sum) ∈ enumerate(new_sums))

    new_mapping_stocks = modified_stocks
    new_mapping_params = modified_params
    # new_mapping_flows = Dict(key[2] => value[2] for (key, value) ∈ modified_flows)
    new_mapping_dyvars = Dict(key[1] => value[1] for (key, value) ∈ modified_dyvars)
    new_mapping_sums = Dict(key[1] => value[1] for (key, value) ∈ modified_sums)



    # REPLACEMENTS:
    # Check dyvars for changed dyvars, sums, stocks, params
    # Check flows for changed stocks (inflow and outflow), dyvars
    # Check sums for changed stocks

    # TODO: Add assert that there only exists one remapping for each var (gonna be somewhere else in this code)
    dyvar_replace_dict::Dict{Symbol, Symbol} = Dict(new_mapping_stocks..., new_mapping_dyvars..., new_mapping_params..., new_mapping_sums...)

    dyvar_replace_pairs::Vector{Pair{Symbol, Symbol}} = [k => v for (k, v) ∈ dyvar_replace_dict]

    for (i, dyvar) ∈ enumerate(new_dyvars)
        old_args::Vector{Symbol} = [x for x ∈ dyvar[2].args] # getting around the Vector[Any]



        replace!(old_args, dyvar_replace_pairs...) # I don't know why I need to do this, but when
        # using regular replace, I get an error.
        # can't replicate the error in the REPL so idk what's happening
        # TODO: Get to the bottom of this
        new_expression = Expr(:call, old_args...)
        new_dyvars[i] = (dyvar[1], new_expression)
    end

    for (i, flow) ∈ enumerate(new_flows)
        new_inflow = replace_symbol(flow[1], new_mapping_stocks)

        flow_definition = flow[2]
        new_variable = replace_symbol(flow_definition.args[2], new_mapping_dyvars)
        new_flow_definition = Expr(:call, flow_definition.args[1], new_variable)
        new_outflow = replace_symbol(flow[3], new_mapping_stocks)

        new_flows[i] = (new_inflow, new_flow_definition, new_outflow)
    end

    for (i, sum) ∈ enumerate(new_sums)
        old_summands = copy(sum[2])
        replace!(old_summands, new_mapping_stocks...)
        # Same replace error as above
        new_sums[i] = (sum[1], old_summands)
    end
    


    new_sfblock = StockAndFlowBlock(new_stocks, new_params, new_dyvars, new_flows, new_sums)
    new_sfarguments = stock_and_flow_syntax_to_arguments(new_sfblock)

    return StockAndFlowF(
        new_sfarguments.stocks,
        new_sfarguments.params,
        map(kv -> kv.first => get(kv.second), new_sfarguments.dyvars),
        new_sfarguments.flows,
        new_sfarguments.sums,
    )

    


end


"""
Buy it, use it, break it, fix it, trash it, change it, mail, upgrade it 
"""
function add_remove_modify(old, added, removed, modified)
    new = Vector()
    for value ∈ old
        if value ∈ removed # so note, removed takes precedence
            # TODO: If removed, delete all links.  Possibly go further and delete all variables linked to it.
            continue
        elseif value ∈ keys(modified)
            push!(new, modified[value])
        else
            push!(new, value)
        end
    end
    return [new ; added]


end




end



















