module Rewrite


export sfrewrite


using StockFlow
using StockFlow.Syntax
import StockFlow.Syntax: parse_dyvar!, parse_flow!, parse_sums!
using MLStyle



struct RewriteModifications
    added::Vector{Dict{Symbol, Vector{Dict}}} # :S => [Dict(:sname => :pop), ...], :LS = [(:lss => 1, :lsv => 1), ...], ...
    removed::Vector{Dict{Symbol, Vector{Any}}} 
    swapped::Dict{Symbol, Symbol} # to be parsed at the end and be added to added/removed, if hasn't been already.

    RewriteModifications() = new(
        Dict(:S => [], :LS => [], :SV => [], :LSV => [], :LV => [], :I => [], :O => [], :F => [], :V => [], :LVV => [], :LPV => [], :P => [], :vop => []),
        Dict(:S => [], :LS => [], :SV => [], :LSV => [], :LV => [], :I => [], :O => [], :F => [], :V => [], :LVV => [], :LPV => [], :P => [], :vop => []),
        Dict{Symbol, Symbol}()
    )
end
"""
Takes a symbol :S, :SV, :V, :P, and returns the corresponding object for a dv link
"""
function dv_link_type(object::Symbol)::Symbol
    @switch object begin
        @case :S => :LV
        @case :SV => :LSV
        @case :V => :LVV
        @case :P => :LPV
        @case _ => error("Unknown symbol found when looking for variable link: $object")
    end
end

function return_link_index(src::Symbol, tgt::Symbol, all_names::Dict{Symbol, Dict}, position::Int, sf::AbstractStockAndFlowF) # tgt is a dv
    link_type, other... = dv_link_type(all_names[src][1])
    srcindex = other[1]

    _, tgtindex, _, _ = all_names[tgt]
    

    match = (srcindex, tgtindex, position) 
    # TODO: REFACTOR TO NOT BE GARBAGE
    @switch link_type begin # probably more efficient to cache rather than find and create a vector every time.
        @case :LV => return only(filter((i, x) -> x == match, enumerate(zip(get_lvs(sf), get_lvv(sf), get_lvvposition(sf)))))[1]
        @case :LSV => return only(filter((i, x) -> x == match, enumerate(zip(get_lsvsv(sf), get_lsvv(sf), get_lsvvposition(sf)))))[1]
        @case :LVV => return only(filter((i, x) -> x == match, enumerate(zip(get_lvsrc(sf), get_lvtgt(sf), get_lvvposition(sf)))))[1] #TODO: Verify I got src and tgt right
        @case :LPV => return only(filter((i, x) -> x == match, enumerate(zip(get_lpvp(sf), get_lpvv(sf), get_lpvvposition(sf)))))[1]
        @case _ => error("Did not recognize link type $link_type")
    end


end

function add_to_dict(line, d_symbol::Symbol, modifications)
    @match line begin
        Expr(:call, op, val::Symbol) => begin
            if op == :+
                push!(modifications.added[d_symbol], val)
            else if op == :-
                push!(modifications.removed[d_symbol], val)
            else
                error("Unknown operator $op")
            end
        end
    end
end

# function add_flow()

# TODO: get the links before getting to this point and shove them in all_names
# EG, :LS = Dict(1 => (3, 4), 2 => (4, 5))
# means we don't need to use return_link_index.
function addswap!(sw, modifications, all_names, sf) # might need to go over this last, so we can add new names to this first.
    @match sw begin
        :($A => $B) => begin # stock => stock, flow => flow, dv => dv, dv (nonflow) => p, etc.
            typeA, otherA... = all_names[A]
            typeB, otherB... = all_names[B]
            @assert (typeA == typeB), "Cannot swap a $typeA for a $typeB !"
            # TODO: Add the ability to swap S, P and SV, implying you can move links from one to the other
            # DV will be more complicated.

            push!(modifications.swapped, (A => B)) # we then grab B using all_names
        end
        Expr(:(=), tgt, Expr(:call, :(=>), old, new)) => begin # sum = old => new, dv = old => new
            @match (old, new) begin
                (Expr(:vect, oldsum...), Expr(:vect, newsum...)) => begin 
                    for (i, (olds, news)) ∈ enumerate(zip(oldsum.args, newsum.args))
                        if olds != news
                            match = (all_names[olds], all_names[tgt])
                            old_link = only(filter((i, x) -> x == match, enumerate(zip(get_lss(sf), get_lssv(sf)))))[1]
                            push!(modifications.added[:LS], tgt => src)
                            push!(modification.deleted[:LS], old_link)
                    end

                end #TODO: THIS
                (Expr(:call, op1, oldoperands...), Expr(:call, op2, newoperands...)) => begin
                    #HAHHAHA THIS MEANS I CAN CHANGE THE OPERATOR TOO
                    # pretty sure changing the operator only requires a set_subpart! call
                    _, index, _, isflow = all_names[tgt] # already know tgt is a :V, and we already know its operator
                    if op1 != op2
                        push!(modifications.added["op"], all_names[tgt][2] => op2) # operator index
                        push!(modifications.removed["op"], all_names[tgt][2])

                        # push!(modifications.swapped, tgt => (tgt, :V, index, op2, isflow)) # things are gonna be really screwed up if you try swapping a flow variable
                        # Is that even possible?
                        # I think so, but the links are gonna be real tangled.
                    end 
                    for (i, (oldop, newop)) ∈ enumerate(zip(oldoperands.args, newoperands.args))
                        if oldop != newop # ok, so we're swapping one operand for another.  We're going to need to delete oldop's link from I, and add newop's to R.
                            oldoptype = all_names[oldop][1] #:V, :P, :SV, :S
                            newoptype = all_names[newop][1] # "

                            # at this point, we should know the name of the tgt, the type of source, and can know the position (if we enumerate the zip)

                            
                            oldoplink = dv_link_type(oldoptype)
                            newoplink = dv_link_type(newoptype)

                            link_index = return_link_index(oldop, tgt, all_names, i, sf)

                            push!(modifications.removed[oldoplink], link_index)

                            push!(modifications.added[newoplink], tgt => (newop, i))

                            

                            # Oh boy have fun doing math to determine what the link index is after adding all the new stuff and subtracting the old.
                            # OR: add the current link to delete, add a representation of the new to add, and deal with it later.



                    end
            end
        end
    end
end



function sfrewrite(sf, block)
    Base.remove_linenums!(block)

    sf_snames = Dict(S => i for (i, S) in enumerate(snames(sf)))
    sf_svnames = Dict(S => i for (i, S) in enumerate(svnames(sf)))
    sf_vnames = Dict(S => i for (i, S) in enumerate(vnames(sf)))
    sf_fnames = Dict(S => i for (i, S) in enumerate(fnames(sf)))
    sf_pnames = Dict(S => i for (i, S) in enumerate(pnames(sf)))

    
    

    name_vector = [collect(keys(names))... for names ∈ [sf_snames, sf_svnames, sf_vnames, sf_fnames, sf_pnames]]

    println(name_vector)
    @assert allunique(name_vector), "Not all names are unique!"
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

    L = StockAndFlowF()
    I = StockAndFlowF()
    R = StockAndFlowF()

    dyvars::Vector{Tuple{Symbol,Expr}} = []
    flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
    sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []

    # So the way this is gonna work:
    # prefix an item with ! to indicate it's being removed from I
    # prefix an item with ⊕ to indicate it's being added to R
    # under swapping:
    # - use X => Y to indicate all instances of X in sf are to be swapped with Y in R
    # - use v1 = A + B => C + B to indicate this particular link is being swapped
    # (this makes it easier to swap links)
    # when all links to an object have been removed, it can be inferred that the object itself can be removed, but this isn't currently done.

    # Once we gather all the information:
    # - Create stockflow L from all links attached to parts indicated to be removed and swapped, without them removed
    # - Create stockflow I from all links attached to parts to be removed or swapped, with the removed and swapped parts removed
    # - Create stockflow R from all links attached to parts to be added or swapped to.
    # Define hom(I, L), hom(I, R)
    # Rewrite
    
    modifications = RewriteModifications()


    current_phase = (_, _) -> ()
    for statement in block.args
        @match statement begin
            QuoteNode(:swapping) => begin
                current_phase = sw -> add_swap!(sw, modifications, all_names) # this one will be different than the others
            end
            QuoteNode(:stocks) => begin
                current_phase = s -> add_to_dict(statement, :S, modifications)

                # current_phase = s -> read_rewrite_line_and_update_dictionaries!(s, strata_snames, type_snames, aggregate_snames, strata_stock_mappings_dict, aggregate_stock_mappings_dict)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> add_to_dict(statement, :P, modifications)
                # current_phase = p -> read_rewrite_line_and_update_dictionaries!(p, strata_pnames, type_pnames, aggregate_pnames, strata_param_mappings_dict, aggregate_param_mappings_dict)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = d -> begin 
                    Base.remove_linenums!(d)
                    @match d begin # TODO: Make this format consistent, maybe
                        Expr(:call, :+, val) => push!(modifications.added[:V], val)
                        Expr(:call, :-, val) => push!(modifications.removed[:V], val)
                        Expr(:(=), Expr(:call, :+, tgt), Expr(:block, R)) => push!(modifications.added[:V], parse_dyvar!(Expr(:(=), tgt, R)))
                        Expr(:(=), Expr(:call, :-, tgt), Expr(:block, R)) => push!(modifications.removed[:V], parse_dyvar!(Expr(:(=), tgt, R)))
                        _ => error("No match on $d for dv")
                   end
                end
            end
            QuoteNode(:flows) => begin
                current_phase = f -> begin
                Base.remove_linenums!(f)
                @match f begin # TODO: Make this format consistent, maybe
                    Expr(:call, :+, val) => push!(modifications.added[:F], val)
                    Expr(:call, :-, val) => push!(modifications.removed[:F], val)
                    Expr(:(=), Expr(:call, :+, tgt), Expr(:block, R)) => push!(modifications.added[:F], parse_dyvar!(Expr(:(=), tgt, R)))
                    Expr(:(=), Expr(:call, :-, tgt), Expr(:block, R)) => push!(modifications.removed[:F], parse_dyvar!(Expr(:(=), tgt, R)))
                    _ => error("No match on $f for f")
               end
            end
            QuoteNode(:sums) => begin
                current_phase = s -> begin
                Base.remove_linenums!(sv)
                @match s begin # TODO: Make this format consistent, maybe
                    Expr(:call, :+, val) => push!(modifications.added[:SV], val)
                    Expr(:call, :-, val) => push!(modifications.removed[:SV], val)
                    Expr(:(=), Expr(:call, :+, tgt), Expr(:block, R)) => push!(modifications.added[:SV], parse_dyvar!(Expr(:(=), tgt, R)))
                    Expr(:(=), Expr(:call, :-, tgt), Expr(:block, R)) => push!(modifications.removed[:SV], parse_dyvar!(Expr(:(=), tgt, R)))
                    _ => error("No match on $sv for sv")
               end
            QuoteNode(kw) =>
                error("Unknown block type for stratify syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end
end


end