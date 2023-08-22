export convertStockFlowToSystemStructure, convertSystemStructureToStockFlow, rebuildStratifiedModelByFlattenSymbols,
extracStocksStructureAndFlatten,extracFlowsStructureAndFlatten,extracSumVStructureAndFlatten,extracVStructureAndFlatten,extracPsStructureAndFlatten,
extracVAndAttrStructureAndFlatten, extracVStructureAndFlatten, args_vname, args, add_prefix!, add_suffix!

flattenTupleNames(sn::Tuple)=Symbol(foldr(string,map(x->string(x), sn)))
flattenTupleNames(sn::Symbol)=sn
flattenTupleNames(sn::SubArray{Symbol})=sn
flattenTupleNames(sn::SubArray)=map(y->Symbol(foldr(string,map(x->string(x), y))),sn)


function extracStocksStructureAndFlatten(p::AbstractStockAndFlowStructure)
    s=[]
    
    for is in 1:ns(p)
        sn=sname(p,is)
        sn=flattenTupleNames(sn)
        
        ifs=inflows(p,is)
        ofs=outflows(p,is)
        vss=vsstock(p,is)
        svss=svsstock(p,is)
        
        ifns=isempty(ifs) ? :F_NONE : fname(p,ifs)
        ofns=isempty(ofs) ? :F_NONE : fname(p,ofs)
        vsns=isempty(vss) ? :V_NONE : vname(p,vss)
        svsns=isempty(svss) ? :SV_NONE : svname(p,svss)
        
        ifns=flattenTupleNames(ifns)
        ofns=flattenTupleNames(ofns)
        vsns=flattenTupleNames(vsns)
        svsns=flattenTupleNames(svsns)

        ss=sn=>(ifns,ofns,vsns,svsns)
        s=vcat(s,ss)
    end
    
    return s
end

function extracStocksStructureAndFlatten(p::AbstractStockAndFlowStructureF)
    s=[]
    
    for is in 1:ns(p)
        sn=sname(p,is)
        sn=flattenTupleNames(sn)
        
        ifs=inflows(p,is)
        ofs=outflows(p,is)
        svss=svsstock(p,is)
        
        ifns=isempty(ifs) ? :F_NONE : fname(p,ifs)
        ofns=isempty(ofs) ? :F_NONE : fname(p,ofs)
        svsns=isempty(svss) ? :SV_NONE : svname(p,svss)
        
        ifns=flattenTupleNames(ifns)
        ofns=flattenTupleNames(ofns)
        svsns=flattenTupleNames(svsns)

        ss=sn=>(ifns,ofns,svsns)
        s=vcat(s,ss)
    end
    
    return s
end

function extracFlowsStructureAndFlatten(p::AbstractStockAndFlowStructure)
    f=[]
    
    if nf(p)>0
        for ifl in 1:nf(p)
            fn=fname(p,ifl)
            vn=vname(p,fv(p,ifl))
            
            fn=flattenTupleNames(fn)
            vn=flattenTupleNames(vn)           
            
            fvs=fn=>vn
            f=vcat(f,fvs)
        end
    end
    
    return f
end

function extracPsStructureAndFlatten(p::AbstractStockAndFlowStructureF)
    pns=[]
    
    if np(p)>0
        for pr in 1:np(p)
            pn=pname(p,pr)            
            pn=flattenTupleNames(pn)        
            pns=vcat(pns,pn)
        end
    end
    
    return pns
end

function extracSumVStructureAndFlatten(p::AbstractStockAndFlowStructure)
    sv=[]
    
    if nsv(p)>0
        for svi in 1:nsv(p)
            svn=svname(p,svi)
            vsvs=vssv(p,svi)
            vsvns=isempty(vsvs) ? :SVV_NONE : vname(p,vsvs)
            
            svn=flattenTupleNames(svn)
            vsvns=flattenTupleNames(vsvns)
            
            svs=svn=>vsvns
            sv=vcat(sv,svs)
        end        
    end
    
    return sv
end

function extracSumVStructureAndFlatten(p::AbstractStockAndFlowStructureF)
    sv=[]
    
    if nsv(p)>0
        for svi in 1:nsv(p)
            svn=svname(p,svi)            
            svn=flattenTupleNames(svn)            
            sv=vcat(sv,svn)
        end        
    end
    
    return sv
end

function args_vname(p::AbstractStockAndFlowStructureF,v)
    srcsv=map(i->(flattenTupleNames(sname(p,i))),stocksv(p,v))
    srcsvv=map(i->(flattenTupleNames(svname(p,i))),svsv(p,v))
    srcpv=map(i->(flattenTupleNames(pname(p,i))),vpsrc(p,v))
    srcvv=map(i->(flattenTupleNames(vname(p,i))),vsrc(p,v))

    return (srcsv,srcsvv,srcpv,srcvv)
end

function args_vnamexxxx(p::AbstractStockAndFlowStructureF,v)
    srcsv=map(i->Symbol(join(sname(p,i))),stocksv(p,v))
    srcsvv=map(i->Symbol(join(svname(p,i))),svsv(p,v))
    srcpv=map(i->Symbol(join(pname(p,i))),vpsrc(p,v))
    srcvv=map(i->Symbol(join(vname(p,i))),vsrc(p,v))

    return (srcsv,srcsvv,srcpv,srcvv)
end

function args(p::AbstractStockAndFlowStructureF,v)
    (srcsv,srcsvv,srcpv,srcvv)=args_vname(p,v)
    return vcat(srcsv,srcsvv,srcpv,srcvv)
end

function args(p::AbstractStockAndFlowF,v)
    (srcsv,srcsvv,srcpv,srcvv)=args_vname(p,v)
       
    lvvp=lvvposition(p,v)
    lvtgtp=lvtgtposition(p,v)
    lsvvp=lsvvposition(p,v)
    lpvvp=lpvvposition(p,v)
    
    # create dictionary of (key=position, value=symbole of source argument)
    position_src=merge(make_dict(lvvp,srcsv),make_dict(lsvvp,srcsvv),make_dict(lpvvp,srcpv),make_dict(lvtgtp,srcvv))    
    ordered_position_src=sort(collect(position_src), by = x->x[1])    
    srcs=map(x->last(x),ordered_position_src)
    
    return srcs
end

extracVStructureAndFlatten(p::AbstractStockAndFlowStructureF) = begin

    vs=[]
    
    if nvb(p)>0
        for v in 1:nvb(p)
            vn = flattenTupleNames(vname(p,v))
            vnp = vn=>args(p,v)
            vs = vcat(vs,vnp)
        end
    end
    return vs
    
end

extracVAndAttrStructureAndFlatten(p::AbstractStockAndFlowF) = begin

    vs=[]
    
    if nvb(p)>0
        for v in 1:nvb(p)
            vn = flattenTupleNames(vname(p,v))
            v_op = foldr(==,vop(p,v)) ? vop(p,v)[1] : error("operators $(vop(p,v)) in the stratified model's auxiliary variable: $(join(vname(p,v))) should be the same!")
            vnp = vn=>(args(p,v)=>v_op)
            vs = vcat(vs,vnp)
        end
    end
    return vs
end


function rebuildStratifiedModelByFlattenSymbols(p::AbstractStockAndFlowStructure)
    s=extracStocksStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    
    return StockAndFlowStructure(s,f,sv)
end

function rebuildStratifiedModelByFlattenSymbols(p::AbstractStockAndFlowF)
    s=extracStocksStructureAndFlatten(p)
    pr=extracPsStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    v=extracVAndAttrStructureAndFlatten(p)
    
    return StockAndFlowF(s,pr,v,f,sv)
end


function convertSystemStructureToStockFlow(p::AbstractStockAndFlowStructure,v)   
    s=extracStocksStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    
    return StockAndFlow(s,f,v,sv)
end

function convertSystemStructureToStockFlow(p::AbstractStockAndFlowStructureF,v)   
    s=extracStocksStructureAndFlatten(p)
    pr=extracPsStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    
    return StockAndFlowF(s,pr,v,f,sv)
end

function convertStockFlowToSystemStructure(p::AbstractStockAndFlowF)
    
    s=extracStocksStructureAndFlatten(p)
    pr=extracPsStructureAndFlatten(p)
    v=extracVStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)

    return StockAndFlowStructureF(s,pr,v,f,sv)
end


"""
Concatenate Symbols.
"""
++(a::Symbol,b::Symbol) = Symbol(string(a, b))

"""
    add_suffix!(sf::AbstractStockAndFlow0, suffix)

Modify a AbstractStockAndFlow0 so named elements end with suffix.
Suffix can be anything which can be cast to a Symbol."""
function add_suffix!(sf::AbstractStockAndFlowStructureF, suffix)
    suffix = Symbol(suffix)
    set_snames!(sf, snames(sf) .++ suffix)
    set_fnames!(sf, fnames(sf) .++ suffix)
    set_svnames!(sf, svnames(sf) .++ suffix)
    set_vnames!(sf, vnames(sf) .++ suffix)
    set_pnames!(sf, pnames(sf) .++ suffix)
    return sf
end


"""
    add_suffix!(sf::AbstractStockAndFlow0, suffix)

Modify a AbstractStockAndFlow0 so named elements end with suffix.
Suffix can be anything which can be cast to a Symbol.
For feet.
"""
function add_suffix!(sf::AbstractStockAndFlow0, suffix)
    suffix = Symbol(suffix)
    set_snames!(sf, snames(sf) .++ suffix)
    set_svnames!(sf, svnames(sf) .++ suffix)
    return sf
end

"""
    add_prefix!(sf::AbstractStockAndFlowStructureF, prefix)

Modify a AbstractStockAndFlowStructureF so named elements begin with prefix
Prefix can be anything which can be cast to a Symbol.
"""
function add_prefix!(sf::AbstractStockAndFlowStructureF, prefix)
    prefix = Symbol(prefix)
    set_snames!(sf, prefix .++ snames(sf))
    set_fnames!(sf, prefix .++ fnames(sf))
    set_svnames!(sf, prefix .++ svnames(sf))
    set_vnames!(sf, prefix .++ vnames(sf))
    set_pnames!(sf, prefix .++ pnames(sf))
    return sf
end

"""
    add_prefix!(sf::AbstractStockAndFlowStructureF, prefix)

Modify a AbstractStockAndFlowStructureF so named elements begin with prefix
Prefix can be anything which can be cast to a Symbol.For feet.
"""
function add_prefix!(sf::AbstractStockAndFlow0, prefix)
    prefix = Symbol(prefix)
    set_snames!(sf, prefix .++ snames(sf))
    set_svnames!(sf, prefix .++ svnames(sf))
    return sf
end


