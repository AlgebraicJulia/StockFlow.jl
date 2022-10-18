export convertStockFlowToSystemStructure, convertSystemStructureToStockFlow, rebuildSystemStructureByFlattenSymbols

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

function rebuildSystemStructureByFlattenSymbols(p::AbstractStockAndFlowStructure)
    s=extracStocksStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    
    return StockAndFlowStructure(s,f,sv)
end


function convertSystemStructureToStockFlow(p::AbstractStockAndFlowStructure,v)   
    s=extracStocksStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    
    return StockAndFlow(s,f,v,sv)
end


function convertStockFlowToSystemStructure(p::AbstractStockAndFlow)
    
    s=extracStocksStructureAndFlatten(p)
    f=extracFlowsStructureAndFlatten(p)
    sv=extracSumVStructureAndFlatten(p)
    return StockAndFlowStructure(s,f,sv)
end



