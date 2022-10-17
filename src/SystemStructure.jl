export convertStockFlowToSystemStructure, convertSystemStructureToStockFlow

function extracStocksStructure(p::AbstractStockAndFlowStructure)
    s=[]
    
    for is in 1:ns(p)
        sn=sname(p,is)
        
        ifs=inflows(p,is)
        ofs=outflows(p,is)
        vss=vsstock(p,is)
        svss=svsstock(p,is)
        
        ifns=isempty(ifs) ? :F_NONE : fname(p,ifs)
        ofns=isempty(ofs) ? :F_NONE : fname(p,ofs)
        vsns=isempty(vss) ? :V_NONE : vname(p,vss)
        svsns=isempty(svss) ? :SV_NONE : svname(p,svss)

        ss=sn=>(ifns,ofns,vsns,svsns)
        s=vcat(s,ss)
    end
    
    return s
end

function extracFlowsStructure(p::AbstractStockAndFlowStructure)
    f=[]
    
    if nf(p)>0
        for ifl in 1:nf(p)
            fn=fname(p,ifl)
            vn=vname(p,fv(p,ifl))
            fvs=fn=>vn
            f=vcat(f,fvs)
        end
    end
    
    return f
end

function extracSumVStructure(p::AbstractStockAndFlowStructure)
    sv=[]
    
    if nsv(p)>0
        for svi in 1:nsv(p)
            svn=svname(p,svi)
            vsvs=vssv(p,svi)
            vsvns=isempty(vsvs) ? :SVV_NONE : vname(p,vsvs)
            svs=svn=>vsvns
            sv=vcat(sv,svs)
        end        
    end
    
    return sv
end

function convertStockFlowToSystemStructure(p::AbstractStockAndFlow)
    
    s=extracStocksStructure(p)
    f=extracFlowsStructure(p)
    sv=extracSumVStructure(p)
    return StockAndFlowStructure(s,f,sv)
end

function convertSystemStructureToStockFlow(p::AbstractStockAndFlowStructure,v)   
    s=extracStocksStructure(p)
    f=extracFlowsStructure(p)
    sv=extracSumVStructure(p)
    
    return StockAndFlow(s,f,v,sv)
end

