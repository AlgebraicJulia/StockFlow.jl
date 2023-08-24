using Pkg;
Pkg.activate(".")

using Test

using StockFlow
using StockFlow.Syntax
using StockFlow.PremadeModels

using Catlab.WiringDiagrams
using Catlab.ACSets
using Catlab.CategoricalAlgebra

l_type = @stock_and_flow begin 
    :stocks
    pop
    
    :parameters
    μ
    δ
    rFstOrder
    rage
    
    :dynamic_variables
    v_aging = pop * rage
    v_fstOrder = pop * rFstOrder
    v_birth = N * μ
    v_death = pop * δ
    
    :flows
    pop => f_aging(v_aging) => pop
    pop => f_fstOrder(v_fstOrder) => pop
    CLOUD => f_birth(v_birth) => pop
    pop => f_death(v_death) => CLOUD
    
    :sums
    N = [pop]
    
end;
l_type_noatts = map(l_type, Name=name->nothing, Op=op->nothing, Position=pos->nothing);


WeightModel = @stock_and_flow begin
    :stocks
    NormalWeight
    OverWeight
    Obese
    
    :parameters
    μ
    δw
    rw
    ro
    δo
    rage
    
    :dynamic_variables
    v_NewBorn = N * μ
    v_DeathNormalWeight = NormalWeight * δw
    v_BecomingOverWeight = NormalWeight * rw
    v_DeathOverWeight = OverWeight * δw
    v_BecomingObese = OverWeight * ro
    v_DeathObese = Obese * δo
    v_idNW = NormalWeight * rage
    v_idOW = OverWeight * rage
    v_idOb = Obese * rage
    
    :flows
    CLOUD => f_NewBorn(v_NewBorn) => NormalWeight
    NormalWeight => f_DeathNormalWeight(v_DeathNormalWeight) => ClOUD
    NormalWeight => f_BecomingOverWeight(v_BecomingOverWeight) => OverWeight
    OverWeight => f_DeathOverWeight(v_DeathOverWeight) => CLOUD
    
    OverWeight => f_BecomingObese(v_BecomingObese) => Obese
    Obese => f_DeathObese(v_DeathObese) => CLOUD
    NormalWeight => f_idNW(v_idNW) => NormalWeight
    OverWeight => f_idOW(v_idOW) => OverWeight
    Obese => f_idOb(v_idOb) => Obese
    
    :sums
    N = [NormalWeight, OverWeight, Obese]
    
end;


ageWeightModel = @stock_and_flow begin
    :stocks
    Child
    Adult
    Senior
    
    :parameters
    μ
    δC
    δA
    δS
    rageCA
    rageAS
    r
    
    :dynamic_variables
    v_NB = N * μ
    v_DeathC = Child * δC
    v_idC = Child * r
    v_agingCA = Child * rageCA
    v_DeathA = Adult * δA
    v_idA = Adult * r
    v_agingAS = Adult * rageAS
    v_DeathS = Senior * δS
    v_idS = Senior * r
    
    :flows
    CLOUD => f_NB(v_NB) => Child
    Child => f_idC(v_idC) => Child
    Child => f_DeathC(v_DeathC) => CLOUD
    Child => f_agingCA(v_agingCA) => Adult
    Adult => f_idA(v_idA) => Adult
    Adult => f_DeathA(v_DeathA) => CLOUD
    Adult => f_agingAS(v_agingAS) => Senior
    Senior => f_idS(v_idS) => Senior
    Senior => f_DeathS(v_DeathS) => CLOUD
    
    :sums
    N = [Child, Adult, Senior]
    



    
end;

begin
    s, = parts(l_type, :S)
    N, = parts(l_type, :SV)
    lsn, = parts(l_type, :LS)
    f_aging, f_fstorder, f_birth, f_death = parts(l_type, :F)
    i_aging, i_fstorder, i_birth = parts(l_type, :I)
    o_aging, o_fstorder, o_death = parts(l_type, :O)
    v_aging, v_fstorder, v_birth, v_death = parts(l_type, :V)
    lv_aging1, lv_fstorder1, lv_death1 = parts(l_type, :LV)
    lsv_birth1, = parts(l_type, :LSV)
    p_μ, p_δ, p_rfstOrder, p_rage = parts(l_type, :P)
    lpv_aging2, lpv_fstorder2, lpv_birth2, lpv_death2 = parts(l_type, :LPV)
end;

typed_WeightModel=ACSetTransformation(WeightModel, l_type_noatts,
  S = [s,s,s],
  SV = [N],
  LS = [lsn,lsn,lsn],   
  F = [f_birth, f_death, f_fstorder, f_death, f_fstorder, f_death, f_aging, f_aging, f_aging],    
  I = [i_birth, i_aging, i_fstorder, i_aging, i_fstorder, i_aging], 
  O = [o_death, o_fstorder, o_aging, o_death, o_fstorder, o_aging, o_death, o_aging],
  V = [v_birth, v_death, v_fstorder, v_death, v_fstorder, v_death, v_aging, v_aging, v_aging],
  LV = [lv_death1, lv_fstorder1, lv_death1, lv_fstorder1, lv_death1, lv_aging1, lv_aging1, lv_aging1],
  LSV = [lsv_birth1],
  P = [p_μ, p_δ, p_rfstOrder, p_rfstOrder, p_δ, p_rage],
  LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_aging2, lpv_aging2, lpv_aging2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_WeightModel);


typed_ageWeightModel=ACSetTransformation(ageWeightModel, l_type_noatts,
  S = [s,s,s],
  SV = [N],
  LS = [lsn,lsn,lsn],   
  F = [f_birth, f_fstorder, f_death, f_aging, f_fstorder, f_death, f_aging, f_fstorder, f_death],    
  I = [i_birth, i_fstorder, i_aging, i_fstorder, i_aging, i_fstorder], 
  O = [o_fstorder, o_death, o_aging, o_fstorder, o_death, o_aging, o_fstorder, o_death],
  V = [v_birth, v_death, v_fstorder, v_aging, v_death, v_fstorder, v_aging, v_death, v_fstorder],
  LV = [lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1],
  LSV = [lsv_birth1],
  P = [p_μ, p_δ, p_δ, p_δ, p_rage, p_rage, p_rfstOrder],
  LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2],
  Name = name -> nothing, Op=op->nothing, Position=pos->nothing
);
@assert is_natural(typed_ageWeightModel);

aged_weight = pullback(typed_WeightModel, typed_ageWeightModel) |> apex |> rebuildStratifiedModelByFlattenSymbols;

#########################################

age_weight_2 = @eval (@stratify (WeightModel, l_type, ageWeightModel) begin
    :stocks
    NormalWeight, OverWeight, Obese => pop <= Child, Adult, Senior

    :flows
    f_NewBorn => f_birth <= f_NB
    f_DeathNormalWeight, f_DeathOverWeight, f_DeathObese => f_death <= f_DeathC, f_DeathA, f_DeathS
    f_idNW, f_idOW, f_idOb => f_aging <= f_agingCA, f_agingAS
    f_BecomingOverWeight, f_BecomingObese => f_fstOrder <= f_idC, f_idA, f_idS

    :dynamic_variables
    v_NewBorn => v_birth <= v_NB
    v_DeathNormalWeight, v_DeathOverWeight, v_DeathObese => v_death <= v_DeathC, v_DeathA, v_DeathS
    v_idNW, v_idOW, v_idOb  => v_aging <= v_agingCA, v_agingAS
    v_BecomingOverWeight, v_BecomingObese => v_fstOrder <= v_idC, v_idA, v_idS

    :parameters
    μ => μ <= μ
    δw, δo => δ <= δC, δA, δS
    rw, ro => rFstOrder <= r
    rage => rage <= rageCA, rageAS

    :sums
    N => N <= N
end)
#########################################


age_weight_3 = @eval (@stratify (WeightModel, l_type, ageWeightModel) begin
    :stocks
    _ => pop <= _

    :flows
    f_NewBorn => f_birth <= f_NB
    _Death => f_death <= _Death
    _id => f_aging <= _aging
    _Becoming => f_fstOrder <= _id

    :dynamic_variables
    v_NewBorn => v_birth <= v_NB
    _Death => v_death <= _Death
    _id  => v_aging <= _aging
    _Becoming => v_fstOrder <= _id

    :parameters
    μ => μ <= μ
    _δ => δ <= _δ
    rage => rage <= rageCA, rageAS
    _ => rFstOrder <= _

    :sums
    _ => N <= _


end) 

    
age_weight_4 = @eval (@stratify (WeightModel, l_type, ageWeightModel) begin


    :flows
    _MATCHESNOTHING => f_birth <= _ALSOMATCHESNOTHING
    _New => f_birth <= _N
    _id => f_aging <= _aging
    _Death => f_death <= _Death
    _ => f_fstOrder <= _
    # everything has matched at this point, so the following are ignored.  Keys must exist, unless matching on _.  It can optionally throw an error with STRICT_MATCHES = true
    _ => f_fstOrder <= _
    f_DeathNormalWeight => f_fstOrder <= _aging

    :parameters
    μ => μ <= μ
    _δ => δ <= _δ
    rage => rage <= rageCA, rageAS
    _ => rFstOrder <= _

    _foo => rage <= _baz # would error if just did foo, as opposed to _foo.


    :dynamic_variables
    v_NewBorn => v_birth <= v_NB
    _Death => v_death <= _Death
    _id  => v_aging <= _aging
    _ => v_fstOrder <= _



end) 



@testset "Pullback computed in standard way is equal to DSL pullbacks" begin
    @test aged_weight == age_weight_2
    @test aged_weight == age_weight_3
    @test aged_weight == age_weight_4
    
    
    

end



