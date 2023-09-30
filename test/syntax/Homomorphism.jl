using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax: NothingFunction
using StockFlow.Syntax.Homomorphism

using StockFlow.PremadeModels


using Catlab.CategoricalAlgebra

@testset "hom macro creates correct homomorphisms" begin
  empty = @stock_and_flow begin end
  empty_hom = ACSetTransformation(empty, empty; :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction )
  @test (@hom empty empty begin end) == empty_hom

  sfA = @stock_and_flow begin; :stocks; A; end;
  sfB = @stock_and_flow begin; :stocks; B; end;
  @test (@hom sfA sfB begin; :stocks; A => B; end;) == ACSetTransformation(sfA, sfB ; :S => [1], :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)

  seir = PremadeModels.seir()

  @test (@hom seir seir begin
    :stocks
    S => S
    E => E
    I => I
    R => R

    :parameters
    μ => μ
    β => β
    tlatent => tlatent
    trecovery => trecovery
    δ => δ
    c => c

    :dynamic_variables
    v_prevalence => v_prevalence
    v_meanInfectiousContactsPerS => v_meanInfectiousContactsPerS
    v_perSIncidenceRate => v_perSIncidenceRate
    v_newIncidence => v_newIncidence
    v_birth => v_birth
    v_inf => v_inf
    v_rec => v_rec
    v_deathS => v_deathS
    v_deathE => v_deathE
    v_deathI => v_deathI
    v_deathR => v_deathR

    :flows
    f_birth => f_birth
    f_incid => f_incid
    f_deathS => f_deathS
    f_inf => f_inf
    f_deathE => f_deathE
    f_rec => f_rec
    f_deathI => f_deathI
    f_deathR => f_deathR

    :sums
    N => N
    NI => NI
    NS => NS
  end) == ACSetTransformation(seir, seir ; :S => (1:4), :F => (1:8), :V => (1:11),
    :SV => (1:3), :P => (1:6), :LS => (1:9), :I => (1:4), :O => (1:7), :LV => (1:7),
    :LSV => (1:3), :LVV => (1:3), :LPV => (1:9), 
    :Op => NothingFunction, :Position => NothingFunction, :Name => NothingFunction)



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

  # l_type_noatts = map(l_type, Name=NothingFunction, Op=NothingFunction, Position=NothingFunction);


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
    NormalWeight => f_DeathNormalWeight(v_DeathNormalWeight) => CLOUD
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

  typed_WeightModel=ACSetTransformation(WeightModel, l_type,
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
  Name=NothingFunction, Op=NothingFunction, Position=NothingFunction
);

(@hom WeightModel l_type begin
  :stocks
  NormalWeight => pop
  OverWeight => pop
  Obese => pop

  :parameters
  μ => μ
  δw => δ
  rw => rFstOrder
  ro => rFstOrder
  δo => δ
  rage => rage

  :dynamic_variables
  v_NewBorn => v_birth
  v_DeathNormalWeight => v_death
  v_BecomingOverWeight => v_fstOrder
  v_DeathOverWeight => v_death
  v_BecomingObese => v_fstOrder
  v_DeathObese => v_death
  v_idNW => v_aging
  v_idOW => v_aging
  v_idOb => v_aging

  :flows
  f_NewBorn => f_birth
  f_DeathNormalWeight => f_death
  f_BecomingOverWeight => f_fstOrder
  f_DeathOverWeight => f_death
  f_BecomingObese => f_fstOrder
  f_DeathObese => f_death
  f_idNW => f_aging
  f_idOW => f_aging
  f_idOb => f_aging

  :sums
  N => N

end) == typed_WeightModel






end