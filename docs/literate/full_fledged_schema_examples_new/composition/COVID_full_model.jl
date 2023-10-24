using GraphViz

using StockFlow
using StockFlow.Syntax

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.WiringDiagrams

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"));

seir = @stock_and_flow begin
    
    :stocks
    S
    E
    IA
    IYU
    IYN
    R
    
    :parameters
    β
    rLatent
    rIncubation
    rDevelop
    rRecovery
    rw
    
    :dynamic_variables
    v_NewIncidence₁ = β * NIC
    v_NewIncidence₂ = v_NewIncidence₁ / N # λ
    v_NewIncidence₃ = v_NewIncidence₂ * S


    v_NewInfectious = E * rLatent
    v_BecomingSymptomatic = IA * rIncubation
    v_SymptomicsNotDevelopingComplications = IYU * rDevelop
    v_NewRecovery = IYN * rRecovery
    v_WaningImmunityR = R * rw
    
    :flows
    S => fNewIncidence(v_NewIncidence₃) => E
    E => fNewInfectious(v_NewInfectious) => IA
    IA => fBecomingSymptomatic(v_BecomingSymptomatic) => IYU
    IYU => fSymptomicsNotDevelopingComplications(v_SymptomicsNotDevelopingComplications) => IYN
    IYN => fNewRecovery(v_NewRecovery) => R
    R => fWaningImmunityR(v_WaningImmunityR) => S
    
    
    
    :sums
    N = [S, E, IA, IYU, IYN, R]
    NIC = [IA, IYU, IYN]
    NI = [IA, IYU, IYN]
    
end


GraphF(seir)

v = @stock_and_flow begin
    :stocks
    S
    E
    VP
    VF
    
    
    :parameters
    rv
    rw
    λ
    eP_complement # 1.0 - eP
    eF_complement # 1.0 - eF


    :dynamic_variables
    v_NewIncidenceVP = VP * eP_complement
    v_NewIncidenceVF = VF * eF_complement

    v_infVP = v_NewIncidenceVP * λ
    v_infVF = v_NewIncidenceVF * λ

    
    :flows
    S => fFirstdoseVaccine(S * rv) => VP
    VP => fWaningImmunityVP(VP * rw) => S
    VP => fSeconddoseVaccine(VP * rv) => VF
    VF => fWaningImmunityVF(VF * rw) => VP
    VP => fNewIncidenceVP(v_infVP) => E
    VF => fNewIncidenceVF(v_infVF) => E

    
    
    :sums
    N = [S, E, VP, VF]
    NIC = []
end

GraphF(v;rd="TB")

ia = @stock_and_flow begin
    :stocks
    IA
    IA2
    IA3
    R
    
    :parameters
    rIncubationIA
    rDevelopIA2
    rRecoveryIA3
    
    :flows
    IA => fDevelopmentOfPersistentAsymptomaticity(IA * rIncubationIA) => IA2
    IA2 => fProgressionIA2_IA3(IA2 * rDevelopIA2) => IA3
    IA3 => fNewrecoveryIA3(IA3 * rRecoveryIA3) => R
    
    :sums
    N = [IA, IA2, IA3, R]
    NIC = [IA, IA2, IA3]
    NI = [IA, IA2, IA3]
end

GraphF(ia)

h = @stock_and_flow begin
    :stocks
     IYU
     HICU
     HNICU
     R
     D
     
     :parameters
     rAdmICU
     rAdmNICU
     rOutICU
     rrH
     rDeathICU
     rDeathNICU
     
     :dynamic_variables
     v_HICUAdmission = IYU * rAdmICU
     v_HNICUAdmission = IYU * rAdmNICU
     v_OutICU = HICU * rOutICU
     v_RecoveryH = HNICU * rrH
     v_NewDeathsHICU = HICU * rDeathICU
     v_NewDeathsHNICU = HNICU * rDeathNICU
     
     :flows
     IYU => f_HICUAdmission(v_HICUAdmission) => HICU
     IYU => f_HNICUAdmission(v_HNICUAdmission) => HNICU
     HICU => f_OutICU(v_OutICU) => HNICU
     HNICU => f_RecoveryH(v_RecoveryH) => R
     HICU => f_NewDeathsHICU(v_NewDeathsHICU) => D
     HNICU => f_NewDeathsHNICU(v_NewDeathsHNICU) => D
     
     :sums
     N = [IYU, HICU, HNICU, R]
     NIC = [IYU]
     NI = [IYU, HICU, HNICU]
     
     
 end

GraphF(h;rd="TB")

footIYUN=foot(:IYU, (:NI,:NIC,:N), (:IYU=>:NI, :IYU=>:NIC, :IYU=>:N))
GraphF(footIYUN;schema="C0")

footRN=foot(:R, :N, :R=>:N)
GraphF(footRN;schema="C0")


footIAN=foot(:IA, (:NI,:NIC,:N), (:IA=>:NI, :IA=>:NIC, :IA=>:N))
GraphF(footIAN;schema="C0")

footSN=foot(:S, :N, :S=>:N)
GraphF(footSN;schema="C0")

footEN=foot(:E, :N, :E=>:N)
GraphF(footEN;schema="C0")

footNIC=foot((),:NIC,())
GraphF(footNIC;schema="C0")

covid = @relation (footSN, footEN, footIAN, footIYUN, footRN, footNIC) begin
    modelA(footSN, footEN, footIAN, footIYUN, footRN, footNIC)
    modelB(footSN, footEN, footNIC)
    modelC(footIAN, footRN)
    modelD(footIYUN, footRN)
end;
display_uwd(covid)

open_modelA=Open(seir, footSN, footEN, footIAN, footIYUN, footRN, footNIC)
open_modelB=Open(v,footSN,footEN, footNIC)
open_modelC=Open(ia,footIAN,footRN)
open_modelD=Open(h,footIYUN,footRN)
# Compose those three models according the UWD-algebra
openCOVID19 = oapply(covid, [open_modelA, open_modelB, open_modelC, open_modelD])
# composed model
COVID19=apex(openCOVID19)

GraphF(COVID19)

L = @stock_and_flow begin
    :stocks
    VP
    VF
    E

    :parameters
    eP_complement
    eF_complement
    λ
    β


    :flows
    VP => fNewIncidenceVP(v_infVP) => E
    VF => fNewIncidenceVF(v_infVF) => E

    :dynamic_variables
    v_NewIncidence₁ = β * NIC
    v_NewIncidence₂ = v_NewIncidence₁ / N # λ

    v_NewIncidenceVP = VP * eP_complement
    v_NewIncidenceVF = VF * eF_complement

    v_infVP = v_NewIncidenceVP * λ
    v_infVF = v_NewIncidenceVF * λ
    

    :sums
    N = [VP, VF, E]
    NIC = []
end;

GraphF(L)

I = @stock_and_flow begin
    :stocks
    VP
    VF
    E

    :parameters
    eP_complement
    eF_complement
    β

    :flows
    VP => fNewIncidenceVP(v_infVP) => E
    VF => fNewIncidenceVF(v_infVF) => E

    :dynamic_variables
    v_NewIncidence₁ = β * NIC
    v_NewIncidence₂ = v_NewIncidence₁ / N # λ
    v_NewIncidenceVP = VP * eP_complement
    v_NewIncidenceVF = VF * eF_complement


    v_infVP = *(v_NewIncidenceVP)
    v_infVF = *(v_NewIncidenceVF)

    :sums
    N = [VP, VF, E]
    NIC = []
end;

GraphF(I)

R = @stock_and_flow begin
    :stocks
    VP
    VF
    E

    :parameters
    eP_complement
    eF_complement
    β

    :flows
    VP => fNewIncidenceVP(v_infVP) => E
    VF => fNewIncidenceVF(v_infVF) => E

    :dynamic_variables
    v_NewIncidence₁ = β * NIC
    v_NewIncidence₂ = v_NewIncidence₁ / N # λ

    v_NewIncidenceVP = VP * eP_complement
    v_NewIncidenceVF = VF * eF_complement


    v_infVP = v_NewIncidenceVP * v_NewIncidence₂
    v_infVF = v_NewIncidenceVF * v_NewIncidence₂


    :sums
    N = [VP, VF, E]
    NIC = []
end;

GraphF(R)

using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = Catlab.CategoricalAlgebra.homomorphism

rule = Rule(hom(I,L), hom(I,R))

Covid19_rewritten = rewrite(rule, COVID19)
GraphF(Covid19_rewritten)

Covid19_rewritten

GraphF(Covid19_rewritten; type="SF", rd="TB")

# define constant parameters
p_COVID19_raw = LVector(
    β=0.8, tLatent=2.9, tIncubation=2.72, tDevelop=6.0, tRecovery=3.5,tw=2*365.0,
    fH=0.002, fICU=0.23, tICU=6.0, tH = 12.0, tOutICU=6.0, fractionIA=0.4
)


p_COVID19 = LVector(
    β=p_COVID19_raw.β, rLatent=1.0/p_COVID19_raw.tLatent, rIncubation=(1.0-p_COVID19_raw.fractionIA)/p_COVID19_raw.tIncubation,
    rDevelop=(1.0-p_COVID19_raw.fH)/p_COVID19_raw.tDevelop, rRecovery=1.0/p_COVID19_raw.tRecovery, rw=1.0/p_COVID19_raw.tw,
    rv=0.01, eP=0.6, eF=0.85, rIncubationIA=p_COVID19_raw.fractionIA/p_COVID19_raw.tIncubation, 
    rDevelopIA2=1.0/p_COVID19_raw.tDevelop, rRecoveryIA3=1.0/p_COVID19_raw.tRecovery,
    rAdmICU=p_COVID19_raw.fH*p_COVID19_raw.fICU/p_COVID19_raw.tDevelop,
    rAdmNICU=p_COVID19_raw.fH*(1.0-p_COVID19_raw.fICU)/p_COVID19_raw.tDevelop,
    rrH=1.0/p_COVID19_raw.tH, rOutICU=1.0/p_COVID19_raw.tOutICU,rDeathICU=0.085, rDeathNICU=0.018,
    eP_complement = 0.4, eF_complement = 0.15
)
# define initial values for stocks
u0_COVID19 = LVector(
    S=38010000.0, E=0.0, IYU=10.0, IA=0.0, R=0.0, HICU=0.0, HNICU=0.0, VP=0.0, VF=0.0,
    D=0.0, IA2=0.0, IA3=0.0, IYN=0.0
)

# results are tested the same as the Anylogic model
prob_COVID19 = ODEProblem(vectorfield(Covid19_rewritten),u0_COVID19,(0.0,100.0),p_COVID19);
sol_COVID19 = solve(prob_COVID19,Tsit5(),abstol=1e-8);
plot(sol_COVID19)

# to have the figures plotted fix to the wider of the cells
HTML("""
<style>
.output_svg div{
  width: 100% !important;
  height: 100% !important;
}
</style>
""")



