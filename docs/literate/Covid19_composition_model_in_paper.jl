# # COVID-19 Composition Model

using StockFlow

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.WiringDiagrams

display_uwd(ex) = to_graphviz(ex, box_labels=:name, junction_labels=:variable, edge_attrs=Dict(:len=>"1"));

# Define functions ϕ of flows in the SEIRH model
fNewIncidence(u,p,t)=p.β*u.S*u.I/p.N
fNewInfectious(u,p,t)=u.E*p.ri
fNewRecovery(u,p,t)=u.I/p.tr * (1.0-p.fH )
fWaningImmunityR(u,p,t)=u.R/p.tw
fHICUAdmission(u,p,t) = u.I/p.tr * p.fH * p.fICU
fHNICUAdmission(u,p,t) = u.I/p.tr * p.fH * (1.0-p.fICU)
fOutICU(u,p,t) = u.HICU/p.tICU
fRecoveryH(u,p,t)= u.HNICU/p.tH

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)
seirh = StockAndFlowp((:S, :E, :I, :R, :HICU, :HNICU), 
   ((:NewIncidence=>fNewIncidence, :S=>:E) => (:S, :I), 
    (:NewInfectious=>fNewInfectious, :E=>:I) => :E,
    (:NewRecovery=>fNewRecovery, :I=>:R) => :I, 
    (:WaningImmunityR=>fWaningImmunityR, :R=>:S) => :R,
    (:HICUAdmission=>fHICUAdmission, :I=>:HICU) => :I, 
    (:HNICUAdmission=>fHNICUAdmission, :I=>:HNICU) => :I, 
    (:OutICU=>fOutICU, :HICU=>:HNICU) => :HICU,
    (:RecoveryH=>fRecoveryH, :HNICU=>:R) => :HNICU))

# Graph(primitive stock-flow model, direction of the diagram - the default value is "LR" from left to right; 
# users could also use "TB" from top to bottom)
# Graph(seirh)

# Define functions ϕ of flows in the Vaccine model
fFirstdoseVaccine(u,p,t) = u.S * p.rv
fSeconddoseVaccine(u,p,t) = u.VP * p.rv
fWaningImmunityVP(u,p,t) = u.VP / p.tw
fWaningImmunityVF(u,p,t) = u.VF / p.tw
fNewIncidenceVP(u,p,t) = p.β*u.VP*u.I*(1.0-p.eP)/p.N
fNewIncidenceVF(u,p,t) = p.β*u.VF*u.I*(1.0-p.eF)/p.N

# Stock and flow diagram of Vaccine model
v = StockAndFlowp((:S, :E, :I, :VP, :VF), 
   ((:FirstdoseVaccine=>fFirstdoseVaccine, :S=>:VP) => :S, 
    (:SeconddoseVaccine=>fSeconddoseVaccine, :VP=>:VF) => :VP, 
    (:WaningImmunityVP=>fWaningImmunityVP, :VP=>:S) => :VP,
    (:WaningImmunityVF=>fWaningImmunityVF, :VF=>:VP) => :VF,
    (:NewIncidenceVP=>fNewIncidenceVP, :VP=>:E) => (:VP, :I),
    (:NewIncidenceVF=>fNewIncidenceVF, :VF=>:E) => (:VF, :I)))

# Graph(v,"TB")

# Define functions ϕ of flows in the Persist Asymptomaticity model
fNewPersistentAsymptomaticity(u,p,t) = u.E * p.ria
fNewRecoveryIA(u,p,t) = u.IA / p.tr

# Stock and flow diagram of Persistent Asymptomaticity Model
ia = StockAndFlowp((:E, :IA, :R), 
   ((:NewPersistentAsymptomaticity=>fNewPersistentAsymptomaticity, :E=>:IA) => :E, 
    (:NewRecoveryIA=>fNewRecoveryIA, :IA=>:R) => :IA))

# Graph(ia)

covid = @relation (S, E, I, R) begin
    modelA(S,E,I,R)
    modelB(S,E,I)
    modelC(E,R)
end;
display_uwd(covid)

# Open three Stock and Flow Diagrams
openseirh = Open(seirh, [:S], [:E], [:I], [:R])
openv = Open(v, [:S], [:E], [:I])
openia = Open(ia, [:E], [:R])
# Compose those three models according the UWD-algebra
openCOVID19 = oapply(covid, [openseirh, openv, openia])
# Generate the composed model (Stock and Flow Diagram)
COVID19 = apex(openCOVID19)

# Graph(COVID19)

# Define constant parameters
p_COVID19 = LVector(
    β=0.8, N=38010001.0, tr=12.22, tw=2*365.0,
    fH=0.002, fICU=0.23, tICU=6.0, tH = 12.0,
    rv=0.01, eP=0.6, eF=0.85, ri=0.207, ria=0.138  
)
# Define initial values for stocks
u0_COVID19 = LVector(
    S=38010000.0, E=0.0, I=1.0, IA=0.0, R=0.0, HICU=0.0, HNICU=0.0, VP=0.0, VF=0.0
)

# Check the dependencies of links of all flows' functions
# checkfls(COVID19, u0_COVID19, p_COVID19)

# Solve the ODEs
prob_COVID19 = ODEProblem(vectorfield(COVID19),u0_COVID19,(0.0,300.0),p_COVID19);
sol_COVID19 = solve(prob_COVID19,Tsit5(),abstol=1e-8);
plot(sol_COVID19)

# Flow of ```FirstdoseVaccine```
ϕFirstdoseVaccine = map(x->fFirstdoseVaccine(sol_COVID19.u[x],p_COVID19,sol_COVID19.t[x]),collect(1:length(sol_COVID19.t)))
