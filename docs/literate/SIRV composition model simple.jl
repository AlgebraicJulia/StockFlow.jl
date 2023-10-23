using StockFlow

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.WiringDiagrams

# define functions ϕ of flows in the SIR model
fNewIncidence(u,p,t)=p.cβ*u.S*u.I/p.N
fNewRecovery(u,p,t)=u.I/p.tr

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)
sir = StockAndFlowp((:S, :I, :R), 
   ((:NewIncidence=>fNewIncidence, :S=>:I)=>(:S,:I),
    (:NewRecovery=>fNewRecovery, :I=>:R)=>:I)
)

Graph(sir)

# define functions ϕ of flows in the SVI model
fNewIncidenceFromV(u,p,t)=p.cβ*u.V*u.I*(1-p.e)/p.N
fNewVaccinated(u,p,t)=u.S*p.rv

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)
svi = StockAndFlowp(
)

Graph(svi)





sirv1=apex()

Graph(sirv1)

uwd_sirv = 
display_uwd(uwd_sirv)

sirv2=oapply(uwd_sirv,Dict()) |> apex

Graph(sirv2)

p_sirv = LVector(
    cβ=0.2, N=1000, tr=12, rv=0.02, e=0.9
)
# define initial values for stocks
u0_sirv = LVector(
    S=990, I=10, R=0, V=0
)

prob_sirv1 = ODEProblem(vectorfield(sirv1),u0_sirv,(0.0,100.0),p_sirv);
sol_sirv1 = solve(prob_sirv1,Tsit5(),abstol=1e-8);
plot(sol_sirv1)

prob_sirv2 = ODEProblem(vectorfield(sirv2),u0_sirv,(0,100.0),p_sirv)
sol_sirv2 = solve(prob_sirv2,Tsit5(),sbstol=1e-8);
plot(sol_sirv2)
