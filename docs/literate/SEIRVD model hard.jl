using StockFlow

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.WiringDiagrams

finc(u,p,t)=p.cβ*u.S*u.I/p.N
finf(u,p,t)=u.E * p.rlatent
frec(u,p,t)=u.I * p.rrec

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)


Graph()

fvac(u,p,t)=u.S * p.rv

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)


fdeath(u,p,t)=u.I * p.rd

# StockAndFlowp(stocks,
#               (flow=>function, upstream=>downstream) => stocks linked)


uwd_seirvd = @relation (S, I) begin

end;
display_uwd(uwd_seirvd)

seirvd=oapply(uwd_seirvd,Dict()) |> apex
Graph(seirvd)

p = LVector(
    cβ=0.2, N=1000, rrec=0.083, rv=0.02, rlatent=0.2, rd=0.0001
)
# define initial values for stocks
u0 = LVector(
    S=990, E=0, I=10, R=0, V=0, D=0
)

prob = ODEProblem(vectorfield(seirvd),u0,(0.0,100.0),p);
sol = solve(prob,Tsit5(),abstol=1e-8);
plot(sol)
