# # SEIR Full Model Measles Chickenpox

using StockFlow
using StockFlow.Syntax

using Catlab
using Catlab.CategoricalAlgebra
using LabelledArrays
using OrdinaryDiffEq
using Plots

using Catlab.Graphics
using Catlab.Programs
using Catlab.Theories
using Catlab.WiringDiagrams

seir = @stock_and_flow begin
    :stocks
    S
    E
    I
    R

    :parameters
    μ
    β
    tlatent
    trecovery
    δ

    :flows
    CLOUD => fbirth(μ * N) => S
    S => fincid(β * S * I / N) => E
    S => fdeathS(S * δ) => CLOUD
    E => finf(E / tlatent) => I
    E => fdeathE(E * δ) => CLOUD
    I => frec(I / trecovery) => R
    I => fdeathI(I * δ) => CLOUD
    R => fdeathR(R * δ) => CLOUD

    :sums
    N = [S, E, I, R]

end



GraphF(seir)

# define parameter values and initial values of stocks
# define constant parameters
p_measles = LVector(
    β=49.598, μ=0.03/12, δ=0.03/12, tlatent=8.0/30, trecovery=5.0/30
)
# define initial values for stocks
u0_measles = LVector(
    S=90000.0-930.0, E=0.0, I=930.0, R=773545.0
)

# solve the ODEs
# The model results are compared with the same model built by Anylogic, and the resules are the same!
prob_measles = ODEProblem(vectorfield(seir),u0_measles,(0.0,120.0),p_measles);
sol_measles = solve(prob_measles,Tsit5(),abstol=1e-8);
plot(sol_measles)

# define parameter values and initial values of stocks
# define constant parameters
p_chickenpox = LVector(
    β=18.0, μ=0.03/12.0, δ=0.03/12.0, tlatent=14.0/30.0, trecovery=5.0/30.0
)
# define initial values for stocks
u0_chickenpox = LVector(
    S=295354.0, E=0.0, I=1000.0, R=567191.0
)

# solve the ODEs
# The model results are compared with the same model built by Anylogic, and the resules are the same!
prob_chickenpox = ODEProblem(vectorfield(seir),u0_chickenpox,(0.0,120.0),p_chickenpox);
sol_chickenpox = solve(prob_chickenpox,Tsit5(),abstol=1e-8);
plot(sol_chickenpox)

sol_chickenpox

# to have the figures plotted fix to the wider of the cells
HTML("""
<style>
.output_svg div{
  width: 100% !important;
  height: 100% !important;
}
</style>
""")



