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
using Catlab.Theories
using Catlab.WiringDiagrams

sir = @stock_and_flow begin

    :stocks
    S
    I 
    R 

    :parameters
    rbirth
    cbeta
    rrecovery
    rdeath

    :dynamic_variables
    v_inf₁ = I / N
    v_inf₂ = v_inf₁ * cbeta
    v_inf₃ =  v_inf₂ * S

    :flows
    CLOUD => f_birth(rbirth * N) => S 
    S => f_inf(v_inf₃) => I 
    I => f_rec(rrecovery * I) => R
    S => f_deathS(S * rdeath) => CLOUD
    I => f_deathI(I * rdeath) => CLOUD
    R => f_deathR(R * rdeath) => CLOUD
    

    :sums
    N = [I, R, S]

end

GraphF(sir)

svi = @stock_and_flow begin
    
    :stocks
    S
    V
    I

    :parameters
    rvaccine
    rdeath
    lambda
    evaccine_complement # 1.0 - evaccine

    :dynamic_variables
    v_vacV = evaccine_complement * V
    v_infV = v_vacV * lambda

    
    :flows
    S => f_vacc(S * rvaccine) => V
    V => f_deathV(V * rdeath) => CLOUD
    V => f_infV(v_infV) => I 

    :sums
    N = [V, I, S]

end

GraphF(svi)

# define the UWD-algebra
sirv_uwd = @relation (S,I) begin
    seir(S,I)
    svi(S,I)
end;
display_uwd(sirv_uwd)

# define a foot of the structured multicospan
footS=foot(:S, :N, :S=>:N)
GraphF(footS;schema="C0")
     

# define a foot of the structured multicospan
footI=foot(:I, :N, :I=>:N)
GraphF(footI;schema="C0")
     

# open sir and svi stock and flow diagram with the feet defined before
open_sir=Open(sir, footS, footI)
open_svi=Open(svi, footS, footI)
# Compose those two models according the UWD-algebra
open_sirv = oapply(sirv_uwd, [open_sir, open_svi])
# the composed stock and flow diagram is the apex of the composed open stock and flow diagram
sirv = apex(open_sirv)

GraphF(sirv)

L = @stock_and_flow begin
    :stocks
    V
    I

    :parameters
    lambda
    cbeta
    evaccine_complement

    :dynamic_variables
    v_inf₁ = I / N
    v_inf₂ = v_inf₁ * cbeta
    
    v_vacV = evaccine_complement * V


    v_infV = v_vacV * lambda


    :flows
    V => f_infV(v_infV) => I

    :sums
    N = [V, I]
    end;
     


GraphF(L)

I = @stock_and_flow begin
    :stocks
    V
    I

    :parameters
    cbeta
    evaccine_complement

    :dynamic_variables
    v_inf₁ = I / N
    v_inf₂ = v_inf₁ * cbeta
    v_vacV = evaccine_complement * V
    v_infV = *(v_vacV)


    :flows
    V => f_infV(v_infV) => I

    :sums
    N = [V, I]
    end;
     


GraphF(I)

R = @stock_and_flow begin
    :stocks
    V
    I

    :parameters
    cbeta
    evaccine_complement

    :dynamic_variables
    v_inf₁ = I / N
    v_inf₂ = v_inf₁ * cbeta

    v_vacV = evaccine_complement * V

    v_infV = v_vacV * v_inf₂

    :flows
    V => f_infV(v_infV) => I

    :sums
    N = [V, I]
    end;
     


GraphF(R)

using AlgebraicRewriting
using AlgebraicRewriting: rewrite
const hom = Catlab.CategoricalAlgebra.homomorphism
rule = Rule(hom(I,L), hom(I,R))

sirv_rewritten = rewrite(rule, sirv)
GraphF(sirv_rewritten)

# define values of constant parameters
p = LVector(
    cbeta=0.1, rbirth=0.001, rdeath=0.001, rrecovery=0.05, # for model sir
    rvaccine=0.01, evaccine=0.3, evaccine_complement = 0.7 # for model svi
)
# define initial values for stocks
u0 = LVector(
    S=990.0, I=10.0, R=0.0, V=0.0
)

prob_sirv = ODEProblem(vectorfield(sirv_rewritten),u0,(0.0,100.0),p);
sol_sirv = solve(prob_sirv,Tsit5(),abstol=1e-8);
plot(sol_sirv)

# to have the figures plotted fix to the wider of the cells
HTML("""
<style>
.output_svg div{
  width: 100% !important;
  height: 100% !important;
}
</style>
""")



