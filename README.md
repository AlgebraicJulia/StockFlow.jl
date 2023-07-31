# AlgebraicStockFlow
[![Tests](https://github.com/AlgebraicJulia/StockFlow.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/AlgebraicJulia/StockFlow.jl/actions/workflows/tests.yml)
<!-- TODO: Set up on codecov.io for repo [![codecov](https://codecov.io/github/AlgebraicJulia/StockFlow.jl/branch/maaster/graph/badge.svg)](https://app.codecov.io/github/AlgebraicJulia/StockFlow.jl) -->


StockFlow.jl is a Julia library for the creation and execution of [stock and flow modeling diagrams](https://en.wikipedia.org/wiki/System_dynamics#Stock_and_flow_diagrams), designed with model composition and stratification in mind through a [categorical approach](https://arxiv.org/abs/2211.01290).

Stock-flow diagrams are used to represent systems where its population can move between different states, such as tracking how many people have been infected or recovered from a disease.  By providing initial values, and values for the parameters, we can solve the differential equation which the stock-flow diagram represents, which gives us a graph for how the population varies over time.

In this particular schema, stock-flow diagrams consist of stocks, flows, dynamic variables, sum variables, parameters, and the links between them.

* Stocks represent accumulations of a population in a particular state, such as 'susceptible', 'infected', 'recovered', etc.
* Flows represent transitions between states.  A flow can go from a stock to a stock, from nothing to a stock, or from a stock to nothing, the latter two representing a population entering or leaving a model.  Flows are dependent on a dynamic variable, which determines its rate.
* Dynamic variables are the equations which determine flow rates.  For instance, if the rate at which people move from uninfected to infected is the uninfected stock S times the parameter c, we can say v\_infectionRate = c * S, and use this as a flow variable.  Dynamic variables can contain stocks, sum variables, parameters and other dynamic variables (but the definition can't be circular; you can't set a dynamic variable equal to itself, for instance).  Currently supports common binary and unary operations such as addition, multiplication, log, exp.
* Sum variables represent the sum of a subpopulation.  The value of a sum variable is the sum of all the stocks which link to it.  Common examples include the sum of the entire population N, the sum of infected individuals NI, or the sum of a particular subpopulation NS.
* Parameters are variables for which the particular values are defined outside the model definition.  You can provide parameters and initial values, then solve the differential equation.

Stock-flow diagrams can be created and manipulated using the `StockFlow` and `StockFlow.Syntax` modules, the latter of which includes the domain specific language for easier creation.

Stock-flow diagrams can be composed and stratified - that is, combined on stocks and sum variables or split into subpopulations.  One can also do algebraic rewriting for more manual, specific fixes, such as preventing a dynamic variable from unnecessarily being computed twice, or substituting one parameter for another.

Examples of the domain specific language, composition, stratification and algebraic rewriting can be seen in the [examples](examples) folder.  In particular, full\_fledged\_schema\_examples\_new contains examples which use the domain specific language.



 ## Example interpretation of a stock and flow diagram using an ODE solver
 
 From the [SEIR Composition Example](examples/full_fledged_schema_examples_new/composition/SEIR_full_model_measles_chickenpox.ipynb):
 
 ```julia
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
     ☁ => fbirth(μ * N) => S # dynamic variables can be defined implicitly or with :dynamic_variables
     S => fincid(β * S * I / N) => E
     S => fdeathS(S * δ) => ☁
     E => finf(E / tlatent) => I
     E => fdeathE(E * δ) => ☁
     I => frec(I / trecovery) => R
     I => fdeathI(I * δ) => ☁
     R => fdeathR(R * δ) => ☁
 
     :sums
     N = [S, E, I, R]
 
 end
 
 # define parameter values and initial values of stocks
 # define constant parameters
 p_measles = LVector(
     β=49.598, μ=0.03/12, δ=0.03/12, tlatent=8.0/30, trecovery=5.0/30
 )
 
 # define initial values for stocks
 u0_measles = LVector(
     S=90000.0-930.0, E=0.0, I=930.0, R=773545.0
 )
 
 prob_measles = ODEProblem(vectorfield(seir),u0_measles,(0.0,120.0),p_measles);
 sol_measles = solve(prob_measles,Tsit5(),abstol=1e-8);
 plot(sol_measles)
 ```
 
 See the full example for additional comments and the chickenpox model.
 ## Note
 
 This library is under active development and does not yet have an official release.
 
