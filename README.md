# AlgebraicStockFlow
[![Tests](https://github.com/AlgebraicJulia/StockFlow.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/AlgebraicJulia/StockFlow.jl/actions/workflows/tests.yml)
<!-- TODO: Set up on codecov.io for repo [![codecov](https://codecov.io/github/AlgebraicJulia/StockFlow.jl/branch/maaster/graph/badge.svg)](https://app.codecov.io/github/AlgebraicJulia/StockFlow.jl) -->

StockFlow.jl is a Julia library for creating stock and flow diagrams using category theory.  Stockflows are used to represent systems where its population can move between different states, such as tracking how many people have been infected or recovered from a disease.  By providing initial values, and values for the parameters, we can solve the differential equation which the stockflow represents, which gives us a graph for how the population varies over time.

In this particular schema, stockflows consist of stocks, flows, dynamic variables, sum variables, parameters, and the links between them.

* Stocks represent accumulations of a population in a particular state, such as 'susceptible', 'infected', 'recovered', etc.
* Flows represent transitions between states.  Flow can go from a stock to a stock, from nothing to a stock, or from a stock to nothing, the latter two representing a population entering or leaving a model.  Flows are dependent on a dynamic variable, which determines its rate.
* Dynamic variables are the equations which determine flow rates.  For instance, if the rate at which people move from uninfected to infected is the uninfected stock S times the parameter c, we can say v\_infectionRate = c * S, and use this as a flow variable.  Dynamic variables can contain stocks, sum variables, parameters and other dynamic variables (but the definition can't be circular; you can't set a dynamic variable equal to itself, for instance).  Currently supports common binary and unary operations such as addition, multiplication, log, exp.
* Sum variables represent the sum of a subpopulation.  The value of a sum variable is the sum of all the stocks which link to it.  Common examples include the sum of the entire population N, the sum of infected variables NI, or the sum of a particular subpopulation NS.
* Parameters are variables for whcih the particular values are defined outside the model definition.  You can provide parameters and initial values, then solve the differential equation.

Stockflows can be created and manipulated using StockFlow and Syntax, the latter of which includes the domain specific language for easier creation.

Stockflows can be composed and stratified - that is, combined on stocks and sum variables or split into subpopulations.  One can also do algebraic rewriting for more manual, specific fixes, such as preventing a dynamic variable from unnecessarily being computed twice, or substituting one parameter for another.

Examples of the domain specific language, composition, stratification and algebraic rewriting can be seen in the examples folder.  In particular, full\_fledged\_schema\_examples\_new contains examples which use the domain specific language.

