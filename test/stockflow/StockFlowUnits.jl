using Test
using StockFlow
using StockFlow.Syntax


@testset "Extracting exponents" begin

  @test extract_exponents(:A) == Dict(:A => 1)
  @test extract_exponents(:(A * B)) == Dict(:A => 1, :B => 1)
  @test extract_exponents(:(A * B * C)) == Dict(:A => 1, :B => 1, :C => 1)

  @test extract_exponents(:(A/B * C)) == Dict(:A => 1, :B => -1, :C => 1)
  @test extract_exponents(:(A^2/B^2 * C^3)) == Dict(:A => 2, :B => -2, :C => 3)

  # (A^1.5 / B^0.5) * 1/C^0.1
  @test extract_exponents(:(A^1.5/B^0.5 / C^0.1)) == Dict(:A => 1.5, :B => -0.5, :C => -0.1)

  @test extract_exponents(:(1/B * C)) == Dict(:B => -1, :C => 1)
end


@testset "Basic Stock and flow unit constructions" begin
  SIR = StockAndFlowU(
    (
      (:S, :people) => (:F_NONE,:inf,:N),
      (:I, :people) =>(:inf,:rec,:N), 
      (:R, :people) => (:rec,:F_NONE,:N)
    ), # stocks

      (
        (:c => Symbol(:(1/time))), 
        (:beta => :NONE), 
        (:tRec => Symbol(:(1/time)))
      ),# parameters
(:v_prevalence=>((:I,:N)=>:/),:v_meanInfectiousContactsPerS=>((:c,:v_prevalence)=>:*),:v_perSIncidenceRate=>((:beta,:v_meanInfectiousContactsPerS)=>:*),
        :v_newInfections=>((:S,:v_perSIncidenceRate)=>:*),:v_newRecovery=>((:I,:tRec)=>:/)),# dynamical variables
(:inf=>:v_newInfections, :rec=>:v_newRecovery),# flows
(:N), # sum dynamical variables
  (:people => [:people => 1.0], Symbol(:(1/time)) => [:time => -1.0], :NONE => Pair{Symbol, Float64}[]), # derived units
  (:people, :time) # units
)

  #= 
    All stocks link to first dunit, params link to 2nd, 3rd, 2nd
    There are 2 links between dunits and units
    people has exponent 1, and 1/time has exponent -1 for time
  =#
  @test (unames(SIR) == [:people, :time]) && 
    (dunames(SIR) == [:people, Symbol(:(1/time)), :NONE]) && 
    get_sdus(SIR) == [1,1,1] && get_pdus(SIR) == [2,3,2] && 
    lunames(SIR) == [(:people, :people), (:time, Symbol(:(1/time)))] && get_exps(SIR) == [1, -1]

  @test SIR == (@stock_and_flow_U begin
    :stocks
    S: people
    I: people
    R: people

    :parameters
    c: 1/time
    beta
    tRec: 1/time

    :dynamic_variables
    v_prevalence = I / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = beta * v_meanInfectiousContactsPerS
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I / tRec

    :flows
    S => inf(v_newInfections) => I
    I => rec(v_newRecovery) => R
    :sums
    N = [S,I,R]
  end)

end