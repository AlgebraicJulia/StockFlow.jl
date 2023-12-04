using Pkg
Pkg.activate(".")

using Test

using StockFlow
using StockFlow.Syntax
using StockFlow.Syntax.Rewrite

using StockFlow.Syntax.Stratification
using Catlab.CategoricalAlgebra
using AlgebraicRewriting


@testset "Trivial rewrite examples act as expected." begin
  empty = StockAndFlowF()
  A = @stock_and_flow begin; :stocks; A; end
  p = @stock_and_flow begin; :parameters; p; end

  Ap = @stock_and_flow begin; :stocks; A; :parameters; p; end

  @test (@rewrite empty begin end) == empty
  @test (@rewrite A begin end) == A

  @test (@rewrite empty begin
    :stocks
    A    
    end) == A

  @test (@rewrite A begin
    :removes
    A
    end) == empty

  @test (@rewrite A begin
    :removes
    B
    end) == A

  @test (@rewrite p begin
    :removes
    p
  end) == empty

  @test (@rewrite Ap begin
    :removes
    p
  end) == A

  @test (@rewrite Ap begin
    :removes
    A
  end) == p

  @test (@rewrite p begin
    :stocks
    A
  end) == Ap

  @test (@rewrite A begin
    :parameters
    p
  end) == Ap

    # @test (@rewrite A begin
    #   :stocks
    #   B
    #   end) == A

    # delete A, but shouldn't
    # We might want to just add a :removes section and if you put it under
    # a :parameters or :stocks header it's just added.

    # @test (@rewrite A begin 
    #   :parameters
    #   A
    #   end) == A
    

end


@testset "Basic dynamic variable manipulation" begin
  # You should never swap both variables in a dynamic variable in one rewrite.


  ABv = @stock_and_flow begin
    :stocks
    A
    B
    :dynamic_variables
    v = A + B
  end

  AAv = @stock_and_flow begin
    :stocks
    A
    :dynamic_variables
    v = A + A
  end

  BBv = @stock_and_flow begin
    :stocks
    B
    :dynamic_variables
    v = B + B
  end

  Av = @stock_and_flow begin
    :stocks
    A
    :dynamic_variables
    v = +(A)
  end





  @test (@rewrite ABv begin
    :swaps
    B => A
    :removes
    B
  end) == AAv

  @test (@rewrite ABv begin
    :swaps
    B => A
    :removes
    B
  end) == AAv

  @test (@rewrite ABv begin
    :redefs
    v := A + A
    :removes
    B
  end) == AAv

  @test (@rewrite ABv begin
    :redefs
    v := +(A)
    :removes
    B
  end) == Av

  @test (@rewrite ABv begin
    :removes
    B
    :redefs
    v := +(A)
  end) == Av

  @test (@rewrite ABv begin
    :removes
    B

  end) == Av


  @test (@rewrite Av begin
    :dynamic_variables
    v = B + B
    :stocks
    B
    :removes
    v
    A
  end) == BBv


end



@testset "Modifying flows" begin

  ABvf = @stock_and_flow begin
    :stocks
    A
    B

    :dynamic_variables
    v = A + B

    :flows
    A => f(v) => B
  end

  ACvf = @stock_and_flow begin
    :stocks
    A
    C

    :dynamic_variables
    v = A + C

    :flows
    A => f(v) => CLOUD
  end

  @test (@rewrite ABvf begin
    :swaps
    B => C
    :stocks
    C
    :removes
    B
  end) == ACvf



end


@testset "Sums" begin
  sv1 = @stock_and_flow begin
    :stocks
    S
    I

    :sums
    N = [S, I]
  end

  sv2 = @stock_and_flow begin
    :stocks
    S

    :sums
    N = [S]
  end


  sv3 = @stock_and_flow begin
    :stocks
    S
    I
    R

    :dynamic_variables
    v1 = S * I
    v2 = R * I
    v3 = R + R

    :sums
    N = [S, I, R]
    NI = [I]
  end

  sv4 = @stock_and_flow begin
    :stocks
    S
    R

    :dynamic_variables
    v1 = *(S)
    v2 = *(R)
    v3 = R + R

    :sums
    N = [S, R]
    NI = []
  end


  @test (@rewrite sv1 begin
    :removes
    I
  end) == sv2

  @test (@rewrite sv3 begin
    :removes
    I
  end) == sv4

  sv4_rewrite1 = (@rewrite sv4 begin
    :stocks
    I

    :redefs
    v1 := S * I
    v2 := R * I
    N := [S, I, R]
    NI := [I]
  end)

  @test is_natural(homomorphism(sv4_rewrite1, sv3)) && is_natural(homomorphism(sv3, sv4_rewrite1))



end


@testset "Substantial examples" begin
  
  age2 = @stock_and_flow begin
    :stocks
    Child
    Adult
    
    :parameters
    c_C
    β
    r
    rAge
    c_A
    
    :dynamic_variables
    v_INC = Child / NC
    v_cINC = c_C * v_INC
    v_cβINC = β * v_cINC
    
    v_infC = Child * v_cβINC
    v_fstC = Child * r
    v_agingC = Child * rAge
    
    
    v_INA = Adult / NA
    v_cINA = c_A * v_INA
    v_cβINA = β * v_cINA
    
    v_infA = Adult * v_cβINA
    v_fstA = Adult * r
    
    :flows
    Child => f_infC(v_infC) => Child
    Child => f_frsC(v_fstC) => Child
    Child => f_aging(v_agingC) => Adult
    Adult => f_infA(v_infA) => Adult
    Adult => f_frsA(v_fstA) => Adult
    
    
    :sums
    NC = [Child]
    NA = [Adult]
  end

  sir = @stock_and_flow begin
    :stocks
    S
    I
    R
    
    :parameters
    c
    β
    rRec
    rAge
    
    :dynamic_variables
    v_prevalence = I / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I * rRec
    v_idS = S * rAge
    v_idI = I * rAge
    v_idR = R * rAge
    
    :flows
    S => f_idS(v_idS) => S
    S => f_inf(v_newInfections) => I
    I => f_idI(v_idI) => I
    I => f_rec(v_newRecovery) => R
    R => f_idR(v_idR) => R
    
    :sums
    N = [S, I, R]
    
    
  end

  s_type = @stock_and_flow begin
    :stocks
    pop
    
    :parameters
    c
    β
    rFstOrder
    rAge
    
    :dynamic_variables
    v_prevalence = pop / N
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_inf = pop * v_perSIncidenceRate
    v_fstOrder = pop * rFstOrder
    v_aging = pop * rAge
    
    :flows
    pop => f_inf(v_inf) => pop
    pop => f_fstOrder(v_fstOrder) => pop
    pop => f_aging(v_aging) => pop

    
    :sums
    N = [pop]
    
    
  end





  aged_sir = @stratify sir s_type age2 begin
    
    :parameters
    c => c <= c_C, c_A
    β => β <= β
    rRec => rFstOrder <= r
    rAge => rAge <= rAge

    :dynamic_variables
    v_prevalence => v_prevalence <= ~v_IN
    v_meanInfectiousContactsPerS => v_meanInfectiousContactsPerS <= ~v_cIN 
    v_perSIncidenceRate => v_perSIncidenceRate <= ~v_cβIN
    v_newInfections => v_inf <= ~v_inf
    v_newRecovery => v_fstOrder <= ~v_fst
    ~id => v_aging <= ~v_aging

    :flows
    f_inf => f_inf <= ~f_inf
    f_rec => f_fstOrder <= ~f_frs
    ~id => f_aging <= f_aging

  end


  aged_sir_rewritten = @rewrite aged_sir begin

    :redefs
    v_meanInfectiousContactsPerSv_cINC := cc_C * v_prevalencev_INC_post
    v_meanInfectiousContactsPerSv_cINA := cc_A * v_prevalencev_INA_post


    :parameters
    fcc
    fca
    fac
    faa

    :dynamic_variables

    v_CCContacts = fcc * v_prevalencev_INC
    v_CAContacts = fca * v_prevalencev_INA
    
    v_ACContacts = fac * v_prevalencev_INC
    v_AAContacts = faa * v_prevalencev_INA
    
    v_prevalencev_INC_post = v_CCContacts + v_CAContacts
    v_prevalencev_INA_post = v_ACContacts + v_AAContacts

  end


  LS = @stock_and_flow begin
    :stocks
    SChild
    IChild
    SAdult
    IAdult
    
    :parameters
    cc_C
    cc_A
    
    :dynamic_variables
    v_prevalencev_INC = IChild / NNC
    v_prevalencev_INA = IAdult / NNA
    v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC
    v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA
    
    
    :sums
    NNC = [SChild, IChild]
    NNA = [SAdult, IAdult]
    
    
  end



  IS = @stock_and_flow begin
      :stocks
      SChild
      IChild
      SAdult
      IAdult
      
      :parameters
      cc_C
      cc_A
      
      :dynamic_variables
      v_prevalencev_INC = IChild / NNC
      v_prevalencev_INA = IAdult / NNA
      v_meanInfectiousContactsPerSv_cINC = *(cc_C)
      v_meanInfectiousContactsPerSv_cINA = *(cc_A)
      
      :sums
      NNC = [SChild, IChild]
      NNA = [SAdult, IAdult]
      
      
  end


  RS = @stock_and_flow begin
      :stocks
      SChild
      IChild
      SAdult
      IAdult
      
      :parameters
      fcc
      fca
      fac
      faa
      cc_C
      cc_A
      
      :dynamic_variables
      v_prevalencev_INC = IChild / NNC
      v_prevalencev_INA = IAdult / NNA
      v_CCContacts = fcc * v_prevalencev_INC
      v_CAContacts = fca * v_prevalencev_INA
      
      v_ACContacts = fac * v_prevalencev_INC
      v_AAContacts = faa * v_prevalencev_INA
      
      v_prevalencev_INC_post = v_CCContacts + v_CAContacts
      v_prevalencev_INA_post = v_ACContacts + v_AAContacts
      v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC_post
      v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA_post
      
      :sums
      NNC = [SChild, IChild]
      NNA = [SAdult, IAdult]
      
      
  end


  aged_sir_rewritten2 = rewrite(Rule(homomorphism(IS, LS), homomorphism(IS, RS)), aged_sir)

  @test is_natural(homomorphism(aged_sir_rewritten, aged_sir_rewritten2)) && is_natural(homomorphism(aged_sir_rewritten2, aged_sir_rewritten))


end

