module PremadeModels
using StockFlow.Syntax


export seir, sis, sir, svi

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
        c
        rAge 
        

        :dynamic_variables
        v_prevalence = NI / NS
        v_meanInfectiousContactsPerS = c * v_prevalence # where c doesn't matter, can just make it 1
        v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
        v_newIncidence = S * v_perSIncidenceRate

        v_birth = μ * N

        v_inf = E / tlatent # may be better to multiply as opposed to divide

        v_rec = I / trecovery # may be better to multiply than divide


        v_deathS = δ * S       
        v_deathE = δ * E
        v_deathI = δ * I
        v_deathR = δ * R
        
        v_idS = rAge * S
        v_idE = rAge * E
        v_idI = rAge * I
        v_idR = rAge * R

        :flows
        CLOUD => f_birth(v_birth) => S
        S => f_incid(v_newIncidence) => E
        S => f_deathS(v_deathS) => CLOUD
        E => f_inf(v_inf) => I
        E => f_deathE(v_deathE) => CLOUD
        I => f_rec(v_rec) => R
        I => f_deathI(v_deathI) => CLOUD
        R => f_deathR(v_deathR) => CLOUD

        S => f_idS(v_idS) => S
        E => f_idE(v_idE) => E
        I => f_idI(v_idI) => I       
        R => f_idR(v_idR) => R

        :sums
        N = [S, E, I, R]
        NI = [I]
        NS = [S, E, I, R]

end

sis = @stock_and_flow begin
    :stocks
    S
    I

    
    :parameters
    μ
    β
    trec # 1 / trecovery.  This corresponds to σ.
    δ
    rAgeS
    rAgeI
    c


    
      
    :dynamic_variables
    v_deathsX = δ * S
    v_births = μ * N

    v_prevalence = NI / NS
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_newIncidence = S * v_perSIncidenceRate

    # v_newInfectious
    v_newRecovery = I * trec
    v_deathsI = I * δ
    v_agingS = S * rAgeS
    v_agingI = I * rAgeI 

    :flows

    S => f_deathsX(v_deathsX) => CLOUD
    CLOUD => f_births(v_births) => S
    S => f_newInfectious(v_newIncidence) => I
    I => f_newRecovery(v_newRecovery) => S

    I => f_deathsI(v_deathsI) => CLOUD
    S => f_agingX(v_agingS) => S
    I => f_agingI(v_agingI) => I


    :sums
    N = [S, I]
    NI = [I]
    NS = [S, I]
  
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
    v_prevalence = NI / NS
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
    NI = [I]
    NS = [S,I,R]
    
    
end



svi = @stock_and_flow begin
    
    :stocks
    S
    V
    I

    :parameters
    rvaccine
    δ
    evaccine
    c
    β
    rAge

    :dynamic_variables

    v_vacc = S * rvaccine 
    v_deathV = δ * V
    
    v_prevalence = NI / NS
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_vaccineInfectionRate = V / evaccine # same thing as multiplying by complement
    v_perSIncidenceVaccinated = v_vaccineInfectionRate * v_perSIncidenceRate

    v_idS = S * rAge
    v_idV = V * rAge
    v_idI = I * rAge
    

    :flows
    S => f_vacc(v_vacc) => V
    V => f_deathV(v_deathV) => CLOUD
    V => f_infV(v_perSIncidenceVaccinated) => I 

    S => f_idS(v_idS) => S
    V => f_idV(v_idV) => V
    I => f_idI(v_idI) => I
    
    :sums
    N = [S, V, I]
    NI = [I]
    NS = [S, V, I]

end








end
