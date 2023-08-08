module PremadeModels
using StockFlow.Syntax


export seir, sis, sir, svi

function seir()
    return deepcopy(seir_model)
end

function sis()
    return deepcopy(sis_model)
end

function sir()
    return deepcopy(sir_model)
end

function svi()
    return deepcopy(svi_model)
end



seir_model = @stock_and_flow begin

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
        

        :dynamic_variables
        v_prevalence = NI / NS
        v_meanInfectiousContactsPerS = c * v_prevalence # where c doesn't matter, can just make it 1
        v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
        v_newIncidence = S * v_perSIncidenceRate

        v_birth = μ * N

        v_inf = E / tlatent

        v_rec = I / trecovery 


        v_deathS = δ * S       
        v_deathE = δ * E
        v_deathI = δ * I
        v_deathR = δ * R
        

        :flows
        CLOUD => f_birth(v_birth) => S
        S => f_incid(v_newIncidence) => E
        S => f_deathS(v_deathS) => CLOUD
        E => f_inf(v_inf) => I
        E => f_deathE(v_deathE) => CLOUD
        I => f_rec(v_rec) => R
        I => f_deathI(v_deathI) => CLOUD
        R => f_deathR(v_deathR) => CLOUD

        :sums
        N = [S, E, I, R]
        NI = [I]
        NS = [S, E, I, R]

end

sis_model = @stock_and_flow begin
    :stocks
    S
    I

    
    :parameters
    μ
    β
    trec # 1 / trecovery.  This corresponds to σ.
    δ
    c

      
    :dynamic_variables
    v_deathsX = δ * S
    v_births = μ * N

    v_prevalence = NI / NS
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_newIncidence = S * v_perSIncidenceRate

    v_newRecovery = I * trec
    v_deathsI = I * δ

    :flows

    S => f_deathsX(v_deathsX) => CLOUD
    CLOUD => f_births(v_births) => S
    S => f_newInfectious(v_newIncidence) => I
    I => f_newRecovery(v_newRecovery) => S

    I => f_deathsI(v_deathsI) => CLOUD

    :sums
    N = [S, I]
    NI = [I]
    NS = [S, I]
  
end


sir_model = @stock_and_flow begin
    :stocks
    S
    I
    R
    
    :parameters
    c
    β
    rRec
    
    :dynamic_variables
    v_prevalence = NI / NS
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_newInfections = S * v_perSIncidenceRate
    v_newRecovery = I * rRec
    
    :flows
    S => f_inf(v_newInfections) => I
    I => f_rec(v_newRecovery) => R
    
    :sums
    N = [S, I, R]
    NI = [I]
    NS = [S,I,R]
    
    
end



svi_model = @stock_and_flow begin
    
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

    :dynamic_variables
    v_vacc = S * rvaccine 
    v_deathV = δ * V
    
    v_prevalence = NI / NS
    v_meanInfectiousContactsPerS = c * v_prevalence
    v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
    v_vaccineInfectionRate = V / evaccine # same thing as multiplying by complement
    v_perSIncidenceVaccinated = v_vaccineInfectionRate * v_perSIncidenceRate


    :flows
    S => f_vacc(v_vacc) => V
    V => f_deathV(v_deathV) => CLOUD
    V => f_infV(v_perSIncidenceVaccinated) => I 

    
    :sums
    N = [S, V, I]
    NI = [I]
    NS = [S, V, I]

end








end
