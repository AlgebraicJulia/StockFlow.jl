# using StockFlow
# using StockFlow.Syntax
# using StockFlow.PremadeModels

# using Catlab.WiringDiagrams

# openid = Open((@stock_and_flow begin
#     :stocks
#     pop

#     :parameters
#     rAge

#     :dynamic_variables
#     v_id = S * rAge

#     :flows
#     S => f_id(v_id) => S
# end), @foot pop => ())

# openseir = Open(PremadeModels.seir(), @feet S => (); E => (); I => (); R => (););

# uwd = @relation (S,E,I,R) begin
#     SEIR(S,E,I,R)
#     S(S)
#     E(E)
#     I(I)
#     R(R)
# end

# GraphF(apex(oapply(uwd, [openseir, openid, openid, openid, openid])))


# s_type = @stock_and_flow begin
#     :stocks
#     pop
    
#     :parameters
#     c
#     β
#     rFstOrder
#     rAge
#     μ
#     δ
    
#     :dynamic_variables

#     v_birth = μ * N
#     v_death = δ * pop

#     v_prevalence = pop / N
#     v_meanInfectiousContactsPerS = c * v_prevalence
#     v_perSIncidenceRate = β * v_meanInfectiousContactsPerS
#     v_inf = pop * v_perSIncidenceRate
#     v_fstOrder = pop * rFstOrder
#     v_aging = pop * rAge
    
#     :flows
#     CLOUD => f_birth(v_birth) => pop
#     pop => f_death(v_death) => CLOUD
#     pop => f_inf(v_inf) => pop
#     pop => f_fstOrder(v_fstOrder) => pop
#     pop => f_aging(v_aging) => pop

    
#     :sums
#     N = [pop]
    
    
# end


# age = @stock_and_flow begin
#     :stocks
#     Child
#     Adult
    
#     :parameters
#     c_C
#     β
#     r
#     rAge
#     c_A
#     μ
#     δc
#     δa
    
#     :dynamic_variables

#     v_birth = μ * N
#     v_cdeath = δc * Child
#     v_adeath = δa * Adult

#     v_INC = Child / NC
#     v_cINC = c_C * v_INC
#     v_cβINC = β * v_cINC
    
#     v_infC = Child * v_cβINC
#     v_fstC = Child * r
#     v_agingC = Child * rAge
    
    
#     v_INA = Adult / NA
#     v_cINA = c_A * v_INA
#     v_cβINA = β * v_cINA
    
#     v_infA = Adult * v_cβINA
#     v_fstA = Adult * r
    
#     :flows
#     CLOUD => f_birth(v_birth) => Child
#     Child => f_cdeath(v_cdeath) => CLOUD
#     Adult => f_adeath(v_adeath) => CLOUD

#     Child => f_infC(v_infC) => Child
#     Child => f_frsC(v_fstC) => Child
#     Child => f_aging(v_agingC) => Adult
#     Adult => f_infA(v_infA) => Adult
#     Adult => f_frsA(v_fstA) => Adult
    
    
#     :sums
#     NC = [Child]
#     NA = [Adult]
# end






# # @stratify (seir, s_type, sex) begin
# #     :stocks
# #     _ => pop <= _

# #     :flows
# #     f_birth => f_birth <= f_birth
# #     f_deathS, f_deathE, f_deathI, f_deathR => f_death <= f_adeath, f_cdeath
# #     f_incid => f_inf <= f_infC, f_infA
# #     f_inf, f_rec => f_fstOrder <= frsC, frsA
# #     f_



# # end

