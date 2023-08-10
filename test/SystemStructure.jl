using Test
using StockFlow
using StockFlow.Syntax

empty = @stock_and_flow begin end

p = @stock_and_flow begin 
    :stocks
    S
    I
    
    :dynamic_variables
    v = S + I
    
    :parameters
    p1
    p2
    
    :flows
    S => f1(v) => CLOUD
    
    :sums
    N = [S]
    NI = [I]
end

p_prefixed = @stock_and_flow begin 
    :stocks
    prefS
    prefI
    
    :dynamic_variables
    prefv = prefS + prefI
    
    :parameters
    prefp1
    prefp2
    
    :flows
    prefS => preff1(prefv) => CLOUD
    
    :sums
    prefN = [prefS]
    prefNI = [prefI]
end

p_suffixed = @stock_and_flow begin 
    :stocks
    Ssuf
    Isuf
    
    :dynamic_variables
    vsuf = Ssuf + Isuf
    
    :parameters
    p1suf
    p2suf
    
    :flows
    Ssuf => f1suf(vsuf) => CLOUD
    
    :sums
    Nsuf = [Ssuf]
    NIsuf = [Isuf]
end

empty_foot = @foot () => ()

footA = @foot S => N
footA_pref = @foot prefS => prefN
footA_suf = @foot Ssuf => Nsuf


@testset "changing names act the same as if the stock flow/foot/open was created with the changed name" begin
    @test empty == add_suffix!(deepcopy(empty), "AAA") == add_prefix!(deepcopy(empty), "BBB")
    @test add_prefix!(deepcopy(p), "pref") == p_prefixed
    @test add_suffix!(deepcopy(p), "suf") == p_suffixed
    
    @test empty_foot == add_suffix!(deepcopy(empty_foot), "CCC") == add_prefix!(deepcopy(empty_foot), "DDD")
    @test add_prefix!(deepcopy(footA), "pref") == footA_pref
    @test add_suffix!(deepcopy(footA), "suf") == footA_suf

end


