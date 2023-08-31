using StockFlow.Syntax.Stratification

using StockFlow.Syntax.Stratification: interpret_stratification_notation
using StockFlow.Syntax: NothingFunction, DSLArgument, unwrap_expression, substitute_symbols

using Catlab.WiringDiagrams
using Catlab.ACSets
using Catlab.CategoricalAlgebra




@testset "Pullback computed in standard way is equal to DSL pullbacks" begin


    l_type = @stock_and_flow begin 
        :stocks
        pop
        
        :parameters
        μ
        δ
        rFstOrder
        rage
        
        :dynamic_variables
        v_aging = pop * rage
        v_fstOrder = pop * rFstOrder
        v_birth = N * μ
        v_death = pop * δ
        
        :flows
        pop => f_aging(v_aging) => pop
        pop => f_fstOrder(v_fstOrder) => pop
        CLOUD => f_birth(v_birth) => pop
        pop => f_death(v_death) => CLOUD
        
        :sums
        N = [pop]
        
    end;
    l_type_noatts = map(l_type, Name=NothingFunction, Op=NothingFunction, Position=NothingFunction);
    
    
    WeightModel = @stock_and_flow begin
        :stocks
        NormalWeight
        OverWeight
        Obese
        
        :parameters
        μ
        δw
        rw
        ro
        δo
        rage
        
        :dynamic_variables
        v_NewBorn = N * μ
        v_DeathNormalWeight = NormalWeight * δw
        v_BecomingOverWeight = NormalWeight * rw
        v_DeathOverWeight = OverWeight * δw
        v_BecomingObese = OverWeight * ro
        v_DeathObese = Obese * δo
        v_idNW = NormalWeight * rage
        v_idOW = OverWeight * rage
        v_idOb = Obese * rage
        
        :flows
        CLOUD => f_NewBorn(v_NewBorn) => NormalWeight
        NormalWeight => f_DeathNormalWeight(v_DeathNormalWeight) => ClOUD
        NormalWeight => f_BecomingOverWeight(v_BecomingOverWeight) => OverWeight
        OverWeight => f_DeathOverWeight(v_DeathOverWeight) => CLOUD
        
        OverWeight => f_BecomingObese(v_BecomingObese) => Obese
        Obese => f_DeathObese(v_DeathObese) => CLOUD
        NormalWeight => f_idNW(v_idNW) => NormalWeight
        OverWeight => f_idOW(v_idOW) => OverWeight
        Obese => f_idOb(v_idOb) => Obese
        
        :sums
        N = [NormalWeight, OverWeight, Obese]
        
    end;
    
    
    ageWeightModel = @stock_and_flow begin
        :stocks
        Child
        Adult
        Senior
        
        :parameters
        μ
        δC
        δA
        δS
        rageCA
        rageAS
        r
        
        :dynamic_variables
        v_NB = N * μ
        v_DeathC = Child * δC
        v_idC = Child * r
        v_agingCA = Child * rageCA
        v_DeathA = Adult * δA
        v_idA = Adult * r
        v_agingAS = Adult * rageAS
        v_DeathS = Senior * δS
        v_idS = Senior * r
        
        :flows
        CLOUD => f_NB(v_NB) => Child
        Child => f_idC(v_idC) => Child
        Child => f_DeathC(v_DeathC) => CLOUD
        Child => f_agingCA(v_agingCA) => Adult
        Adult => f_idA(v_idA) => Adult
        Adult => f_DeathA(v_DeathA) => CLOUD
        Adult => f_agingAS(v_agingAS) => Senior
        Senior => f_idS(v_idS) => Senior
        Senior => f_DeathS(v_DeathS) => CLOUD
        
        :sums
        N = [Child, Adult, Senior]
            
    end;
    
    begin
        s, = parts(l_type, :S)
        N, = parts(l_type, :SV)
        lsn, = parts(l_type, :LS)
        f_aging, f_fstorder, f_birth, f_death = parts(l_type, :F)
        i_aging, i_fstorder, i_birth = parts(l_type, :I)
        o_aging, o_fstorder, o_death = parts(l_type, :O)
        v_aging, v_fstorder, v_birth, v_death = parts(l_type, :V)
        lv_aging1, lv_fstorder1, lv_death1 = parts(l_type, :LV)
        lsv_birth1, = parts(l_type, :LSV)
        p_μ, p_δ, p_rfstOrder, p_rage = parts(l_type, :P)
        lpv_aging2, lpv_fstorder2, lpv_birth2, lpv_death2 = parts(l_type, :LPV)
    end;
    
    typed_WeightModel=ACSetTransformation(WeightModel, l_type_noatts,
      S = [s,s,s],
      SV = [N],
      LS = [lsn,lsn,lsn],   
      F = [f_birth, f_death, f_fstorder, f_death, f_fstorder, f_death, f_aging, f_aging, f_aging],    
      I = [i_birth, i_aging, i_fstorder, i_aging, i_fstorder, i_aging], 
      O = [o_death, o_fstorder, o_aging, o_death, o_fstorder, o_aging, o_death, o_aging],
      V = [v_birth, v_death, v_fstorder, v_death, v_fstorder, v_death, v_aging, v_aging, v_aging],
      LV = [lv_death1, lv_fstorder1, lv_death1, lv_fstorder1, lv_death1, lv_aging1, lv_aging1, lv_aging1],
      LSV = [lsv_birth1],
      P = [p_μ, p_δ, p_rfstOrder, p_rfstOrder, p_δ, p_rage],
      LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_fstorder2, lpv_death2, lpv_aging2, lpv_aging2, lpv_aging2],
      Name=NothingFunction, Op=NothingFunction, Position=NothingFunction
    );
    @assert is_natural(typed_WeightModel);
    
    
    
    typed_ageWeightModel=ACSetTransformation(ageWeightModel, l_type_noatts,
      S = [s,s,s],
      SV = [N],
      LS = [lsn,lsn,lsn],   
      F = [f_birth, f_fstorder, f_death, f_aging, f_fstorder, f_death, f_aging, f_fstorder, f_death],    
      I = [i_birth, i_fstorder, i_aging, i_fstorder, i_aging, i_fstorder], 
      O = [o_fstorder, o_death, o_aging, o_fstorder, o_death, o_aging, o_fstorder, o_death],
      V = [v_birth, v_death, v_fstorder, v_aging, v_death, v_fstorder, v_aging, v_death, v_fstorder],
      LV = [lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1, lv_aging1, lv_death1, lv_fstorder1],
      LSV = [lsv_birth1],
      P = [p_μ, p_δ, p_δ, p_δ, p_rage, p_rage, p_rfstOrder],
      LPV = [lpv_birth2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2, lpv_aging2, lpv_death2, lpv_fstorder2],
      Name =NothingFunction, Op=NothingFunction, Position=NothingFunction
    );
    @assert is_natural(typed_ageWeightModel);
    
    aged_weight = pullback(typed_WeightModel, typed_ageWeightModel) |> apex |> rebuildStratifiedModelByFlattenSymbols;
    
    # #########################################
    
    age_weight_2 = sfstratify(WeightModel, l_type, ageWeightModel, quote 
        :stocks
        NormalWeight, OverWeight, Obese => pop <= Child, Adult, Senior
    
        :flows
        f_NewBorn => f_birth <= f_NB
        f_DeathNormalWeight, f_DeathOverWeight, f_DeathObese => f_death <= f_DeathC, f_DeathA, f_DeathS
        f_idNW, f_idOW, f_idOb => f_aging <= f_agingCA, f_agingAS
        f_BecomingOverWeight, f_BecomingObese => f_fstOrder <= f_idC, f_idA, f_idS
    
        :dynamic_variables
        v_NewBorn => v_birth <= v_NB
        v_DeathNormalWeight, v_DeathOverWeight, v_DeathObese => v_death <= v_DeathC, v_DeathA, v_DeathS
        v_idNW, v_idOW, v_idOb  => v_aging <= v_agingCA, v_agingAS
        v_BecomingOverWeight, v_BecomingObese => v_fstOrder <= v_idC, v_idA, v_idS
    
        :parameters
        μ => μ <= μ
        δw, δo => δ <= δC, δA, δS
        rw, ro => rFstOrder <= r
        rage => rage <= rageCA, rageAS
    
        :sums
        N => N <= N
        
    end)
    #########################################
    
    age_weight_3 =  sfstratify(WeightModel, l_type, ageWeightModel, quote
    
        :flows
        f_NewBorn => f_birth <= f_NB
        ~Death => f_death <= ~Death
        ~id => f_aging <= ~aging
        ~Becoming => f_fstOrder <= ~id
    
        :dynamic_variables
        v_NewBorn => v_birth <= v_NB
        ~Death => v_death <= ~Death
        ~id  => v_aging <= ~aging
        ~Becoming => v_fstOrder <= ~id
    
        :parameters
        μ => μ <= μ
        ~δ => δ <= ~δ
        rage => rage <= rageCA, rageAS
        _ => rFstOrder <= _
    
    end) 
    
    age_weight_4 =  sfstratify(WeightModel, l_type, ageWeightModel, quote
    
        :flows
        ~NO_MATCHES => f_birth <= ~NO_MATCHES
        f_NewBorn => f_birth <= f_NB
        ~Death => f_death <= ~Death
        ~id => f_aging <= ~aging
        ~Becoming => f_fstOrder <= ~id
        ~Becoming => f_aging <= ~id # Everything already matched; ignored
        _ => f_aging <= _ # also ignored
    
        :dynamic_variables
        v_NewBorn => v_birth <= v_NB
        ~Death => v_death <= ~Death
        ~id  => v_aging <= ~aging
        _ => v_fstOrder <= _
    
        :parameters
        μ => μ <= μ
        ~δ => δ <= ~δ
        rage => rage <= rageCA, rageAS
        _ => rFstOrder <= _
    
    end) 




    @test aged_weight == age_weight_2
    @test aged_weight == age_weight_3
    @test aged_weight == age_weight_4
end

@testset "Ensuring interpret_stratification_notation correctly reads lines" begin # This should be all valid cases.  There's always going to be at least one value on both sides.
    # Note the orders.  The lists produced go left to right.  A1, A2 => B <= C1, C2 results in [A1 => B, A2 => B], [C1 => B. C2 => B]


    @test interpret_stratification_notation(:(A => B <= C)) == ([DSLArgument(:A, :B, Set{Symbol}())], [DSLArgument(:C, :B, Set{Symbol}())])
    
    @test interpret_stratification_notation(:(A1, A2 => B <= C)) == (
        [DSLArgument(:A1, :B, Set{Symbol}()), DSLArgument(:A2, :B, Set{Symbol}())],
        [DSLArgument(:C, :B, Set{Symbol}())]
    )
    @test interpret_stratification_notation(:(A => B <= C1, C2)) == (
        [DSLArgument(:A, :B, Set{Symbol}())],
        [DSLArgument(:C1, :B, Set{Symbol}()), DSLArgument(:C2, :B, Set{Symbol}())],
    )
    @test interpret_stratification_notation(:(_ => B <= _)) == (
        [DSLArgument(:_, :B, Set{Symbol}())],
        [DSLArgument(:_, :B, Set{Symbol}())],
    )
    @test interpret_stratification_notation(:(~A => B <= ~C)) == (
        [DSLArgument(:A, :B, Set{Symbol}([:~]))],
        [DSLArgument(:C, :B, Set{Symbol}([:~]))],
    )
    @test interpret_stratification_notation(:(~A1, A2 => B <= ~C)) == (
        [DSLArgument(:A1, :B, Set{Symbol}([:~])), DSLArgument(:A2, :B, Set{Symbol}())],
        [DSLArgument(:C, :B, Set{Symbol}([:~]))],
    )

    @test interpret_stratification_notation(:(~_ => B <= ~_, C)) == ( # Weird case.  Matches everything with _ as a substring.
        [DSLArgument(:_, :B, Set{Symbol}([:~]))],
        [DSLArgument(:_, :B, Set{Symbol}([:~])), DSLArgument(:C, :B, Set{Symbol}())]
    )

end


@testset "Unwrapping expressions works correctly" begin
    @test unwrap_expression(:S) == (:S, Set{Symbol}())
    @test unwrap_expression(:(~S)) == (:S, Set{Symbol}([:~]))
    @test unwrap_expression(:(~_)) == (:_, Set{Symbol}([:~]))
end


# function substitute_symbols(s::Dict{Symbol, Int}, t::Dict{Symbol, Int}, m::Vector{DSLArgument} ; use_flags::Bool=true)::Dict{Int, Int}

@testset "Testing substituting symbols" begin # underscore matching occurs at the very end, after this step.
    s1 = Dict(:A => 1)
    t1 = Dict(:B => 2)
    m1₁ = [DSLArgument(:A, :B, Set{Symbol}())]
    m1₂ = [DSLArgument(:A, :B, Set{Symbol}([:~]))]
 
    @test substitute_symbols(s1, t1, m1₁) == Dict(1 => 2) # A=>B -> 1=>2
    @test substitute_symbols(s1, t1, m1₂) == Dict(1 => 2) # A=>B -> 1=>2
    @test substitute_symbols(s1, t1, m1₂, use_flags=false) == Dict(1 => 2) # A=>B -> 1=>2

    s2 = Dict(:A1 => 10, :A2 => 20, :A3 => 30) # Unfortunately, cannot do substring matches starting with numbers, since it would require a symbol starting with a numbre.  Might need to add something for this...
    t2 = Dict(:B1 => 1, :B2 => 2)
    m2₁ = [DSLArgument(:A, :B1, Set{Symbol}([:~]))]

    @test substitute_symbols(s2, t2, m2₁) == Dict(10 => 1, 20 => 1, 30 => 1) #~A=>B -> 10=>1, 20=>1, 30=>1
    # @test substitute_symbols(s2, t2, m2₁, use_flags=false) == Dict()  # deliberately throws an error

    s3 = Dict{Symbol, Int}()
    t3 = Dict{Symbol, Int}()
    m3 = Vector{DSLArgument}()

    @test substitute_symbols(s3, t3, m3) == Dict()
    @test substitute_symbols(s3, t3, m3, use_flags=false) == Dict()

    s4 = Dict(:A1 => 1, :A2 => 2, :AB3 => 3, :AB4 => 4, :A5 => 5)
    t4 = Dict(:B1 => 1, :B2 => 2, :B3 => 3)
    m4 = [DSLArgument(:A1, :B1, Set{Symbol}()), DSLArgument(:B, :B2, Set{Symbol}([:~])), DSLArgument(:A, :B3, Set{Symbol}([:~]))]

    # always goes with first match.  A1 is taken, B matches AB3 and AB4, then A matches A2 and A5
    @test substitute_symbols(s4, t4, m4) == Dict(1 => 1, 3 => 2, 4 => 2, 2 => 3, 5 => 3)
end


