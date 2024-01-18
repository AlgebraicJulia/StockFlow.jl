using Test
using StockFlow

@testset "StockFlow DSL" begin
    include("Syntax.jl")
end

@testset "Attribute names prefixes and suffixes" begin
    include("SystemStructure.jl")
end

@testset "Causal Loop F" begin
    include("CausalLoop.jl")
end

@testset "Visualization" begin
    include("Visualization.jl")
end