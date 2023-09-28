using Test
using StockFlow

@testset "StockFlow DSL" begin
  include("Syntax.jl")
end

@testset "Attribute names prefixes and suffixes" begin
  include("SystemStructure.jl")
end