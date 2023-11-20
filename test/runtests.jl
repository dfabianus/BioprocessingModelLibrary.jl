using BioprocessingModelLibrary
using Test

@testset "BioprocessingModelLibrary.jl" begin
    # Write your tests here.
end

@testset "Refolding" begin
    include("test_refolding_models.jl")
end