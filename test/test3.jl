using Test

@testset "Testing solution to Exercise 3" begin

@testset "Running ex3.jl" begin
   include("../ex3.jl")
end;

@testset "Testing that variables exist" begin
   @test @isdefined(response_1b)
   @test @isdefined(response_1d)
   @test @isdefined(response_1e)
end;

@testset "Testing that variables are not missing" begin
   @test !ismissing(response_1b)
   @test !ismissing(response_1d)
   @test !ismissing(response_1e)
end;

@testset "Add your tests here" begin
   @test 1 == 1
end;

end; # Exercise 3
