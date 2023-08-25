using Test

@testset "Testing solution to Exercise 2" begin

@testset "Running ex2.jl" begin
   include("../ex2.jl")
end;

@testset "Testing that variables exist" begin
   @test @isdefined(response_2b)
   @test @isdefined(response_2d)
   @test @isdefined(response_2e)
end;

@testset "Testing that variables are not missing" begin
   @test !ismissing(response_2b)
   @test !ismissing(response_2d)
   @test !ismissing(response_2e)
end;

@testset "Add your tests here" begin
   @test 1 == 1
end;

end; # Exercise 2
