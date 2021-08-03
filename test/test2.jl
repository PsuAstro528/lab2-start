using Test

@testset "Testing solution to Exercise 2" begin

@testset "Running ex2.jl" begin
   include("../ex2.jl")
end;

@testset "Testing that variables exist" begin
   #@test @isdefined(response_1p)
end;

end; # Exercise 2
