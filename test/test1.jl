using Test

@testset "Testing solution to Exercise 1" begin

@testset "Running ex1.jl" begin
   include("../ex1.jl")
end;

@testset "Testing that variables exist" begin 
   #@test @isdefined student_year 
end;

end; # Exercise 1

