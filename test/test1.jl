using Test

@testset "Testing solution to Exercise 1" begin

@testset "Running ex1.jl" begin
   include("../ex1.jl")
end;

@testset "Testing that variables exist" begin
   @test @isdefined(response_1a)
   @test @isdefined(response_1b)
   @test @isdefined(response_1c)
   @test @isdefined(response_1d)
   @test @isdefined(response_1e)
   @test @isdefined(response_1f)
   @test @isdefined(response_1g)
   @test @isdefined(response_1h)
   @test @isdefined(response_1i)
   @test @isdefined(response_1j)
   @test @isdefined(response_1k)
   @test @isdefined(response_1l)
   #@test @isdefined(response_1m)
   @test @isdefined(response_1n)
   @test @isdefined(response_1o)
   @test @isdefined(response_1p)
end;

@testset "Testing that variables are not missing" begin
   @test !ismissing(response_1a)
   @test !ismissing(response_1b)
   @test !ismissing(response_1c)
   @test !ismissing(response_1d)
   @test !ismissing(response_1e)
   @test !ismissing(response_1f)
   @test !ismissing(response_1g)
   @test !ismissing(response_1h)
   @test !ismissing(response_1i)
   @test !ismissing(response_1j)
   @test !ismissing(response_1k)
   @test !ismissing(response_1l)
   #@test !ismissing(response_1m)
   @test !ismissing(response_1n)
   @test !ismissing(response_1o)
   @test !ismissing(response_1p)
end;

@testset "Testing that functions' structure" begin
   @test my_function_1_arg_is_good_to_go
   @test my_function_2_args_is_good_to_go
end

@testset "Testing numerical answers" begin
   maxN_4gb = floor(Int,sqrt(4*2^30/8))
	@test response_1c == maxN_4gb
	maxN_1tb = floor(Int,sqrt(1024*2^30/8))
   @test response_1f == maxN_1tb
	my_flops = 1.5e9
	my_est_1d = (2//3)*maxN_4gb^3 / my_flops
	@test my_est_1d/10 <= response_1d <= 10*my_est_1d
	my_est_1g = (2//3)*maxN_1tb^3 / my_flops
	@test my_est_1g/10 <= response_1g <= 10*my_est_1g
	my_est_1i = (2//3)*100^3 / my_flops
	@test my_est_1i/10 <= response_1i <= 10*my_est_1i
end;


end; # Exercise 1
