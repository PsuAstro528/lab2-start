### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 0854551f-fc6d-4c51-bf85-a9c7030e588b
using BenchmarkTools

# ╔═╡ cb099d3d-5c09-4c56-9c59-83adab60f651
using LinearAlgebra

# ╔═╡ a1122848-347f-408b-99c2-a7a514073864
using Plots

# ╔═╡ 4df26174-9aa3-46d0-879b-aec9a771714b
using LaTeXStrings

# ╔═╡ ea0a7ca2-503f-4d61-a3d4-42503f322782
begin
	using PlutoUI, PlutoTeachingTools
	eval(Meta.parse(code_for_check_type_funcs))
end

# ╔═╡ 4ac5cc73-d8d6-43d5-81b2-944d559fd2ca
md"""
# Astro 528 Lab 2, Exercise 1

## Benchmarking Code


Julia provides several tools for measuring code performance. Perhaps the simplest way is using the [`@time`](https://docs.julialang.org/en/v1.0/base/base/#Base.@time) or [`@elapsed`](https://docs.julialang.org/en/v1.0/base/base/#Base.@elapsed) macros, such as
"""

# ╔═╡ 637f6a84-ad01-43a7-899b-b7867fbe3b3d
@elapsed randn(1000)

# ╔═╡ c23045a0-56b2-4e1e-b5d3-8e248fd1bffd
md"""
The `@time` macro prints the time, but returns the value of the following expression. (Pluto doesn't normally show the output printed to the terminal inside the notebook.  You can either find it in the window where you're running the Pluto server or you can use the `with_terminal()` function provided by PlutoUI.jl to view the output.)  The `@elapsed` macro discards the following expressions return value and returns the elapsed time evaluating the expression.
"""

# ╔═╡ 6dcedd62-c606-43ba-b2f6-e018aefcb035
with_terminal() do 
	@time rand(1000)
end

# ╔═╡ e3bfd6bc-261d-4e69-9ce6-e78d959c13be
md"""There are even more sophisticated macros in the [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) package which provides `@btime` and `@belapsed` that provide outputs similar to `@time` and `@elapsed`, but take longer than @time, since it runs the code multiple times in an attempt to give more accurate results.  It also provides a `@benchmkark` macro that is quite flexible and provides more detailed information.  
"""

# ╔═╡ 2a5eb6c7-9b8c-4923-8fe4-e7211c6c820f
md"Let's define a function, `my_function_0_args` that takes a single arguement, the number of samples, and calls a mathematical function with zero arguments (e.g., `rand`) to be benchmarked for a given problem size.  "

# ╔═╡ b3e64508-319f-4506-8512-211e30be4bee
my_function_0_args(N::Integer) = rand(N)

# ╔═╡ a0a24490-3136-45e4-8f55-44448d8366d1
md"Next, we'll specificy what problem sizes you'd like to benchmark as a vector of integers named `num_list`.  (At some point you may want to add larger problem sizes, but beware that that will result in a delay while the rest of the notebook updates.)"

# ╔═╡ 69c39140-ba66-4345-869c-a5499d5d4376
num_list = [1,2,4,8,16,32,64,128,256,512]

# ╔═╡ 2e0b53c8-4107-485f-9f33-d921b6fc0c05
md"I've provided a function `benchmark_my_funciton` at the bottom of the notebook that will help us compactly compute benchmarks for `my_function_0_args`."

# ╔═╡ 928fe18c-2b34-427f-afdb-6273a0ce135e
md"Now, we'll generate several random datasets to use for testing function with different size arrays here, so that we don't have to keep regenerating datasets over and over."

# ╔═╡ fd7c88ef-02fb-408f-a551-45a93d873ded
begin
	x_list = rand.(num_list)
	y_list = rand.(num_list)
end;

# ╔═╡ ba4b35c8-0fe8-4a25-9817-d13c0cd98c6f
md"We'll start by benchmarking the square root function when applied to an array."

# ╔═╡ bb4e3e7e-d8da-4c53-b861-0e66237aae3c
sqrt_broadcasted(x) = sqrt.(x);

# ╔═╡ 6e1e699f-4f4a-4a21-946f-e763eb106f37
md"Now try benchmarking one or two univariate functions of your own and comparing them to `sqrt`.  Write a function `my_function_1_arg` that takes an array and applies a function of your choice to the input array."

# ╔═╡ 39ed7918-8926-4866-8b25-c61dbbd35991
my_function_1_arg(x::Array) = missing

# ╔═╡ a403c0f3-fc4a-4bc1-87c4-8a2d5d15e1a0
begin
	my_function_1_arg_is_good_to_go = false
	if !@isdefined(my_function_1_arg)
		func_not_defined(:my_function_1_arg)
	elseif length(methods(my_function_1_arg,[Array,])) <1
		PlutoTeachingTools.warning_box(md"`my_function_1_arg` should take an  `Array` as its arguement")
	elseif ismissing(my_function_1_arg([1,2,3]))
		still_missing()
	elseif size(my_function_1_arg([1,2,3]))!=(3,)
		almost(md"The size of the `my_function_1_arg`'s output doesn't match the size of the input.")
	else
		my_function_1_arg_is_good_to_go = true
		correct()
	end
end

# ╔═╡ 5d234d6c-be3d-4d0f-bd58-774b6786db54
md"Now, let's try benchmarking functions that take two variables."

# ╔═╡ 11df2b04-2bef-497e-a064-cbcd159aecc7
add_broadcasted(x,y) = x.+y;

# ╔═╡ 8320739f-ec14-4a49-b374-a0baa198f646
mul_broadcasted(x,y) = x.*y;

# ╔═╡ e352d85a-b8ce-4355-8b0b-9c28465dd006
div_broadcasted(x,y) = x./y;

# ╔═╡ 0ce6409c-3aa4-4446-ae13-81c4e743d322
md"Create a function `my_function_2_args` that takes two arrays and applies computes a function of your choice to them."

# ╔═╡ f13619d1-6ece-42be-965d-ab6bb9e9b8cd
my_function_2_args(x::Array, y::Array) = missing

# ╔═╡ 340808ea-e99b-4de7-ad55-b2d812ff0f4d
begin
	my_function_2_args_is_good_to_go = false
	if !@isdefined(my_function_2_args)
		func_not_defined(:my_function_2_args)
	elseif length(methods(my_function_2_args,[Array,Array])) <1
		PlutoTeachingTools.warning_box(md"`my_function_2_args` should take two `Array`'s as arguemetns")
	elseif ismissing(my_function_2_args([1,2,3],[4,5,6]))
		still_missing()
	elseif size(my_function_2_args([1,2,3],[4,5,6]))!=(3,)
		almost(md"The size of the `my_function_2_args`'s output doesn't match the size of the input.")
	else
		my_function_2_args_is_good_to_go = true
		correct()
	end
end

# ╔═╡ 3621793a-427b-40d2-b28d-7bb9c6f3a28c
md"""
Now, it's your turn.  Try updating `my_function_1_arg` and `my_function_2_args` to compute a few different mathematical functions.  For example, try a couple of trig functions, and a logarithm. 
"""

# ╔═╡ 10499798-1ba6-4a2f-827b-aecc4a1f8346
md"""
a.  How much longer did it take to compute a trig function than simple arithmetic?
How much longer did it take to compute a logarithm than simple arithmetic?
"""

# ╔═╡ ffde1960-8792-4092-9a0b-42b2726bb2da
response_1a = missing

# ╔═╡ 10f54651-4b05-4806-9a0c-cb0b6afa245b
display_msg_if_fail(check_type_isa(:response_1a,response_1a,Markdown.MD)) 

# ╔═╡ 2291567f-71a4-4a53-a80c-aab58f29ddf8
md"""
b.  Did the number of evaluations per second vary significantly depending on the number of elements in the array?  
How large of an array was necessary before the performance reached its asymptote?
"""

# ╔═╡ 05e772d4-a2ad-4f16-866c-b8aa3ab7a550
response_1b = missing

# ╔═╡ eb1ff48b-7dbc-4b37-a18d-77eadc568f5a
display_msg_if_fail(check_type_isa(:response_1b,response_1b,Markdown.MD)) 

# ╔═╡ c888afb4-59da-4d7b-a050-df84f3735202
md"""
## Memory Requirements for Linear Algebra Problems^[1]

Consider a modern laptop with 4 GB ($=4*2^{30}$) of usable memory. Assume it uses 8 bytes of memory to store each floating point number (i.e., Float64, double precision, real\*8). 
"""

# ╔═╡ b9157cf9-99da-46ee-8ba0-ed6a95153f74
md"""
### Theory
c. What is the number of rows in the largest square matrix that the above computer could fit into its available memory at one time? 
(Feel free to add some cells with code to compute the answer.)
"""

# ╔═╡ b0873452-b9b5-4bc5-bb76-5f299a6d366e
response_1c = missing

# ╔═╡ cc584d2a-2a3b-4898-b753-1a7275712ad2
display_msg_if_fail(check_type_isa(:response_1c,response_1c,Integer)) 

# ╔═╡ e69c7809-ec6c-4f3a-b869-100381d40bf9
if !@isdefined(response_1c)  || ismissing(response_1c)
	nothing
elseif response_1c == floor(Int,sqrt(4*2^30/8))
	correct()
else
	almost(md"Try again")
end

# ╔═╡ ea68d6eb-a542-44bf-b0e7-815e2372bc33
md"""
d. Estimate how long (in seconds) would it take to solve the maximum size linear system that would fit into memory at once, if we use LU factorization to solve a linear system.  You may assume the computation is maximally efficient, the computer reaches peak performance and the LU decomposition requires $(2/3)*N^3$ flops, where $N$ refers to the number of rows in the square array being factorized.
Use an approximation for the number of floating point operations per second based on your results above."""

# ╔═╡ 484b47d6-65de-497a-b78f-6996c8787de8
response_1d = missing

# ╔═╡ 8e2045c2-c5de-4895-9719-efb9de113dec
display_msg_if_fail(check_type_isa(:response_1d,response_1d,Real)) 

# ╔═╡ 2da17c61-a8d7-4dc2-8209-69fe809e9d4f
begin
	local N = floor(Int,sqrt(4*2^30/8))
	local num_ops = (2//3)*N^3
	local my_flops = 1.5e9
	local my_est = num_ops/my_flops
	if !@isdefined(response_1d)  || ismissing(response_1d)
		nothing
	elseif response_1d < my_est/10
		almost(md"Are you sure?  That seems low to me.")
		correct()
	elseif response_1d > my_est*10
		almost(md"Are you sure?  That seems high to me.")
	else
		correct(md"That seems like a plaussible runtime to me.")	
	end
end

# ╔═╡ 4e8e86ed-92b0-4aa5-945b-ef6193cb2559
md"""
e. Does memory or compute time limit the size of system that can be practically solved with LU decomposition for this modern laptop?"""

# ╔═╡ 1d4275fd-3208-477d-b26e-be23940035ba
response_1e = missing

# ╔═╡ c64bd8d8-6801-46a7-876b-c75b83300416
display_msg_if_fail(check_type_isa(:response_1e,response_1e,Markdown.MD)) 

# ╔═╡ 5cc54a66-461c-4bb3-8e29-2079adeff04f
md"""f. Now consider a high-end server with 1TB of RAM (such as ICS-ACI's high-memory nodes).  
How many rows are in the largest square matrix that would fit into its memory at once?"""

# ╔═╡ 30efc3ea-a2a2-48f8-833f-6b51845e47dc
response_1f = missing

# ╔═╡ 33bc47b6-f48b-4e24-bb54-fa794f04e166
display_msg_if_fail(check_type_isa(:response_1f,response_1f,Integer)) 

# ╔═╡ 3fab2e27-766c-4464-9b6c-e5235a80183a
if !@isdefined(response_1f)  || ismissing(response_1f)
	nothing
elseif response_1f == floor(Int,sqrt(1024*2^30/8))
	correct()
else
	almost(md"Try again")
end

# ╔═╡ 9af89fd2-a5f6-4807-beb6-d3fa9fa89b69
md"""1g. How long do you estimate it would take (assuming performance similar to the system you're using)?"""

# ╔═╡ 7e42beb2-f525-40c4-9332-5c0b4bd27e9c
response_1g = missing

# ╔═╡ 4d3b7121-8409-492f-8de5-7cc7859f193f
display_msg_if_fail(check_type_isa(:response_1g,response_1g,Real)) 

# ╔═╡ a847804b-886e-4046-a7b8-52f3fbfc9376
begin
	local N = floor(Int,sqrt(1024*2^30/8))
	local num_ops = (2//3)*N^3
	local my_flops = 1.5e9
	local my_est = num_ops/my_flops
	if !@isdefined(response_1g)  || ismissing(response_1g)
		nothing
	elseif response_1g < my_est/10
		almost(md"Are you sure?  That seems low to me.")
		correct()
	elseif response_1g > my_est*10
		almost(md"Are you sure?  That seems high to me.")
	else
		correct(md"That seems like a plaussible runtime to me.")	
	end
end

# ╔═╡ 702843dc-83a9-4fde-b7fa-989e1e5e6c88
md"1h. Does memory or run-time limit the largest matrix that can facotrized on a high-end server?  Why?"

# ╔═╡ 59e0467c-a8a8-43fe-a735-146e001a14fb
response_1h = missing

# ╔═╡ 64c04240-24a4-4d1e-b03d-3ca8dd1fdd38
display_msg_if_fail(check_type_isa(:response_1h,response_1h,Markdown.MD)) 

# ╔═╡ 607ebb54-065f-44e9-a66c-686dae2dcd54
md"""
### In practice

i. Following your work above, estimate how long (in seconds) it would take to solve a linear system with $N=100$ via LU factorization.  
"""

# ╔═╡ dc2116e4-5995-4b81-8d19-3656da4c34dc
response_1i = missing

# ╔═╡ 5dc61f62-57f0-4c1e-9373-d6d00931d12d
display_msg_if_fail(check_type_isa(:response_1i,response_1i,Real)) 

# ╔═╡ 62dd7c5b-21a3-4c53-99a7-b6f7e57e39f2
begin
	local N = 100
	local num_ops = (2//3)*N^3
	local my_flops = 1.5e9
	local my_est = num_ops/my_flops
	if !@isdefined(response_1i)  || ismissing(response_1i)
		nothing
	elseif response_1i < my_est/10
		almost(md"Are you sure?  That seems low to me.")
		correct()
	elseif response_1i > my_est*10
		almost(md"Are you sure?  That seems high to me.")
	else
		correct(md"That seems like a plaussible runtime to me.")	
	end
end

# ╔═╡ 261c65e6-c93e-4253-b828-a92ef90ed797
md"""h.  Now, we'll benchmark how long it actually takes to solve a linear system with $N=100$ via [LU factorization](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html#LinearAlgebra.lu) and the ["left division operator" (`\`)](https://docs.julialang.org/en/v1/base/math/#Base.:\\-Tuple{Any,%20Any}) using the following function and `@time` (_not_ `@btime`).  Importantly, we're going to repeat this a few times.  
"""

# ╔═╡ 9a8c7c33-0e71-42fa-b2a3-e6bbd17c2b81
N = 100  # Set problem size

# ╔═╡ ca8ca136-d01a-4909-90df-331b82f1b8f5
begin  # Create problem data
	A = rand(N,N)
	x = rand(N)
	y = A*x
end;

# ╔═╡ 23f7567b-3493-407e-a251-58d5dbecaed1
function solve_Ax_via_LU_factorization(A::Matrix, y::Vector)
   	local F = lu(A)   # [See manual](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/index.html#LinearAlgebra.lu)
   	local x = F \ y
end

# ╔═╡ 1d3f387e-caed-4b60-a74a-8576875b9270
with_terminal() do 
	@time x1 = solve_Ax_via_LU_factorization(A,x)
	@time x2 = solve_Ax_via_LU_factorization(A,x)
	@time x3 = solve_Ax_via_LU_factorization(A,x)
end

# ╔═╡ afeda1b5-75ae-4c82-a5eb-d6853083ce14
md"""
j.  Is there any noticeable difference in your three results?  If so, what do you think explains the difference?"""

# ╔═╡ 06da82e6-256e-4a46-b8d5-090e6595940c
response_1j = missing

# ╔═╡ 324e7a1d-cb0b-4c01-893f-672ecffcd2fa
display_msg_if_fail(check_type_isa(:response_1j,response_1j,Markdown.MD)) 

# ╔═╡ 5330262e-74a1-4c0e-8ba9-c458deeb8f6a
md"""Now, try try benchmarking the same code using the `@benchmark` macro.  (This will take several seconds.)"""

# ╔═╡ ef6e093f-613a-47ed-bcb1-f593f41cf73a
@benchmark solve_Ax_via_LU_factorization($A,$x) seconds=5

# ╔═╡ f43b83a3-2560-4036-b6f4-20cc3860efb0
md"""
k.  Is there a significant difference between the minimum and maximum time required?  If so, what do you think is the biggest effect in explaining the difference?  Which output do you think is most relevant for your typical scientific applications?  """

# ╔═╡ d985b122-1ea7-4dc7-b724-45eb77bcf146
response_1k = missing

# ╔═╡ 03a78cda-38ba-4c74-b777-1a51dacc223c
display_msg_if_fail(check_type_isa(:response_1k,response_1k,Markdown.MD)) 

# ╔═╡ e86ed6db-6e4a-4652-9e1b-7a6a5c168927
md"""l.  How does your result compare to what you estimated analytically in part i?"""

# ╔═╡ a8fe454d-34c4-4619-ace4-610975a8df5f
response_1l = missing

# ╔═╡ 50902545-34bb-4e77-88b7-afc4dac439ec
display_msg_if_fail(check_type_isa(:response_1l,response_1l,Markdown.MD)) 

# ╔═╡ 7c98f340-2719-4230-a7ae-e69db51e66cb
md"""
## Scaling with Problem Size
"""

# ╔═╡ 616fabd6-0ea7-4524-a6a1-c4c60d4f87d0
md"""
Next, lets consider several relatively small problem sizes (e.g., no larger than a few hundred) and benchmark `solve_Ax_via_LU_factorization` to see how the performance scales as we increase $N$.  We'll plot the results on a log-log scale and compare to a simplistic analytic model.
"""

# ╔═╡ cdfbad83-736b-44ac-a792-3bbc37fb1076
N_list = 2 .^(1:9)

# ╔═╡ 84539bf8-d5e4-4f75-b5fb-1197ca6221f2
begin
	time_for_one_fma = @belapsed 17.0*π+4.0
	model_time_list = time_for_one_fma * 2//3 .* N_list.^3 
end;

# ╔═╡ fb07fab4-89e4-4e59-bfb5-2e3fd5e26f37
begin
	time_list = zeros(length(N_list))
	for (i,N) in enumerate(N_list)
	  A = rand(N,N)
	  x = rand(N)
	  y = A*x
	  time_list[i] = @belapsed solve_Ax_via_LU_factorization($A,$y)
	end
	time_list
end;

# ╔═╡ 3715f7c6-48c7-4706-a7b7-1757c2891240
begin
	plt = plot()
	plot!(plt,log10.(N_list),log10.(model_time_list), xlabel=L"\log_{10} N", ylabel = L"\log_{10}(\mathrm{Time}/s)", label="Model", legend=:bottomright) 
	scatter!(plt,log10.(N_list),log10.(time_list), label="Actual") 
	plt
end

# ╔═╡ 20cd5f33-d865-470a-9ef4-0cd18a5bba00
md"""
n. How does the actual performance compare the analytic model?  What is your guess for the cause of any deviations?"""

# ╔═╡ ac0374ce-9a02-4fdb-a1e9-5e6fac9192c2
response_1n = missing

# ╔═╡ ba2c71eb-ff66-43a1-8056-982ccfb41510
display_msg_if_fail(check_type_isa(:response_1n,response_1n,Markdown.MD)) 

# ╔═╡ 9ad38265-3d12-4df5-b251-bf7f19fe8947
md"o. For real life problems, what other considerations are likely to limit performance?"

# ╔═╡ c06d4df7-b8ff-493e-b933-eee61baaf651
response_1o = missing

# ╔═╡ 732cf450-0330-43fa-aa0b-6c0fb1b02fd7
display_msg_if_fail(check_type_isa(:response_1o,response_1o,Markdown.MD)) 

# ╔═╡ 19617de2-d198-43d6-a3da-71b7cdfc0be1
md"p. How could one practically solve even larger linear systems?"

# ╔═╡ 7c0a4bda-cd4b-4bef-a2f6-84efe76670aa
response_1p = missing

# ╔═╡ ec7d8268-7f97-4ba3-a0ab-717506002fb6
display_msg_if_fail(check_type_isa(:response_1p,response_1p,Markdown.MD)) 

# ╔═╡ fb3f1372-a7fc-4698-a8ce-d96354520a63
md"""[^1]: Acknowledgment:  The questions in this subsection are based on Oliveira & Stewarts Writing Scientific Software, Chapter 5, Problem #6.
"""

# ╔═╡ 07653065-2ef3-4a63-a25b-1b308c22aff5
md"## Helper code"

# ╔═╡ 5b88a5a8-425b-4ecf-a26a-08722d12ef95
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1;

# ╔═╡ 9b2cdb77-9330-4a8b-84b5-deed94c19662
begin
"""
	   benchmark_my_funciton(f,n_list)
	   benchmark_my_funciton(f,x_list)
	   benchmark_my_funciton(f,x_list, y_list)
	
	Benchmarks a user-proved function.  
	User-provided function may take the number of samples, one array or two arrays.
	Returns NamedTuple with two lists (num_list, times_list) containing the number of samples and the runtime.
"""
	function benchmark_my_funciton end
	
	function benchmark_my_funciton(f::Function, num_list::Vector{T} ) where { T<:Integer }
		times_list = zeros(length(num_list))
		for (i,n) in enumerate(num_list)
			times_list[i] = @belapsed $f($n)
		end
		return (;num_list, times_list)
	end
	
	function benchmark_my_funciton(f::Function, x_list::Vector{A} ) where { A<:AbstractArray } 
		times_list = zeros(length(x_list))
		for (i,x) in enumerate(x_list)
			times_list[i] = @belapsed $f($x)
		end
		return (;num_list, times_list)
	end
	
	function benchmark_my_funciton(f::Function, x_list::Vector{A}, y_list::Vector{A} ) where { A<:AbstractArray } 
		@assert length(x_list) == length(y_list)
		times_list = zeros(length(x_list))
		for i in 1:length(x_list)
			x = x_list[i]
			y = y_list[i]
			times_list[i] = @belapsed $f($x,$y)
		end
		return (;num_list, times_list)
	end
end

# ╔═╡ 6b2d4a5c-80bc-4620-9e8a-d8684302a9f2
benchmarks_0 = benchmark_my_funciton(my_function_0_args, num_list)

# ╔═╡ 17c29c66-91f8-45e2-b1f2-337dd51a6e03
begin
		plt0 = scatter(benchmarks_0.num_list, log10.(benchmarks_0.num_list./benchmarks_0.times_list), xscale=:log10, label=:none)
		xlabel!(plt0, "Size of Array")
		ylabel!(plt0, "log_10 (Evals/s)")
		title!(plt0,"Benchmarks for my_function_0_args")
		plt0
end

# ╔═╡ 4e1f8fa4-b846-4cbd-a5de-b9e137ec04f9
benchmarks_sqrt = benchmark_my_funciton(sqrt_broadcasted, x_list)

# ╔═╡ 185fac9d-ea9e-460a-8f90-d5dad528ed4b
if my_function_1_arg_is_good_to_go
	benchmarks_1 = benchmark_my_funciton(my_function_1_arg, x_list)
end

# ╔═╡ 79716270-4570-41cb-9746-394eead121ee
begin
		plt1 = plot()
		scatter!(plt1,benchmarks_sqrt.num_list, log10.(benchmarks_sqrt.num_list./benchmarks_sqrt.times_list), xscale=:log10, label="sqrt", legend=:topleft)
		if my_function_1_arg_is_good_to_go
			scatter!(plt1,benchmarks_1.num_list, log10.(benchmarks_1.num_list./benchmarks_1.times_list), label="my_function_1_args")
	end
		xlabel!(plt1, "Size of Array")
		ylabel!(plt1, "log_10 (Evals/s)")
		title!(plt1,"Benchmarks for univariate functions")
		plt1
end

# ╔═╡ d5322db3-501f-467f-a3f4-297258bc2570
benchmarks_add = benchmark_my_funciton(add_broadcasted, x_list, y_list)

# ╔═╡ fffd69c6-47d6-411c-ba57-9c52e659f598
benchmarks_mul = benchmark_my_funciton(mul_broadcasted, x_list, y_list)

# ╔═╡ 4a45fa07-f4dd-4a76-8bd1-8627060a6fc5
benchmarks_div = benchmark_my_funciton(div_broadcasted, x_list, y_list)

# ╔═╡ 930a527f-9eef-43bf-970d-3a17285a233c
if my_function_2_args_is_good_to_go
	benchmarks_2 = benchmark_my_funciton(my_function_2_args, x_list, y_list)
end

# ╔═╡ d73bc561-1063-4020-bba8-d89464d31254
begin
		plt2 = plot()
		scatter!(plt2,benchmarks_add.num_list, log10.(benchmarks_add.num_list./benchmarks_add.times_list), xscale=:log10, label="add Float64s", legend=:topleft)
		scatter!(plt2,benchmarks_mul.num_list, log10.(benchmarks_mul.num_list./benchmarks_mul.times_list), label="multiply Float64s")
		scatter!(plt2,benchmarks_div.num_list, log10.(benchmarks_div.num_list./benchmarks_div.times_list), label="divide Float64s")
		if my_function_2_args_is_good_to_go
			scatter!(plt2,benchmarks_2.num_list, log10.(benchmarks_2.num_list./benchmarks_2.times_list), label="my_function_2_args")
		end

		xlabel!(plt2, "Size of Array")
		ylabel!(plt2, "log_10 (Runtime/s)")
		title!(plt2,"Benchmarks for functions of 2 variables")
		plt2
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.1.1"
LaTeXStrings = "~1.2.1"
Plots = "~1.19.3"
PlutoTeachingTools = "~0.1.2"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "c31ebabde28d102b602bada60ce8922c266d205b"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random", "StaticArrays"]
git-tree-sha1 = "ed268efe58512df8c7e224d2e170afd76dd6a417"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.13.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "344f143fa0ec67e47917848795ab19c6a455f32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.32.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f473cdf6e2eb360c576f9822e7c765dd9d26dbc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "eaf96e05a880f3db5ded5a5a8a7817ecba3c7392"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "44e3b40da000eab4ccb1aecdc4801c040026aeb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.13"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "94bf17e83a0e4b20c8d77f6af8ffe8cc3b386c0a"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.1"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "1e72752052a3893d0f7103fbac728b60b934f5a5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.19.4"

[[PlutoTeachingTools]]
deps = ["LaTeXStrings", "Markdown", "PlutoUI", "Random"]
git-tree-sha1 = "265980831960aabe7e1f5ae47c898a8459588ee7"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.1.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "885838778bb6f0136f8317757d7803e0d81201e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.9"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─4ac5cc73-d8d6-43d5-81b2-944d559fd2ca
# ╠═637f6a84-ad01-43a7-899b-b7867fbe3b3d
# ╟─c23045a0-56b2-4e1e-b5d3-8e248fd1bffd
# ╠═6dcedd62-c606-43ba-b2f6-e018aefcb035
# ╟─e3bfd6bc-261d-4e69-9ce6-e78d959c13be
# ╠═0854551f-fc6d-4c51-bf85-a9c7030e588b
# ╟─2a5eb6c7-9b8c-4923-8fe4-e7211c6c820f
# ╠═b3e64508-319f-4506-8512-211e30be4bee
# ╟─a0a24490-3136-45e4-8f55-44448d8366d1
# ╠═69c39140-ba66-4345-869c-a5499d5d4376
# ╟─2e0b53c8-4107-485f-9f33-d921b6fc0c05
# ╠═6b2d4a5c-80bc-4620-9e8a-d8684302a9f2
# ╟─17c29c66-91f8-45e2-b1f2-337dd51a6e03
# ╟─928fe18c-2b34-427f-afdb-6273a0ce135e
# ╠═fd7c88ef-02fb-408f-a551-45a93d873ded
# ╟─ba4b35c8-0fe8-4a25-9817-d13c0cd98c6f
# ╠═bb4e3e7e-d8da-4c53-b861-0e66237aae3c
# ╟─4e1f8fa4-b846-4cbd-a5de-b9e137ec04f9
# ╟─6e1e699f-4f4a-4a21-946f-e763eb106f37
# ╠═39ed7918-8926-4866-8b25-c61dbbd35991
# ╟─a403c0f3-fc4a-4bc1-87c4-8a2d5d15e1a0
# ╟─185fac9d-ea9e-460a-8f90-d5dad528ed4b
# ╟─79716270-4570-41cb-9746-394eead121ee
# ╟─5d234d6c-be3d-4d0f-bd58-774b6786db54
# ╠═11df2b04-2bef-497e-a064-cbcd159aecc7
# ╠═d5322db3-501f-467f-a3f4-297258bc2570
# ╠═8320739f-ec14-4a49-b374-a0baa198f646
# ╠═fffd69c6-47d6-411c-ba57-9c52e659f598
# ╠═e352d85a-b8ce-4355-8b0b-9c28465dd006
# ╠═4a45fa07-f4dd-4a76-8bd1-8627060a6fc5
# ╟─0ce6409c-3aa4-4446-ae13-81c4e743d322
# ╠═f13619d1-6ece-42be-965d-ab6bb9e9b8cd
# ╟─340808ea-e99b-4de7-ad55-b2d812ff0f4d
# ╠═930a527f-9eef-43bf-970d-3a17285a233c
# ╟─d73bc561-1063-4020-bba8-d89464d31254
# ╟─3621793a-427b-40d2-b28d-7bb9c6f3a28c
# ╟─10499798-1ba6-4a2f-827b-aecc4a1f8346
# ╠═ffde1960-8792-4092-9a0b-42b2726bb2da
# ╟─10f54651-4b05-4806-9a0c-cb0b6afa245b
# ╟─2291567f-71a4-4a53-a80c-aab58f29ddf8
# ╠═05e772d4-a2ad-4f16-866c-b8aa3ab7a550
# ╟─eb1ff48b-7dbc-4b37-a18d-77eadc568f5a
# ╟─c888afb4-59da-4d7b-a050-df84f3735202
# ╟─b9157cf9-99da-46ee-8ba0-ed6a95153f74
# ╠═b0873452-b9b5-4bc5-bb76-5f299a6d366e
# ╟─cc584d2a-2a3b-4898-b753-1a7275712ad2
# ╟─e69c7809-ec6c-4f3a-b869-100381d40bf9
# ╟─ea68d6eb-a542-44bf-b0e7-815e2372bc33
# ╠═484b47d6-65de-497a-b78f-6996c8787de8
# ╟─8e2045c2-c5de-4895-9719-efb9de113dec
# ╟─2da17c61-a8d7-4dc2-8209-69fe809e9d4f
# ╟─4e8e86ed-92b0-4aa5-945b-ef6193cb2559
# ╠═1d4275fd-3208-477d-b26e-be23940035ba
# ╟─c64bd8d8-6801-46a7-876b-c75b83300416
# ╟─5cc54a66-461c-4bb3-8e29-2079adeff04f
# ╠═30efc3ea-a2a2-48f8-833f-6b51845e47dc
# ╟─33bc47b6-f48b-4e24-bb54-fa794f04e166
# ╟─3fab2e27-766c-4464-9b6c-e5235a80183a
# ╟─9af89fd2-a5f6-4807-beb6-d3fa9fa89b69
# ╠═7e42beb2-f525-40c4-9332-5c0b4bd27e9c
# ╟─4d3b7121-8409-492f-8de5-7cc7859f193f
# ╟─a847804b-886e-4046-a7b8-52f3fbfc9376
# ╟─702843dc-83a9-4fde-b7fa-989e1e5e6c88
# ╠═59e0467c-a8a8-43fe-a735-146e001a14fb
# ╟─64c04240-24a4-4d1e-b03d-3ca8dd1fdd38
# ╟─607ebb54-065f-44e9-a66c-686dae2dcd54
# ╠═dc2116e4-5995-4b81-8d19-3656da4c34dc
# ╟─5dc61f62-57f0-4c1e-9373-d6d00931d12d
# ╟─62dd7c5b-21a3-4c53-99a7-b6f7e57e39f2
# ╟─261c65e6-c93e-4253-b828-a92ef90ed797
# ╟─9a8c7c33-0e71-42fa-b2a3-e6bbd17c2b81
# ╠═ca8ca136-d01a-4909-90df-331b82f1b8f5
# ╠═cb099d3d-5c09-4c56-9c59-83adab60f651
# ╠═23f7567b-3493-407e-a251-58d5dbecaed1
# ╠═1d3f387e-caed-4b60-a74a-8576875b9270
# ╟─afeda1b5-75ae-4c82-a5eb-d6853083ce14
# ╠═06da82e6-256e-4a46-b8d5-090e6595940c
# ╟─324e7a1d-cb0b-4c01-893f-672ecffcd2fa
# ╟─5330262e-74a1-4c0e-8ba9-c458deeb8f6a
# ╠═ef6e093f-613a-47ed-bcb1-f593f41cf73a
# ╟─f43b83a3-2560-4036-b6f4-20cc3860efb0
# ╟─d985b122-1ea7-4dc7-b724-45eb77bcf146
# ╟─03a78cda-38ba-4c74-b777-1a51dacc223c
# ╟─e86ed6db-6e4a-4652-9e1b-7a6a5c168927
# ╟─a8fe454d-34c4-4619-ace4-610975a8df5f
# ╟─50902545-34bb-4e77-88b7-afc4dac439ec
# ╟─7c98f340-2719-4230-a7ae-e69db51e66cb
# ╟─616fabd6-0ea7-4524-a6a1-c4c60d4f87d0
# ╠═cdfbad83-736b-44ac-a792-3bbc37fb1076
# ╟─84539bf8-d5e4-4f75-b5fb-1197ca6221f2
# ╟─fb07fab4-89e4-4e59-bfb5-2e3fd5e26f37
# ╟─3715f7c6-48c7-4706-a7b7-1757c2891240
# ╟─20cd5f33-d865-470a-9ef4-0cd18a5bba00
# ╠═ac0374ce-9a02-4fdb-a1e9-5e6fac9192c2
# ╟─ba2c71eb-ff66-43a1-8056-982ccfb41510
# ╠═9ad38265-3d12-4df5-b251-bf7f19fe8947
# ╠═c06d4df7-b8ff-493e-b933-eee61baaf651
# ╠═732cf450-0330-43fa-aa0b-6c0fb1b02fd7
# ╟─19617de2-d198-43d6-a3da-71b7cdfc0be1
# ╟─7c0a4bda-cd4b-4bef-a2f6-84efe76670aa
# ╠═ec7d8268-7f97-4ba3-a0ab-717506002fb6
# ╟─fb3f1372-a7fc-4698-a8ce-d96354520a63
# ╟─07653065-2ef3-4a63-a25b-1b308c22aff5
# ╠═a1122848-347f-408b-99c2-a7a514073864
# ╠═4df26174-9aa3-46d0-879b-aec9a771714b
# ╟─ea0a7ca2-503f-4d61-a3d4-42503f322782
# ╠═5b88a5a8-425b-4ecf-a26a-08722d12ef95
# ╠═9b2cdb77-9330-4a8b-84b5-deed94c19662
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
