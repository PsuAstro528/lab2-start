### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ af508570-b20f-4dd3-a995-36c79fc41823
begin
	using PlutoUI, PlutoTeachingTools, PlutoTest
	using BenchmarkTools
	#eval(Meta.parse(code_for_check_type_funcs))
end

# ╔═╡ 27667e0a-8ebc-4397-8ac3-33a0f19f6987
md"""
#  Astro 528, Lab 2, Exercise 3
## Unit Tests & Assertions
### (as demonstrated on Solving Kepler's equation)
"""


# ╔═╡ bedbe9b9-03be-4537-b938-89d2857d0cba
md"""
[Kepler's Equation](https://en.wikipedia.org/wiki/Kepler%27s_equation) relates the *mean anomaly* ($M$) and *eccentric anomaly* ($E$), two angles that can be used to specify where a body is along its orbit in the two-body (aka, Kepler) problem.  
```math
M = E - e \sin(E),
```
where $e$ is the orbital eccentricity.  The mean anomaly increases linearly with time, but the eccentric anomaly does not.
"""

# ╔═╡ 8116b382-c927-4563-8fde-dc034dd96ab9
md"1a.  Update the function `calc_mean_anom` below to compute the mean anomaly for a given eccentric anomaly and eccentricity."

# ╔═╡ c9eef6f8-4583-429b-a3ae-09f1ec8e5ecf
"""   calc_mean_anom(ecc_anom, e) 
Calculate mean anomaly for a given eccentric anomaly and eccentricity.
"""
function calc_mean_anom(ecc_anom::Real, e::Real) 
	missing
end;

# ╔═╡ 05cf7bc1-cc04-4679-9534-05834f097371
if !@isdefined(calc_mean_anom)
   func_not_defined(:calc_mean_anom)
else
	if !(length(methods(calc_mean_anom,[Float64,Float64])) >= 1)
		   keep_working(md"Your calc_mean_anom can't take two Float64's as arguments.")
	elseif ismissing(calc_mean_anom(1.0,0.5))
		still_missing()
	elseif calc_mean_anom(1.0,0.5) ≈ 0.5792645075960517
		correct()
	else
		keep_working()	
	end
end

# ╔═╡ 3ebe8069-11d3-4885-aec8-72e2f3f5f906
md"""Solving the Kepler Equation for $E$ given for a given $e$ and $E$ is critical to be able to determine where a body is at a given time.  Since it is a transendental equation (i.e., there is no closed form algebraic solution), it is solved itteratively.  Over the centuries, there have been numerous studies of how to solve the Kepler equation efficiently.  Below, I've coded up an implementation for you to help you get stared.  For this exercise, you don't need to understand the details of how the algorithm works.  Instead, we will focus on how we can use principles of modern software development to improve these functions via assertions and unit tests. 
"""

# ╔═╡ f6be4fa8-351a-4c1c-b389-b79f09db2a4b
md"## Starter code to solve Kepler's Equation"

# ╔═╡ 6421c10b-9429-45ef-8ffb-b5117fae9e58
"""
   ecc_anom_init_guess_danby(mean_anomaly, eccentricity)

Returns initial guess for the eccentric anomaly for use by itterative solvers of Kepler's equation for bound orbits.  

Based on "The Solution of Kepler's Equations - Part Three"
Danby, J. M. A. (1987) Journal: Celestial Mechanics, Volume 40, Issue 3-4, pp. 303-312 (1987CeMec..40..303D)
"""
function ecc_anom_init_guess_danby(M::Real, ecc::Real)
	@assert -2π<= M <= 2π
	@assert 0 <= ecc < 1.0
    if  M < zero(M)
		M += 2π
	end
    E = (M<π) ? M + 0.85*ecc : M - 0.85*ecc
end;

# ╔═╡ fbb8c035-da0d-45d0-856c-0668f2ef954f
"""
   update_ecc_anom_laguerre(eccentric_anomaly_guess, mean_anomaly, eccentricity)

Update the current guess for solution to Kepler's equation
  
Based on "An Improved Algorithm due to Laguerre for the Solution of Kepler's Equation"
   Conway, B. A.  (1986) Celestial Mechanics, Volume 39, Issue 2, pp.199-211 (1986CeMec..39..199C)
"""
function update_ecc_anom_laguerre(E::Real, M::Real, ecc::Real)
  es = ecc*sin(E)
  ec = ecc*cos(E)
  F = (E-es)-M
  Fp = one(M)-ec
  Fpp = es
  n = 5
  root = sqrt(abs((n-1)*((n-1)*Fp*Fp-n*F*Fpp)))
  denom = Fp>zero(E) ? Fp+root : Fp-root
  return E-n*F/denom
end;

# ╔═╡ d48ca14f-2b62-4f35-8c8a-07aa3563b579
"""
   calc_ecc_anom_itterative_laguerre( mean_anomaly, eccentricity )

Estimates eccentric anomaly for given mean_anomaly and eccentricity.
Optional parameter `tol` specifies tolerance (default 1e-8)
"""
function calc_ecc_anom(mean_anom::Real, ecc::Real; tol::Real = 1.0e-8)
  	@assert 0 <= ecc < 1.0
	@assert 1e-16 <= tol < 1
  	M = rem2pi(mean_anom,RoundNearest)
    E = ecc_anom_init_guess_danby(M,ecc)
	local E_old
    max_its_laguerre = 200
    for i in 1:max_its_laguerre
       E_old = E
       E = update_ecc_anom_laguerre(E_old, M, ecc)
       if abs(E-E_old) < tol break end
    end
    return E
end;

# ╔═╡ 76c08fb3-89e1-4897-90f1-2f73ba298010
md"""

### Assertions

Sometimes a programmer calls a function with arguments that either don't make sense or represent a case that the function was not originally designed to handle properly. The worst possible function behavior in such a case is returning an incorrect result without any warning that something bad has happened. Returning an error at the end is better, but can make it difficult to figure out the problem. Generally, the earlier the problem is spotted, the easier it will be to fix the problem. Therefore, good developers often include assertions to verify that the function arguments are acceptable.  

For example, in `ecc_anom_init_guess_danby` above, we included an assertion that the eccentricity was positive-semidefinite and less than or equal to unity.  

1b. What other preconditions should be met for the inputs to the functions above?    What is about `calc_mean_anom` from 1a?
"""

# ╔═╡ 693371c4-8e35-4a8e-9fa2-8c0e441515ac
response_1b = missing  # md"YOUR RESPONSE"

# ╔═╡ 33f00289-8f32-4512-968e-d29fcd6968ff
if !@isdefined(response_1b)
   var_not_defined(:response_1b)
elseif ismissing(response_1b)
	still_missing()
else
	nothing
end

# ╔═╡ a3e0ad59-288a-423d-a7e9-ecadb055e0dc
md"1c. Update the code above to include at least one additional assertion."

# ╔═╡ e50297d5-8599-48dd-935e-0c3975e4e379
begin
	num_evals = 100
	mean_anoms_for_benchmarks = 2π*rand(num_evals)
	eccs_for_benchmarks = rand(num_evals)
end

# ╔═╡ c0a54587-86bf-4f6d-9e5c-574badc06865
md"1d.  Adding assertions creates extra work.  It's good to think about whether an assertion will result in a significant performance hit.  Benchmark your code before an after adding the assertions.  How does the typical run time compare?  What are the implications for whether it makes sense to leaves the assertions in a production code?"

# ╔═╡ 61db6108-8be6-4544-b1be-6e448d99748f
@benchmark calc_ecc_anom.($mean_anoms_for_benchmarks,$eccs_for_benchmarks)

# ╔═╡ a394967a-16cb-42fe-a1fc-797e108439d3
response_1d = missing

# ╔═╡ 6569d30f-0c9d-4533-a421-086fcf1b5f62
if !@isdefined(response_1d)
   var_not_defined(:response_1d)
elseif ismissing(response_1d)
	still_missing()
else
	nothing
end

# ╔═╡ 90471f40-c548-494f-84cb-871ee0f3f5f9
md"""
## Unit Tests
Units tests check that the post-conditions are met for at least some certain test inputs.  I'll demonstrate a couple of unit tests below.
"""

# ╔═╡ e1a70885-ac8b-4b4c-9080-92c3178f5a03
@test calc_ecc_anom(0,0.1) ≈ 0 atol = 1e-8

# ╔═╡ a6034884-89a1-43d7-b30e-f389ae2e05d3
@test calc_ecc_anom(Float64(π),0.5) ≈ π atol = 1e-8

# ╔═╡ 48aed012-2da1-4a06-9404-a34cdd4b3eee
md"""
Note that testing equality or inequality is straightforward for integers, but dangerous for floating point numbers.  If a floating point number is equal to an integer, then it's generally ok, but I think it's better to always be cautious about testing floating  point numbers for equality.  Instead, you can test that two numbers are approximately equal using $\simeq$.  When testing that two numbers are approximately equal you need to specify a tolerance.  `atol` refers to an absolute tolerance, while `rtol` refers to a relative or fractional tolerance. For further information on using these, see the [Julia manual](https://docs.julialang.org/en/v1/stdlib/Test/index.html#Basic-Unit-Tests-1).  (Technically, we're using `PlutoTest.@test`, which in implemented in terms of `Test.@test`.
"""

# ╔═╡ 9a89ccd7-f690-4dc1-9736-779aa7845ab2
tip(md"In Pluto, Jupyter, VSCode and Atom (and probably many other modern IDEs), you can get unicode characters like ≈ by typing `\approx<tab>`.")

# ╔═╡ b17e6d66-9a19-48a8-a6ba-c8d5145387f3
md"""
1d. What other unit tests could be useful for diagnosing any errors or non-robust behavior?  Write at least three new unit tests that help to check whether the above code is accurate.  After writing your first one, try to think of how to make the next ones more useful than just doing very similar tests over and over again.
Try think of a corner case where a non-robust algorithm might not be accurate.  
Explain your testing plan below.
"""

# ╔═╡ 4d6d13f7-47b0-4e8a-af69-3aa0a8cb8d2d
response_1e = missing

# ╔═╡ c5b8a22d-d47b-4f65-b90e-544e8b6c3a88
if !@isdefined(response_1e)
   var_not_defined(:response_1e)
elseif ismissing(response_1e)
	still_missing()
else
	nothing
end

# ╔═╡ 4ec0d903-4e3f-42d6-8853-86c848bb92ce
md"""
Check what happens when your tests run.  Do your functions pass all of them?  If not, correct the function (or the tests if necessary) and rerun the tests. """

# ╔═╡ 869e40bc-0679-4b20-acff-f8db6637a887
md"""
### Testing the assertions.
In this case, the assertions are probably pretty simple.  But sometimes, the assertions can be complicated enough that you'll need to test that they're working as intended.  When an assert statement is followed by an expression that evaluates to false, then it ["throws an exception"](https://docs.julialang.org/en/v1.0/manual/control-flow/#Exception-Handling-1).  We want to make sure that our code is throwing an exception when we pass our function invalid arguments.  I'll demonstrate with a test of passing an invalid eccentricity.
"""

# ╔═╡ dccdcf37-3354-4717-8cf3-b391bbaf6e9a
@test_throws AssertionError calc_ecc_anom(rand()*2π,-0.5) 

# ╔═╡ 94891c38-8e36-43f5-b0a6-9b8b38145ef9
md"""1f. Add tests to make sure that the assertions you added in 1c trigger as you intended them to."""


# ╔═╡ a9af0fa1-99af-447d-86ad-cb926f7b9de1
# INSERT YOUR TEST(S) HERE

# ╔═╡ 2ca10bfe-549b-4910-95a3-e68a292df0d6
md"""
### Continous Integration Testing.
Often, a well-intentioned programmer introduces a bug, but doesn't notice until long after the bug was written.  One way to reduce the risk of such bugs is to have an comprehensive set of unit tests that are applied _automatically_ each time a developer commits a change.  If some new code causes a test to fail, we want to know that promptly, so it can be fixed and before it causes scientists to lose time running the buggy code or trying to interpret results of a buggy code.

The '.github/workflows/test.yaml' file provided in this repository already provides instructions for GitHub.com to automatically run tests each time you commit changes and push them to GitHub.  The tests for this notebook are in `tests/test3.jl`.  

1g.  Add the tests that you wrote above to `tests/test3.jl`, so that they become part of your repository's _continuous integration_ testing.  Commit those changes to your local repository and push them to github.  Go to your repository on GitHub, click "Actions", then click "Test notebooks", then look at the line for most recent commit (at the top of the list).  Is there a green check or a red x?  If a red x, then click on the commit message for the most recent commit, and then click "test (v1.6, x86, ubuntu-latest)" to see the messages associated with that test run.  See if you can figure out why a test failed.  
"""

# ╔═╡ a93f0326-78ca-407d-aec2-5de318c002ca
tip(md"If a test does fail, then it may be useful to run the tests on Roar or your local machine, so you get faster feedback on whether you're tests are performing as expected.  From your local repository's directory, you can run `julia --project=test test/test3.jl ` to run the tests.  Once you get those working, then commit the changes to your local repo, push to GitHub, and see if the new version of your code passes all its test.")

# ╔═╡ b760fedd-41ea-4784-845f-ede0163c0d12
md"## Helper Code"

# ╔═╡ 2e893623-2d05-48ac-b7c6-0ba167dc7419
ChooseDisplayMode()

# ╔═╡ bfdd8ecf-5f05-4056-a9d8-f3404774ff52
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.1.1"
PlutoTeachingTools = "~0.1.2"
PlutoTest = "~0.1.0"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Statistics", "UUIDs"]
git-tree-sha1 = "c31ebabde28d102b602bada60ce8922c266d205b"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.1.1"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[HypertextLiteral]]
git-tree-sha1 = "1e3ccdc7a6f7b577623028e0095479f4727d8ec1"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.8.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "94bf17e83a0e4b20c8d77f6af8ffe8cc3b386c0a"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.1"

[[PlutoTeachingTools]]
deps = ["LaTeXStrings", "Markdown", "PlutoUI", "Random"]
git-tree-sha1 = "265980831960aabe7e1f5ae47c898a8459588ee7"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.1.3"

[[PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "3479836b31a31c29a7bac1f09d95f9c843ce1ade"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.1.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╟─27667e0a-8ebc-4397-8ac3-33a0f19f6987
# ╟─bedbe9b9-03be-4537-b938-89d2857d0cba
# ╟─8116b382-c927-4563-8fde-dc034dd96ab9
# ╠═c9eef6f8-4583-429b-a3ae-09f1ec8e5ecf
# ╟─05cf7bc1-cc04-4679-9534-05834f097371
# ╟─3ebe8069-11d3-4885-aec8-72e2f3f5f906
# ╟─f6be4fa8-351a-4c1c-b389-b79f09db2a4b
# ╠═6421c10b-9429-45ef-8ffb-b5117fae9e58
# ╠═fbb8c035-da0d-45d0-856c-0668f2ef954f
# ╠═d48ca14f-2b62-4f35-8c8a-07aa3563b579
# ╟─76c08fb3-89e1-4897-90f1-2f73ba298010
# ╠═693371c4-8e35-4a8e-9fa2-8c0e441515ac
# ╟─33f00289-8f32-4512-968e-d29fcd6968ff
# ╟─a3e0ad59-288a-423d-a7e9-ecadb055e0dc
# ╠═e50297d5-8599-48dd-935e-0c3975e4e379
# ╟─c0a54587-86bf-4f6d-9e5c-574badc06865
# ╠═61db6108-8be6-4544-b1be-6e448d99748f
# ╠═a394967a-16cb-42fe-a1fc-797e108439d3
# ╟─6569d30f-0c9d-4533-a421-086fcf1b5f62
# ╟─90471f40-c548-494f-84cb-871ee0f3f5f9
# ╠═e1a70885-ac8b-4b4c-9080-92c3178f5a03
# ╠═a6034884-89a1-43d7-b30e-f389ae2e05d3
# ╟─48aed012-2da1-4a06-9404-a34cdd4b3eee
# ╟─9a89ccd7-f690-4dc1-9736-779aa7845ab2
# ╟─b17e6d66-9a19-48a8-a6ba-c8d5145387f3
# ╠═4d6d13f7-47b0-4e8a-af69-3aa0a8cb8d2d
# ╟─c5b8a22d-d47b-4f65-b90e-544e8b6c3a88
# ╟─4ec0d903-4e3f-42d6-8853-86c848bb92ce
# ╟─869e40bc-0679-4b20-acff-f8db6637a887
# ╠═dccdcf37-3354-4717-8cf3-b391bbaf6e9a
# ╟─94891c38-8e36-43f5-b0a6-9b8b38145ef9
# ╠═a9af0fa1-99af-447d-86ad-cb926f7b9de1
# ╟─2ca10bfe-549b-4910-95a3-e68a292df0d6
# ╟─a93f0326-78ca-407d-aec2-5de318c002ca
# ╟─b760fedd-41ea-4784-845f-ede0163c0d12
# ╠═2e893623-2d05-48ac-b7c6-0ba167dc7419
# ╠═af508570-b20f-4dd3-a995-36c79fc41823
# ╠═bfdd8ecf-5f05-4056-a9d8-f3404774ff52
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
