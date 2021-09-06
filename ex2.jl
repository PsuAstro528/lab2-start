### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ c841608b-12c3-4220-aff7-843176ff052f
begin
	# Load packages that we'll be using for this exercise
	begin 
		using DifferentialEquations
		using LinearAlgebra
		using Plots
		using Printf
		using LaTeXStrings
		using PlutoUI, PlutoTeachingTools
	end
	eval(Meta.parse(code_for_check_type_funcs))
end

# ╔═╡ 308b9284-956b-4a9f-93f9-8397f34bed94
md"""
# Astro 528, Lab 2, Exercise 2
"""

# ╔═╡ e7cdc134-e45b-4d3a-b614-7153a237659e
ChooseDisplayMode()

# ╔═╡ 83e6e5a6-1347-46a5-aa06-61535cdc1120
md"""
## Numerical stability of N-body integrators

According to our reading, one of the top priorities for scientific computing is numerical stability.  It can be hard to appreciate why numerical stability is important.  In order to illustrate the importance of numerical stability, we will consider a classic problem, integrating the trajectory of a star and planet interacting under Newtonian gravity.  

This can be expressed as integrating a set of $3N$ second-order ordinary differential equations (ODEs).
```math
\frac{d^2 x_i}{dt^2} = \sum_{j=1}^N \frac{G m_j (x_j-x_i)}{\left|x_j-x_i\right|^3}
```
where $x_i$ is the position vector of body $i$, $m_i$ is the mass of body $i$, $G$ is the gravitational constant, and $N$ is the number of bodies in the system.  
Typically, this is rewritten as a set of $6N$ first-order ODEs.

```math
\begin{eqnarray}
\frac{dx_i}{dt} & = & v_i \\
\frac{dv_i}{dt} & = & \sum_{j=1}^N \frac{G m_j (x_j-x_i)}{\left|x_j-x_i\right|^3}
\end{eqnarray}
```
where the $v_i$'s are the velocity vectors.  
"""

# ╔═╡ 90139e48-c737-4117-b8ec-730d40e3541f
md"""
In this exercise, you will compare the results of integrating these equations using multiple different algorithms and parameters.  We'll start with simple integration algorithms that one could easily code on their own, but then progress to more sophisticated integration algorithms that would require a fair bit of work to implement oneself.  In order to allow you to focus on the main point of this exercise, we'll use Julia's [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).  The package has [good documentation](http://docs.juliadiffeq.org/stable/), but I'll provide nearly all the code you need to do the calculations.  Your job is to read through the code enough that you can tell what it's doing, to think about what you expect to happen, and then compare your expectations to the actual results.  If you read over the code, it may help you start to get comfortable with julia programming, but you don't need to worry about every little bit of syntax for now.  

First, we need to write a function that computes the derivatives of each of the positions and velocities.  
"""

# ╔═╡ de3ffc8c-1837-4df7-a977-363d5bd20117
"""
   gravity_as_first_order_system(du, u, mass, t)

Computes velocities and accelerations due to Newtonian gravity
for a system of N massive particles.

# Inputs:  
- u: array of 3*N positions followed by 3*N velocities
  [x_1, y_1, z_1, x_2, y_2, z_2, ...., vx_1, vy_1, vz_1, ..., vz_N]
- mass: array of N masses
- t: time

# Outputs:
- du: Aray of 3*N velocities followed by 3*N accelerations
  [vx_1, vy_1, vz_1, vx_2, vy_2, vz_2, ...., ax_1, ay_1, az_1, ..., az_N]

# Assumptions:
- Gravitational constant = 1.  
Therefore, if masses are in solar masses and distances are in AU, 
then the time unit is such that one year = 2pi.

"""
function gravity_as_first_order_system(du, u, mass, t)
    @assert length(du) == length(u) == 6 * length(mass)
    N = length(mass)
    @assert 2 <= N      # require at least two massive bodies
    for k in 1:(3*N) 
       du[k] = u[3N+k]  # derivative of position is simply velocity
    end
    du[3N+1:end] .= 0.0   # initialize accelerations to zero
    for i in 0:(N-2)      # loop over each pairs of bodies
        for j in (i+1):(N-1)
            dx = u[3*i+1] - u[3*j+1]  # displacements
            dy = u[3*i+2] - u[3*j+2]
            dz = u[3*i+3] - u[3*j+3]
            # calculate distance once per pair
            d = sqrt(dx^2+dy^2+dz^2)  
            # derivatives of velocities are accelerations 
            # acceleration on body i due to body j
            du[3N+3*i+1] -= mass[j+1] * dx / d^3 
            du[3N+3*i+2] -= mass[j+1] * dy / d^3 
            du[3N+3*i+3] -= mass[j+1] * dz / d^3 
            # reuse computation of d for acceleration on body j due to body i
            du[3N+3*j+1] += mass[i+1] * dx / d^3 
            du[3N+3*j+2] += mass[i+1] * dy / d^3 
            du[3N+3*j+3] += mass[i+1] * dz / d^3 
        end
    end
end

# ╔═╡ 05f2bcdf-9ec0-43b9-9d94-5a8592e80ce1
md"""
We could have simplified the algorithm slightly by having both `i` and `j` loop over each planet.  Instead, the limits of each loop are set so that we only consider each pair once.  The most expensive part of calculating the acceleration is the `sqrt` function.  Thus, we speed up the calculation by about a factor of two. 
"""

# ╔═╡ 2fee29a5-ecdf-4516-af26-876191b0271f
md"""
Next, we will set the initial conditions for our integration.  For testing purposes, we want to create a system where we have a good understanding of the correct behavior.  We'll setup just two bodies, assign them masses analogous to the Sun and a Jupiter-mass planet with a semi-major axis of 1 AU, and place them on the $x$-axis, moving in the $x$-$y$ plane, so that they should follow a circular orbit
"""

# ╔═╡ f4470846-8d69-47e0-aa55-df3a12dedd10
begin
	# Set initial conditions
	m = [1.0, 0.001]               # Masses: Sun and Jupiter
	init_separation = 1.0          # separation: 1 AU
	# Place star and planet on a nearly circular orbit
	init_velocity = sqrt(sum(m)/init_separation) # Uniform circular motion
	year = 2pi*sqrt(init_separation^3/sum(m))    # Kepler's third law
	r_init = init_separation .* [-m[2]/sum(m), 0, 0, m[1]/sum(m), 0, 0 ]
	v_init = init_velocity   .* [0, -m[2]/sum(m), 0, 0, m[1]/sum(m), 0 ]
	# DifferentialEquations.jl wants the initial conditions as a single vector
	u_init = vcat(r_init, v_init);  # concatenate positions & velocities
end

# ╔═╡ c7721744-061e-49e1-9470-0fd13bb960e9
md"""
Rather than writing our own code to solve differential equations, we'll make use of the [DifferentialEquations.jl package](http://docs.juliadiffeq.org/stable/index.html).  It's big, powerful and complex.  Here we'll use just a few of it's features.  Below I illustrate how to create an "ODEProblem" object by providing the function to calculate derivatives (`gravity_as_first_order_system`), the initial conditions (`u_init`), a "Tuple" (here used to store an ordered-pair) with the start and stop time, and a list of fixed parameters (here the masses of the bodies, `m`).  
"""

# ╔═╡ 6bd037f3-c3c9-4bf2-a7dd-525056881519
time_span_short = (0.0,4*year)       # Set how long to integrate for

# ╔═╡ e6c08b6a-9bd1-4efb-b61e-7e63d900b617
# Setup the ODE problem, but don't actually solve it yet
prob_a = ODEProblem(gravity_as_first_order_system,u_init,time_span_short,m)

# ╔═╡ 2095146e-4ec3-495d-ac57-45e800e8f2df
md"""
You can see from the output, that this returned a variable that is a [composite type](https://docs.julialang.org/en/v1/manual/types/index.html#Composite-Types-1), containing the information that we provided explicitly, as well as some additional variables that it set to default values for us.  For example, even though, we didn't explicitly specify that the `u` array would contain `Float64`'s, it figured that out and will make use of that information to compile code that is optimized for that type.

Before you execute the cells below to integrate the equations numerically, please think about what the trajectory of the planet _should_ look like if the equations were integrated perfectly.  

First, we'll integrate the differential equations using the simplest possible algorithm, [Euler's method](https://en.wikipedia.org/wiki/Euler_method).  
At each step of the integration, it sets $u(t_{i+1}) = u(t_i) + \Delta t \frac{du}{dt}(t_i,u_i)$.  

What do you expect the numerical solution will look like?
"""

# ╔═╡ 48835468-052e-47e2-abab-b65688482d0f
md"""
Ok, let's integrate this system.  To keep things simple and make issues readily apparent, we'll manually set a fixed time step of one 36th of an orbit.  We'll use the `@time` macro to measure how long it take and how much memory is allocated.  
"""

# ╔═╡ 732424a0-833b-44f2-a8cf-a964fcbec8d0
with_terminal() do 
	@time sol_tmp = solve(prob_a, Euler(), dt=year/36);
end

# ╔═╡ 6854f09c-f0b2-41e3-bc1b-e0bcbb420ad1
md"""
Wow, that took a lot of time and memory for such a short calculating.  I've been singing the praises of Julia as a very efficient language for scientific computing.  What happened here?
Let's try doing that that again (in the cell below, so you can compare).
"""

# ╔═╡ 3dafb51e-80e3-4f80-8d38-881bcccfc171
with_terminal() do 
	prob_a2 = ODEProblem(gravity_as_first_order_system,u_init,time_span_short,m)
	@time sol_tmp = solve(prob_a2, Euler(), dt=year/36);
end

# ╔═╡ d6848180-96d3-4ff8-ad27-ab9f7f8398e7
md"""
a.  How does the run time and memory compare to the first time you executed the cell?  
What's changed?  """

# ╔═╡ de5585dd-2c16-46ce-aaef-daa42fea08ea
response_2a = missing # INSERT your response as a String or Markdown 

# ╔═╡ 1660f55e-6834-48f1-a9c0-f3108faadd93
display_msg_if_fail(check_type_isa(:response_2a,response_2a,[AbstractString,Markdown.MD])) 

# ╔═╡ 0e8e5bd2-575e-4b1d-a303-076af7b231b1
md"""
Next, we'll plot the trajectory of the planet in the $x$-$y$ plane.  We could use the normal syntax to plot using Julia's [Plots.jl](http://docs.juliaplots.org/latest/) package.  However, the DifferentialEquations.jl package kindly provides a "recipie" that makes it even easier to plot the results of calling `DifferentialEquations.solve`.  We can call `plot` with just one arguments, the result of calling solve, and it would plot each of the 12 variables as a function of time.  I've added an optional argument to label the $x$-axis.
"""

# ╔═╡ 11bf2fbc-e1b1-4d91-a4d0-cf4b3ab5d34b
sol_a = solve(prob_a, Euler(), dt=year/36);

# ╔═╡ 2e3d18cb-57f1-456a-8fd0-e79bef3a4381
plot(sol_a, xlabel="Time")

# ╔═╡ 767ddda0-ffd4-4e0d-a730-fc755ecacec8
md"""
That's not the most intuitive way to visualize the results.  Below, I'll create a scatter plot and use `vars=(10,11)` to specify that we want the 10th and 11th variables (x and y position of body 2), so we can see the trajectory of the planet more clearly.  I've added some additional [optional arguements](http://docs.juliaplots.org/latest/attributes/) to make the plot a little prettier.
"""

# ╔═╡ 1d08ed26-aee4-4760-ac3b-fa1166c4189e
scatter(sol_a,vars=(10,11), markersize = 0, xlabel="x (AU)", ylabel="y (AU)",
        title="Trajectory (4 years, Δt=1/36yr)", legend = false)

# ╔═╡ 647d70ed-7e81-4b63-8428-99b5c8545325
md"""
b.  Is that what you expected?  What's wrong with this trajectory?  What could we do to numerically solve the system more accurately?  
"""

# ╔═╡ 0d77e7f0-3e67-4beb-a014-b4953a7feb95
response_2b = missing; # md"INSERT RESPONSE"

# ╔═╡ f1db3d6b-7216-412b-a596-d3492a3fdac6
display_msg_if_fail(check_type_isa(:response_2b,response_2b,[AbstractString,Markdown.MD])) 

# ╔═╡ 2f12108d-2dbd-48a7-828a-5382df5b0579
md"""
## Reducing the time-step
The Euler method has one parameter (`dt`) the ammount of simulation time that the system will be advanced in each step of the integrator.  Try adjusting the timestep in the cell below and observe how the trajectory changes.  
"""

# ╔═╡ 41a23efe-6827-4044-8d18-6759030105f8
begin
	prob_c = ODEProblem(gravity_as_first_order_system,u_init,time_span_short,m)
	sol_c = solve(prob_c, Euler(), dt=year/72)
end;

# ╔═╡ 33931a2d-a393-42dd-822f-cfb5af743de0
scatter(sol_c,vars=(10,11),markersize = 0, xlabel="x (AU)", ylabel="y (AU)", 
	        title="Trajectory (4 years, Δt=1/72yr)", legend = false)

# ╔═╡ e88f5024-ebb0-436e-a9e0-682e5402f9a6
md"""
c.  Would reducing `dt` be a practical way to get accurate results with the Euler integrator?
"""

# ╔═╡ c13d6d51-421b-42f7-94cb-f3691fd319b1
response_2c = missing; # md"INSERT RESPONSE"

# ╔═╡ 70371442-6d13-4f9d-b609-0934151dc83b
display_msg_if_fail(check_type_isa(:response_2c,response_2c,[AbstractString,Markdown.MD])) 

# ╔═╡ dfd79c5e-314a-4feb-a6d2-9b27bd6d924c
md"""
## Increasing the order of the integrator

Euler's method is first-order.  We could try to improve on our result by using a second-order integration, variously known as modified Euler's method or Heun's method.  It first computes $\tilde{u}(t_{i+1}) = u(t_i) + \Delta t \frac{du}{dt}(t_i,u_i)$ and then improves the estimate with
$$u(t_{i+1}) = u(t_i) + \frac{\Delta t}{2} \left[ \frac{du}{dt}(t_i,u_i) + \frac{d\tilde{u}}{dt}(t_{i+1},u_i) \right].$$   
d.  What do expect the trajectory integrated with the improved Euler's method will look like?
"""

# ╔═╡ 0f284d05-d5e0-4490-9f99-e9de1163566f
response_2d = missing; # md"INSERT RESPONSE"

# ╔═╡ ba74243b-4cee-4b4e-bf21-9ece6923a3ef
display_msg_if_fail(check_type_isa(:response_2d,response_2d,[AbstractString,Markdown.MD])) 

# ╔═╡ 27352445-00c7-4c3f-9f41-b6127bfb1de2
begin
	prob_d = ODEProblem(gravity_as_first_order_system,u_init,time_span_short,m)
	sol_d = solve(prob_d,Heun(),dt=year/36)
end;

# ╔═╡ 2525e649-cf6b-45dd-9532-4fec2aba2a38
scatter(sol_d,vars=(10,11),markersize = 0, xlabel="x (AU)", ylabel="y (AU)", title="Trajectory (4 years, Δt=1/36yr, 2nd order)", legend = false)

# ╔═╡ ce8aab4e-dcfc-43b0-b129-9335cfd1e5f5
md"""
e.  Did the results match your expectation?  If not, explain.  
How is the result better?  Based on your results, would you be comfortable using the algorithm above?   Why or why not?
"""

# ╔═╡ 54f4c47c-4cf3-4a9e-8f32-2f7bb0c89881
response_2e = missing; # md"INSERT RESPONSE"

# ╔═╡ 0c79007e-9721-45eb-9c2d-5a603dd7efaa
display_msg_if_fail(check_type_isa(:response_2e,response_2e,[AbstractString,Markdown.MD])) 

# ╔═╡ 55e5ce08-e95f-40f6-bc01-9dc712aab257
md"""
Now, let's try using the same algorithm, but for a longer integration.  Instead of a simulation of 4 orbits, we'll try integrating for a thousand orbits.  What do you expect the trajectory will look like? 
"""

# ╔═╡ 83dd34a6-3ced-4d25-a125-2848fe384d47
begin
	time_span_long = (0.0,100*year) 
	prob_e = ODEProblem(gravity_as_first_order_system,u_init,time_span_long,m)
	sol_e = solve(prob_e,Heun(),dt=year/36,saveat=year)
end;

# ╔═╡ 43f85b14-a159-4946-a27c-22da4d50c846
scatter(sol_e,vars=(10,11),markersize = 0, xlabel="x (AU)", ylabel="y (AU)", title="Trajectory  (1000 years, Δt=1/36yr, 2nd order)", legend = false)

# ╔═╡ 588f0dfc-ee79-4617-ac24-62cae9e682b7
md"""
f.  Did the results match your expectation?  If not, explain.  
Based on these results, would you be comfortable using the algorithm above?   Why or why not?  How could we test whether an algorithm is providing acceptable results?
"""

# ╔═╡ a44cbf9f-ad8c-4943-a65f-f82c2d61093e
response_2f = missing; # md"INSERT RESPONSE"

# ╔═╡ c518c2a0-c053-4451-aec1-8e89128600f1
display_msg_if_fail(check_type_isa(:response_2f,response_2f,[AbstractString,Markdown.MD])) 

# ╔═╡ 90175c63-d4a9-4f00-943f-df63f9d37575
md"""
## Quantifying  accuracey
In order to provide a more quantiative assessment of the integration accuracey, we can make use of our knowledge about the physics of this system.  Conservation laws dictate that the energy and angluar momentum should be conserved.  For a circular orbit, the star-planet separation should remain constant, and the orbital phase should increase linearly.  Below, I've written some helper functions that will allow you to more easily visualize how well each integration performs.

"""

# ╔═╡ cbf61a08-2780-47d2-9200-37d84f79bc6c
"Calculate the phase of the orbit of the pl_id'th planet relative to x-axis (internal)"
function calc_angle(u::Vector; pl_id::Integer = 1  )
    @assert length(u) >= 6  
    @assert 1 <= pl_id <= length(u)//6 
    dx = u[3*pl_id+1]-u[1]
    dy = u[3*pl_id+2]-u[2]
    atan(dy,dx)
end

# ╔═╡ 77e9fb29-eb09-40a1-a2c7-de90a546031d
"""
   calc_phase_error(sol; phase_init = 0.0, year = 2pi)

Calculate a vector of the difference in orbital phase of the 
second body at each time relative to a linearly increasing phase.

Inputs:
   - sol is the result of DifferentialEquations.solve
   - year: time required for the phase to increase by 2pi
Outputs:
   - Vector of differences in phase
Assumptions:
   - sol.u contains an array of positions [x_1, y_1, z_1, x_2, y_2, z_2, ...]
   - sol.t contains an array of simulation times
"""
function calc_phase_error(sol::ODESolution; year = 2pi, pl_id::Integer = 1)
	@assert 1 <= pl_id <= length(sol.u)//6 
    dtheta_dt = 2pi/year  
    phase_init = calc_angle(sol.u[1], pl_id=pl_id)  # Phase at first stored time
    delta = mod.(calc_angle.(sol.u,pl_id=pl_id).-(sol.t.-sol.t[1]).*dtheta_dt.-phase_init,2pi)
    # mod returns returns values in [0,2pi), but we want between [-pi,pi)
    delta[delta .< -pi] .+= 2pi
    delta[delta .>  pi] .-= 2pi
    delta
end

# ╔═╡ 494c1b1f-6f18-4928-a721-9956904ed260
"Calculate the separation between first and second bodies"
function calc_separation(u::Vector; pl_id::Integer = 1 )
    @assert length(u) >= 6
	@assert 1 <= pl_id <= length(u)//6 
    dx = u[3*pl_id+1]-u[1]
    dy = u[3*pl_id+2]-u[2]    
    dz = u[3*pl_id+3]-u[3]
    sqrt(dx^2+dy^2+dz^2)
end

# ╔═╡ 32dae7aa-e65a-4add-8b81-04c5eb30499a
"""
   calc_energy(u::Vector, mass::Vector)

Calculate the energy of the system
Inputs:
- u: array of 3*N positions followed by 3*N velocities
  [x_1, y_1, z_1, x_2, y_2, z_2, ...., vx_1, vy_1, vz_1, ..., vz_N]
- mass: array of N masses
Assumes:
- Gravitational constants = 1
"""
function calc_energy(u::Vector, m::Vector)
    @assert length(u) == 6*length(m)
    @assert length(m) == 2  
    # Assumes 2 bodies each with 3 coordinates, so velocities begin at 7
    kinetic = 0.5*m[1]*sum(u[7:9].^2) + m[2]*sum(u[10:12].^2) 
    d = calc_separation(u)
    potential = -m[1]*m[2]/d
    kinetic + potential
end

# ╔═╡ 5bbe0d8f-62fc-4375-8111-f8774f0b2da8
"""
   calc_angular_momentum(u::Vector,m::Vector)

Calculate the angular momentum of the system"
Inputs:
- u: array of 3*N positions followed by 3*N velocities
  [x_1, y_1, z_1, x_2, y_2, z_2, ...., vx_1, vy_1, vz_1, ..., vz_N]
- mass: array of N masses
Assumes:
- Gravitational constants = 1
"""
function calc_angular_momentum(u::Vector,m::Vector)
    @assert length(u) == 6*length(m)
    N = length(m)
    @assert length(m) == 2  
    # Reshape u into a 2-d array that's 3x4
    # u[:,i] is the ith 3-vector. 
    # positions in i=[1,2] velocities in i[3,4]
    L  = m[1] * cross(reshape(u,3,4)[:,1], reshape(u,3,4)[:,3])
    L += m[2] * cross(reshape(u,3,4)[:,2], reshape(u,3,4)[:,4])
    L[3]
    #L

end

# ╔═╡ ffd49b82-5bc5-4f01-bff6-12102282ee51
md"""
The first few times we did an integration and made plots, we did it one step at a time.  Going forward, you're going to try integrating this system many times changing the algorithm and a few parameters.  In order to make that more efficient (both for the computer and for you), we'll package all those stepts into a function, `make_test_plots_v1`.
"""

# ╔═╡ 4ab79682-ac15-4b51-a630-464bdf5f9108
md"""
Now, we can easily integrate systems and inspect the results using just one line of code for each algorithm that we test.  For example:
"""

# ╔═╡ bb4b371f-c5f4-4d56-9ffb-71cafc3a2ef7
md"""
g.  Based on these plots, are you happy with the accuracy of the integration?  
"""

# ╔═╡ 20a77293-c0c1-445a-a9a5-239028238a2d
response_2g =  missing; # md"INSERT RESPONSE"

# ╔═╡ a5507391-1d49-4709-bc9d-bcf79f2c1ca3
display_msg_if_fail(check_type_isa(:response_2g,response_2g,[AbstractString,Markdown.MD])) 

# ╔═╡ 090cfcd4-ee1a-4707-a523-26efa2bc2e9c
md"""
Try experimenting with alternative integration algorithms by replacing `alg=Heun()` in the cell above with [other integration algorithms](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html#Full-List-of-Methods-1) such as `Midpoint()`, `RK4()`, `Tsit5()`, `DP8()`.  (For higher-order algorithms, you may want to reduce the number of steps per orbits.)
"""

# ╔═╡ 74f1ea84-e247-4d25-a232-fb486589f604
md"""
h.  Did any of the integrators you tried perform acceptably?  If so, which?  
"""

# ╔═╡ ab3da7b3-c190-490c-aba9-54aaadcda44a
response_2h =  missing; # md"INSERT RESPONSE"

# ╔═╡ af16d7e3-1b71-4c7d-a126-48c32122fa7e
display_msg_if_fail(check_type_isa(:response_2h,response_2h,[AbstractString,Markdown.MD])) 

# ╔═╡ bbeb1677-2f15-46fd-9c0f-7d399b610f5b
md"""
## Choosing Appropriate Algorithms 
Next, we will rewrite the problem in a slightly different way that allows the DifferentialEquations.jl package to make use of the nice mathematical properties of a Hamiltonian system.  There are special mathematical properties of the N-body problem.
There are [specialized integration algorithms](http://docs.juliadiffeq.org/stable/solvers/dynamical_solve.html#Specialized-OrdinaryDiffEq.jl-Integrators-1) that can be applied when the derivative of the positions is proportional to the velocities and the derivative of the velocities does not depend on the velocities.  To make use of these algorithms, I've provided new functions that calculate the derivatives of the positions ("drift") separately from calculating the derivatives of the velocities ("kick").   
"""

# ╔═╡ 0c4b3da0-af77-42b7-acea-a8613f907235
"""
   gravity_drift(du, v, u, mass, t)

Sets derivative of positions equal to velocities for a system of N massive particles.

# Inputs:  
- v: array of 3*N velocities
- u: array of 3*N positions 
  [x_1, y_1, z_1, x_2, y_2, z_2, ...., z_N]
- mass: array of N masses
- t: time

# Outputs:
- du: Aray of 3*N position derivatives

"""
function gravity_drift(du, v, u, mass, t)
    @assert length(du) == length(v) == length(u) == 3 * length(mass)
    N = length(mass)
    for k in 1:(3*N) 
       du[k] = v[k]     # derivative of positions is simply velocity
    end
end

# ╔═╡ 29ce2040-35c1-474b-96b9-d2d3455beee7
"""
   gravity_kick(dv, v, u, mass, t)

Sets derivative of velocities equal to acceleration due to Newtonian gravity for a system of N massive particles.

# Inputs:  
- v: array of 3*N velocities
- u: array of 3*N positions 
  [x_1, y_1, z_1, x_2, y_2, z_2, ...., z_N]
- mass: array of N masses
- t: time

# Outputs:
- dv: Aray of 3*N velocity derivatives

# Assumptions:
- Gravitational constant = 1.  
Therefore, if masses are in solar masses and distances are in AU, 
then the time unit is such that one year = 2pi.

"""
function gravity_kick(dv, v, u, mass, t)
    @assert length(dv) == length(v) == length(u) == 3 * length(mass)
    N = length(mass)
    @assert 2 <= N      # require at least two massive bodies
    dv .= 0.0             # initialize accelerations to zero
    for i in 0:(N-2)      # loop over each pairs of bodies
        for j in (i+1):(N-1)
            dx = u[3*i+1] - u[3*j+1]  # displacements
            dy = u[3*i+2] - u[3*j+2]
            dz = u[3*i+3] - u[3*j+3]
            # calculate distance once per pair
            d = sqrt(dx^2+dy^2+dz^2)  
            # derivatives of velocities are accelerations 
            # acceleration on body i due to body j
            dv[3*i+1] -= mass[j+1] * dx / d^3 
            dv[3*i+2] -= mass[j+1] * dy / d^3 
            dv[3*i+3] -= mass[j+1] * dz / d^3 
            # acceleration on body j due to body i
            dv[3*j+1] += mass[i+1] * dx / d^3 
            dv[3*j+2] += mass[i+1] * dy / d^3 
            dv[3*j+3] += mass[i+1] * dz / d^3 
        end
    end
    dv
end

# ╔═╡ 4f623022-2d89-46c5-901a-9efdcfaed825
md"""
It turns out that DifferentialEquations stores the results in a slightly different format, so we need to make new versions of the function `calc_phase_error`, `calc_energy`, and `calc_angular_momentum`.  
"""

# ╔═╡ ac70198f-f581-4068-89f2-079b11dceead
"""
   calc_phase_error_v2(sol; phase_init = 0.0, year = 2pi)

Calculate a vector of the difference in orbital phase of the 
second body at each time relative to a linearly increasing phase.

Inputs:
   - sol is the result of DifferentialEquations.solve
   - year: time required for the phase to increase by 2pi
Outputs:
   - Vector of differences in phase
Assumptions:
   - Each sol.u[i].x[1] contains an array of positions at time i [x_1, y_1, z_1, x_2, y_2, z_2, ...]
   - sol.t contains an array of simulation times
"""
function calc_phase_error_v2(sol::ODESolution; year = 2pi)
    dtheta_dt = 2pi/year  
    phase_init = calc_angle(sol.u[1].x[1])
    delta = mod.(map(i->calc_angle(sol.u[i].x[1]),1:length(sol.u)).-(sol.t.-sol.t[1]).*dtheta_dt.-phase_init,2pi)
    # mod returns returns values in [0,2pi), but we want between [-pi,pi)
    delta[delta .< -pi] .+= 2pi
    delta[delta .>  pi] .-= 2pi
    delta
end

# ╔═╡ ba8f50a7-8eff-473d-960e-b5255fc08bb4
"Calculate the energy of the system"
function calc_energy(u::AbstractVector, v::AbstractVector, m::AbstractVector)
    @assert length(u) == length(v) == 3*length(m)
    @assert length(m) == 2 
    kinetic = 0.5*m[1]*sum(v[1:3].^2) + m[2]*sum(v[4:6].^2) 
    d = calc_separation(u)
    potential = -m[1]*m[2]/d
    kinetic + potential
end

# ╔═╡ 298476f3-66b5-460e-858e-1d6ecec77dec
"Calculate the angular momentum of the system"
function calc_angular_momentum(u::AbstractVector, v::AbstractVector, m::AbstractVector)
    @assert length(u) == length(v) == 3*length(m)
    N = length(m)
    @assert length(m) == 2 
    # cross works on 3-vectors, so use view(v,range) to provide view that looks like a 3-vector
    L  = m[1] * cross(view(u,1:3), view(v,1:3))
    L += m[2] * cross(view(u,4:6), view(v,4:6))
    L[3]
end

# ╔═╡ dc413bad-9d28-487f-855a-e68295abee75
"""
   make_test_plots_v1(opts)

Integrate a two-body system and plot the change in energy, angular modmentum, 
radial separation of the planet and deviation of its phase from linear growth.

# Optional names arguements: (default value)
- alg: algorithm to use (DP8(); see http://docs.juliadiffeq.org/stable/solvers/ode_solve.html#Full-List-of-Methods-1)
- duration: duration to integrate in orbits (100)
- steps_per_orbit: number of time steps per orbit for Euler algorithm (ignored by other integrators) (1000)
- save_every_n_orbits: how often to store results for plotting (1)
- init_separation: initial separation (1)
- mass: masses of bodies ([1, 0.001])

"""
function make_test_plots_v1(; alg=DP8(), duration=100, steps_per_orbit=1000, 
        maxiters=1_000_000, save_every_n_orbits=1, init_separation= 1, mass = [1.0, 0.001], plt_title = "")
    @assert length(mass) == 2  # Phase and separation error only make sense for 1 planet
    # Setup initial conditions
    init_velocity = sqrt(sum(mass)/init_separation) # Uniform circular motion
    year = 2pi*sqrt(init_separation^3/sum(mass))    # Kepler's third law
    r_init = init_separation .* [-mass[2]/sum(mass), 0, 0, mass[1]/sum(mass), 0, 0 ]
    v_init = init_velocity   .* [0, -mass[2]/sum(mass), 0, 0, mass[1]/sum(mass), 0 ]
    # DifferentialEquations.jl wants the initial conditions as a single vector
    u_init = vcat(r_init, v_init);  # concatenate positions & velocities
    @assert length(u_init) == 6 * length(mass)
    year = 2pi*sqrt(init_separation^3/sum(mass))    # Kepler's third law
    time_span = (0.0,0.01*year) 
    prob = ODEProblem(gravity_as_first_order_system,u_init,time_span,mass)
    if alg==Euler() # Euler requires a specified time step dt
        # First do a very short integration to make sure code is compiled before timing
        sol = solve(prob,alg,dt=year/steps_per_orbit,saveat=year*save_every_n_orbits,maxiters=maxiters);
        time_span = (0.0,duration*year) 
        prob = ODEProblem(gravity_as_first_order_system,u_init,time_span,mass)
        # Now do the requested integration and time how long it takes
        @time sol = solve(prob,alg,dt=year/steps_per_orbit,saveat=year*save_every_n_orbits,maxiters=maxiters, force_dtmin=false);
    else # Other algorithms heuristically pick a timestep
        # First do a very short integration to make sure code is compiled before timing
        sol = solve(prob,alg,saveat=year*save_every_n_orbits);
        time_span = (0.0,duration*year) 
        prob = ODEProblem(gravity_as_first_order_system,u_init,time_span,mass)
        # Now do the requested integration and time how long it takes
        @time sol = solve(prob,alg,saveat=year*save_every_n_orbits,maxiters=maxiters);
    end
    separation_init = calc_separation(u_init)
    Lz_init = calc_angular_momentum(u_init,mass)
    E_init = calc_energy(u_init,mass)
    # Make plots
    plot_angle    = scatter(sol.t,calc_phase_error(sol,year=year), xlabel = "Time", ylabel = "Phase error", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2,sigdigits=1):duration*year)
    plot_distance = scatter(sol.t,calc_separation.(sol.u).-separation_init, xlabel = "Time", ylabel = "Separation", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2,sigdigits=1):duration*year)
    plot_energy   = scatter(sol.t,map(x->calc_energy(x,mass).-E_init,sol.u), xlabel = "Time", ylabel = "Energy Error", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2,sigdigits=1):duration*year)
    plot_Lz       = scatter(sol.t,map(x->calc_angular_momentum(x,mass).-Lz_init,sol.u), xlabel = "Time", ylabel = "L_z Error", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2,sigdigits=1):duration*year, yformatter=((x)->Printf.@sprintf "%0.1e" x))
    if length(plt_title) >= 1 
		plt = plot( plot_energy, plot_Lz, plot_distance, plot_angle, layout = (2,2), plot_title=plt_title )
	else
		plt = plot( plot_energy, plot_Lz, plot_distance, plot_angle, layout = (2,2) )
	end
	
	plt
end

# ╔═╡ f14d96ef-d885-403a-8f2c-4023abb8db71
make_test_plots_v1(alg=Euler(), duration=100, plt_title="Euler, 1st order")

# ╔═╡ 590ad389-9037-4a23-ab35-5bea96e32b11
make_test_plots_v1(alg=Heun(), steps_per_orbit=1000, duration=1000, plt_title="Heun, 2nd order")

# ╔═╡ cbf0d4f7-a10c-4499-9bd5-052a21e3c984
make_test_plots_v1(alg=Midpoint(), steps_per_orbit=1000, duration=1000, plt_title="Midpoint, 2nd order")

# ╔═╡ c49f5f6d-1ceb-40a3-8032-a74b465936a7
make_test_plots_v1(alg=RK4(), steps_per_orbit=10, duration=1000, plt_title="Runge-Kutta 4th order")

# ╔═╡ b4d95d86-254c-4982-bd15-5ff263ff750d
make_test_plots_v1(alg=Tsit5(), steps_per_orbit=1000, duration=1000,plt_title="Tsit5")

# ╔═╡ 8f68ba9a-ab88-4472-b19b-9f30b3f562e6
make_test_plots_v1(alg=DP8(), steps_per_orbit=36, duration=1000, plt_title="Dormand-Prince, ~8th order")

# ╔═╡ efcc2856-fb3e-45c1-9f70-8e6f4cb979ec
make_test_plots_v1(alg=Vern8(), steps_per_orbit=1000, duration=1000, plt_title="Verner's Runge-Kutta,  ~8th order")

# ╔═╡ 8e59c119-1850-424f-9a78-1d6852054c26
md"""
Now, we'll make a new function `make_test_plots_v2` that allows us to perform similar tests, but using the symplectic integration algorithms that are applicable to our problem.
"""

# ╔═╡ 31c7f905-7291-4f9f-9cc3-2ab5d155f538
"""
   make_test_plots_v2(opts)

Integrate a two-body system and plot the change in energy, L_z, radial separation of the planet and deviation of phase from linear growth.
Setup for algorithms that uses knowledge that this is a Hamiltonian system

# Optional names arguements: (default value)
- alg: algorithm to use (KahanLi6(); see http://docs.juliadiffeq.org/stable/solvers/dynamical_solve.html#Symplectic-Integrators-1)
- duration: duration to integrate in orbits (100)
- steps_per_orbit: number of time steps per orbit (36)
- save_every_n_orbits: how often to store results for plotting (1)
- init_separation: initial separation (1)
- mass: masses of bodies ([1, 0.001])            
"""
function make_test_plots_v2(; duration=100, alg=KahanLi6(), steps_per_orbit=36, 
        save_every_n_orbits=1, init_separation = 1, mass=[1.0,0.001], plt_title="")
    @assert length(mass) == 2  # Phase and separation error only make sense for 1 planet
    # Setup initial conditions
    init_velocity = sqrt(sum(mass)/init_separation) # Uniform circular motion
    year = 2pi*sqrt(init_separation^3/sum(mass))    # Kepler's third law
    r_init = init_separation .* [-mass[2]/sum(mass), 0, 0, mass[1]/sum(mass), 0, 0 ]
    v_init = init_velocity   .* [0, -mass[2]/sum(mass), 0, 0, mass[1]/sum(mass), 0 ]
    @assert length(r_init) == length(v_init) == 3 * length(mass)
    # First do a very short integration to make sure code is compiled before timing
    time_span = (0.0,0.01*year) 
    prob = DynamicalODEProblem(gravity_kick,gravity_drift,v_init,r_init,time_span,m)
    sol = solve(prob,alg,dt=year/steps_per_orbit,saveat=year*save_every_n_orbits);
    # Now do the requested integration and time how long it takes
    time_span = (0.0,duration*year) 
    prob = DynamicalODEProblem(gravity_kick,gravity_drift,v_init,r_init,time_span,m)
    walltime = @elapsed sol = solve(prob,alg,dt=year/steps_per_orbit,saveat=year*save_every_n_orbits);
    # Make plots
    angles = calc_phase_error_v2(sol,year=year)
    plot_angle     = scatter(sol.t,angles, xlabel = "Time", ylabel = "Phase error", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2, sigdigits=1):duration*year, yformatter=((x)->Printf.@sprintf "%0.1e" x))
    separation_init = calc_separation(r_init)
    separations = map(i->calc_separation(sol.u[i].x[2]),1:length(sol.u)) .- separation_init

    plot_separations = scatter(sol.t,separations, xlabel = "Time", ylabel = "Separation", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2, sigdigits=1):duration*year, yformatter=((x)->Printf.@sprintf "%0.1e" x))
	annotate!(duration*year/2,0,(Printf.@sprintf "%.2g sec" walltime))
    E_init = calc_energy(r_init,v_init,mass)
    energies = map(i->calc_energy(sol.u[i].x[2],sol.u[i].x[1],mass),1:length(sol.u)) .- E_init
    plot_energy   = scatter(sol.t,energies, xlabel = "Time", ylabel = "Energy Error", markersize=0, legend = false,
      grid=:no, xticks=0:round(duration*year/2, sigdigits=1):duration*year, yformatter=((x)->Printf.@sprintf "%0.1e" x ))
	annotate!(duration*year/2,0,(Printf.@sprintf "ΔE=%0.1e" energies[end]))
    Lz_init = calc_angular_momentum(r_init,v_init,mass)
    Lzs =  map(i->calc_angular_momentum(sol.u[i].x[2],sol.u[i].x[1],mass),1:length(sol.u)) .- Lz_init
    Lz_range = minimum(Lzs):maximum(Lzs)
    plot_Lz       = scatter(sol.t,Lzs, xlabel = "Time", ylabel = "L_z Error", markersize=0, legend = false, grid=:no, xticks=0:round(duration*year/2, sigdigits=1):duration*year, yformatter=((x)->Printf.@sprintf "%0.1e" x))
	annotate!(duration*year/2,0,(Printf.@sprintf "ΔLz=%0.1e" Lzs[end]))
	alg_str = string(alg)
	println("# Algorithm: ", alg_str)
	println("# Runtime ", walltime, " seconds.")
    println("# Final errors: E: ", energies[end], " Lz: ", Lzs[end], " r: ",separations[end], " θ: ", angles[end])
    plot( plot_energy, plot_Lz, plot_separations, plot_angle, layout = (2,2), fmt = :png )
	if length(plt_title) >= 1 
		plt = plot( plot_energy, plot_Lz, plot_separations, plot_angle, layout = (2,2), plot_title=plt_title, fmt = :png )
	else
		plt = plot( plot_energy, plot_Lz, plot_separations, plot_angle, layout = (2,2), fmt = :png )
	end

end


# ╔═╡ 7e473a6a-e2ac-4c67-894a-6965deac31db
md"""
First try a lower-order algorithm such as [Leapfrog integration](https://en.wikipedia.org/wiki/Leapfrog_integration).  The `VerletLeapfrog()` integrator is second-order, but symplectic and time-reversible.  How do you expect the results using this integrator to differ from the results with previous integrators?
"""

# ╔═╡ 830717e6-9064-4485-80aa-abbae0f2ae58
make_test_plots_v2(alg=VerletLeapfrog(),steps_per_orbit=1000,duration=1000, plt_title="Verlet Leapfrog, 2nd order")

# ╔═╡ 06638ff4-3960-4783-99e8-33494de25de1
md"""
i.  How does the accuracy of the results with the Verlet Leapfrog integrator compare to previous results?  Based on these results, would you be comfortable using the algorithm above?   Why or why not?
"""

# ╔═╡ b463565b-a8b0-4fc9-83e0-ae688dc1643e
response_2i =  missing; # md"INSERT RESPONSE"

# ╔═╡ ccd4bd59-c30f-488c-a5b5-f88aa61ef849
display_msg_if_fail(check_type_isa(:response_2i,response_2i,[AbstractString,Markdown.MD])) 

# ╔═╡ acaffdff-67a1-43f5-ada6-b3d3a9cb0875
md"""
## Choosing efficient algorithms

j.  Try several of the [symplectic integrators](https://diffeq.sciml.ai/stable/solvers/dynamical_solve/#Symplectic-Integrators), such as `McAte3()`, `CalvoSanz4()`, `McAte5()`, `KahanLi6()`, `KahanLi8()`.  The number refers to the order of the integrator.  Compare the accuracy of the results with symplectic integrators of different orders.  You may want to extend the duration of the integrations by replacing `num_orbits_symplectic` with a large value.
Also compare the wall time required (see the duration printed in the separation panel).  
"""

# ╔═╡ a287218b-6152-4af7-bce7-8a8559165ae4
response_2j =  missing; # md"INSERT RESPONSE"

# ╔═╡ 9d8908a1-b164-4b8f-a555-47314094b5b6
display_msg_if_fail(check_type_isa(:response_2j,response_2j,[AbstractString,Markdown.MD])) 

# ╔═╡ a3decfbf-befc-474b-abdb-daada98148ee
num_orbits_symplectic = 10_000 

# ╔═╡ 96fbe112-7f1b-462f-9699-abc93e9e4054
make_test_plots_v2(alg=McAte3(),steps_per_orbit=30,duration=num_orbits_symplectic, plt_title="McAte3")

# ╔═╡ 8b39624e-247f-47dc-8eb3-48ec8239b9a4
make_test_plots_v2(alg=CalvoSanz4(),steps_per_orbit=30,duration=num_orbits_symplectic, plt_title="CalvoSanz4")

# ╔═╡ e22eb346-9383-471a-97db-918d9008be7a
make_test_plots_v2(alg=McAte5(),steps_per_orbit=30,duration=num_orbits_symplectic, plt_title="McAte5")

# ╔═╡ 4d199c96-3de6-4078-951d-501a04431e80
make_test_plots_v2(alg=KahanLi6(),steps_per_orbit=30,duration=num_orbits_symplectic, plt_title="KahanLi6")

# ╔═╡ 2d1349d8-b390-40bf-bd58-6f1c72d715b0
make_test_plots_v2(alg=KahanLi8(),steps_per_orbit=30,duration=num_orbits_symplectic, plt_title="KahanLi8")

# ╔═╡ 128a01f8-e597-4e5f-9654-609865522f96
md"""
For an $n$th-order integrator the error term due to _truncation_ is generally of order $\Delta~t^{n+1}$.  That means that you could likely use a smaller number of `steps_per_orbit` with a higher-order algorithm to achieve a similar numerical precision.  
"""

# ╔═╡ 62844569-a5b2-4257-9f2b-1fe8257082dd
md"""
k.  Pick an accuracy target (e.g., $\Delta~$ or phase error at the end of the simulation).  Then tinker with the the number of steps per orbit, so that the results with different integrators achieve similar accuracy.  (No need to be super precise about this.  E.g., you might aim for the exponent to be the same, but not worry about the leading digit.)  Then compare the time required to achieve the target accuracy with different integrators.  Which algorithms and time-steps would you recommend for performing an N-body integration of the solar system?
"""

# ╔═╡ 3a55e7e6-b091-402f-a960-087a21edb6ae
response_2k =  missing; # md"INSERT RESPONSE"

# ╔═╡ 28b61d7f-53dc-4570-b76d-5a22973cb1c9
display_msg_if_fail(check_type_isa(:response_2j,response_2j,[AbstractString,Markdown.MD])) 

# ╔═╡ cb29423a-b1cf-43f2-8b4e-a6fe74ae2c71
md"""

l.  How long would it take to integrate this system for $10^8$ orbits?  
"""

# ╔═╡ 82fdb87f-6d17-407a-9fdf-4e1385a60a91
response_2l =  missing; # md"INSERT RESPONSE"

# ╔═╡ 675d7d77-8efb-4519-872f-50003c555610
display_msg_if_fail(check_type_isa(:response_2l,response_2l,[AbstractString,Markdown.MD])) 

# ╔═╡ 2204f9b6-4f84-47a1-a7db-1a4295b772d6
md"""
## Lessons learned

m.  Based on this lesson, what properties of an n-body integrator are most important for obtaining accurate scientific results when studying the long-term orbital of a planetary system?  
"""

# ╔═╡ 3f7ea976-3c75-4b78-9b90-803ae8202e89
response_2m = missing; # md"INSERT RESPONSE"

# ╔═╡ 137014c4-7da5-4bed-937f-0b33ace84be7
display_msg_if_fail(check_type_isa(:response_2m,response_2m,[AbstractString,Markdown.MD])) 

# ╔═╡ d365f5d7-79c9-4792-a1b1-41358bd8b29e
md"## Helper Code"

# ╔═╡ bd845b9c-1339-45a3-a8eb-d12cc5227bd3
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
DifferentialEquations = "~6.18.0"
LaTeXStrings = "~1.2.1"
Plots = "~1.19.3"
PlutoTeachingTools = "~0.1.2"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "a71d224f61475b93c9e196e83c17c6ac4dedacfa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.18"

[[ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "0f7998147ff3d112fad027c894b6b6bebf867154"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.3"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "d17071d7fc9a98ca2d958cd217e62a17c5eeebed"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.16.10"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "fe34902ac0c3a35d016617ab7032742865756d7d"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.7.1"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "9723e1c07c1727082e169ca50789644a552fb023"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.2.3"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f53ca8d41e4753c41cdafa6ec5f7ce914b34be54"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.13"

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

[[Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "1060c5023d2ac8210c73078cb7c0c567101d201c"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.2"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dc7dedc2c2aa9faf59a55c622760a25cbefbe941"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.31.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

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

[[DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "6eba402e968317b834c28cd47499dd1b572dd093"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.31.1"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DataStructures", "DocStringExtensions", "FastBroadcast", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "9d312bb0b7c8ace440a71c64330cf1bea0ade0c8"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.70.0"

[[DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "RecipesBase", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "0972ca167952dc426b5438fc188b846b7a66a1f3"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.16.1"

[[DiffEqFinancial]]
deps = ["DiffEqBase", "DiffEqNoiseProcess", "LinearAlgebra", "Markdown", "RandomNumbers"]
git-tree-sha1 = "db08e0def560f204167c58fd0637298e13f58f73"
uuid = "5a0ffddc-d203-54b0-88ba-2c03c0fc2e67"
version = "2.4.0"

[[DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "LightGraphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "51441835fed64ebc8ed52e43a9734bfb58d5ecd7"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "7.0.0"

[[DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[DiffEqPhysics]]
deps = ["DiffEqBase", "DiffEqCallbacks", "ForwardDiff", "LinearAlgebra", "Printf", "Random", "RecipesBase", "RecursiveArrayTools", "Reexport", "StaticArrays"]
git-tree-sha1 = "8f23c6f36f6a6eb2cbd6950e28ec7c4b99d0e4c9"
uuid = "055956cb-9e8b-5191-98cc-73ae4a59e68a"
version = "3.9.0"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "214c3fcac57755cfda163d91c58893a8723f93e9"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.0.2"

[[DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqFinancial", "DiffEqJump", "DiffEqNoiseProcess", "DiffEqPhysics", "DimensionalPlotRecipes", "LinearAlgebra", "MultiScaleArrays", "OrdinaryDiffEq", "ParameterizedFunctions", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "ececc535bd2aa55a520131d955639288704e3851"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "6.18.0"

[[DimensionalPlotRecipes]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "af883a26bbe6e3f5f778cb4e1b81578b534c32a6"
uuid = "c619ae07-58cd-5f6d-b883-8f17bd6a98f9"
version = "1.2.0"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3889f646423ce91dd1055a76317e9a1d3a23fff1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.11"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics", "Test"]
git-tree-sha1 = "6cdd99d0b7b555f96f7cb05aa82067ee79e7aef4"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.2"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "e9d82a6f35d199d3821c069932115e19ca2a2b3d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.3.19"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "ad435656c49da7615152b856c0f9abe75b0b5dc9"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.8.4"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

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

[[FastBroadcast]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "26be48918640ce002f5833e8fc537b2ba7ed0234"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.8"

[[FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8c8eac2af06ce35973c3eadb4ab3243076a408e7"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.1"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

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

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "e2af66012e08966366a43251e1fd421522908be6"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.18"

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

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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
git-tree-sha1 = "15ff9a14b9e1218958d3530cc288cf31465d9ae2"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.3.13"

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
git-tree-sha1 = "c6a1fff2fd4b1da29d3dccaffb1e1001244d844e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.12"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3395d4d4aeb3c9d31f5929d32760d8baeee88aaf"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.5.0+0"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1a8c6237e78b714e901e406c096fc8a65528af7d"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.1"

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

[[JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Pkg", "Tokenize"]
git-tree-sha1 = "03d48b801c3d497c301a137ad32be7a70f1a64cb"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "0.15.2"

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

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "41fc666d11a346e55f7fb70318e7078bfc0ae7cb"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.3"

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

[[LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoopVectorization]]
deps = ["ArrayInterface", "DocStringExtensions", "IfElse", "LinearAlgebra", "OffsetArrays", "Polyester", "Requires", "SLEEFPirates", "Static", "StrideArraysCore", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "1578cc856a165170a9572602e885c8511d81da1e"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.56"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "6a8a2a625ab0dea913aba95c11370589e0239ff0"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.6"

[[ManualMemory]]
git-tree-sha1 = "71c64ebe61a12bad0911f8fc4f91df8a448c604c"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.4"

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

[[ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "InteractiveUtils", "JuliaFormatter", "LabelledArrays", "Latexify", "Libdl", "LightGraphs", "LinearAlgebra", "MacroTools", "NaNMath", "NonlinearSolve", "RecursiveArrayTools", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SafeTestsets", "SciMLBase", "Serialization", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "b6225318d687dbde588c6d0cd0b40161826763b1"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "5.26.0"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MultiScaleArrays]]
deps = ["DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "SparseDiffTools", "Statistics", "StochasticDiffEq", "TreeViews"]
git-tree-sha1 = "258f3be6770fe77be8870727ba9803e236c685b8"
uuid = "f9640e96-87f6-5992-9c3b-0743c6a49ffa"
version = "1.8.1"

[[MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "45c9940cec79dedcdccc73cc6dd09ea8b8ab142c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.3.18"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50608f411a1e178e0129eab4110bd56efd08816f"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.0"

[[NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "ef18e47df4f3917af35be5e5d7f5d97e8a83b0ec"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.8"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "4f825c6da64aebaa22cc058ecfceed1ab9af1c7e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.3"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "3afbf5398ff3e51427ed620e5ffbe96c7fdba67c"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.4.0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "NLsolve", "Polyester", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "379ed10814cba5087978bce9df31161a496620d3"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "5.60.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport", "SciMLBase"]
git-tree-sha1 = "d290c172dae21d73ae6a19a8381abbb69ef0a624"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.10.0"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "c8abc88faa3f7a3950832ac5d6e690881590d6dc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.0"

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
git-tree-sha1 = "1bbbb5670223d48e124b388dee62477480e23234"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.19.3"

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

[[PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[Polyester]]
deps = ["ArrayInterface", "IfElse", "ManualMemory", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities", "VectorizationBase"]
git-tree-sha1 = "4b692c8ce1912bae5cd3b90ba22d1b54eb581195"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.3.7"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "56ead4aaafc41d83694e17b0dd89d3e929d01a14"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.0"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "0426474f50756b3b47b08075604a41b460c45d17"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.16.1"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization"]
git-tree-sha1 = "2e1a88c083ebe8ba69bc0b0084d4b4ba4aa35ae0"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.1.13"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "5975a4f824533fa4240f40d86f1060b9fc80d7cc"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "bfdf9532c33db35d2ce9df4828330f0e92344a52"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.25"

[[SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "f0bf114650476709dd04e690ab2e36d88368955e"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.18.2"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "d5640fc570fb1b6c54512f0bd3853866bd298b3e"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.0"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

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

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "LightGraphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "bdaa461dbfebb839b740f4582e8728f0f92ca375"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.15.0"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "508822dca004bf62e210609148511ad03ce8f1d8"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "5b2f81eeb66bcfe379947c500aae773c85c31033"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.8"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2f6792d523d7448bbe2fec99eca9218f06cc746d"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.8"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "3df66a4a9ba477bea5cb10a3ec732bb48a2fc27d"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.6.4"

[[StochasticDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "d9e996e95ad3c601c24d81245a7550cebcfedf85"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.36.0"

[[StrideArraysCore]]
deps = ["ArrayInterface", "ManualMemory", "Requires", "ThreadingUtilities", "VectorizationBase"]
git-tree-sha1 = "e1c37dd3022ba6aaf536541dd607e8d5fb534377"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.1.17"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "75412a0ce4cd7995d7445ba958dd11de03fd2ce5"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.5.3"

[[Sundials_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "013ff4504fc1d475aa80c63b455b6b3a58767db2"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.0+1"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[SymbolicUtils]]
deps = ["AbstractTrees", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs"]
git-tree-sha1 = "17fecd52e9a82ca52bfc9d28a2d31b33458a6ce0"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.13.1"

[[Symbolics]]
deps = ["ConstructionBase", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "dae26a27018d0cad7efd585a9a0012c6a2752a88"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "1.4.2"

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
git-tree-sha1 = "8ed4a3ea724dac32670b062be3ef1c1de6773ae8"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.4.4"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "03013c6ae7f1824131b2ae2fc1d49793b51e8394"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.6"

[[TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "209a8326c4f955e2442c07b56029e88bb48299c7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.12"

[[Tokenize]]
git-tree-sha1 = "eee92eda3cc8e104b7e56ff4c1fcf0d78ca37c89"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.18"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[VectorizationBase]]
deps = ["ArrayInterface", "Hwloc", "IfElse", "Libdl", "LinearAlgebra", "Static"]
git-tree-sha1 = "a4bc1b406dcab1bc482ce647e6d3d53640defee3"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.20.25"

[[VertexSafeGraphs]]
deps = ["LightGraphs"]
git-tree-sha1 = "b9b450c99a3ca1cc1c6836f560d8d887bcbe356e"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.1.2"

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

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"

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
# ╟─308b9284-956b-4a9f-93f9-8397f34bed94
# ╟─e7cdc134-e45b-4d3a-b614-7153a237659e
# ╟─83e6e5a6-1347-46a5-aa06-61535cdc1120
# ╟─90139e48-c737-4117-b8ec-730d40e3541f
# ╠═de3ffc8c-1837-4df7-a977-363d5bd20117
# ╟─05f2bcdf-9ec0-43b9-9d94-5a8592e80ce1
# ╟─2fee29a5-ecdf-4516-af26-876191b0271f
# ╠═f4470846-8d69-47e0-aa55-df3a12dedd10
# ╟─c7721744-061e-49e1-9470-0fd13bb960e9
# ╠═6bd037f3-c3c9-4bf2-a7dd-525056881519
# ╠═e6c08b6a-9bd1-4efb-b61e-7e63d900b617
# ╟─2095146e-4ec3-495d-ac57-45e800e8f2df
# ╟─48835468-052e-47e2-abab-b65688482d0f
# ╠═732424a0-833b-44f2-a8cf-a964fcbec8d0
# ╟─6854f09c-f0b2-41e3-bc1b-e0bcbb420ad1
# ╠═3dafb51e-80e3-4f80-8d38-881bcccfc171
# ╟─d6848180-96d3-4ff8-ad27-ab9f7f8398e7
# ╠═de5585dd-2c16-46ce-aaef-daa42fea08ea
# ╟─1660f55e-6834-48f1-a9c0-f3108faadd93
# ╟─0e8e5bd2-575e-4b1d-a303-076af7b231b1
# ╠═11bf2fbc-e1b1-4d91-a4d0-cf4b3ab5d34b
# ╠═2e3d18cb-57f1-456a-8fd0-e79bef3a4381
# ╟─767ddda0-ffd4-4e0d-a730-fc755ecacec8
# ╠═1d08ed26-aee4-4760-ac3b-fa1166c4189e
# ╟─647d70ed-7e81-4b63-8428-99b5c8545325
# ╠═0d77e7f0-3e67-4beb-a014-b4953a7feb95
# ╟─f1db3d6b-7216-412b-a596-d3492a3fdac6
# ╟─2f12108d-2dbd-48a7-828a-5382df5b0579
# ╠═41a23efe-6827-4044-8d18-6759030105f8
# ╠═33931a2d-a393-42dd-822f-cfb5af743de0
# ╟─e88f5024-ebb0-436e-a9e0-682e5402f9a6
# ╠═c13d6d51-421b-42f7-94cb-f3691fd319b1
# ╟─70371442-6d13-4f9d-b609-0934151dc83b
# ╟─dfd79c5e-314a-4feb-a6d2-9b27bd6d924c
# ╠═0f284d05-d5e0-4490-9f99-e9de1163566f
# ╟─ba74243b-4cee-4b4e-bf21-9ece6923a3ef
# ╠═27352445-00c7-4c3f-9f41-b6127bfb1de2
# ╟─2525e649-cf6b-45dd-9532-4fec2aba2a38
# ╟─ce8aab4e-dcfc-43b0-b129-9335cfd1e5f5
# ╠═54f4c47c-4cf3-4a9e-8f32-2f7bb0c89881
# ╟─0c79007e-9721-45eb-9c2d-5a603dd7efaa
# ╟─55e5ce08-e95f-40f6-bc01-9dc712aab257
# ╠═83dd34a6-3ced-4d25-a125-2848fe384d47
# ╟─43f85b14-a159-4946-a27c-22da4d50c846
# ╟─588f0dfc-ee79-4617-ac24-62cae9e682b7
# ╠═a44cbf9f-ad8c-4943-a65f-f82c2d61093e
# ╟─c518c2a0-c053-4451-aec1-8e89128600f1
# ╟─90175c63-d4a9-4f00-943f-df63f9d37575
# ╠═cbf61a08-2780-47d2-9200-37d84f79bc6c
# ╠═77e9fb29-eb09-40a1-a2c7-de90a546031d
# ╠═494c1b1f-6f18-4928-a721-9956904ed260
# ╠═32dae7aa-e65a-4add-8b81-04c5eb30499a
# ╠═5bbe0d8f-62fc-4375-8111-f8774f0b2da8
# ╟─ffd49b82-5bc5-4f01-bff6-12102282ee51
# ╠═dc413bad-9d28-487f-855a-e68295abee75
# ╟─4ab79682-ac15-4b51-a630-464bdf5f9108
# ╠═f14d96ef-d885-403a-8f2c-4023abb8db71
# ╠═590ad389-9037-4a23-ab35-5bea96e32b11
# ╟─bb4b371f-c5f4-4d56-9ffb-71cafc3a2ef7
# ╠═20a77293-c0c1-445a-a9a5-239028238a2d
# ╟─a5507391-1d49-4709-bc9d-bcf79f2c1ca3
# ╟─090cfcd4-ee1a-4707-a523-26efa2bc2e9c
# ╠═cbf0d4f7-a10c-4499-9bd5-052a21e3c984
# ╠═c49f5f6d-1ceb-40a3-8032-a74b465936a7
# ╠═b4d95d86-254c-4982-bd15-5ff263ff750d
# ╠═8f68ba9a-ab88-4472-b19b-9f30b3f562e6
# ╠═efcc2856-fb3e-45c1-9f70-8e6f4cb979ec
# ╟─74f1ea84-e247-4d25-a232-fb486589f604
# ╠═ab3da7b3-c190-490c-aba9-54aaadcda44a
# ╟─af16d7e3-1b71-4c7d-a126-48c32122fa7e
# ╟─bbeb1677-2f15-46fd-9c0f-7d399b610f5b
# ╠═0c4b3da0-af77-42b7-acea-a8613f907235
# ╠═29ce2040-35c1-474b-96b9-d2d3455beee7
# ╠═4f623022-2d89-46c5-901a-9efdcfaed825
# ╠═ac70198f-f581-4068-89f2-079b11dceead
# ╠═ba8f50a7-8eff-473d-960e-b5255fc08bb4
# ╠═298476f3-66b5-460e-858e-1d6ecec77dec
# ╟─8e59c119-1850-424f-9a78-1d6852054c26
# ╠═31c7f905-7291-4f9f-9cc3-2ab5d155f538
# ╟─7e473a6a-e2ac-4c67-894a-6965deac31db
# ╟─830717e6-9064-4485-80aa-abbae0f2ae58
# ╟─06638ff4-3960-4783-99e8-33494de25de1
# ╠═b463565b-a8b0-4fc9-83e0-ae688dc1643e
# ╟─ccd4bd59-c30f-488c-a5b5-f88aa61ef849
# ╠═acaffdff-67a1-43f5-ada6-b3d3a9cb0875
# ╠═a287218b-6152-4af7-bce7-8a8559165ae4
# ╟─9d8908a1-b164-4b8f-a555-47314094b5b6
# ╠═a3decfbf-befc-474b-abdb-daada98148ee
# ╠═96fbe112-7f1b-462f-9699-abc93e9e4054
# ╠═8b39624e-247f-47dc-8eb3-48ec8239b9a4
# ╠═e22eb346-9383-471a-97db-918d9008be7a
# ╠═4d199c96-3de6-4078-951d-501a04431e80
# ╠═2d1349d8-b390-40bf-bd58-6f1c72d715b0
# ╟─128a01f8-e597-4e5f-9654-609865522f96
# ╟─62844569-a5b2-4257-9f2b-1fe8257082dd
# ╠═3a55e7e6-b091-402f-a960-087a21edb6ae
# ╟─28b61d7f-53dc-4570-b76d-5a22973cb1c9
# ╟─cb29423a-b1cf-43f2-8b4e-a6fe74ae2c71
# ╠═82fdb87f-6d17-407a-9fdf-4e1385a60a91
# ╟─675d7d77-8efb-4519-872f-50003c555610
# ╟─2204f9b6-4f84-47a1-a7db-1a4295b772d6
# ╠═3f7ea976-3c75-4b78-9b90-803ae8202e89
# ╟─137014c4-7da5-4bed-937f-0b33ace84be7
# ╟─d365f5d7-79c9-4792-a1b1-41358bd8b29e
# ╟─c841608b-12c3-4220-aff7-843176ff052f
# ╠═bd845b9c-1339-45a3-a8eb-d12cc5227bd3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
