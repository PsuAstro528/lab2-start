println(basename(pwd()))
if !occursin(r"lab2-start",pwd())
   println("# Writing PlutoDeployment to avoid exporting ex3.jl everytime.")
   settings = 
   """
   [Export]
   exclude = ["ex3.jl","ex*backup*.jl"] 

   [SliderServer]
   exclude = ["ex3.jl","ex*backup*.jl"]
   """
   write("../PlutoDeployment.toml",settings)   
end

