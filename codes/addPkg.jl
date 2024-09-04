using Pkg

Pkg.activate(".")

Pkg.add("FiniteMPS")
#Pkg.add("FiniteLattices")
#Pkg.add(url = "git@github.com:Qiaoyi-Li/FiniteLattices.jl.git")
Pkg.add(url = "https://github.com/Qiaoyi-Li/FiniteLattices.jl.git")
Pkg.add("MKL")
Pkg.add("JLD2")
Pkg.add("CairoMakie")
Pkg.resolve()
Pkg.gc()