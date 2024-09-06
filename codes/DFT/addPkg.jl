using Pkg

Pkg.activate(".")


Pkg.add("LinearAlgebra")
Pkg.add("ColorSchemes")
Pkg.add("CairoMakie")
Pkg.add("FiniteDifferences")
Pkg.add("LaTeXStrings")
Pkg.add("Interpolations")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add(url = "https://github.com/Qiaoyi-Li/FiniteLattices.jl.git")
Pkg.resolve()
Pkg.gc()