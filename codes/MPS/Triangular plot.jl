using FiniteLattices,CairoMakie,LinearAlgebra,ColorSchemes

include("supplementary/supplementary.jl")

Lx = 12
Ly = 4

Latt  = InfTria(Lx,Ly)

onsite = sin.(pi .* range(0,1,size(Latt)) .- pi/2)

fig,ax = plot_Triangular(Latt;
site = true,
sitecolor = get(colorschemes[:bwr],onsite,:extrema),
sitesize = 20 .* onsite,
sitealpha = abs.(onsite))

#= fig,ax = plot_Triangular(Latt;
site = true,
) =#

display(fig)

#save("lattices/$(Lx)x$(Ly)Triangular.pdf",fig)

