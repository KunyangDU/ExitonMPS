using CSV,DataFrames,CairoMakie,JLD2
include("../../analysis/analysis.jl")

#Latt = XCYYKag(4,4;shift = Tuple(Tuple.(eachcol(unitRbasis2*Nbpos2))))
posname = "data/POSCAR"
unitRbasis2,LattConst = RBasis2(posname)
Nbpos2 = NbPosition2(posname)

Latt = XCmultiYYKag(6,6;norbs = 3,
atomshift = Tuple(Tuple.(eachcol(unitRbasis2*Nbpos2))))

p₀=1
S²,Sz = 1,1
RspaceCC = CSV.read("data/RCC/Rspace CpCoeff-$(p₀)-$(S²) $(Sz).dat",DataFrame)[:,1]

fig = Figure()
ax = Axis(fig[1,1],
width = 400,height = 400,
autolimitaspect = 1,
title = "p₀=$(p₀), S²=$(S²), Sz=$(Sz)")

XCYYKag!(ax,Latt)
XCYYKag!(ax,Latt;bond=false,
site = true,sitelabel = false,
sitesize = 150 * RspaceCC,sitecolor = repeat(:red,size(Latt)))

R₀ = collect(FiniteLattices.coordinate(Latt,get_index(Latt,p₀,map(x -> div(x,2),get_cellsize(Latt)))))

scatter!(ax,R₀...,marker = :star5,color = :blue,markersize = 15)
resize_to_layout!(fig)
display(fig)

save("figures/Rspace merge-$(p₀)-$(S²) $(Sz).png",fig)


