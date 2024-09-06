using CSV,DataFrames,CairoMakie,JLD2
include("../../analysis/analysis.jl")

#Latt = XCYYKag(4,4;shift = Tuple(Tuple.(eachcol(unitRbasis2*Nbpos2))))
posname = "data/POSCAR"
unitRbasis2,LattConst = RBasis2(posname)
Nbpos2 = NbPosition2(posname)

Latt = XCmultiYYKag(6,6;norbs = 1,
atomshift = Tuple(Tuple.(eachcol(unitRbasis2*Nbpos2))))
fig = Figure()
ax = Axis(fig[1,1],
width = 500,height = 500,
autolimitaspect = 1)

XCYYKag!(ax,Latt,site = true,sitelabel =false,
sitesize = 10*ones(size(Latt)))
resize_to_layout!(fig)
display(fig)

save("figures/lattice $(get_cellsize(Latt)).png",fig)

