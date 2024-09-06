using CairoMakie,ColorSchemes,LinearAlgebra
include("../src/analysis.jl")

wnband = CSV.read("data/wannier TB/band.dat", DataFrame)

klabel = importLabel("data/KLABELS")

fig = Figure()
ax = Axis(fig[1,1],
xticks = map(x -> klabel[1:4,x],(2,1)),
width = 300,height = 250,
xlabel = L"k\ /\ (1/\mathrm{A_0})",
ylabel = L"E(k)\ /\ \mathrm{eV}",
title = "Wannier wnband")

band!(ax,wnband;color = :black)

xlims!(ax,extrema(wnband."kr" )...)

resize_to_layout!(fig)
display(fig)

save("data/Wannier wnband.png",fig)
