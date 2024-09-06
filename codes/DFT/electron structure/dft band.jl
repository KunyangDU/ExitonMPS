using CairoMakie,ColorSchemes
include("../src/analysis.jl")

bandname = "data/BAND_soc_afm.dat"
labelname = "data/KLABELS"

Nk,Nband = importBandSize(bandname)
dftband = importDFTBand(bandname)
klabel = importLabel(labelname)

fig = Figure()
ax = Axis(fig[1,1],
xticks = map(x -> klabel[1:4,x],(2,1)),
width = 300,height = 250,
xlabel = L"k\ /\ (1/\mathrm{A_0})",
ylabel = L"E(k)\ /\ \mathrm{eV}",
title = "DFT band")

band!(ax,dftband;color = :black)

xlims!(ax,klabel[[1,4],2]...)
ylims!(ax,-3,2)
resize_to_layout!(fig)
display(fig)
save("data/DFT band.png",fig)

