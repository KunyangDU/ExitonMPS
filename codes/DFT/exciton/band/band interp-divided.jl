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


level = 1

for ii in 2:size(dftband)[2]
    lines!(ax,dftband[:,1],dftband[:,ii],color = :black)
end

energydf = CSV.read("data/interp energy-$(level).dat", DataFrame)
fatbanddf = CSV.read("data/interp fatband-$(level).dat", DataFrame)

for jj in 2:size(energydf,2)
    scatter!(ax,energydf[:,"kr"],energydf[:,jj],markersize = 2 .* fatbanddf[:,jj],color = :red)
end

xlims!(ax,klabel[[1,4],2]...)
ylims!(ax,-2,2)
resize_to_layout!(fig)
display(fig)
save("data/interp fatband-$(level).png",fig)

