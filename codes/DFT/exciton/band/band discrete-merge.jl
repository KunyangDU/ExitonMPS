using CairoMakie,ColorSchemes
include("../src/analysis.jl")

bandname = "data/BAND_soc_afm.dat"
labelname = "data/KLABELS"

Nk,Nband = importBandSize(bandname)
dftband = importDFTBand(bandname)
klabel = importLabel(labelname)

fig = Figure()

for level in 1:4

    figloc = FigureLocation(level,2,2)

    ax = Axis(fig[figloc...],
    xticks = map(x -> klabel[1:4,x],(2,1)),
    width = 250,height = 200,
    xlabel = L"k\ /\ (1/\mathrm{A_0})",
    ylabel = L"E(k)\ /\ \mathrm{eV}",
    title = "level = $(level)")

    for ii in 2:size(dftband)[2]
        lines!(ax,dftband[:,1],dftband[:,ii],color = :black)
    end

    energydf = CSV.read("data/discrete energy-$(level).dat", DataFrame)
    fatbanddf = CSV.read("data/discrete fatband-$(level).dat", DataFrame)

    for jj in 2:size(energydf,2)
        scatter!(ax,energydf[:,"kr"],energydf[:,jj],markersize = 2 .* fatbanddf[:,jj];
        color = :blue)
    end

    xlims!(ax,klabel[[1,4],2]...)
    ylims!(ax,-1.5,1.5)

    FigureDecoration(ax,figloc,2,2)

end


resize_to_layout!(fig)
display(fig)
save("figures/FATBAND discrete band merge.png",fig)
