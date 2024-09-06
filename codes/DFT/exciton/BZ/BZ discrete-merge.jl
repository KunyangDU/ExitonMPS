using CairoMakie,ColorSchemes
include("../src/analysis.jl")

fig = Figure()

for level in 1:4

    path = "data/BSE-$(level).dat"
    df = importFatBand(path)

    figloc = FigureLocation(level,2,2)

    ax = Axis(fig[figloc...],
    width = 200,height = 200,
    autolimitaspect  = 1,
    title = "level = $(level)",
    xlabel = L"k_x",
    ylabel = L"k_y")
    FBZboundary!(ax,FBZBD)
    osite = originksite(df)
    value = cpcoeff(df,osite,Tuple(173:178),Tuple(179:184))

    cpdict = datadict(value,osite)
    cpdict,osite = complete(cpdict,osite)

    cood = coordinate(osite)

    x = cood[1,:]
    y = cood[2,:]
    sctcolor = get(colorschemes[:berlin],[cpdict[osite[:,oii]] for oii in 1:size(osite)[2]],:extrema)
    for ii in 1:size(cood)[2]
        scatter!(ax,x[ii],y[ii],color = sctcolor[ii],markersize = 15)
    end
   
    FigureDecoration(ax,figloc,2,2)

    if level == 4 
        Colorbar(fig[:,3],limits = extrema([cpdict[osite[:,oii]] for oii in 1:size(osite)[2]]), colormap = :berlin,
        label = L"\sum_{\alpha=v,c} \zeta_{\alpha,\mathbf{k}}^{\mathbf{q}}")
    end

end

resize_to_layout!(fig)
display(fig)
save("figures/FATBAND discrete BZ merge.png",fig)

