using CairoMakie,ColorSchemes
include("../../analysis/analysis.jl")


level = 1

path = "data/BSE-$(level).dat"
df = importFatBand(path)

fig = Figure()

ax = Axis(fig[1,1],
width = 300,
height = 300,
autolimitaspect  = 1,
title = "AFM, q = 0, level=1\n band index = (177,178)x(179,180)",
xlabel = L"k_x",
ylabel = L"k_y")
FBZboundary!(ax,FBZBD,showbasis=true)
osite = originksite(df)
value = cpcoeff(df,osite,Tuple(174:178),Tuple(179:182))

cpdict = datadict(value,osite)
cpdict,osite = complete(cpdict,osite)

cood = coordinate(osite)

x = cood[1,:]
y = cood[2,:]
sctcolor = get(colorschemes[:berlin],[cpdict[osite[:,oii]] for oii in 1:size(osite)[2]],:extrema)
for ii in 1:size(cood)[2]
    scatter!(ax,x[ii],y[ii],color = sctcolor[ii],markersize = 15)
end

Colorbar(fig[1,2],limits = extrema([cpdict[osite[:,oii]] for oii in 1:size(osite)[2]]), colormap = :berlin,
label = L"\sum_{\alpha=v,c} \zeta_{\alpha,\mathbf{k}}^{\mathbf{q}}")
resize_to_layout!(fig)
display(fig)

save("figures/FATBAND discrete BZ divided-$(level).png",fig)

