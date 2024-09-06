using GeoStats
import CairoMakie as Mke


include("../../analysis/analysis.jl")


level = 2

path = "data/BSE-$(level).dat"
df = importFatBand(path)

fig = Figure()

ax = Axis(fig[1,1],
width = 300,
height = 300,
autolimitaspect  = 1,
title = "FATBAND, level = $(level)",
xlabel = L"k_x",
ylabel = L"k_y")

osite = originsite(df)
value = cpcoeff(df,osite,Tuple(174:178),Tuple(179:182))
cood = coordinate(osite)

# attribute table
table = (; Z=value)
# coordinates for each row
coord = Tuple.(eachcol(cood))
# georeference data
geotable = georef(table, coord)

# interpolation domain
grid = CartesianGrid((200,200),(-1,-1),(0.01,0.01))

# choose an interpolation model
model = Kriging(GaussianVariogram(range=0.5))
# perform interpolation over grid
interp = geotable |> Interpolate(grid, model)

# visualize the solution
#viz(interp.geometry, color = interp.Z)
Mke.plot!(ax,interp.geometry,color = interp.Z)
FBZboundary!(ax,FBZBD)

Mke.resize_to_layout!(fig)
display(fig)

save("figures/FATBAND intep BZ Kriging-$(level).png",fig)



