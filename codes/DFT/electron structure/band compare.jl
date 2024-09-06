# compare DFT and Wannier band structure

using CairoMakie,ColorSchemes
include("../src/analysis.jl")

# DFT band

filename = "data/BAND_soc_afm.dat"
Nk,Nband = Tuple(parse.(Int64,split(CSV.read(filename, DataFrame;header = false,comment=" ",delim='=')[2,1]," ")[end-1:end]))
dfo = CSV.read(filename, DataFrame;header = false,comment="#",delim='=')

datam = hcat(map(x -> parse.(Float64,x),filter(x -> x != [""],map(x -> replace.(split(x,"    ")," " => ""),dfo[:,1])))...)
kr,band = map(x -> datam[x,:],(1,2))
bandm = hcat(band[1:Nk,1],zeros(Nk,Nband-1))


for ii in 1:Nband-1
    !isequal(kr[Nk*(ii-1)+1:Nk*ii],kr[Nk*(ii+1):-1:Nk*ii+1]) && @show ii
    bandm[:,ii+1] = let 
        if iseven(ii)
            band[Nk*ii+1:Nk*(ii+1)]
        else
            band[Nk*(ii+1):-1:Nk*ii+1]
        end
    end
end

df = DataFrame(hcat(kr[1:Nk],bandm),vcat(["kr"],string.(1:Nband)))

klabeldf = importLabel("data/KLABELS")

fig = Figure()
ax = Axis(fig[1,1],
xticks = map(x -> klabeldf[1:4,x],(2,1)),
width = 300,height = 250,
xlabel = L"k\ /\ (1/\mathrm{A_0})",
ylabel = L"E(k)\ /\ \mathrm{eV}",
title = "DFT band")

lgd = []

for ii in 2:size(df)[2]
    dft = lines!(ax,df[:,1],df[:,ii],color = :black,linewidth = 0.8)
    ii == 2 && push!(lgd,dft)
end

# Wannier band

WannierEshift = -3.03
N = 36

band = CSV.read("data/wannier TB/band.dat", DataFrame)

for ii in 1:N
    wn = lines!(ax,band."kr" ,band[:,string(ii)] .+ WannierEshift,color = :red)
    ii == 1 && push!(lgd,wn)
end

Legend(fig[1,2],lgd,["DFT","Wannier"])

xlims!(ax,extrema(band."kr" )...)
ylims!(ax,-2,2)
resize_to_layout!(fig)
display(fig)

save("figures/band compare.png",fig)



