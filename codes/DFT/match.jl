#= fatband = Dict()
for vb in vbs
    fatband[vb] = cpcoeff(ofatband,osite,vb,cbs)
end

for cb in cbs
    fatband[cb] = cpcoeff(ofatband,osite,vbs,cb)
end =#

#= pathsite = []
for site in eachcol(osite)
    ((site[1] == site[2] && site[1] >= 0) || (site[1] >= 0 && site[2] == 0) || ((site .- [1/2,0],site .- [1/3,1/3]) |> x -> approx(abs(dot(x[1],x[2])),norm(x[1])*norm(x[2])) )) && push!(pathsite,site)
end

kpath = basism(KBASIS2)*hcat(pathsite...)
kr = pathlength(kpath)
pathsite =#

using CairoMakie,CSV,DataFrames
include("src/analysis.jl")

fatbandname =  "data/BSE-1.dat"
bandname = "data/BAND_soc_afm.dat"
labelname = "data/KLABELS"

Nk,Nband = importBandSize(bandname)
dftband = importDFTBand(bandname)

ofatband = importFatBand(fatbandname)
osite = originksite(ofatband)
maxvb = MAXVB
broaden = 6
vbs = maxvb-broaden+1:maxvb
cbs = maxvb+1:maxvb+broaden
totalband = vcat(vbs,cbs)


fatbandkr = []
tempenergy = Dict()
tempbroaden = Dict()

for (iisite,site) in enumerate(eachcol(basism(KBASIS2)*osite)),(iipt,point) in enumerate(eachcol(kpath))
    if approx(norm(site .- point),0)
        push!(fatbandkr,kr[iipt])

        tempenergy[length(fatbandkr)] = zeros(length(totalband))

        for (iibi,bandindex) in enumerate(totalband)
            tempenergy[length(fatbandkr)][iibi] = dftband[iipt,string(bandindex)]
        end

        tempbroaden[length(fatbandkr)] = zeros(length(totalband))

        for (iivb,vb) in enumerate(vbs)
            tempbroaden[length(fatbandkr)][iivb] = cpcoeff(ofatband,osite[:,iisite][:,1:1],vb,cbs)[1]
        end
        
        for (iicb,cb) in enumerate(cbs)
            tempbroaden[length(fatbandkr)][broaden + iicb] = cpcoeff(ofatband,osite[:,iisite][:,1:1],vbs,cb)[1]
        end
    else
        continue
    end
end

dfenergy,dfbroaden = map(x -> DataFrame(hcat(fatbandkr,hcat([x[ii] for ii in eachindex(fatbandkr)]...)'),vcat("kr",string.(totalband))),(tempenergy,tempbroaden))
