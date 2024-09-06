using CairoMakie,CSV,DataFrames,Interpolations,DataInterpolations
include("../src/analysis.jl")

broaden = 6
vbs,cbs,totalband = indexrange(broaden;maxvb = MAXVB)

Ninterp = 50
bandname = "data/BAND_soc_afm.dat"
dftband = importDFTBand(bandname)

kr = collect(range(extrema(dftband[1:div(importBandSize(bandname)[1],3),"kr"])...,Ninterp))
tempenergy =  ExtractBand(BandInterp(dftband,totalband),kr,totalband)
energydf = CoordinateDf(tempenergy,kr,"kr",totalband)

for level in 1:4
    fatbandname =  "data/BSE-$(level).dat"
    ofatband = importFatBand(fatbandname)

    pathsite = FatbandSite(ofatband)
    pathr = MeasureDistance(pathsite)

    # periodize
    discbroaden = CoordinateDf(
    let 
        SelectFatband(ofatband,pathsite,vbs,cbs) |> x -> vcat(x,x[1,:]')
    end,vcat(pathr,kr[end]),"kr",totalband)

    fatbanddf = CoordinateDf(
    let 
        ExtractBand(FatbandInterp(discbroaden,totalband),kr,totalband)
    end,kr,"kr",totalband)

    CSV.write("data/interp energy-$(level).dat", energydf)
    CSV.write("data/interp fatband-$(level).dat", fatbanddf)
end

