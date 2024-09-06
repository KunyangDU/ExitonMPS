using CairoMakie,CSV,DataFrames,Interpolations
include("../src/analysis.jl")

broaden = 6
vbs,cbs,totalband = indexrange(broaden;maxvb = MAXVB)

bandname = "data/BAND_soc_afm.dat"
dftband = importDFTBand(bandname)
dftenergy = BandInterp(dftband,totalband)

for level in 1:4
    fatbandname =  "data/BSE-$(level).dat"
    ofatband = importFatBand(fatbandname)

    pathsite = FatbandSite(ofatband)
    pathr = MeasureDistance(pathsite)

    energydf = CoordinateDf(
    let
        ExtractBand(dftenergy,pathr,totalband) |> x -> vcat(x,x[1,:]')
    end,vcat(pathr,dftband."kr"[div(importBandSize(bandname)[1],3)]),"kr",totalband)
    
    fatbanddf = CoordinateDf(
    let
        SelectFatband(ofatband,pathsite,vbs,cbs) |> x -> vcat(x,x[1,:]')
    end,vcat(pathr,dftband."kr"[div(importBandSize(bandname)[1],3)]),"kr",totalband)

    CSV.write("data/discrete energy-$(level).dat", energydf)
    CSV.write("data/discrete fatband-$(level).dat", fatbanddf)
end

