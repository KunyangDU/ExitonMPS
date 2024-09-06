using CairoMakie,ColorSchemes,LinearAlgebra
include("../src/analysis.jl")

bandname = "data/wannier TB/dz2 dxy dx2-y2.dat"
wncoeff = importWannierBand(bandname)
N = 36 # size of WN orbitals

kpath = basism(KBASIS2)*vrange(KPATH;eachstep = 10)
kr = pathlength(kpath)

band = DataFrame(:kr => kr)
for ii in 1:N
    band[!,string(ii)] = zeros(length(kr))
end

for (ii,k) in enumerate(eachcol(kpath))
    progress = ii/size(kpath,2)
    @show progress

    H = zeros(N,N) .* 1im
    for ii in 1:size(wncoeff)[1]
        H[wncoeff[ii,"orbital1"],wncoeff[ii,"orbital2"]] = H[wncoeff[ii,"orbital1"],wncoeff[ii,"orbital2"]] + (wncoeff[ii,"tRe"] + im * wncoeff[ii,"tIm"]) * exp(im * dot(vcat(k,0.0),basism(RBASIS3)*collect(wncoeff[ii,1:3])))
    end

    egv = eigvals(H)
    !approx(real.(egv),abs.(egv)) && @error "non Real!"
    band[ii,2:end] = abs.(egv)'
end

CSV.write("data/wannier TB/band.dat", band)
