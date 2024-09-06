using CairoMakie,ColorSchemes,LinearAlgebra,JLD2
include("../src/analysis.jl")

bandname = "data/wannier TB/dz2 dxy dx2-y2.dat"
wncoeff = importWannierBand(bandname)
N = 36 # size of WN orbitals
level = 1
fatbandname =  "data/BSE-$(level).dat"

ofatband = importFatBand(fatbandname)

kpoints = coordinate(originsite(ofatband))

vecs = DataFrame(:kpoints => collect.(eachcol(kpoints)))
for ii in 2:N+1
    vecs[:,string(ii)] = [1im .* zeros(N) for _ in 1:size(kpoints,2)]
end


for (ii,k) in enumerate(collect.(eachcol(kpoints)))
    progress = ii/size(kpoints,2)
    @show progress

    H = zeros(N,N) .* 1im
    for ii in 1:size(wncoeff)[1]
        H[wncoeff[ii,"orbital1"],wncoeff[ii,"orbital2"]] = H[wncoeff[ii,"orbital1"],wncoeff[ii,"orbital2"]] + (wncoeff[ii,"tRe"] + im * wncoeff[ii,"tIm"]) * exp(im * dot(vcat(k,0.0),basism(RBASIS3)*collect(wncoeff[ii,1:3])))
    end

    egv = eigvecs(H)
    vecs[ii,2:end] = collect.(eachcol(egv))
end

save("data/wannier TB/eigenstates.jld2","eigenstates",vecs)


