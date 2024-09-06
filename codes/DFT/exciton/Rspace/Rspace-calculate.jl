using CSV,DataFrames,CairoMakie,JLD2
include("../analysis/analysis.jl")

#Latt = XCYYKag(4,4;shift = Tuple(Tuple.(eachcol(unitRbasis2*Nbpos2))))
posname = "data/POSCAR"
unitRbasis2,LattConst = RBasis2(posname)
Nbpos2 = NbPosition2(posname)

Latt = XCmultiYYKag(6,6;norbs = 3,
atomshift = Tuple(Tuple.(eachcol(unitRbasis2*diagm(LattConst)*Nbpos2))),
scale = LattConst[1])

eigenstates = load("data/wannier TB/eigenstates.jld2","eigenstates")
broaden = 6
level = 1
fatbandname =  "data/BSE-$(level).dat"
ofatband = importFatBand(fatbandname)
vbs,cbs,totalband = indexrange(broaden;maxvb = MAXVB)
tbvbs,tbcbs,tbtotalband = indexrange(broaden;maxvb = TBMAXVB)
dft2TB = datadict(tbtotalband,totalband)


pvec = 1:div(size(eigenstates,2)-1,2)

spin = [(1,1),(1,0),(1,-1),(0,0)]
for (iispin,(S²,Sz)) in enumerate(spin)
    @show (S²,Sz)
    RspaceCC = zeros(div(size(eigenstates,2)-1,2),size(Latt))

    for (ip,p₀) in enumerate(pvec)
        R₀ = collect(FiniteLattices.coordinate(Latt,get_index(Latt,p₀,map(x -> div(x,2),get_cellsize(Latt)))))
        for ii in 1:size(Latt)
            progress = round(100*((ii/size(Latt)+ip-1)/length(pvec) + iispin-1)/(length(spin)),digits = 3)
            @show progress
            R = collect(FiniteLattices.coordinate(Latt,ii))
            p = Latt[ii][1]
            RspaceCC[ip,ii] = let
                RspaceCpCoeff(ofatband,eigenstates,
                vbs,cbs,dft2TB,
                R₀,p₀,R,p,S²,Sz)
            end
        end
        CSV.write("data/RCC/Rspace CpCoeff-$(p₀)-$(S²) $(Sz).dat",DataFrame(RspaceCpCoeff=RspaceCC[ip,:]))
        #CSV.write("data/test-$(p₀).dat",DataFrame(RspaceCpCoeff=RspaceCC[ip,:]))
    end
    RSCCdf = CoordinateDf(RspaceCC,collect(pvec),"p",string.(1:size(Latt)))
    CSV.write("data/RCC/Rspace CpCoeff-$(S²) $(Sz).dat",RSCCdf)
    #CSV.write("data/test.dat",RSCCdf)
end



