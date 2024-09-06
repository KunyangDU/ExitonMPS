
global A₀      = (6.80130,6.80130,12.76793)
global RBASIS3 = [(1/2,sqrt(3)/2,0.0),(1/2,-sqrt(3)/2,0.0),(0.0,0.0,1.0)] |> x -> [A₀ .* x[ii] for ii in 1:3]
global KBASIS2 = kbasis2(RBASIS3)
global FBZBD   = hcat(collect.([(1/3,1/3),(2/3,-1/3),(1/3,-2/3),(-1/3,-1/3),(-2/3,1/3),(-1/3,2/3),(1/3,1/3)])...)
global INTERPKBASIS2 = [KBASIS2[1],KBASIS2[2] .- KBASIS2[1]]
global KPATH   = [(0.0,0.0),(0.5,0.0),(1/3,1/3),(0.0,0.0)]
global MAXVB   = 178
global TBMAXVB = 18

function complete(cpdict::Dict,osite::Matrix)

    cpdict[[-0.5,-0.0]] = cpdict[[0.5,-0.0]]
    cpdict[[-0.0,0.5]] = cpdict[[-0.0,-0.5]]
    cpdict[[0.5,-0.5]] = cpdict[[-0.5,0.5]]
    
    for completevec in [[-0.5,-0.0],[-0.0,0.5],[0.5,-0.5]]
        osite = hcat(osite,completevec)
    end

    return cpdict,osite

end

function checkpath(site::Vector)
    return ((site[1] == site[2] && site[1] >= 0) || (site[1] >= 0 && site[2] == 0) || ((site .- [1/2,0],site .- [1/3,1/3]) |> x -> approx(abs(dot(x[1],x[2])),norm(x[1])*norm(x[2])) ))
end

function RBasis3(filename::String)
    poscar = CSV.read(filename,DataFrame)
    rbasis3 = hcat(map(x -> parse.(Float64,split(x)),poscar[2:4,1])...)
    unitRbasis3 = hcat(map(x -> x ./ norm(x),eachcol(rbasis3))...)
    LattConst = norm.(collect.(eachcol(rbasis3)))
    return unitRbasis3,LattConst
end

function RBasis2(filename::String)
    unitRbasis3,LattConst = RBasis3(filename)
    return unitRbasis3[1:2,1:2],LattConst[1:2]
end

function NbPosition2(filename::String)
    poscar = CSV.read(filename,DataFrame)
    atompos3 = hcat(map(x -> parse.(Float64,split(x)[1:end-1]),poscar[8:end,1])...)
    Nbpos2 = atompos3[1:2,1:6]
    return Nbpos2
end

