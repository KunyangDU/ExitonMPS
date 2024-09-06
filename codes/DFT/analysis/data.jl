function singlecpcoeff(df::DataFrame,site::Matrix,VB::Number,CB::Number;key = "ACC")
    return [df[(df."h1" .== site[1,ii]) .& (df."h2" .== site[2,ii]) .& (df."h3" .== 0.0).& (df."VB" .== VB).& (df."CB" .== CB),key][1]  for ii in 1:size(site)[2]]
end


function cpcoeff(df::DataFrame,site::Matrix,VB,CB;key = "ACC")
    return sum([singlecpcoeff(df,site,vb,cb;key=key) for vb in VB, cb in CB])
end


function datadict(value::Vector,site::Matrix)
    return Dict([(site[:,ii],value[ii]) for ii in 1:size(site)[2]])
end

function datadict(value::Vector,site::Vector)
    return Dict([(site[ii],value[ii]) for ii in eachindex(site)])
end

function FatBand(filename::String;broaden::Int64 = 2,maxvb::Int64 = MAXVB)
    ofatband = importFatBand(filename)

    vbs = maxvb-broaden+1:maxvb
    cbs = maxvb+1:maxvb+broaden
    fatband = Dict()

    for vb in vbs
        fatband[vb] = cpcoeff(ofatband,originksite(ofatband),vb,cbs)
    end

    for cb in cbs
        fatband[cb] = cpcoeff(ofatband,originksite(ofatband),vbs,cb)
    end

    return fatband

end


function FatbandSite(fatband::DataFrame,)
    pathsite = []
    for site in map(collect,eachcol(originsite(fatband)))
        if checkpath(site)
            push!(pathsite,site)
        end
    end

    return hcat(pathsite...)
end

function SelectFatband(fatband::DataFrame,selectsites::Matrix,vbs::Vector,cbs::Vector)

    tempbroaden = zeros(size(selectsites,2),sum(map(length,(vbs,cbs))) )
    
    for (iivb,vb) in enumerate(vbs)
        tempbroaden[:,iivb] = cpcoeff(fatband,selectsites,vb,cbs)
    end

    for (iicb,cb) in enumerate(cbs)
        tempbroaden[:,length(vbs)+iicb] = cpcoeff(fatband,selectsites,vbs,cb)
    end
    
    return tempbroaden
    
end

function ExtractBand(itps,selectr::Vector,totalband::Vector)
    
    temp = zeros(length(selectr),length(totalband))
    for ii in eachindex(selectr)
        temp[ii,:] =  [itps[jj](selectr[ii]) for jj in eachindex(totalband)]
    end

    return temp
end

function CoordinateDf(datam::Matrix,x::Vector,xindex::String,yindexs::Vector)
    return DataFrame(hcat(x,datam),string.(vcat(xindex,yindexs)))
end

function RspaceCpCoeff(fatband::DataFrame,eigenstates::DataFrame,
    vbs::Vector,cbs::Vector,dft2TB::Dict,
    R1::Vector,p1::Int64,R2::Vector,p2::Int64,
    S²::Int64,Sz::Int64)

    oη = 1im .* zeros(size(eigenstates,1))
    for vb in vbs,cb in cbs
        As = map(x -> dot(x,[1,1im]),collect.(eachcol(hcat(map(key -> cpcoeff(fatband,originsite(ofatband),vb,cb;key),["RCC","ICC"])...)')))
        vbBs = map(p -> map(x -> x[p], eigenstates[:,string(dft2TB[vb])]),(p1,div(size(eigenstates,2)-1,2) + p1))
        cbBs = map(p -> map(x -> x[p], eigenstates[:,string(dft2TB[vb])]),(p2,div(size(eigenstates,2)-1,2) + p2))

        if S² == 1
            Sz == 1 && (@. oη +=  As*vbBs[1]*conj(cbBs[1]))
            Sz == 0 && (@. oη +=  As*(vbBs[1]*conj(cbBs[2]) + vbBs[2]*conj(cbBs[1]))/sqrt(2))
            Sz == -1 && (@. oη +=  As*vbBs[2]*conj(cbBs[2]))
        elseif S² == 0
            Sz != 0 && error("Sz is not 0")
            @. oη +=  As*(vbBs[1]*conj(cbBs[2]) - vbBs[2]*conj(cbBs[1]))/sqrt(2)
        else
            error("S² not in {0,1}")
        end
        
    end
    oη = dot(map(x -> exp(1im * dot(x,R1 .- R2)),eigenstates."kpoints"),oη)
    η = abs(oη)

    return η
    
end

function RspaceCpCoefftest(fatband::DataFrame,eigenstates::DataFrame,
    vbs::Vector,cbs::Vector,dft2TB::Dict,
    R1::Vector,p1::Int64,R2::Vector,p2::Int64,
    S²::Int64,Sz::Int64)

    oη = (1+0.0im) .* ones(size(eigenstates,1))
    oη = dot(map(x -> exp(1im * dot(x,R1 .- R2)),eigenstates."kpoints"),oη)
    η = abs.(oη)

    return η
    
end
