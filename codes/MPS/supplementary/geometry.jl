function approxeq(testd::Number,ruled::Number;kwargs...)
    precision = get(kwargs,:precision,1e-5)

    if abs(testd-ruled)<precision
        return true
    else
        return false
    end
end

function isboundary(Latt::CompositeLattice,pair::Tuple)
    Lx,Ly = Latt.subLatts[1].sites[end]
    ind = Dict([(ii,get_index(Latt,ii)) for ii in 1:3])
    for jj in 1:2
        if pair[jj] in ind[1] && mod(pair[jj]-14,3*Ly)==0 && jj == 1 && pair[2] in ind[2] ||  pair[jj] in ind[3] && mod(pair[jj]-8,3*Ly)==0 && jj==2 && pair[1] in ind[2]
            pair == (14,20) && println(pair)
            return true
        end
    end

    return false
end

function get_index(Latt::AbstractLattice,orbital)

    orbital_ind = []

    for ii in 1:size(Latt)
        Latt[ii][1] in orbital && push!(orbital_ind,ii)
    end

    return orbital_ind

end

function get_cellsize(Latt::CompositeLattice)
    Lx = maximum([Latt.subLatts[1].sites[ii][1] for ii in 1:div(size(Latt),3)])
    Ly = maximum([Latt.subLatts[1].sites[ii][2] for ii in 1:div(size(Latt),3)])
    return Lx,Ly
end

function get_cellsize(Latt::SimpleLattice)
    Lx = maximum([Latt.sites[ii][1] for ii in 1:size(Latt)])
    Ly = maximum([Latt.sites[ii][2] for ii in 1:size(Latt)])
    return Lx,Ly
end