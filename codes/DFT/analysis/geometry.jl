


function coordinate(a::Matrix;basis = KBASIS2)
    return basism(basis)*a
end

function originsite(df::DataFrame;dims::Int64 = 2)
    return hcat(unique([collect(df[ii,1:dims]) for ii in 1:size(df)[1]])...)
end

function kbasis3(basis::Vector)
    basis = collect.(basis)
    V = dot(basis[1],cross(basis[2],basis[3]))
    b1 = cross(basis[1],basis[2])*2*pi/V
    b2 = cross(basis[2],basis[3])*2*pi/V
    b3 = cross(basis[3],basis[1])*2*pi/V
    return Tuple.([b1,b2,b3])
end

function kbasis2(basis::Vector)
    kbasis = kbasis3(basis)
    kbasis2 = []
    for kvec in kbasis
        kvec[1] != kvec[2] && push!(kbasis2,kvec[1:2])
    end
    return kbasis2
end


function FBZboundary!(ax::Axis,BDpoint::Vector;
    BDlinewidth::Number = 2.0,
    BDcolor::Symbol = :black,
    showbasis::Bool = false,
    arrowsize::Number = 0.2,
    arrowwidth::Number = 2.0,
    arrowcolor::Symbol = :blue,
    )
    coord = [basism(KBASIS2)*collect(vec) for vec in BDpoint]
    x = vcat([coord[ii][1] for ii in eachindex(coord)],coord[1][1])
    y = vcat([coord[ii][2] for ii in eachindex(coord)],coord[1][2])

    lines!(ax,x,y,linewidth = BDlinewidth,color = BDcolor)

    if showbasis
        arrow0!(ax,0,0,KBASIS2[1]...;arrowsize = arrowsize,color = arrowcolor,linewidth = arrowwidth)
        arrow0!(ax,0,0,KBASIS2[2]...;arrowsize = arrowsize,color = arrowcolor,linewidth = arrowwidth)
    end

end


function vrange(beginvec::Union{Vector,Tuple},endvec::Union{Vector,Tuple};step::Int64 = 100)
    return hcat([collect(beginvec .+ (endvec .- beginvec) .* t)  for t in range(0,1,step)]...)
end

function vrange(ipath::Vector;eachstep::Number = 100)
    finalpath = collect(ipath[1])

    for ii in eachindex(ipath)[1:end-1]
        finalpath = hcat(finalpath,vrange(ipath[ii],ipath[ii+1];step=eachstep+1)[:,2:end])
    end

    return finalpath
end

function pathlength(finalpath::Matrix)
    return cumsum(norm.(eachcol(hcat([0.0,0.0],diff(finalpath,dims = 2)))))
end

function ddist(tempdist,x::Number,domain::Tuple;h::Number=1e-5) 
    x == extrema(domain)[1] && return (tempdist(x+h) - tempdist(x))/(h)
    x == extrema(domain)[2] && return (tempdist(x) - tempdist(x-h))/(h)
    return (tempdist(x+h) - tempdist(x-h))/(2h)
end

function dist(itp,x::Number,target::Vector;basis::Matrix = basism(KBASIS2))
    return norm(map(y -> itp[y](x),collect(eachindex(itp))) .- basis*target)
end


function get_cellsize(Latt::CompositeLattice)
    return map(x -> maximum([Latt.subLatts[1].sites[ii][x] for ii in 1:div(size(Latt),length(Latt.subLatts))]),(1,2))
end
