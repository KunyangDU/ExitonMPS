function arrow0!(ax::Axis,x, y, u, v; arrowsize=0.386, color=:black, transparency=1,linewidth = 1.2)
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = arrowsize*nuv*v4, arrowsize*nuv*v5
    lines!(ax,[x,x+u], [y,y+v], color=(color,transparency),linewidth = linewidth)
    lines!(ax,[x+u,x+u-v5[1]], [y+v,y+v-v5[2]], color=(color,transparency),linewidth = linewidth)
    lines!(ax,[x+u,x+u-v4[1]], [y+v,y+v-v4[2]], color=(color,transparency),linewidth = linewidth)
end

function FBZboundary!(ax::Axis,BDpoint::Matrix = FBZBD;
    BDlinewidth::Number = 2.0,
    BDcolor::Symbol = :black,
    showbasis::Bool = false,
    arrowsize::Number = 0.2,
    arrowwidth::Number = 2.0,
    arrowcolor::Symbol = :blue,
    )
    coord = coordinate(BDpoint)
    x = coord[1,:]
    y = coord[2,:]

    lines!(ax,x,y,linewidth = BDlinewidth,color = BDcolor)

    if showbasis
        arrow0!(ax,0,0,KBASIS2[1]...;arrowsize = arrowsize,color = arrowcolor,linewidth = arrowwidth)
        arrow0!(ax,0,0,KBASIS2[2]...;arrowsize = arrowsize,color = arrowcolor,linewidth = arrowwidth)
    end

end


function site!(ax::Axis,site::Vector)
    cood = coordinate(site)
    x = [cood[ii][1] for ii in eachindex(cood)]
    y = [cood[ii][2] for ii in eachindex(cood)]
    scatter!(ax,x,y)
end

function kpath!(ax::Axis,path::Vector = KPATH;
    basis::Vector = KBASIS2,color::Symbol = :red,linewidth::Number=2.0)
    
    lines!(ax,map(x -> (basism(basis)*hcat(collect.(path)...))[x,:],(1,2))...;
    color = color,linewidth = linewidth)

end

function band!(ax::Axis,df::DataFrame;kwargs...)
    for ii in 2:size(df)[2]
        lines!(ax,df."kr",df[:,ii];kwargs...)
    end
end

function FigureLocation(level::Int64,row::Int64,column::Int64)
    return convert(Int64,ceil(level/row)),renormalize(level,column)
end

function FigureDecoration(ax::Axis,figloc::Tuple,row::Int64,column::Int64)
    figloc[1] != row && hidexdecorations!(ax,grid = false)
    figloc[2] != 1 && hideydecorations!(ax,grid = false)
end











