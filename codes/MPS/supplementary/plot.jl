function plot_XCKagome(Latt::CompositeLattice;
    figsize = (100,95) .* get_cellsize(Latt) ,
    bond = true,
    tplevel = (1),
    site = false,
    sitelabel = true,
    sitesize = 16 .* ones(size(Latt)),
    sitecolor = repeat(:grey,size(Latt)),
    sitealpha = ones(size(Latt)),
    )
    Lx,Ly = get_cellsize(Latt) 

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
    xticks = 1:2*Lx)


    if bond
        for level in tplevel
            # NN bond 
            for (i, j) in neighbor(Latt;level = level)

                    x = map([i, j]) do i
                        coordinate(Latt, i)[1]
                    end
                    y = map([i, j]) do i
                        coordinate(Latt, i)[2]
                    end

                    if norm(coordinate(Latt,i) .- coordinate(Latt,j)) > sqrt(3)/2 - 1e-5 
                        y[findmin(y)[2]] = y[findmin(y)[2]] + Ly*sqrt(3)/2
                    end

                    lines!(ax, x, y;
                        linewidth=2,
                        color=RGBf(0.5, 0.5, 0.5)
                    )
            end
        end
    end

    if site
        for i in 1:size(Latt)
                x, y = coordinate(Latt, i)

                CairoMakie.scatter!(ax, x, y;
                    markersize=sitesize[i],
                    color=sitecolor[i],
                    alpha = sitealpha[i])
        
                sitelabel && text!(ax, x*1.01, y*1.01; text = "$i")
        end
    end

    return fig,ax
end    

function plot_InfKagome(Latt::CompositeLattice;
    figsize = (100,95) .* get_cellsize(Latt) ,
    bond = true,
    tplevel = (1),
    site = false,
    sitelabel = true,
    sitesize = 16 .* ones(size(Latt)),
    sitecolor = repeat(:grey,size(Latt)),
    sitealpha = ones(size(Latt)),
    )
    Lx,Ly = get_cellsize(Latt) 

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
    xticks = 1:2*Lx)


    if bond
        for level in tplevel
            # NN bond 
            for (i, j) in neighbor(Latt;level = level)

                    x = map([i, j]) do i
                        coordinate(Latt, i)[1]
                    end
                    y = map([i, j]) do i
                        coordinate(Latt, i)[2]
                    end

                    if norm(coordinate(Latt,i) .- coordinate(Latt,j)) > sqrt(3)/2 - 1e-5 
                        y[findmin(y)[2]] = y[findmin(y)[2]] + Ly*sqrt(3)/2
                        x[findmin(x)[2]] = x[findmin(x)[2]] + Ly*1/2
                    end

                    lines!(ax, x, y;
                        linewidth=2,
                        color=RGBf(0.5, 0.5, 0.5)
                    )
            end
        end
    end

    if site
        for i in 1:size(Latt)
                x, y = coordinate(Latt, i)

                CairoMakie.scatter!(ax, x, y;
                    markersize=sitesize[i],
                    color=sitecolor[i],
                    alpha = sitealpha[i])
        
                sitelabel && text!(ax, x*1.01, y*1.01; text = "$i")
        end
    end

    return fig,ax
end    


function plot_Triangular(Latt::SimpleLattice;
    figsize = (100,95) .* get_cellsize(Latt) ,
    bond = true,
    tplevel = (1),
    site = false,
    sitelabel = true,
    sitesize = 16 .* ones(size(Latt)),
    sitecolor = repeat(:grey,size(Latt)),
    sitealpha = ones(size(Latt)),
    )

    Lx,Ly = get_cellsize(Latt) 

    fig = Figure(size=figsize)
    ax = Axis(fig[1, 1],
    xticks = 1:2*Lx)

    if bond
        for level in tplevel
            # NN bond 
            for (i, j) in neighbor(Latt;level = level)

                    x = map([i, j]) do i
                        coordinate(Latt, i)[1]
                    end
                    y = map([i, j]) do i
                        coordinate(Latt, i)[2]
                    end

                    if norm(coordinate(Latt,i) .- coordinate(Latt,j)) > sqrt(3) - 1e-5 
                        y[findmin(y)[2]] = y[findmin(y)[2]] + Ly*sqrt(3)/2
                        x[findmin(x)[2]] = x[findmin(x)[2]] + Ly*1/2
                    end

                    lines!(ax, x, y;
                        linewidth=2,
                        color=RGBf(0.5, 0.5, 0.5)
                    )
            end
        end
    end

    if site
        for i in 1:size(Latt)
                x, y = coordinate(Latt, i)

                CairoMakie.scatter!(ax, x, y;
                    markersize=sitesize[i],
                    color=sitecolor[i],
                    alpha = sitealpha[i])
        
                sitelabel && text!(ax, x, y + 0.1; text = "$i")
        end
    end

    return fig,ax
end 
