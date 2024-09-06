using CairoMakie,ColorSchemes,Interpolations
include("../src/analysis.jl")

fig = Figure()
limit = 1.0 .* (0,20)
for level in 1:4
    path = string("data/BSE-",level,".dat")
    df = importFatBand(path)

    osite = originksite(df)
    value = cpcoeff(df,osite,Tuple(173:178),Tuple(179:184))
    cpdict = datadict(value,osite)

    cpdict,osite = complete(cpdict,osite)

    cood = coordinate(osite)

    interpdata = zeros(3,3)

    for (ix,x) in enumerate(0:0.25:0.5),(iy,y) in enumerate(0:0.25:0.5)
        tempsite = basism(INTERPKBASIS2)*[x,y]
        for (icd,cd) in enumerate([cood[:,ii] for ii in 1:size(cood)[2]])
            if approx(cd,tempsite)
                interpdata[ix,iy] = cpdict[osite[:,icd]]
                break
            end
        end
    end

    x,y,z = ietp_cubic2(0:0.25:0.5,0:0.25:0.5,interpdata;
    xlims = (0,0.5),ylims = (0,0.5),interp_num = (10,10))


    itpsite = hcat([[xx,yy] for xx in x, yy in y]...)
    itpvalue = vec(z)

    figloc = convert(Int64,ceil(level/2)),renormalize(level,2)



    ax = Axis(fig[figloc...],
    width = 200,height = 200,
    autolimitaspect  = 1,
    title = "level = $(level)",
    xlabel = L"k_x",
    ylabel = L"k_y")
    FBZboundary!(ax,[tp ./ pi for tp in FBZBD],showbasis=false)

    sctcolor = get(colorschemes[:hot],itpvalue,limit)
    for n in 0:2
        itpksite = [-1/2 -sqrt(3)/2;sqrt(3)/2 -1/2]^n*coordinate(itpsite;basis = INTERPKBASIS2)
        scatter!(ax,map(x -> itpksite[x,:] ./ pi,[1,2])...,color = sctcolor)
    end

    figloc[1] != 2 && hidexdecorations!(ax,grid = false)
    figloc[2] != 1 && hideydecorations!(ax,grid = false)

end
Colorbar(fig[:,3],limits = limit, colormap = :hot,
label = L"\sum_{\alpha=v,c} \zeta_{\alpha,\mathbf{k}}^{\mathbf{q}}")
resize_to_layout!(fig)
display(fig)

save("figures/FATBAND intep BZ merge.png",fig)

