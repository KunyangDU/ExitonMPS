

function importFatBand(filename::String)
    dfo = CSV.read(filename, DataFrame;header = false)
    return unique(DataFrame(hcat(map(x -> map(x -> parse(Float64, x), split(replace(x,"+i*" => " "))),dfo[:,"Column1"])...)[[1,2,3,7,8,9,10,6],:]',["h1","h2","h3","VB","CB","RCC","ICC","ACC"]))
end

function importDFTBand(filename::String)
    dfo = CSV.read(filename, DataFrame;header = false,comment="#",delim='=')
    return DataFrame(hcat(map(x -> parse.(Float64,x),filter(x -> x != [""],map(x -> replace.(split(x,"    ")," " => ""),dfo[:,1])))...)',["kr","band"])
end

function importLabel(filename::String;replacedict = Dict([("GAMMA", "Γ")]))
    df = CSV.read(filename,DataFrame;delim="=")
    splited = map(x -> split(x," "),df[:,1])
    modisplited = []

    for (isp,sp) in enumerate(splited)
        tempsp = []
        for ssp in sp
            ssp == "" && continue
            push!(tempsp,ssp)
        end
        length(tempsp) != 2 && continue
        push!(modisplited,tempsp)
    end

    klinesdf = DataFrame(hcat(collect(map(x -> (let
        if x == 1
            [get(replacedict,modisplited[ii][x],modisplited[ii][x]) for ii in eachindex(modisplited)]
        else
            [parse(Float64,modisplited[ii][x]) for ii in eachindex(modisplited)]
        end
    end),
    (1,2)))...),["KLABELS","KSTICKS"])
    
    return klinesdf
end

function importWannierBand(filename::String)
    dfo = CSV.read(filename, DataFrame;header = false,delim="   ")[:,2:end]
    return hcat(map(x3 -> DataFrame(hcat(map( x2 -> map(x1 -> parse(x3[2],x1), dfo[:,x2]),x3[1])...),x3[3]),[(1:5,Int64,["n1","n2","n3","orbital1","orbital2"]),(6:7,Float64,["tRe","tIm"])])...)
end

function importDFTBand(filename::String)
    dfo = CSV.read(filename, DataFrame;header = false,comment="#",delim='=')

    kr,band = map(x -> hcat(map(x -> parse.(Float64,x),filter(x -> x != [""],map(x -> replace.(split(x,"    ")," " => ""),dfo[:,1])))...)[x,:],(1,2))
    Nk,Nband = importBandSize(filename)
    bandm = hcat(band[1:Nk,1],zeros(Nk,Nband-1))

    for ii in 1:Nband-1
        !isequal(kr[Nk*(ii-1)+1:Nk*ii],kr[Nk*(ii+1):-1:Nk*ii+1]) && @show ii
        bandm[:,ii+1] = let 
            if iseven(ii)
                band[Nk*ii+1:Nk*(ii+1)]
            else
                band[Nk*(ii+1):-1:Nk*ii+1]
            end
        end
    end

    return DataFrame(hcat(kr[1:Nk],bandm),vcat(["kr"],string.(1:Nband)))
end

function importBandSize(filename::String)
    Nk,Nband = Tuple(parse.(Int64,split(CSV.read(filename, DataFrame;header = false,comment=" ",delim='=')[2,1]," ")[end-1:end]))
    return Nk,Nband
end

