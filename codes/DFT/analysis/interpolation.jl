function ietp_cubic2(x::Union{StepRangeLen,Array},y::Union{StepRangeLen,Array},z::Matrix;
    xlims=extrema(x),ylims=extrema(y),
    interp_num=(100,100))
    
    # row for x, column for y

    !isequal(length.((x,y)),size(z)) && @error "size not compitable!"

    rx,ry = (x,y) |> x -> Tuple([range(extrema(x)...,length(x)) for x in x])

    itp = cubic_spline_interpolation((rx, ry), z;bc=Line(OnGrid()), extrapolation_bc=Throw())

    x2 = range(xlims..., length=interp_num[1])
    y2 = range(ylims..., length=interp_num[2])
    z2 = [itp(x,y) for y in y2, x in x2]

    return (x2,y2,z2)
end


function BandInterp(band::DataFrame,totalband::Vector)
    return map(x -> linear_interpolation(band[:,"kr"],band[:,string(x)]),totalband)
end

function FatbandInterp(broaden::DataFrame,totalband::Vector)
    return map(x -> DataInterpolations.CubicSpline(broaden[:,string(x)], broaden[:,"kr"]),totalband)
end
