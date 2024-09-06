function newton(f, f_prime, x0; tol=1e-5, max_iter=1000)
    x = x0
    for i in 1:max_iter
        fx = f(x)
        fpx = f_prime(x)

        if abs(fx) < tol
            println("Converged to zero point in $i iterations.")
            return x_new
        end
        
        if abs(fpx) < tol
            println("Derivative too small; stopped to avoid division by zero. f(x) = $(fx)")
            return x
        end
        x_new = x - fx / fpx
        x = x_new
    end
    println("Maximum number of iterations reached without convergence.")
    return x
end

function MeasureDistance(target::Vector;
    ipath::Matrix = basism(KBASIS2)*vrange(KPATH),
    x0::Number = 0)
    
    ir = pathlength(ipath)
    itp = map(x -> linear_interpolation(ir,ipath[x,:]),[1,2])
    tempdist(x) = dist(itp,x,target)
    tempddist(x) = ddist(tempdist,x,extrema(ir))
    return newton(tempdist,tempddist,x0)
end

function MeasureDistance(targets::Matrix;
    ipath::Matrix = basism(KBASIS2)*vrange(KPATH),
    x0::Number = 0)
    return map(x -> MeasureDistance(collect(x);ipath=ipath,x0=x0),eachcol(targets))
end

function renormalize(index::Int64,unit::Int64)
    modifiedindex = mod(index,unit)
    if modifiedindex == 0
        return unit
    else
        return modifiedindex
    end
end

function approx(a::Number,b::Number;tol::Number = 1e-5)
    abs(a .- b) < tol && return true
    return false
end

function approxzero(M::Matrix;tol = 1e-5)
    R,C = size(M)
    for rr in 1:R,cc in 1:C
        if abs(M[rr,cc]) > tol
            @show rr,cc
            return false
        end
    end
    return true
end

function basism(basis::Vector)
    return hcat(collect.(basis)...)
end

function approx(a::Vector,b::Vector;tol::Number = 1e-5)
    norm(a .- b) < tol && return true
    return false
end

function indexrange(broaden::Int64;maxvb::Int64 = MAXVB)
    vbs = maxvb-broaden+1:maxvb
    cbs = maxvb+1:maxvb+broaden
    totalband = vcat(vbs,cbs)
    return map(collect,(vbs,cbs,totalband))
end



