function Base.repeat(sign::Symbol,n::Int64)
    return [sign for _ in 1:n]
end

function Base.repeat(tpsign::NTuple{N,Symbol},n::Int64) where N
    return [tpsign for _ in 1:n]
end