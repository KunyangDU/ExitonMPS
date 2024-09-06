function XCKagome(Lx::Int64,Ly::Int64;kwargs...)
    Latt_A = XCTria(Lx, Ly)
    Latt_B = XCTria(Lx, Ly)
    Latt_C = XCTria(Lx, Ly)
    shift = get(kwargs,:shift,((0.75, sqrt(3)/4), (0.5, 0.0), (0.25,sqrt(3)/4)))

    path = get(kwargs,:path,Snake!)
    
    Latt = CompositeLattice(Latt_A, Latt_B, Latt_C, shift) |> path
    
    return Latt
end

function YCKagome(Lx::Int64,Ly::Int64;kwargs...)
    Latt_A = YCTria(Lx, Ly)
    Latt_B = YCTria(Lx, Ly)
    Latt_C = YCTria(Lx, Ly)
    shift = get(kwargs,:shift,((0.75, sqrt(3)/4), (0.5, 0.0), (0.25,sqrt(3)/4)))

    path = get(kwargs,:path,Snake!)
    
    Latt = CompositeLattice(Latt_A, Latt_B, Latt_C, shift) |> path
    
    return Latt
end

function InfTria(L::Int64, W::Int64, θ::Real = 0.0;)
    @assert L ≥ W && iseven(W)
    e = ((1.0, 0.0), (1/2, sqrt(3)/2))
    sites = [(x, y) for x in 1:L for y in 1:W]

    if iszero(θ)
         BC = PeriodicBoundaryCondition((0, W))
    else
         BC = TwistBoundaryCondition((0, W), θ)
    end
    return TriangularLattice(e, sites, BC)
end

function InfKagome(Lx::Int64,Ly::Int64;kwargs...)
    Latt_A = InfTria(Lx, Ly)
    Latt_B = InfTria(Lx, Ly)
    Latt_C = InfTria(Lx, Ly)
    #shift = get(kwargs,:shift,((0.75, sqrt(3)/4), (0.5, 0.0), (0.25,sqrt(3)/4)))
    shift = get(kwargs,:shift,((0.0,0.0), (0.5, 0.0), (0.25,sqrt(3)/4)))

    path = get(kwargs,:path,Snake!)
    
    Latt = CompositeLattice(Latt_A, Latt_B, Latt_C, shift) |> path
    
    return Latt
end