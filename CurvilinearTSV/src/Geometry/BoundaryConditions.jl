# Boundary condition and port specification types for the unified TSVSolver API

abstract type AbstractBoundaryCondition end

struct DirichletBC <: AbstractBoundaryCondition
    conductor_id::Int
    voltage::Float64
end

struct NeumannBC <: AbstractBoundaryCondition
    boundary_id::Int
    flux::Float64
end

struct DielectricInterfaceBC <: AbstractBoundaryCondition
    interface_id::Int
    eps_inner::ComplexF64
    eps_outer::ComplexF64
    surface_charge::Float64
end

struct DepletionFreeBoundaryBC <: AbstractBoundaryCondition
    tsv_id::Int
    N_A::Float64
    phi_bulk::Float64
    psi_inversion_clamp::Float64
end

struct PortSpec
    port_id::Int
    signal_tsv_ids::Vector{Int}
    ground_tsv_ids::Vector{Int}
    Z0::Float64
end
PortSpec(port_id::Int, sig_id::Int, gnd_ids::Vector{Int}; Z0::Float64=50.0) =
    PortSpec(port_id, [sig_id], gnd_ids, Z0)
