using TestItems
using Test
using CurvilinearTSV
using CurvilinearTSV.Geometry
using CurvilinearTSV.Backends
using StaticArrays
using SpecialFunctions
using Krylov
using ArnoldiMethod

@testset "CurvilinearTSV Regression Suite" begin
    @testset "Backend dispatch" begin
        cpu_backend = CPUBackend()
        @test cpu_backend isa CPUBackend
        
        @test CPUBackend() isa AbstractBackend
        @test CUDABackend() isa AbstractBackend
        @test CuTileBackend() isa AbstractBackend
    end

    @testset "2D geometry - circle" begin
        circle = Circle2D(1.0)
        @test circle.R == 1.0
        @test circle.center == SVector(0.0, 0.0)
        
        pts = sample_circle_2d(circle, 8)
        @test length(pts) == 8
    end

    @testset "2D geometry - ellipse" begin
        ellipse = Ellipse2D(2.0, 1.0)
        @test ellipse.a == 2.0
        @test ellipse.b == 1.0
        
        pts = sample_ellipse_2d(ellipse, 8)
        @test length(pts) == 8
    end

    @testset "2D geometry - teardrop" begin
        teardrop = Teardrop2D(1.0, 0.5)
        @test teardrop.R == 1.0
        @test teardrop.h == 0.5
        
        pts = sample_teardrop_2d(teardrop, 8)
        @test length(pts) == 8
    end

    @testset "2D geometry - Fourier star" begin
        coefficients = SVector(0.1, 0.05, -0.02)
        star = FourierStar2D(1.0, coefficients)
        @test star.R == 1.0
        @test length(star.coefficients) == 3
        
        # For small perturbations (|coefficients| << R), the contour remains star-shaped (unique solution)
        # This follows from convexity and monotonicity of the radial function r(θ)
        rθs = [star.R + sum(star.coefficients[i] * cos(i * θ) for i in 1:length(star.coefficients)) for θ in LinRange(0, 2π, 8)]
        @test all(r > 0 for r in rθs)
    end

    @testset "3D geometry - sphere" begin
        sphere = Sphere3D(1.0)
        @test sphere.R == 1.0
        @test sphere.center == SVector(0.0, 0.0, 0.0)
        
        pts = sample_sphere_3d(sphere, 27)
        @test length(pts) == 27
    end

    @testset "Backend availability" begin
        @test backend_available(CPUBackend()) == true
    end
end
