# Test 5th order interpolation with 5th order polynomial
function p5(x::Real)
    Float32(-0.001x^5 + 0.03x^4 - 0.3x^3 + 1.5x^2 - 2x + 5)
end

function dpdx(x::Real)
    Float32(-0.005x^4 + 0.12x^3 - 0.9x^2 + 3x - 2)
end

# Second order polynomial to test boundaries
function p2(x::Real)
    Float32(-0.5x^2 + 2x - 3)
end

@testset "interpolations" begin

    @testset "Extrapolation" begin   
        x_stagger = Float32.(1:10) .- 0.5f0
        x_stagger = reshape(x_stagger,(10,1,1))

        x = xup(x_stagger,false)
        @test x[:,1,1] == Float32.(1:10)

        x = xdn(x_stagger,false)
        @test x[:,1,1] == Float32.(0:9)
    end

    @testset "5th order exact interpolation" begin   
        x = Float32.(1:10) .- 0.5f0
        y = p5.(x)

        # Reshape to 3D array, check against inner parts

        # x direction
        p_x = reshape(y,(10,1,1))
        @test xup(p_x,false)[3:end-3,1,1] ≈ p5.(x .+ 0.5)[3:end-3]
        @test xdn(p_x,false)[4:end-2,1,1] ≈ p5.(x .- 0.5)[4:end-2]

        # y direction
        p_y = reshape(y,(1,10,1))
        @test yup(p_y,false)[1,3:end-3,1] ≈ p5.(x .+ 0.5)[3:end-3]
        @test ydn(p_y,false)[1,4:end-2,1] ≈ p5.(x .- 0.5)[4:end-2]

        # z direction
        p_z = reshape(y,(1,1,10))
        @test zup(p_z,false)[1,1,3:end-3] ≈ p5.(x .+ 0.5)[3:end-3]
        @test zdn(p_z,false)[1,1,4:end-2] ≈ p5.(x .- 0.5)[4:end-2]
    end
end

@testset "derivatives" begin

    @testset "6th order exact derivative" begin   
        x = Float32.(1:10) .- 0.5f0
        dx = ones(Float32,10)
        y = p5.(x)

        # Reshape to 3D array, check against inner parts

        # x direction
        p_x = reshape(y,(10,1,1))
        @test dxup(p_x,dx,false)[3:end-3,1,1] ≈ dpdx.(x .+ 0.5)[3:end-3]
        @test dxdn(p_x,dx,false)[4:end-2,1,1] ≈ dpdx.(x .- 0.5)[4:end-2]

        # y direction
        p_y = reshape(y,(1,10,1))
        @test dyup(p_y,dx,false)[1,3:end-3,1] ≈ dpdx.(x .+ 0.5)[3:end-3]
        @test dydn(p_y,dx,false)[1,4:end-2,1] ≈ dpdx.(x .- 0.5)[4:end-2]

        # z direction
        p_z = reshape(y,(1,1,10))
        @test dzup(p_z,dx,false)[1,1,3:end-3] ≈ dpdx.(x .+ 0.5)[3:end-3]
        @test dzdn(p_z,dx,false)[1,1,4:end-2] ≈ dpdx.(x .- 0.5)[4:end-2]
    end

end