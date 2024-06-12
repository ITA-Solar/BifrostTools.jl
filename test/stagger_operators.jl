@testset "interpolation" begin
    ####
    # Not all combinations of slices/columns are going to pass these tests.
    # If there is a rapid change in the variable, the 5th order Bifrost interpolation
    # and cubic spline interpolation are going to give different results.  
    ####
    @testset "z direction" begin
        # Load a column
        bz_stagger = get_var(xp,xp.snaps,"bz",
            destagger=false,slicex=[10],slicey=[10])

        # Destagger by the 5th order Bifrost Interpolation
        bz = zup(bz_stagger)
        
        # Use Intepolations cubic spline to interpolate
        x = 1:xp.mesh.mz
        x_new = x .+ 0.5
        itp = cubic_spline_interpolation(x,bz_stagger[1,1,:],extrapolation_bc=Line())
        bz_itp = itp(x_new)

        @test isapprox(bz_itp,bz[1,1,:],atol=1e-3)
    end
    
    @testset "y direction" begin
        # load a strange narrow 3x48x3 cube
        by_stagger = get_var(xp,xp.snaps,"by",
            destagger=true,slicex=1:3,slicez=1:3)

        # Destagger by the 5th order Bifrost Interpolation
        by = yup(by_stagger,true)
        
        # Use Intepolations cubic spline to interpolate
        x = 1:xp.mesh.my
        x_new = x .+ 0.5

        xx = 1:3
        yy = 1:3

        itp = cubic_spline_interpolation((xx,x,yy),by_stagger,extrapolation_bc=Line())
        by_itp = itp[xx,x_new,yy]

        @test isapprox(by_itp,by,atol=1e-3)
    end

    @testset "x direction" begin   
        x_stagger = Float32.(1:10) .- 0.5f0
        x_stagger = reshape(x_stagger,(10,1,1))

        x = xup(x_stagger,false)
        @test x[:,1,1] == Float32.(1:10)

        x = xdn(x_stagger,false)
        @test x[:,1,1] == Float32.(0:9)
    end

end