@testset "unit conversions" begin
    
    params = read_params(expname,xp.snaps[1],expdir)

    rho = get_var(xp,xp.snaps,"r")
    rho_si = get_var(xp,xp.snaps,"r",units="si")
    rho_cgs = get_var(xp,xp.snaps,"r",units="cgs")

    u_r = params["u_r"]

    @test rho_cgs ≈ rho .* parse(Float32,u_r)
    @test rho_si ≈ 1_000*rho_cgs

end