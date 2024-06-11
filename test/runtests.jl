using Test
using BifrostTools

expdir = "../data/sp.n064/"
expname = "en48"

xp = BifrostExperiment(expname,expdir)
params = read_params(expname,xp.snaps[1],expdir)

@testset "experiment and mesh" begin
    
    @test xp.snaps == [1]
    @test xp.snapsize == (xp.mesh.mx, xp.mesh.my, xp.mesh.mz)

end

rho = get_var(xp,xp.snaps,"r")
@testset "reading variables" begin
    tmp_array = Float32[
        4.517059f-9, 4.6568345f-9, 4.8888467f-9, 5.1624958f-9, 5.4335088f-9, 
        5.5271556f-9, 5.5458482f-9, 5.5393974f-9, 5.546665f-9, 5.586004f-9, 
        5.703038f-9, 5.9570238f-9, 6.1816396f-9, 6.4527703f-9, 6.774606f-9, 
        7.119954f-9, 7.509514f-9, 7.911374f-9, 8.388659f-9, 8.831276f-9, 
        9.330005f-9, 9.9333075f-9, 1.07750155f-8, 1.1793163f-8, 1.3059435f-8, 
        1.4809497f-8, 1.7210718f-8, 2.0800567f-8, 2.6005944f-8, 3.4854107f-8, 
        5.2834523f-8, 9.635602f-8, 1.9742674f-7, 4.3695874f-7, 9.646677f-7, 
        2.2240474f-6, 4.748619f-6, 1.0463934f-5, 2.5050658f-5, 6.4452746f-5, 
        0.00016707876, 0.0004727985, 0.0013061843, 0.003768678, 0.012945924, 
        0.05415829, 0.22439207, 0.884576, 2.0215719, 2.5673914, 4.0931287, 
        6.104089, 8.702029, 12.321905, 17.355167, 23.90907, 32.994244, 
        45.009403, 61.117073, 82.6865, 111.65955, 150.39403, 202.42265, 277.7372
    ]

    @test rho[10,10,:] == tmp_array
end


@testset "unit conversions" begin
    rho_si = get_var(xp,xp.snaps,"r",units="si")
    rho_cgs = get_var(xp,xp.snaps,"r",units="cgs")

    u_r = params["u_r"]

    @test rho_cgs ≈ rho .* parse(Float32,u_r)
    @test rho_si ≈ 1_000*rho_cgs
end

@testset "derivatives and interpolations" begin
    nothing
end
