@testset "experiment and mesh" begin
    
    @test xp.snaps == [1]
    @test xp.snapsize == (xp.mesh.mx, xp.mesh.my, xp.mesh.mz)

end