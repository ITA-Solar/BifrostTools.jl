@testset "experiment and mesh" begin
    
    @test xp.snaps == [1]
    # (nx, ny, nz)
    @test xp.snapsize == BifrostTools.get_snapsize(xp.mesh) == (48, 48, 64)
    # 8 primary, 2 aux
    @test xp.num_primary_vars == 8

end