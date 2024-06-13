@testset "utility functions" begin
    
    @test size(BifrostTools.squeeze(ones(1,1,10,1,1,10,1))) == (10,10)

    basename,snapname = BifrostTools.get_basename(expname,xp.snaps[1],expdir)

    @test splitpath(snapname)[end] == "en48_001"
    @test splitpath(basename)[end] == "en48"
end