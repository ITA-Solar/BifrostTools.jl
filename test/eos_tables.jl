@testset "calculate electron density" begin
    
    tmp_array = Float32[
        3.049999f14, 3.0051822f14, 2.9655743f14, 2.9357018f14, 2.9122838f14, 
        2.8965733f14, 2.8885182f14, 2.885963f14, 2.8863922f14, 2.891092f14, 
        2.904146f14, 2.929381f14, 2.9665926f14, 3.0098403f14, 3.0582764f14, 
        3.110285f14, 3.1600604f14, 3.207982f14, 3.2572512f14, 3.3144383f14
    ]


    ne = get_electron_density(xp,xp.snaps[1];
        slicex=[10],slicey=5:24,slicez=[10],verbose=false)

    @test ne[1,:,1] == tmp_array

end