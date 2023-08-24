
struct BifrostExperiment
    mesh            ::BifrostMesh
    expname         ::String
    expdir          ::String
    snaps           ::Vector{Int64}
    snapsize        ::Tuple{Int64, Int64, Int64}
    num_snaps       ::Int64
    num_primary_vars::Int64

    function BifrostExperiment(
        expname::String="none",
        expdir::String=pwd(),
        )
        if expname=="none"
            expname = basename(expdir)
        end

        filenames = readdir(expdir)
        # Find mesh-file
        mesh_match = false
        mesh_file = ""
        for filename in filenames
            match_result = match(r"^" * expname * r".*\.mesh$", filename)
            if match_result != nothing
                if mesh_match
                    error("Found multiple mesh files.")
                else
                    mesh_file *= match_result.match
                    mesh_match = true 
                end
            end
        end
        if mesh_match 
            mesh = BifrostMesh(string(expdir, "/", mesh_file))
        else
            error("Did not find mesh file with expname '$expname' in $expdir")
        end
        
        # Find number of snaps
        snaps = get_snap_numbers(expdir, expname; filenames=filenames,
                                    findall=true)

        # Get some snap-independent parameters
        params_file = string(expdir, "/", expname, "_", string(snaps[1], pad=3),
                             ".idl"
                             )
        params = br_read_params(params_file)
        snapsize, num_primary_vars = get_snapsize_and_numvars(params)
        
        new(mesh, expname, expdir, snaps, snapsize, length(snaps),
            num_primary_vars
            )
    end
end     


