
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
        expdir::String=pwd()
        ;
        mesh_file=nothing
        )
        if expname=="none"
            expname = splitpath(expdir)[end]
        end

        filenames = readdir(expdir)
        # Find mesh-file
        if mesh_file == nothing
            mesh_file = ""
            mesh_match = false
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
        else
            mesh_match = true
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
        params_file = string(
            joinpath(expdir,expname),
            "_", lpad(snaps[1],3,"0"),
            ".idl"
        )

        params = read_params(params_file)
        snapsize, num_vars = get_snapsize_and_numvars(params)
   
        new(mesh, expname, expdir, sort(snaps), snapsize, length(snaps),
            num_vars[1]
            )
    end
end     


