
struct EOSTables
    tabparamsf::String
    tabparamsf_root::String
    params::Dict{String,String}
    nRhoBin::Int32
    RhoAxis::Vector{Float32}
    nEiBin::Int32
    EiAxis::Vector{Float32}
    RhoEi_recl::Int32
    RhoEiRadTable_recl::Int
    nTgBin::Int32
    TgAxis::Vector{Float32}
    nNeBin::Int32
    NeAxis::Vector{Float32}
    NeTgRadTable_recl::Int
    nRadBins::Int32
    function EOSTables(tabparams::String)

        tabparamsf = normpath(tabparams)
        tabparamsf_root = dirname(tabparamsf)

        p = read_params(tabparams)

        RhoMin = log(parse(Float64, p["RhoMin"]))
        lnRhor = log(parse(Float64, p["RhoMax"])) - RhoMin
        nRhoBin = parse(Int, p["nRhoBin"])
        lnRho = [RhoMin + Float32(i - 1) / Float32(nRhoBin - 1) * lnRhor for i = 1:nRhoBin]

        EiMin = log(parse(Float64, p["EiMin"]))
        lnEir = log(parse(Float64, p["EiMax"])) - EiMin
        nEiBin = parse(Int,p["nEiBin"])
        lnEi = [EiMin + Float32(i - 1) / Float32(nEiBin - 1) * lnEir for i = 1:nEiBin]

        nRadBins = parse(Int, p["nRadBins"])
        
        RhoEi_recl = nEiBin * nRhoBin * 4
        RhoEiRadTable_recl = nEiBin * nRhoBin * nRadBins


        nTBin = -1
        lnTg = [-1]

        if haskey(p, "TMin") && haskey(p, "TMax")

            TMin = log(parse(Float64, p["TMin"]))
            lnTgr = log(parse(Float64, p["TMax"])) - TMin
            nTBin = parse(Int, p["nTBin"])
            lnTg = [TMin + Float32(i - 1) / Float32(nTBin - 1) * lnTgr for i = 1:nTBin]

        end

        nNeBin = -1
        lnNe = [-1]
        NeTgRadTable_recl = -1

        if haskey(p, "NeMin") && haskey(p, "NeMax")

            NeMin = log(parse(Float64, p["NeMin"]))
            lnNer = log(parse(Float64, p["NeMax"])) - NeMin
            nNeBin = parse(Int, p["nNeBin"])
            lnNe = [NeMin + Float32(i - 1) / Float32(nNeBin - 1) * lnNer for i = 1:nNeBin]
            NeTgRadTable_recl = nNeBin * nTBin * nRadBins * 2

        end

        new(
            tabparamsf,
            tabparamsf_root,
            p,
            nRhoBin,
            lnRho,
            nEiBin,
            lnEi,
            RhoEi_recl,
            RhoEiRadTable_recl,
            nTBin,
            lnTg,
            nNeBin,
            lnNe,
            NeTgRadTable_recl,
            nRadBins,
        )
    end
end


function get_eostable(t::EOSTables)
    f = FortranFile(
        joinpath(t.tabparamsf_root, t.params["EOSTableFile"]),
        "r",
        access="direct",
        recl=t.RhoEi_recl * 4,
    )
    var = read(f, rec=1, (Float32, (t.nEiBin, t.nRhoBin, 4)))
    return var
end


# --- interpolate from eos

function eos_interpolate(eos::EOSTables, nvar::Int)
    
    lnRho = log(parse(Float32,eos.params["RhoMax"]) / parse(Float32,eos.params["RhoMin"]))
    dlnRho = lnRho / (parse(Float32,eos.params["nRhoBin"]) - 1)

    lnEi = log(parse(Float32,eos.params["EiMax"]) / parse(Float32,eos.params["EiMin"]))
    dlnEi = lnEi / (parse(Float32,eos.params["nEiBin"]) - 1)

    eia = eos.EiAxis[1]:dlnEi:eos.EiAxis[end]
    rhoa = eos.RhoAxis[1]:dlnRho:eos.RhoAxis[end]

    tab = get_eostable(eos)

    return cubic_spline_interpolation((eia, rhoa), tab[:, :, nvar], extrapolation_bc=Line())
end
