
# Converts from simulation units to si units
const si_params = Dict(
        # pressure      g/cm/s   --> kg/m/s
        "u_p"  => params["u_p"]*0.1,
        # gas dentisy   g/cm^3   --> kg/m^3
        "u_r"  => params["u_r"]*1e3,
        # momentum      g/cm^2/s --> kg/m^2/s
        "u_m"  => params["u_r"]*params["u_u"]*10,
        # b field
        "u_b"  => params["u_b"]*1e-4
        )