function unknown_itp_bc_error()
    throw(ArgumentError("Uknown itp_bc. Try 'flat' or ['flat', 'periodic']"*
    "or 'Flat()' or '(Flat(), Periodic())'."))
end

function wrong_itp_bc_length_error()
    throw(ArgumentError("itp_bc is the wrong length. Should be 1 or equal to the"*
    "number of dimensions (including time)"))
end

function unknown_itp_type_error()
    throw(ArgumentError("Uknown itp_type. Try 'linear' or 'cubic-splines'."))
end

function could_not_parse_variables_error()
    throw(ArgumentError("Could not parse `variables`. Use strings or symbols, e.g"*
    "' bx by ' or '[:bx, :by]'."))
end

function nonuniform_nonlinear_not_supported_error()
    throw(ArgumentError("""Non-linear interpolation on non-uniform grid:
This is currently not implemented. However, it is possible. For example, by
transforming the coordinates to the uniform index-space before interpolation.
"""
    ))
end
