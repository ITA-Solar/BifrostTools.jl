
function br_squeeze(A::AbstractArray)
  keepdims = Tuple(i for i in size(A) if i != 1)
  return reshape(A, keepdims)
end

function br_heatmap_xz(A::AbstractArray, M::BifrostMesh; kwargs...)
  return heatmap(M.x, M.z, br_squeeze(A)', yflip=true; kwargs...)
end

function br_heatmap_yz(A::AbstractArray, M::BifrostMesh; kwargs...)
  return heatmap(M.y, M.z, br_squeeze(A)', yflip=true; kwargs...)
end

function br_heatmap_xy(A::AbstractArray, M::BifrostMesh; kwargs...)
  return heatmap(M.x, M.y, br_squeeze(A); kwargs...)
end
