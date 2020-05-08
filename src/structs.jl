
struct Plane{T}
    normal::T
    w::Float64
end

mutable struct Polygon{T}
    vertices::Array{T,1}
    plane::Plane
end

mutable struct Node{T}
    plane::Union{Nothing,Plane{T}}
    front::Union{Nothing,Node{T}}
    back::Union{Nothing,Node{T}}
    polygons::Array{Polygon{T},1}
end

struct Solid{T}
    vertices::Array{T,1}
    indices::Array{Int,1}
end
