mutable struct Vertex
    pos::Array{Float64}
end

Base.copy(s::Vertex) = Vertex(copy(s.pos))

mutable struct Plane
    normal::Array{Float64}
    w::Float64
end

mutable struct Polygon
    vertices::Array{Vertex,1}
    plane::Plane
end

mutable struct Node
    plane::Union{Nothing,Plane}
    front::Union{Nothing,Node}
    back::Union{Nothing,Node}
    polygons::Array{Polygon,1}
end

mutable struct Solid
    vertices::Array{Vertex,1}
    indices::Array{Int64,1}
end
