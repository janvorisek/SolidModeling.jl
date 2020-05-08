module SolidModeling

using LinearAlgebra

const PlaneEpsilon = 1e-8

include("structs.jl")

function fromPolygons(polygons::Vector{Polygon{T}}) where T
    csg = Solid(Vector{T}(), Vector{Int}())

    p = 0
    for i = 0:(length(polygons) - 1)
        poly = polygons[i + 1]
        for j = 2:(length(poly.vertices) - 1)
            push!(csg.vertices, poly.vertices[1],
                                poly.vertices[j],
                                poly.vertices[j + 1])
            push!(csg.indices, p, p + 1, p + 2)

            p = p + 3
        end
    end
    return csg
end

function toPolygons(model::Solid{T}) where T
    list = Array{Polygon,1}()
    for i = 0:3:(length(model.indices) - 1)

        triangle = Array{T,1}()
        for j = 0:2
            v = model.vertices[model.indices[i + j + 1] + 1]
            push!(triangle, copy(v))
        end
        push!(list, Polygon{T}(triangle, fromPoints(triangle[1], triangle[2], triangle[3])))
    end
    return list
end

"""
    bunion(first::Solid, second::Solid)

Return new solid after union of `A` and `B`.

    +-------+            +-------+
    |       |            |       |
    |   A   |            |       |
    |    +--+----+   =   |       +----+
    +----+--+    |       +----+       |
         |   B   |            |       |
         |       |            |       |
         +-------+            +-------+

# Examples
```julia
c1 = cube(0.5, 0.5, 0.5, 1.0, 1.0, 1.0)
c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)
r = bunion(c1, c2) # union of solids c1 and c2

v = volume(r) # volume of union c1 c2
```
"""
function bunion(first::Solid{T}, second::Solid{T}) where T
    a = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())
    b = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())

    build(a, toPolygons(first))
    build(b, toPolygons(second))

    clipTo(a, b)
    clipTo(b, a)
    invert(b)
    clipTo(b, a)
    invert(b)

    build(a, allPolygons(b))
    return fromPolygons(allPolygons(a))
end

"""
    bsubtract(first::Solid, second::Solid)

Return new solid after performing `B`-`A`.

    +-------+            +-------+
    |       |            |       |
    |   A   |            |       |
    |    +--+----+   =   |    +--+
    +----+--+    |       +----+
         |   B   |
         |       |
         +-------+

# Examples
```julia
c1 = cube(0.5, 0.5, 0.5, 1.0, 1.0, 1.0)
c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)
r = bsubtract(c1, c2) # subtraction of c2 from c1

v = volume(r) # volume of c1-c2
```
"""
function bsubtract(first::Solid{T}, second::Solid{T}) where T
    a = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())
    b = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())

    build(a, toPolygons(first))
    build(b, toPolygons(second))

    invert(a)
    clipTo(a, b)
    clipTo(b, a)
    invert(b)
    clipTo(b, a)
    invert(b)
    build(a, allPolygons(b))
    invert(a)
    return fromPolygons(allPolygons(a))
end

"""
bintersect(first::Solid, second::Solid)

Return new solid after computing intersection of `A` and `B`.

    +-------+
    |       |
    |   A   |
    |    +--+----+   =   +--+
    +----+--+    |       +--+
         |   B   |
         |       |
         +-------+

# Examples
```julia
c1 = cube(0.5, 0.5, 0.5, 1.0, 1.0, 1.0)
c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)
r = bintersect(c1, c2) # intersection of c1 and c2

v = volume(r)
```
"""
function bintersect(first::Solid{T}, second::Solid{T}) where T
    a = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())
    b = Node{T}(nothing, nothing, nothing, Array{Polygon{T},1}())

    build(a, toPolygons(first))
    build(b, toPolygons(second))

    invert(a)
    clipTo(b, a)
    invert(b)
    clipTo(a, b)
    clipTo(b, a)
    build(a, allPolygons(b))
    invert(a)
    return fromPolygons(allPolygons(a))
end

function interpolate(vertex, other, t)
    return vertex .+ (other .- vertex).*t
end

function fromPoints(a, b, c)
    cross(b - a, c - a)
    n = normalize(cross(b - a, c - a))
    return Plane(n, dot(n, a))
end

function splitPolygon(plane::Plane{T}, polygon, coplanarFront, coplanarBack, front, back) where T
    COPLANAR = 0x00
    FRONT = 0x01
    BACK = 0x02
    SPANNING = 0x03

    # Classify each point as well as the entire polygon into one of the above four classes.
    polygonType = 0x00
    types = Vector{UInt8}(undef, length(polygon.vertices))

    for i in eachindex(polygon.vertices)
        vertex = polygon.vertices[i]
        t = dot(plane.normal, vertex) - plane.w
        type = if t < -PlaneEpsilon
            BACK
        elseif t > PlaneEpsilon
            FRONT
        else
            COPLANAR
        end
        polygonType |= type
        types[i] = type
    end

    if polygonType == COPLANAR
        push!(dot(plane.normal, polygon.plane.normal) > 0 ? coplanarFront : coplanarBack, polygon)
    elseif polygonType == FRONT
        push!(front, polygon)
    elseif polygonType == BACK
        push!(back, polygon)
    else # SPANNING
        f = Array{T,1}()
        b = Array{T,1}()

        for i = 0:(length(polygon.vertices) - 1)
            j = (i + 1) % length(polygon.vertices)
            ti = types[i + 1]
            tj = types[j + 1]
            vi = polygon.vertices[i + 1]
            vj = polygon.vertices[j + 1]

            if ti != BACK push!(f, copy(vi)) end
            if ti != FRONT push!(b, copy(vi)) end
            if (ti | tj) == SPANNING
                t = (plane.w - dot(plane.normal, vi)) / dot(plane.normal, vj .- vi)
                v = interpolate(vi, vj, t)
                push!(f, v)
                push!(b, v)
            end
        end

        if (length(f) >= 3) push!(front, Polygon{T}(f, fromPoints(f[1], f[2], f[3]))) end
        if (length(b) >= 3) push!(back, Polygon{T}(b, fromPoints(b[1], b[2], b[3]))) end
    end
end

function invert(node::Node)
    for poly in node.polygons
        reverse!(poly.vertices)
        poly.plane = Plane(-1 .* poly.plane.normal, -1 .* poly.plane.w)
    end

    if node.plane !== nothing
        node.plane = Plane(-1 .* node.plane.normal, -1 .* node.plane.w)
    end

    if node.front !== nothing invert(node.front) end
    if node.back !== nothing invert(node.back) end

    # swap front and back in place
    node.front, node.back = node.back, node.front
end

function clipPolygons(node::Node{T}, polygons::Array{Polygon{T},1}) where T
    polysize = length(polygons)

    if node.plane === nothing return copy(polygons) end

    front = Array{Polygon{T},1}()
    back = Array{Polygon{T},1}()

    for polygon in polygons
        splitPolygon(node.plane, polygon, front, back, front, back)
    end

    if node.front !== nothing
        front = clipPolygons(node.front, front)
    end

    if node.back !== nothing
        back = clipPolygons(node.back, back)
    else
        back = Array{Polygon{T},1}()
    end

    return vcat(front, back)
end

function clipTo(node::Node, bsp)
    node.polygons = clipPolygons(bsp, node.polygons)

    if node.front !== nothing clipTo(node.front, bsp) end
    if node.back !== nothing clipTo(node.back, bsp) end
end

function allPolygons(node::Node)
    polygons = deepcopy(node.polygons)

    if node.front !== nothing polygons = vcat(polygons, allPolygons(node.front)) end
    if node.back !== nothing polygons = vcat(polygons, allPolygons(node.back)) end

    return polygons
end

function build(node::Node{T}, polygons) where T

    iszero(length(polygons)) && return

    if node.plane === nothing node.plane = deepcopy(polygons[1].plane) end

    front = Array{Polygon{T},1}()
    back = Array{Polygon{T},1}()

    for polygon in polygons
        splitPolygon(node.plane, polygon, node.polygons, node.polygons, front, back)
    end

    if length(front) > 0
        if node.front === nothing node.front = Node{T}(nothing, nothing, nothing, Array{Polygon,1}()) end
        build(node.front, front)
    end

    if length(back) > 0
        if node.back === nothing node.back = Node{T}(nothing, nothing, nothing, Array{Polygon,1}()) end
        build(node.back, back)
    end
end

function cube(center::Vector{Float64}, lx::Float64, ly::Float64, lz::Float64)::Solid
    return cube(center[1] - lx / 2, center[2] - ly / 2, center[3] - lz / 2, center[1] + lx / 2, center[2] + ly / 2, center[3] + lz / 2)
end

"""
    cube(xMin::Float64, yMin::Float64, zMin::Float64, xMax::Float64, yMax::Float64, zMax::Float64)::Solid
    cube(center::Vector{Float64}, lx::Float64, ly::Float64, lz::Float64)::Solid

Creates a cube by specifying its bounding box `xMin`, `yMin`, `zMin` and `xMax`, `yMax`, `zMax`.

Or create a cube using a center point located in `center` and dimensions `lx`, `ly`, and `lz`.

# Examples
```julia-repl
julia> cube(0, 0, 0, 1, 1, 1)
Main.SolidModeling.Solid(...)

julia> cube([0.5, 0.5, 0.5], 1, 1, 1)
Main.SolidModeling.Solid(...)
```
"""
function cube(xMin::Float64, yMin::Float64, zMin::Float64, xMax::Float64, yMax::Float64, zMax::Float64)::Solid
    v1 = [xMin, yMin, zMax]
    v2 = [xMin, yMax, zMax]
    v3 = [xMax, yMax, zMax]
    v4 = [xMax, yMin, zMax]
    v5 = [xMin, yMin, zMin]
    v6 = [xMin, yMax, zMin]
    v7 = [xMax, yMax, zMin]
    v8 = [xMax, yMin, zMin]

    f1p = Polygon([v1, v4, v3, v2], fromPoints(v1, v4, v3))
    f2p = Polygon([v7, v8, v5, v6], fromPoints(v7, v8, v5))
    f3p = Polygon([v1, v2, v6, v5], fromPoints(v1, v2, v6))
    f4p = Polygon([v2, v3, v7, v6], fromPoints(v2, v3, v7))
    f5p = Polygon([v3, v4, v8, v7], fromPoints(v3, v4, v8))
    f6p = Polygon([v4, v1, v5, v8], fromPoints(v4, v1, v5))

    polys = [f1p, f2p, f3p, f4p, f5p, f6p]

    return fromPolygons(polys)
end

function signedVolumeOfTriangle(p1, p2, p3)
   	v321 = p3[1] * p2[2] * p1[3]
   	v231 = p2[1] * p3[2] * p1[3]
   	v312 = p3[1] * p1[2] * p2[3]
   	v132 = p1[1] * p3[2] * p2[3]
    v213 = p2[1] * p1[2] * p3[3]
    v123 = p1[1] * p2[2] * p3[3]

   	return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
end

"""
    volume(c::Solid)::Float64

Return the volume of **Solid** `c`.

# Examples
```julia-repl
julia> volume(cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
0.9999999999999999
```
"""
function volume(c::Solid)::Float64
    volume = 0.0

    for i = 0:3:(length(c.indices) - 1)
        dv = signedVolumeOfTriangle(c.vertices[i + 1], c.vertices[i + 2], c.vertices[i + 3])
        volume += dv
    end

    return volume
end

export cube
export bunion
export bsubtract
export bintersect
export volume

end # module
