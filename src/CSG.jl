module CSG

using Vec
using LinearAlgebra

const PlaneEpsilon = 1e-8

include("structs.jl")

function fromPolygons(polygons)
    csg = Solid(Array{Vertex,1}(), Array{Int64,1}())

    p::Int64 = 0
    for i = 0:(length(polygons) - 1)
        poly = polygons[i + 1]
        for j = 2:(length(poly.vertices) - 1)
            push!(csg.vertices, poly.vertices[1])
            push!(csg.indices, p)
            push!(csg.vertices, poly.vertices[j])
            push!(csg.indices, p + 1)
            push!(csg.vertices, poly.vertices[j + 1])
            push!(csg.indices, p + 2)

            p = p + 3
        end
    end
    return csg
end

function toPolygons(model)
    list = Array{Polygon,1}();
    for i = 0:3:(length(model.indices) - 1)

        triangle = Array{Vertex,1}()
        for j = 0:2
            v = model.vertices[model.indices[i + j + 1] + 1]
            push!(triangle, copy(v));
        end
        push!(list, Polygon(triangle, fromPoints(triangle[1].pos, triangle[2].pos, triangle[3].pos)))
    end
    return list;
end

function union(first, second)
    a = Node(nothing, nothing, nothing, Array{Polygon,1}())
    b = Node(nothing, nothing, nothing, Array{Polygon,1}())

    build(a, toPolygons(first));
    build(b, toPolygons(second));

    clipTo(a, b)
    clipTo(b, a)
    invert(b)
    clipTo(b, a)
    invert(b)

    build(a, allPolygons(b))
    return fromPolygons(allPolygons(a))
end

function subtract(first, second)
    a = Node(nothing, nothing, nothing, Array{Polygon,1}())
    b = Node(nothing, nothing, nothing, Array{Polygon,1}())

    build(a, toPolygons(first));
    build(b, toPolygons(second));

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

function intersect(first, second)
    a = Node(nothing, nothing, nothing, Array{Polygon,1}())
    b = Node(nothing, nothing, nothing, Array{Polygon,1}())

    build(a, toPolygons(first));
    build(b, toPolygons(second));

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
    return Vertex(lerp(vertex.pos, other.pos, t))
end

function fromPoints(a, b, c)
    cross(b - a, c - a)
    n = normalize(cross(b - a, c - a))
    return Plane(n, dot(n, a))
end

function splitPolygon(plane, polygon, coplanarFront, coplanarBack, front, back)
    COPLANAR = 0
    FRONT = 1
    BACK = 2
    SPANNING = 3

    # Classify each point as well as the entire polygon into one of the above four classes.
    polygonType = 0
    types = Array{Int64,1}()

    for vertex in polygon.vertices
        t = dot(plane.normal, vertex.pos) - plane.w
        type::Int64 = (t < -PlaneEpsilon) ? BACK : ((t > PlaneEpsilon) ? FRONT : COPLANAR)
        polygonType |= type
        push!(types, type)
    end

    if polygonType === COPLANAR
        push!(dot(plane.normal, polygon.plane.normal) > 0 ? coplanarFront : coplanarBack, polygon)
    elseif polygonType === FRONT
        push!(front, polygon)
    elseif polygonType === BACK
        push!(back, polygon)
    elseif polygonType === SPANNING
        f = Array{Vertex,1}()
        b = Array{Vertex,1}()

        for i = 0:(length(polygon.vertices) - 1)
            j = (i + 1) % length(polygon.vertices)
            ti = types[i + 1]
            tj = types[j + 1]
            vi = polygon.vertices[i + 1]
            vj = polygon.vertices[j + 1]

            if ti != BACK push!(f, copy(vi)) end
            if ti != FRONT push!(b, copy(vi)) end
            if (ti | tj) == SPANNING
                t = (plane.w - dot(plane.normal, vi.pos)) / dot(plane.normal, vj.pos - vi.pos)
                v = interpolate(vi, vj, t)
                push!(f, v)
                push!(b, v)
            end
        end

        if (length(f) >= 3) push!(front, Polygon(f, fromPoints(f[1].pos, f[2].pos, f[3].pos))) end
        if (length(b) >= 3) push!(back, Polygon(b, fromPoints(b[1].pos, b[2].pos, b[3].pos))) end
    end
end



function invert(node::Node)
    for poly in node.polygons
        reverse!(poly.vertices)
        poly.plane.normal = -1 * poly.plane.normal
        poly.plane.w = -1 * poly.plane.w
    end

    if node.plane !== nothing
        node.plane.normal = -1 * node.plane.normal
        node.plane.w = -1 * node.plane.w
    end

    if node.front !== nothing invert(node.front) end
    if node.back !== nothing invert(node.back) end

    # swap front and back in place
    node.front, node.back = node.back, node.front
end

function clipPolygons(node::Node, polygons::Array{Polygon,1})
    polysize = length(polygons)

    if node.plane === nothing return copy(polygons) end

    front = Array{Polygon,1}()
    back = Array{Polygon,1}()

    for polygon in polygons
        splitPolygon(node.plane, polygon, front, back, front, back);
    end

    if node.front !== nothing
        front = clipPolygons(node.front, front)
    end

    if node.back !== nothing
        back = clipPolygons(node.back, back)
    else
        back = Array{Polygon,1}()
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

function build(node::Node, polygons)
    if length(polygons) === 0 return end

    if node.plane === nothing node.plane = deepcopy(polygons[1].plane) end

    front = Array{Polygon,1}()
    back = Array{Polygon,1}()

    for polygon in polygons
        splitPolygon(node.plane, polygon, node.polygons, node.polygons, front, back)
    end

    if length(front) > 0
        if node.front === nothing node.front = Node(nothing, nothing, nothing, Array{Polygon,1}()) end
        build(node.front, front)
    end

    if length(back) > 0
        if node.back === nothing node.back = Node(nothing, nothing, nothing, Array{Polygon,1}()) end
        build(node.back, back)
    end
end

function cube(xMin::Float64, yMin::Float64, zMin::Float64, xMax::Float64, yMax::Float64, zMax::Float64)::Solid
    v1 = Vertex(VecE3(xMin, yMin, zMax))
    v2 = Vertex(VecE3(xMin, yMax, zMax))
    v3 = Vertex(VecE3(xMax, yMax, zMax))
    v4 = Vertex(VecE3(xMax, yMin, zMax))
    v5 = Vertex(VecE3(xMin, yMin, zMin))
    v6 = Vertex(VecE3(xMin, yMax, zMin))
    v7 = Vertex(VecE3(xMax, yMax, zMin))
    v8 = Vertex(VecE3(xMax, yMin, zMin))

    f1p = Polygon([v1, v4, v3, v2], fromPoints(v1.pos, v4.pos, v3.pos))
    f2p = Polygon([v7, v8, v5, v6], fromPoints(v7.pos, v8.pos, v5.pos))
    f3p = Polygon([v1, v2, v6, v5], fromPoints(v1.pos, v2.pos, v6.pos))
    f4p = Polygon([v2, v3, v7, v6], fromPoints(v2.pos, v3.pos, v7.pos))
    f5p = Polygon([v3, v4, v8, v7], fromPoints(v3.pos, v4.pos, v8.pos))
    f6p = Polygon([v4, v1, v5, v8], fromPoints(v4.pos, v1.pos, v5.pos))

    polys = [f1p, f2p, f3p, f4p, f5p, f6p];

    return fromPolygons(polys);
end

function signedVolumeOfTriangle(p1::VecE3, p2::VecE3, p3::VecE3)
   	v321 = p3.x * p2.y * p1.z
   	v231 = p2.x * p3.y * p1.z
   	v312 = p3.x * p1.y * p2.z
   	v132 = p1.x * p3.y * p2.z
    v213 = p2.x * p1.y * p3.z
    v123 = p1.x * p2.y * p3.z

   	return (1.0 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)
end

function volume(c::Solid)::Float64
    volume = 0.0;

    for i = 0:3:(length(c.indices) - 1)
        dv = signedVolumeOfTriangle(c.vertices[i + 1].pos, c.vertices[i + 2].pos, c.vertices[i + 3].pos)
        volume += dv
    end

    return volume
end

export cube;
export union;
export subtract;
export intersect;
export volume;

end # module
