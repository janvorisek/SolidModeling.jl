using SolidModeling;
using Test;

# Test union of two cubes
function cubeUnion()
    c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0);
    c = bunion(c1, c2);

    return volume(c)
end

function cubeUnionReversed()
    c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0);
    c = bunion(c2, c1);

    return volume(c)
end

# Test subtraction of two cubes
function cubeSubtract()
    c1 = cube([0.5, 0.5, 0.5], 1.0, 1.0, 1.0)
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)
    r = bsubtract(c1, c2)
    @show r

    v = volume(r)

    return v
end

# Test intersection of two cubes
function cubeIntersect()
    c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)
    @show bintersect(c1, c2)
    v = volume(bintersect(c1, c2))

    return v
end

@testset begin
    @test isapprox(cubeUnion(), 4.25)                   # Test against known volume
    @test isapprox(cubeUnion(), cubeUnionReversed())    # Test if ordering doesnt matter
    @test isapprox(cubeSubtract(), (1-0.125))
    @test isapprox(cubeIntersect(), 0.125)
end
