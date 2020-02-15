using CSG;
using Test;

# Test union of two cubes
function cubeUnion()
    c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0);
    c = CSG.union(c1, c2);

    return volume(c)
end


function cubeUnionReversed()
    c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0);
    c = CSG.union(c2, c1);

    return volume(c)
end

@testset begin
    @test isapprox(cubeUnion(), 4.25)                   # Test against known volume
    @test isapprox(cubeUnion(), cubeUnionReversed())    # Test if ordering doesnt matter
end


