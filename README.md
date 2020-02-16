# CSG.jl

**Constructive Solid Geometry (CSG)** is a modeling technique that uses boolean operations like union and intersection to combine 3D solids. This library implements CSG operations on using BSP trees.

It's ported from the excellent javascript library [csg.js](https://evanw.github.io/csg.js/docs/) with some added features (like volume calculation).

## Installation

You can install this package using the following Julia command:

```julia
Pkg.add("CSG")
```

The package can then be loaded and used in a Julia script or a Jupyter Notebook by:

```julia
using CSG
```

## Usage

You can perform 3 basic operations on two solids - `intersection`, `subtraction` and `union`.

### Union of two cubes

```julia
# 1x1x1 cube with center point in [0.5, 0.5, 0.5]
c1 = cube(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

# 1.5x1.5x1.5 cube with center point in [1.25, 1.25, 1.25]
c2 = cube(0.5, 0.5, 0.5, 2.0, 2.0, 2.0)

# calculate union of the two cubes
c = CSG.union(c1, c2)

# calculate volume of the union if needed
vol = volume(c) 
```

### Intersection of two cubes

```julia
# calculate intersection of the two cubes from above
c = CSG.intersect(c1, c2)

# calculate volume of the union if needed
vol = volume(c) 
```

### Subtraction of two cubes

```julia
# calculate intersection of the two cubes from above (c1-c2)
c = CSG.subtract(c1, c2)

# calculate volume of the union if needed
vol = volume(c) 
```

To see more about CSG operations, see the [csg.js docs](https://evanw.github.io/csg.js/docs/).

## Authors

* **Jan Vorisek** <[**jan@vorisek.me**](mailto:jan@vorisek.me)>

Original javascript code written by:

* [**Evan Wallace**](https://github.com/evanw) - [csg.js](https://evanw.github.io/csg.js/docs/)

## License

This project is licensed under the [MIT License](LICENSE.md).
