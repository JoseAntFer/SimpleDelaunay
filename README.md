# SimpleDelaunay
**SimpleDelaunay** is a single-file header-only C++ library to compute 2D, 3D and 4D Delaunay Triangulations.

Multidimensional Delaunay triangulation generators can be found in well known libraries such as [qhull](http://www.qhull.org/) or [CGAL](https://doc.cgal.org/latest/Triangulation/index.html). However, these libraries can be either cumbersome to use or they are undesired dependencies to have. In constrast, **SimpleDelaunay** is designed to be a plug&play and dependency-free solution.

It generates the triangulation by computing the lower half convex hull of the input point cloud with an added dimension corresponding to the norm of the points coordinates. The convex hull is computed using a variant of the famous [quickhull algorithm](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.117.405).

## A Quick Example
```cpp
#include <vector>
#include "SimpleDelaunay.hpp"

int main()
{
    // Input: vector with the point coordinates in 2 dimensions
    std::vector<double> points = { 0.75816742, 0.24371858,
                                   0.12689219, 0.06034812,
                                   0.88915805, 0.24796188,
                                   0.91783859, 0.69470075,
                                   0.23865371, 0.70646204 };
                                   
    // Call SimpleDelaunay for the specific number of dimensions (2)
    std::vector<int> connectivity = SimpleDelaunay::compute<2>(points);

    return 0;
}
```
`SimpleDelaunay::compute<N>` can also be called with custom indices for the points and for arbitrary (non `std::vector`) containers as long as the input data is contiguous in memory (see examples).

## License
**SimpleDelaunay** is an open source project under the MIT license.
