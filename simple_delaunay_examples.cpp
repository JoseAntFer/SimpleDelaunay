/*
MIT License

Copyright (c) 2019 Jose Antonio Fernandez Fernandez

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include <vector>
#include "SimpleDelaunay.hpp"

void onlyPointsExample();
void pointsAndIndicesExample();
void nonVectorExample();

/*
    The same code applies for 2, 3 and 4 dimensions. The only
    difference is the template argument (2, 3 or 4, respectively)
    in the SimpleDelaunay::compute<N> function.
*/

void onlyPointsExample()
{
    // Input: vector with the point coordinates in 3 dimensions
    std::vector<double> points = { 0.75816742, 0.24371858, 0.92870883,
                                   0.12689219, 0.06034812, 0.53746581,
                                   0.88915805, 0.24796188, 0.27345906,
                                   0.91783859, 0.69470075, 0.28810121,
                                   0.23865371, 0.70646204, 0.07248404,
                                   0.55374917, 0.52939551, 0.81973793,
                                   0.81546199, 0.37344253, 0.96196011,
                                   0.61142366, 0.11634092, 0.10327177,
                                   0.92806606, 0.04172719, 0.33958352,
                                   0.62684985, 0.6717684 , 0.81939159 };

    // Call Simple Delaunay for the specific number of dimensions (<3>)
    std::vector<int> connectivity = SimpleDelaunay::compute<3>(points);

    /* The solution looks like:
        
        2, 3, 9, 6, 
        2, 3, 8, 6, 
        2, 5, 9, 6, 
        1, 2, 5, 0, 
        1, 2, 8, 0, 
        2, 6, 8, 0, 
        2, 5, 6, 0, 
        2, 4, 7, 5, 
        2, 3, 9, 5, 
        2, 3, 4, 5, 
        3, 4, 9, 5, 
        1, 4, 7, 5, 
        1, 4, 9, 5, 
        1, 2, 7, 5, 
        1, 7, 8, 2, 
        3, 4, 7, 2
    */
}

void pointsAndIndicesExample()
{
    // Input: vector with the point coordinates in 3 dimensions
    std::vector<double> points = { 0.75816742, 0.24371858, 0.92870883,
                                   0.12689219, 0.06034812, 0.53746581,
                                   0.88915805, 0.24796188, 0.27345906,
                                   0.91783859, 0.69470075, 0.28810121,
                                   0.23865371, 0.70646204, 0.07248404,
                                   0.55374917, 0.52939551, 0.81973793,
                                   0.81546199, 0.37344253, 0.96196011,
                                   0.61142366, 0.11634092, 0.10327177,
                                   0.92806606, 0.04172719, 0.33958352,
                                   0.62684985, 0.6717684 , 0.81939159 };

    // Input: vector with user defined indices different than [0...npoints - 1]
    std::vector<int> indices = {67, 23, 85, 35, 7, 643, 1, 34, 567, 8642};

    // Call Simple Delaunay for the specific number of dimensions (<3>)
    std::vector<int> connectivity = SimpleDelaunay::compute<3>(points, indices);

    /* The solution looks like:
        
        85, 643, 8642, 1, 
        23, 85, 643, 67, 
        23, 85, 567, 67, 
        1, 85, 643, 67,
        1, 85, 567, 67, 
        35, 85, 8642, 1,
        35, 85, 567, 1, 
        35, 85, 8642, 643,
        23, 34, 85, 643, 
        7, 34, 85, 643, 
        7, 35, 85, 643, 
        7, 35, 8642, 643, 
        7, 23, 34, 643, 
        7, 23, 8642, 643, 
        23, 34, 567, 85,
        7, 34, 35, 85
    */
}

void nonVectorExample()
{
    /*
        Simple Delaunay accepts any memory-contiguous point and
        index data. For non std::vector arguments, pass a pointer
        to the first position of the memory block and the number of
        memory items (as plain old C).
    */

    // Array here as non-vector container
    std::array<double, 10*3> points = { 0.75816742, 0.24371858, 0.92870883,
                                        0.12689219, 0.06034812, 0.53746581,
                                        0.88915805, 0.24796188, 0.27345906,
                                        0.91783859, 0.69470075, 0.28810121,
                                        0.23865371, 0.70646204, 0.07248404,
                                        0.55374917, 0.52939551, 0.81973793,
                                        0.81546199, 0.37344253, 0.96196011,
                                        0.61142366, 0.11634092, 0.10327177,
                                        0.92806606, 0.04172719, 0.33958352,
                                        0.62684985, 0.6717684 , 0.81939159 };

    // Array here as non-vector container
    std::array<int, 10> indices = { 67, 23, 85, 35, 7, 643, 1, 34, 567, 8642 };

    // Call Simple Delaunay with pointers and sizes
    std::vector<int> connectivity = SimpleDelaunay::compute<3>(&points[0], 10*3, &indices[0], 10);

    /* The solution looks like:
        - same as above -
    */
}
