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
#ifndef SIMPLE_DELAUNAY
#define SIMPLE_DELAUNAY
#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

namespace SimpleDelaunay
{
// Global parameters
static double TOL = 1e-15;

/* =========== Auxiliar vector utilities =========== */
// Tiny templated vector class
template<int N>
class Vec
{
private:
    std::array<double, N> coeffs;

public:
    // Constructors & Destructor
    Vec() {};

    template<typename... Args>
    Vec(Args... args) : coeffs{ double(args)... }
    {
        static_assert(sizeof...(Args) == N, "Invalid number of constructor arguments.");
    }
    Vec(std::array<double, N> &coeffs_in) : coeffs(coeffs_in) {}

    // Simple operations
    double & operator[](const int pos) { return this->coeffs[pos]; }
    const double & operator[](const int pos) const { return this->coeffs[pos]; }
    Vec<N> operator*(const double val) 
    {
        Vec<N> new_vector = Vec<N>();
        for (int i = 0; i < N; i++)
            new_vector[i] = coeffs[i]*val;
        return new_vector;
    }
    void operator*=(const double val)
    {
        for (int i = 0; i < N; i++)
            coeffs[i] *= val;
    }
    Vec<N> operator/(const double val)
    {
        return (*this)*(1.0/val);
    }
    void operator/=(const double val)
    {
        for (int i = 0; i < N; i++)
            coeffs[i] /= val;
    }
    Vec<N> operator+(const Vec<N> &other)
    {
        Vec<N> new_vector = Vec<N>();
        for (int i = 0; i < N; i++)
            new_vector[i] = coeffs[i] + other[i];
        return new_vector;
    }
    void operator+=(const Vec<N> &other)
    {
        for (int i = 0; i < N; i++)
            coeffs[i] += other[i];
    }
    Vec<N> operator-(const Vec<N> &other)
    {
        Vec<N> new_vector = Vec<N>();
        for (int i = 0; i < N; i++)
            new_vector[i] = coeffs[i] - other[i];
        return new_vector;
    }
    void operator-=(const Vec<N> &other)
    {
        for (int i = 0; i < N; i++)
            coeffs[i] -= other[i];
    }
    double squaredNorm()
    {
        double norm = 0.0;
        for (int i = 0; i < N; i++)
            norm += coeffs[i] * coeffs[i];
        return norm;
    }
    double norm()
    {
        double norm_sq = this->squaredNorm();
        assert(norm_sq >= 0.0);
        return std::sqrt(norm_sq);
    }
    Vec<N> normalized()
    {
        double norm = this->norm();
        assert(norm > TOL);
        return (*this) / norm;
    }
    double dot(const Vec<N> &other)
    {
        double result = 0.0;
        for (int i = 0; i < N; i++)
            result += coeffs[i] * other[i];
        return result;
    }

    // Static methods
    static Vec<N> Zeros()
    {
        Vec<N> ret;
        for (int i = 0; i < N; i++) {
            ret[i] = 0.0;
        }
        return ret;
    }

    // Output
    void print()
    {
        std::cout << "[ ";
        for (int i = 0; i < N - 1; i++)
            std::cout << coeffs[i] << ", ";
        std::cout << coeffs[N - 1] << " ]\n";
    }
};

// (unrolled) vector operations for 2, 3, 4 and 5 dimensions
//// Normals
inline Vec<3> nonUnitaryNormal3(Vec<3> &p0, Vec<3> &p1, Vec<3> &p2)
{
    Vec<3> a = p1 - p0;
    Vec<3> b = p2 - p0;

    double n0 = a[1] * b[2] - a[2] * b[1];
    double n1 = -a[0] * b[2] + a[2] * b[0];
    double n2 = a[0] * b[1] - a[1] * b[0];

    return Vec<3>(n0, n1, n2); // .normalized();
}
inline Vec<4> nonUnitaryNormal4(Vec<4> &p0, Vec<4> &p1, Vec<4> &p2, Vec<4> &p3)
{
    Vec<4> a = p1 - p0;
    Vec<4> b = p2 - p0;
    Vec<4> c = p3 - p0;

    double x0  = b[3] * c[2];
    double x1  = a[2] * c[3];
    double x2  = a[3] * b[2];
    double x3  = b[2] * c[3];
    double x4  = a[2] * b[3];
    double x5  = a[3] * c[2];
    double x6  = a[0] * c[1];
    double x7  = a[1] * b[0];
    double x8  = b[1] * c[0];
    double x9  = a[0] * b[1];
    double x10 = a[1] * c[0];
    double x11 = b[0] * c[1];

    double n0 =  a[1] * x0  - a[1] * x3 + b[1] * x1  - b[1] * x5 + c[1] * x2 - c[1] * x4;
    double n1 = -a[0] * x0  + a[0] * x3 - b[0] * x1  + b[0] * x5 - c[0] * x2 + c[0] * x4;
    double n2 = -a[3] * x11 + a[3] * x8 - b[3] * x10 + b[3] * x6 + c[3] * x7 - c[3] * x9;
    double n3 =  a[2] * x11 - a[2] * x8 + b[2] * x10 - b[2] * x6 - c[2] * x7 + c[2] * x9;

    return Vec<4>(n0, n1, n2, n3); // .normalized();
}
inline Vec<5> nonUnitaryNormal5(Vec<5> &p0, Vec<5> &p1, Vec<5> &p2, Vec<5> &p3, Vec<5> &p4)
{
    Vec<5> a = p1 - p0;
    Vec<5> b = p2 - p0;
    Vec<5> c = p3 - p0;
    Vec<5> d = p4 - p0;

    double x0  = a[1] * b[2];
    double x1  = c[3] * d[4];
    double x2  = a[1] * b[3];
    double x3  = c[4] * d[2];
    double x4  = a[1] * b[4];
    double x5  = c[2] * d[3];
    double x6  = c[4] * d[3];
    double x7  = a[2] * b[1];
    double x8  = a[2] * b[3];
    double x9  = c[1] * d[4];
    double x10 = a[2] * b[4];
    double x11 = c[3] * d[1];
    double x12 = c[2] * d[4];
    double x13 = a[3] * b[1];
    double x14 = c[4] * d[1];
    double x15 = a[3] * b[2];
    double x16 = a[3] * b[4];
    double x17 = c[1] * d[2];
    double x18 = c[3] * d[2];
    double x19 = a[4] * b[1];
    double x20 = c[1] * d[3];
    double x21 = a[4] * b[2];
    double x22 = c[2] * d[1];
    double x23 = a[4] * b[3];
    double x24 = a[0] * b[2];
    double x25 = a[0] * b[3];
    double x26 = a[0] * b[4];
    double x27 = a[2] * b[0];
    double x28 = c[4] * d[0];
    double x29 = c[0] * d[3];
    double x30 = a[3] * b[0];
    double x31 = c[0] * d[4];
    double x32 = c[2] * d[0];
    double x33 = a[4] * b[0];
    double x34 = c[3] * d[0];
    double x35 = c[0] * d[2];
    double x36 = a[0] * b[1];
    double x37 = a[1] * b[0];
    double x38 = c[0] * d[1];
    double x39 = c[1] * d[0];

    double n0 =  x0 * x1  - x0 * x6  - x1  * x7  + x10 * x11 - x10 * x20 - x11 * x21 + x12 * x13 - x12 * x2  - x13 * x3 + x14 * x15  - x14 * x8  - x15 * x9  + x16 * x17 - x16 * x22 - x17 * x23 + x18 * x19 - x18 * x4  - x19 * x5  + x2  * x3 + x20 * x21 + x22 * x23 + x4  * x5  + x6  * x7 + x8  * x9;
    double n1 = -x1 * x24 + x1 * x27 + x10 * x29 - x10 * x34 + x12 * x25 - x12 * x30 - x15 * x28 + x15 * x31 + x16 * x32 - x16 * x35 + x18 * x26 - x18 * x33 - x21 * x29 + x21 * x34 - x23 * x32 + x23 * x35 + x24 * x6  - x25 * x3  - x26 * x5 - x27 * x6  + x28 * x8  + x3  * x30 - x31 * x8 + x33 * x5;
    double n2 =  x1 * x36 - x1 * x37 - x11 * x26 + x11 * x33 + x13 * x28 - x13 * x31 + x14 * x25 - x14 * x30 + x16 * x38 - x16 * x39 + x19 * x29 - x19 * x34 - x2  * x28 + x2  * x31 + x20 * x26 - x20 * x33 - x23 * x38 + x23 * x39 - x25 * x9 - x29 * x4  + x30 * x9  + x34 * x4  - x36 * x6 + x37 * x6;
    double n3 =  x0 * x28 - x0 * x31 - x10 * x38 + x10 * x39 - x12 * x36 + x12 * x37 - x14 * x24 + x14 * x27 - x17 * x26 + x17 * x33 + x19 * x32 - x19 * x35 + x21 * x38 - x21 * x39 + x22 * x26 - x22 * x33 + x24 * x9  - x27 * x9  - x28 * x7 + x3  * x36 - x3  * x37 + x31 * x7  - x32 * x4 + x35 * x4;
    double n4 =  x0 * x29 - x0 * x34 + x11 * x24 - x11 * x27 - x13 * x32 + x13 * x35 - x15 * x38 + x15 * x39 + x17 * x25 - x17 * x30 - x18 * x36 + x18 * x37 + x2  * x32 - x2  * x35 - x20 * x24 + x20 * x27 - x22 * x25 + x22 * x30 - x29 * x7 + x34 * x7  + x36 * x5  - x37 * x5  + x38 * x8 - x39 * x8;

    return Vec<5>(n0, n1, n2, n3, n4); // .normalized();
}
//// Simplex volumes (projected to the last dimension = 0.0)
inline double projectedVolume2(Vec<3> &p0, Vec<3> &p1, Vec<3> &p2)
{
    Vec<3> a = p1 - p0;
    Vec<3> b = p2 - p0;

    double det = a[0] * b[1] - a[1] * b[0];
    return std::abs(det/2.0);
}
inline double projectedVolume3(Vec<4> &p0, Vec<4> &p1, Vec<4> &p2, Vec<4> &p3)
{
    Vec<4> a = p1 - p0;
    Vec<4> b = p2 - p0;
    Vec<4> c = p3 - p0;

    double det = a[0] * b[1] * c[2] - a[0] * b[2] * c[1] - a[1] * b[0] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1] - a[2] * b[1] * c[0];
    return std::abs(det/6.0);
}
inline double projectedVolume4(Vec<5> &p0, Vec<5> &p1, Vec<5> &p2, Vec<5> &p3, Vec<5> &p4)
{
    Vec<5> a = p1 - p0;
    Vec<5> b = p2 - p0;
    Vec<5> c = p3 - p0;
    Vec<5> d = p4 - p0;

    double x0  = a[0] * b[1];
    double x1  = c[2] * d[3];
    double x2  = a[0] * b[2];
    double x3  = c[3] * d[1];
    double x4  = a[0] * b[3];
    double x5  = c[1] * d[2];
    double x6  = c[3] * d[2];
    double x7  = a[1] * b[0];
    double x8  = a[1] * b[2];
    double x9  = c[0] * d[3];
    double x10 = a[1] * b[3];
    double x11 = c[2] * d[0];
    double x12 = c[1] * d[3];
    double x13 = a[2] * b[0];
    double x14 = c[3] * d[0];
    double x15 = a[2] * b[1];
    double x16 = a[2] * b[3];
    double x17 = c[0] * d[1];
    double x18 = c[2] * d[1];
    double x19 = a[3] * b[0];
    double x20 = c[0] * d[2];
    double x21 = a[3] * b[1];
    double x22 = c[1] * d[0];
    double x23 = a[3] * b[2];

    double det = x0 * x1 - x0 * x6 - x1 * x7 + x10 * x11 - x10 * x20 - x11 * x21 + x12 * x13 - x12 * x2 - x13 * x3 + x14 * x15 - x14 * x8 - x15 * x9 + x16 * x17 - x16 * x22 - x17 * x23 + x18 * x19 - x18 * x4 - x19 * x5 + x2 * x3 + x20 * x21 + x22 * x23 + x4 * x5 + x6 * x7 + x8 * x9;
    return std::abs(det/24.0);
}
/* =========== Auxiliar vector utilities =========== */


/* =========== Delaunay Types =========== */
template<int N>
struct Node
{
    Vec<N> p;
    int idx;
    bool active = true;
    Node(Vec<N> &p_in, int idx_in) : p(p_in), idx(idx_in), active(true) {}
};

template<int N>
class FrontEdge
{
public:
    std::array<Node<N>*, N - 1> nodes;

    FrontEdge() {};
    FrontEdge(std::array<Node<N>*, N - 1> &nodes_in) 
    {
        // Store the indices sorted
        int outer_min = -99999;
        for (int outer = 0; outer < N - 1; outer++) {
            int inner_min = 99999;
            int inner_min_node_i;
            for (int inner = 0; inner < N - 1; inner++) {
                if ((nodes_in[inner]->idx > outer_min) && (nodes_in[inner]->idx < inner_min)) {
                    inner_min = nodes_in[inner]->idx;
                    inner_min_node_i = inner;
                }
            }
            this->nodes[outer] = nodes_in[inner_min_node_i];
            outer_min = inner_min;
        }
    };
    ~FrontEdge() {};

    bool operator==(FrontEdge<N> &other)
    {
        for (int i = 0; i < N - 1; i++) {
            if (nodes[i]->idx != other.nodes[i]->idx) {
                return false;
            }
        }
        return true;
    }
};

template<int N>
class FrontFace
{
public:
    // Members
    std::array<Node<N>*, N> nodes;
    Vec<N> normal;
    bool active = true;

    // Methods
    //// Constructors and Destructor
    FrontFace(std::array<Node<3>*, 3> &nodes, Vec<3> &inside_point)
        : nodes(nodes)
    {
        this->normal = nonUnitaryNormal3(nodes[0]->p, nodes[1]->p, nodes[2]->p);
        this->correctNormalOrientation(inside_point);
    }
    FrontFace(std::array<Node<4>*, 4> &nodes, Vec<4> &inside_point)
        : nodes(nodes)
    {
        this->normal = nonUnitaryNormal4(nodes[0]->p, nodes[1]->p, nodes[2]->p, nodes[3]->p);
        this->correctNormalOrientation(inside_point);
    }
    FrontFace(std::array<Node<5>*, 5> &nodes, Vec<5> &inside_point)
        : nodes(nodes)
    {
        this->normal = nonUnitaryNormal5(nodes[0]->p, nodes[1]->p, nodes[2]->p, nodes[3]->p, nodes[4]->p);
        this->correctNormalOrientation(inside_point);
    }
    FrontFace(std::array<Node<N>*, N> &nodes, Vec<N> &normal, bool dummy)
        : nodes(nodes)
    {
        this->normal = normal;
    }
    ~FrontFace() {}

    //// Utilities
    double sdf(Vec<N> point)
    {
        return (point - nodes[0]->p).dot(normal);
    }
    void correctNormalOrientation(Vec<N> &inside_point)
    {
        if (this->hasVision(inside_point)) {
            normal = normal * (-1);
        }
    }
    bool hasVision(Vec<N> &point) 
    {
        return this->sdf(point) > TOL;
    }

    //// Templated cases
    void appendEdges(std::vector<FrontEdge<3>> &edges)
    {
        std::array<Node<3>*, 2> edge1 = { nodes[0], nodes[1] };
        std::array<Node<3>*, 2> edge2 = { nodes[1], nodes[2] };
        std::array<Node<3>*, 2> edge3 = { nodes[2], nodes[0] };
        edges.push_back(FrontEdge<N>(edge1));
        edges.push_back(FrontEdge<N>(edge2));
        edges.push_back(FrontEdge<N>(edge3));
    }
    void appendEdges(std::vector<FrontEdge<4>> &edges)
    {
        std::array<Node<4>*, 3> edge1 = { nodes[0], nodes[1], nodes[2] };
        std::array<Node<4>*, 3> edge2 = { nodes[0], nodes[1], nodes[3] };
        std::array<Node<4>*, 3> edge3 = { nodes[0], nodes[2], nodes[3] };
        std::array<Node<4>*, 3> edge4 = { nodes[1], nodes[2], nodes[3] };
        edges.push_back(FrontEdge<N>(edge1));
        edges.push_back(FrontEdge<N>(edge2));
        edges.push_back(FrontEdge<N>(edge3));
        edges.push_back(FrontEdge<N>(edge4));
    }
    void appendEdges(std::vector<FrontEdge<5>> &edges)
    {
        std::array<Node<5>*, 4> edge1 = { nodes[0], nodes[1], nodes[2], nodes[3] };
        std::array<Node<5>*, 4> edge2 = { nodes[0], nodes[1], nodes[2], nodes[4] };
        std::array<Node<5>*, 4> edge3 = { nodes[0], nodes[1], nodes[3], nodes[4] };
        std::array<Node<5>*, 4> edge4 = { nodes[0], nodes[2], nodes[3], nodes[4] };
        std::array<Node<5>*, 4> edge5 = { nodes[1], nodes[2], nodes[3], nodes[4] };
        edges.push_back(FrontEdge<N>(edge1));
        edges.push_back(FrontEdge<N>(edge2));
        edges.push_back(FrontEdge<N>(edge3));
        edges.push_back(FrontEdge<N>(edge4));
        edges.push_back(FrontEdge<N>(edge5));
    }
    double getProjectedVolume(FrontFace<3> *myself)  // Pass myself for template specialization
    {
        return projectedVolume2(nodes[0]->p, nodes[1]->p, nodes[2]->p);
    }
    double getProjectedVolume(FrontFace<4> *myself)
    {
        return projectedVolume3(nodes[0]->p, nodes[1]->p, nodes[2]->p, nodes[3]->p);
    }
    double getProjectedVolume(FrontFace<5> *myself)
    {
        return projectedVolume4(nodes[0]->p, nodes[1]->p, nodes[2]->p, nodes[3]->p, nodes[4]->p);
    }
};
/* =========== end Delaunay Types =========== */


/* =========== Auxiliar Functions =========== */
/*
Classify the edges in a KDTree and find the non-repeated
ones in the leaves.
*/
template<int N>
inline void findNonRepeatedEdges(std::vector<FrontEdge<N>*> &edges, std::vector<FrontEdge<N>*> &out, int depth)
{
    const unsigned int cap = 25;
    const unsigned int max_depth = 15;

    // Leaf
    int n_edges = (int)edges.size();
    if ((n_edges < cap) || (depth > max_depth)) {
        std::vector<bool> candidate(n_edges, true);
        for (int edge_i = 0; edge_i < n_edges; edge_i++) {
            if (!candidate[edge_i]) {
                continue;
            }
            candidate[edge_i] = false;

            bool repeated = false;
            for (int edge_j = 0; edge_j < n_edges; edge_j++) {
                if (!candidate[edge_j]) {
                    continue;
                }
                if (edge_i != edge_j) {
                    if (*edges[edge_i] == *edges[edge_j]) {
                        repeated = true;
                        candidate[edge_j] = false;
                        break;  // Not sure why it works with this on
                    }
                }
            }

            if (!repeated) {
                out.push_back(edges[edge_i]);
            }
        }
    }

    // Keep splitting
    else {
        int dim = depth % (N - 1);
        int random_pivot = ((n_edges - N) * (n_edges + dim)) % n_edges;
        int pivot = edges[random_pivot]->nodes[dim]->idx;

        std::vector<FrontEdge<N>*> left_edges, rght_edges;
        left_edges.reserve(edges.size());
        rght_edges.reserve(edges.size());

        for (FrontEdge<N> *edge_i : edges) {
            if (edge_i->nodes[dim]->idx < pivot) {
                left_edges.push_back(edge_i);
            }
            else {
                rght_edges.push_back(edge_i);
            }
        }

        // Spawn children
        findNonRepeatedEdges<N>(left_edges, out, depth + 1);
        findNonRepeatedEdges<N>(rght_edges, out, depth + 1);
    }
}
/* =========== End Auxiliar Functions =========== */

/* =========== Main Function =========== */
/*
Delaunay triangulation entry point. Takes pointers to first element in points
and indices arrays so it support any contiguous memory data type 
(ie. vectors, arrays...)
*/
template<int n>
inline std::vector<int> compute(
    const double *points, const unsigned int points_size,
    const int *indices, const unsigned int indices_size)
{
    // Function variables
    const int N = n + 1;
    int n_nodes = points_size / n;
    Vec<N> internal_point;
    std::vector<Node<N>*> nodes(n_nodes);
    std::vector<FrontFace<N>*> active_fronts;
    std::vector<FrontFace<N>*> final_fronts;
    std::vector<int> out_connectivity;

    // Input checking
    {
        if (n_nodes < N) {
            std::cout << "Delaunay Error: Too few points for a Delaunay triangulation." << std::endl;
            abort();
        }

        if (n*n_nodes != points_size) {
            std::cout << "Delaunay Error: Number of points coordinates do not match with number of dimensions." << std::endl;
            abort();
        }

        if (n_nodes != indices_size) {
            std::cout << "Delaunay Error: Number of points coordinates do not match number of indices." << std::endl;
            abort();
        }
    }

    // Delaunay Nodes creation
    {
        // Node creation
        for (int node_i = 0; node_i < n_nodes; node_i++) {

            // Create the point in ndim + 1 dimension
            double height = 0.0;
            std::array<double, N> coords;
            for (int coord_i = 0; coord_i < n; coord_i++) {
                coords[coord_i] = points[n*node_i + coord_i];
                height += coords[coord_i] * coords[coord_i];
            }
            coords[n] = height;

            // Node initialization
            Vec<N> point(coords);
            nodes[node_i] = new Node<N>(point, indices[node_i]);
        }
    }
    
    // Generate Initial Front
    {
        double height = 1e15;
        std::array<Node<N>*, N> initial_nodes;

        // Initial front
        internal_point = Vec<N>::Zeros();
        for (int node_i = 0; node_i < N; node_i++) {
            initial_nodes[node_i] = new Node<N>(nodes[node_i]->p, -(node_i + 1));
            initial_nodes[node_i]->p[N - 1] = height;
            internal_point += initial_nodes[node_i]->p;
        }
        internal_point /= (double)N;

        // Initial front normal
        Vec<N> start_normal = Vec<N>::Zeros();
        start_normal[N - 1] = -1.0;
        active_fronts.push_back(new FrontFace<N>(initial_nodes, start_normal, false));

        // Add initial nodes to the nodes list for later clean up
        for (Node<N> *node : initial_nodes) {
            node->active = false;
            nodes.push_back(node);
        }
    }

    // Main loop
    {
        while (active_fronts.size() > 0) {
            // Pop back
            FrontFace<N> *front = active_fronts.back();
            active_fronts.pop_back();

            // Early exit: Already not visible front
            if (!front->active) {
                delete front;
                continue;
            }

            // Find furthest active outside node
            Node<N> *furthest_node = nullptr;
            {
                double max_dist = 0.0;
                for (Node<N>* node : nodes) {
                    if (node->active) {
                        double dist = front->sdf(node->p);
                        if (dist > max_dist) {
                            max_dist = dist;
                            furthest_node = node;
                        }
                    }
                }
            }

            // Final front: FrontFace cannot see more points
            if (furthest_node == nullptr) {
                final_fronts.push_back(front);
            }

            // Expand Front
            else {
                // Front and node used
                furthest_node->active = false;
                //front->active = false;

                // Find visible fronts
                std::vector<FrontFace<N>*> visible_fronts;
                std::vector<FrontEdge<N>> all_edges;

                //// Current front
                visible_fronts.push_back(front);
                front->appendEdges(all_edges);

                //// Rest of the fronts in the list
                for (FrontFace<N> *front_i : active_fronts) {
                    if (front_i->active) {
                        if (front_i->hasVision(furthest_node->p)) {
                            visible_fronts.push_back(front_i);
                            front_i->appendEdges(all_edges);
                            front_i->active = false;  // For later removal
                        }
                    }
                }

                // Find outer edges
                std::vector<FrontEdge<N>*> outer_edges;
                
                //// Create a vector of pointers to feed the classification
                std::vector<FrontEdge<N>*> ptr_all_edges;
                ptr_all_edges.reserve(all_edges.size());
                for (FrontEdge<N> &edge : all_edges) {
                    ptr_all_edges.push_back(&edge);
                }
                //// Classification acceleration
                findNonRepeatedEdges<N>(ptr_all_edges, outer_edges, 0);

                // New fronts
                for (FrontEdge<N> *edge : outer_edges) {
                    std::array<Node<N>*, N> nodes;
                    for (int node_i = 0; node_i < N - 1; node_i++) {
                        nodes[node_i] = edge->nodes[node_i];
                    }
                    nodes[N - 1] = furthest_node;
                    FrontFace<N> *new_front = new FrontFace<N>(nodes, internal_point);
                    active_fronts.push_back(new_front);
                }

                // This front is now internal
                delete front;
            }
        }
    }

    // Connectivity array
    {
        out_connectivity.reserve(N*final_fronts.size());
        for (FrontFace<N> *front : final_fronts) {

            // Not valid if connected to the starting simplex
            bool valid = true;
            for (Node<N> *node : front->nodes) {
                if (node->idx < 0) {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                continue;
            }

            // Not valid if zero volume simplex
            if (front->getProjectedVolume(front) < TOL) {
                continue;
            }

            // Valid
            for (Node<N> *node : front->nodes) {
                out_connectivity.push_back(node->idx);
            }
        }
        out_connectivity.shrink_to_fit();
    }

    // Cleaning
    {
        for (FrontFace<N> *front : final_fronts) {
            delete front;
        }
        for (Node<N> *node : nodes) {
            delete node;
        }
    }

    return out_connectivity;
}
/* =========== End Main Function =========== */

/* =========== Alternative entry points =========== */
// Common entry points when data is in std::vector.
template<int n>
inline std::vector<int> compute(const std::vector<double> &points)
{
    // Generate range as indices
    std::vector<int> indices(points.size() / n);
    for (int node_i = 0; node_i < (int)indices.size(); node_i++) {
        indices[node_i] = node_i;
    }

    return compute<n>(&points[0], points.size(), &indices[0], indices.size());
}

template<int n>
inline std::vector<int> compute(const std::vector<double> &points, const std::vector<int> &indices)
{
    return compute<n>(&points[0], points.size(), &indices[0], indices.size());
}


/* =========== End Alternative entry points =========== */

} // namespace SimpleDelaunay

#endif
