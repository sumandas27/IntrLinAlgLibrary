#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>

//--- Standard Template Library ---//
#include <algorithm>
#include <functional>
#include <array>

//TODO: Need to write << overloaders for Vector and Matrix.
//TODO: Find a way to initialize matrices and vectors with variable argument list.

//----------------------------------------------------------------------//
//NON-LINEAR ALGEBRA FUNCTIONS:

/* Initializes the IntLinAlg Library.
 */
void IntrLinAlgLibrary_init() {
    std::cout << std::setprecision(3);
}

/* The range of tolerance for double equality.
 */
constexpr double epsilon() {
    return 0.00001;
}

/* Checks equality between two doubles if their difference lies within a tolerance range.
 */
bool is_equal(double val1, double val2) {
    return abs(val1 - val2) < epsilon();
}

/* Converts an angle from degrees to radians.
 */
constexpr double deg_to_rad(double degrees) {
    return degrees * M_PI / 180;
}

//----------------------------------------------------------------------//
//VECTOR STRUCT AND METHODS:

/* A vector is an array of real numbers.
 * Vectors are represented as column vectors as row vectors are rarely used.
 * @param D The dimension of the vector.
 */
template <unsigned int D>
struct Vector {
    std::array<double, D> components;
    Vector(const std::array<double, D>& _components);

    double&       operator[](unsigned int index);
    double const& operator[](unsigned int index) const;

    void print() const;

    template <unsigned int X>
    friend std::ostream& operator<<(std::ostream& os, const Vector<X>& v);
};

/* Constructs a vector.
 * @param _components The input structure containing the components of the vector.
 */
template <unsigned int D>
Vector<D>::Vector(const std::array<double, D>& _components) {
    assert (D != 0);
    components = _components;
}

/* Vector indexing is 1-based. This overload is primarly intended for the left side of an assignment.
 * @param index The index of the vector wanted to be changed/accessed.
 * @returns The component at the argument index.
 */
template <unsigned int D>
double& Vector<D>::operator[](unsigned int index) {
    assert (index >= 1 && index <= D);
    return components[index - 1];
}

/* Vector indexing is 1-based. This overload is exclusively used for the right side of an assignment.
 * @param index The index of the vector wanted to be accessed.
 * @returns The component at the argument index.
 */
template <unsigned int D>
double const& Vector<D>::operator[](unsigned int index) const {
    assert (index >= 1 && index <= D);
    return components[index - 1];
}

/* Prints the vector to the terminal. Each component is rounded to the nearest thousandth.
 */
template <unsigned int D>
void Vector<D>::print() const {
    for (const double& component : components)
        std::cout << "{\t" << (component < epsilon() ? 0 : component) << "\t}\n";
    std::cout << "\n";
}

/* Alternative to the print() function, allows the vector to be printed directly to the standard output.
 * Given a vector "v": "std::cout << v;"
 * @param v The vector to be printed to the standard output.
 */
template <unsigned int X>
std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
    os << "yo awsss cvector b " << X;
    return os;
}

//----------------------------------------------------------------------//
//MATRIX STRUCT AND METHODS:

/* A matrix is an array of arrays of real numbers.
 * @param R The number of rows in the matrix.
 * @param C The number of columns in the matrix.
 */
template <unsigned int R, unsigned int C>
struct Matrix {
    std::array<double, R * C> entries;
    Matrix(const std::array<double, R * C>& _entries);

    class Proxy {

    public:
        Proxy(double* _rowPtr) : rowPtr(_rowPtr) { };

        double&       operator[](unsigned int col);
        double const& operator[](unsigned int col) const;

    private:
        double* rowPtr;
    };

    Proxy operator[](unsigned int row);

    void print() const;

    template <unsigned int X, unsigned int Y>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<X, Y>& m);
};

/* Constructs a matrix.
 * @param _entries The input structure containing the entries of the matrix.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C>::Matrix(const std::array<double, R * C>& _entries) {
    assert (R != 0 && C != 0);
    entries = _entries;
}

/* Accesses the argument column of the proxy's row. Matrix indexing is 1-based. 
 * This overload is intended for the left hand side of an assignment.
 * @param col The column index of the proxy row wanted to be changed/accessed.
 * @returns The entry at the argument column of the proxy row.
 */
template <unsigned int R, unsigned int C>
double& Matrix<R, C>::Proxy::operator[](unsigned int col) {
    assert (col >= 1 && col <= C);
    return rowPtr[col - 1];
}

/* Accesses the argument column of the proxy's row. Matrix indexing is 1-based. 
 * This overload is exclusively for the right hand side of an assignment.
 * @param col The column index of the proxy row wanted to be changed/accessed.
 * @returns The entry at the argument column of the proxy row.
 */
template <unsigned int R, unsigned int C>
double const& Matrix<R, C>::Proxy::operator[](unsigned int col) const {
    assert (col >= 1 && col <= C);
    return rowPtr[col - 1];
}

/* Allows access to the entries of the matrix. Doubly overloaded subscript operators require use of an intermediate proxy class.
 * @param row The row of the entry wanted to be accessed.
 * @returns A proxy class containing a pointer to the first entry of the row.
 * The proxy's subscript overloaders can then be used to access specific entries in the row.
 */
template <unsigned int R, unsigned int C>
typename Matrix<R, C>::Proxy Matrix<R, C>::operator[](unsigned int row) {
    assert (row >= 1 && row <= R);
    return Proxy(&entries[(row - 1) * C]);
}

/* Prints the matrix to the terminal. Each entry is rouded to the nearest thousandths.
 */
template <unsigned int R, unsigned int C>
void Matrix<R, C>::print() const {
    for (int i = 0; i < R; i++) {
        std::cout << "{\t";
        for (int j = 0; j < C; j++)
            std::cout << (abs(entries[i * C + j]) < epsilon() ? 0 : entries[i * C + j]) << "\t";
        std::cout << "}\n";
    }
    std::cout << "\n";
}

/* Alternative to the print() function, allows the matrix to be printed directly to the standard output.
 * Given a matrix "m": "std::cout << m;"
 * @param m The matrix to be printed to the standard output.
 */
template <unsigned int X, unsigned int Y>
std::ostream& operator<<(std::ostream& os, const Matrix<X, Y>& m) {
    os << "yo wassup[ matrx";
    return os;
}

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

/* Two vectors are equal if all corresponding components are equal.
 * @returns true if the vector arguments v1 and v2 are equal, false if otherwise.
 */
template <unsigned int D>
bool operator==(const Vector<D>& lhs, const Vector<D>& rhs) {
    return std::equal(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), is_equal);
}

/* Two matrices are equal if all corresponding entries are equal.
 * @returns true if the matrix arguments m1 and m2 are equal, false if otherwise.
 */
template <unsigned int R, unsigned int C>
bool operator==(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    return std::equal(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), is_equal);
}

/* The sum of two vectors is a vector of the same size with corresponding components added.
 * @returns A vector that is the sum of two argument vectors.
 */
template <unsigned int D>
Vector<D> operator+(const Vector<D>& lhs, const Vector<D> rhs) {
    std::array<double, D> sum{};
    std::transform(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), sum.begin(), std::plus<double>());
    return Vector<D>(sum);
}

/* The sum of two matrices is a matrix of the same size with corresponding entries added.
 * @returns A matrix that is the sum of two argument matrices.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> operator+(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    std::array<double, R * C> sum{};
    std::transform(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), sum.begin(), std::plus<double>());
    return Matrix<R, C>(sum);
}

/* The difference of two vectors is a vector of the same size with corresponding components subtracted.
 * @returns A vector that is the difference of two argument vectors.
 */
template <unsigned int D>
Vector<D> operator-(const Vector<D>& lhs, const Vector<D>& rhs) {
    std::array<double, D> diff{};
    std::transform(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), diff.begin(), std::minus<double>());
    return Vector<D>(diff);
}

//Vector operator-(const Vector& v1, const Vector& v2);
/* The difference of two matrices is a matrix of the same size with corresponding entries subtracted.
 * @returns A matrix that is the difference of two argument matrices.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> operator-(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    std::array<double, R * C> diff{};
    std::transform(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), diff.begin(), std::minus<double>());
    return Matrix<R, C>(diff);
}

/* The product of a scalar and a vector is a vector of the same size with all its components multiplied by the scalar.
 * @returns A vector that is the product of a scalar and a vector.
 */
template <unsigned int D>
Vector<D> operator*(double scalar, const Vector<D>& v) {
    std::array<double, D> product{};
    std::transform(v.components.begin(), v.components.end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Vector<D>(product);
}

/* Scalar-vector multiplication is commutative.
 * @returns A vector that is the product of a vector and a scalar.
 */
template <unsigned int D>
Vector<D> operator*(const Vector<D>& v, double scalar) {
    return scalar * v;
}

/* The product of a scalar and a matrix is a matrix of the same size with all its entries multiplied by the scalar.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> operator*(double scalar, const Matrix<R, C>& m) {
    std::array<double, R * C> product{};
    std::transform(m.entries.begin(), m.entries.end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Matrix<R, C>(product);
}

/* Scalar-matrix multiplication is commutative.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> operator*(const Matrix<R, C>& m, double scalar) {
    return scalar * m;
}

/* A matrix-vector product is the linear combination of the vector's components and the matrix's column vectors.
 * @returns The matrix-vector product of the argument matrix and vector.
 */
template <unsigned int R, unsigned int C>
Vector<R> operator*(const Matrix<R, C>& m, const Vector<C>& v) {
    std::array<double, R> product{};
    for (unsigned int i = 0; i < R; i++) {
        double entry = 0.0;
        for (unsigned int j = 0; j < C; j++)
            entry += v.components[j] * m.entries[i * C + j];
        product[i] = entry;
    }
    return Vector<R>(product);
}

/* A zero vector is a vector where all components are zero.
 * @param D The dimension of the zero matrix.
 * @returns A zero vector of the argument size.
 */
template <unsigned int D>
Vector<D> zero_vector() {
    std::array<double, D> zero{};
    return Vector<D>(zero);
}

/* A zero matrix is a matrix where all entries are zero.
 * @param R The number of rows in the zero matrix.
 * @param C The number of columns in the zero matrix.
 * @returns A zero matrix of the argument size.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> zero_matrix() {
    std::array<double, R * C> zero{};
    return Matrix<R, C>(zero);
}

/* A standard vector is a zero vector with one component being a one instead of a zero.
 * @param D The dimension of the standard vector.
 * @param one_component The location of the "one" component.
 * @returns A standard vector of the argument dimension with the 1 in the argument location.
 */
template <unsigned int D>
Vector<D> standard_vector(unsigned int one_component) {
    assert (one_component >= 1 && one_component <= D);

    std::array<double, D> standardVector{};
    standardVector[one_component - 1] = 1.0;
    return Vector<D>(standardVector);
}

/* An identity matrix is a square zero matrix with diagonal entries being a one instead of a zero.
 * @param S The size: the number of rows and columns of the identity matrix.
 * @returns An identity matrix of the argument size.
 */
template <unsigned int S>
Matrix<S, S> identity_matrix() {
    std::array<double, S * S> identityMatrix{};
    for (unsigned int i = 0; i < S; i++)
        identityMatrix[i * S + i] = 1.0;
    return Matrix<S, S>(identityMatrix);
}

/* A rotation matrix is a 2x2 matrix that rotates an R^2 vector by some amount of degrees counter-clockwise.
 * Given rotation matrix A and some R^2 vector x, Ax = x' where x' is the rotated R^2 vector.
 * @param degrees The angle (in degrees) of the rotation matrix.
 * @returns The 2x2 rotation matrix of the argument angle in degrees.
 */
Matrix<2, 2> rotation_matrix(double degrees) {
    std::array<double, 4> rotationMatrix{
         cos(deg_to_rad(degrees)),
        -sin(deg_to_rad(degrees)),
         sin(deg_to_rad(degrees)),
         cos(deg_to_rad(degrees))
    };
    return Matrix<2, 2>(rotationMatrix);
}

/* The transpose of an nxm matrix is an mxn matrix where (i,j)-entries are transformed to (j,i)-entries.
 * @param m The matrix whose transpose is to be returned.
 * @returns The transpose of the argument matrix.
 */
template <unsigned int R, unsigned int C>
Matrix<C, R> transpose(const Matrix<R, C>& m) {
    std::array<double, C * R> transpose{};
    for (unsigned int i = 0; i < C; i++)
    for (unsigned int j = 0; j < R; j++)
        transpose[i * R + j] = m.entries[j * C + i];
    return Matrix<C, R>(transpose);
}

/* An elementary row operation where two rows are exchanged in a matrix: row1 <--> row2
 * @param m The matrix to be modified.
 * @param row1 The first row to be exchanged.
 * @param row2 The second row to be exchanged.
 */
template <unsigned int R, unsigned int C>
void ERO_row_swap(Matrix<R, C>& m, unsigned int row1, unsigned int row2) {
    assert (row1 >= 1 && row1 <= R);
    assert (row2 >= 1 && row2 <= R);

    if (row1 == row2)
        return;

    for (unsigned int col = 0; col < C; col++)
        std::swap(m.entries[(row1 - 1) * C + col], m.entries[(row2 - 1) * C + col]);
}

/* An elementary row operation where a row is multiplied by a constant in a matrix: scalar * row --> row
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the row by.
 * @param row The row to be multiplied.
 */
//void ERO_scalar_multiplication(Matrix& m, double scalar, unsigned int row);
template <unsigned int R, unsigned int C>
void ERO_scalar_multiplication(Matrix<R, C>& m, double scalar, unsigned int row) {
    assert (row >= 1 && row <= R);

    unsigned int beg = (row - 1) * C;
    unsigned int end = (row - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, m.entries.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
}

/* An elementary row where a multiple of one row is added to another row in a matrix: scalar * scaledRow + outputRow --> outputRow
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the scaled row by.
 * @param scaledRow The row to be scaled by.
 * @param outputRow The output row to add and to copy the results into.
 */
template <unsigned int R, unsigned int C>
void ERO_row_sum(Matrix<R, C>& m, double scalar, unsigned int rowToScale, unsigned int outputRow) {
    assert (rowToScale >= 1 && rowToScale <= R);
    assert (outputRow >= 1 && outputRow <= R);

    std::array<double, C> scaledRow{};
    unsigned int beg = (rowToScale - 1) * C;
    unsigned int end = (rowToScale - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));

    beg = (outputRow - 1) * C;
    end = (outputRow - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), m.entries.begin() + beg, std::plus<double>());
}

/* The row-echelon form (ref) of a matrix is a matrix with the same solution set that follows 2 restrictions:
 *   1. Every nonzero row lies above all zero rows
 *   2. The leading entry of a nonzero row is in a column to the right of every leading entry of a nonzero row above
 * 
 * Matrices may have an infinite amount of row-echelon form, this function returns the one calculated by the forward pass
 * of the Gaussian Elimination.
 * @param m The matrix whose row-echelon form is to be returned.
 * @returns The row-echelon form of the argument matrix.
 */
template <unsigned int R, unsigned int C>
Matrix<R, C> ref(Matrix<R, C> m) {
    unsigned int pivotRow = 0;
    for (unsigned int i = 0; i < C; i++) {
        if (pivotRow >= R - 1)
            break;

        double nonzeroFound = 0.0; /* 0 means a nonzero entry has not been found, else nonzero entry is set to this variable */
        for (unsigned int j = pivotRow; j < R; j++) {
            if (is_equal(m.entries[j * C + i], 0.0))
                continue;

            if (nonzeroFound == 0.0) {
                nonzeroFound = m.entries[j * C + i];
                ERO_row_swap(m, pivotRow + 1, j + 1);
                pivotRow++;
            }
            else {
                double scalar = -m.entries[j * C + i] / nonzeroFound;
                ERO_row_sum(m, scalar, pivotRow, j + 1);
            }
        }
    }

    return m;
}

template <unsigned int R, unsigned int C>
Matrix<R, C> rref(Matrix<R, C> m) {
    return m; //temp
}