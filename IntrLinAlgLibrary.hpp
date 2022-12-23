#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>

#include <string>
#include <cmath>
#include <numeric>

#include <algorithm>
#include <functional>
#include <array>
#include <vector>

/* README
 * ------------------------------------------------------------------------------------------------------------------
 * Link to the Intro to Linear Algebra textbook my course uses:
 * https://cloudflare-ipfs.com/ipfs/bafykbzaced62pai3pki6og6a5jimvjis3qpwo3yvzinlvj6gvjqzzrcqyhkre?filename=Lawrence%20E.%20Spence%2C%20Arnold%20J.%20Insel%2C%20Stephen%20H.%20Friedberg%20-%20Elementary%20Linear%20Algebra%20%282nd%20Edition%29-Prentice%20Hall%20%282007%29.pdf
 *
 * Write 'using ila::Vector, ila::Matrix;' to type 'Vector' and 'Matrix' instead of 'ila::Vector' and 'ila::Matrix'.
 * Write 'using namespace ila;' to remove the prefix entirely. There may be potential naming conflicts with other libraries.
 * 
 * This library uses templates for rows and columns, so all matrix and vector sizes must be known at compile time.
 * Matrix and vector sizes cannot change after being initialized.
 * 
 * Creating a vector:
 * ila::Vector<4> myVector(1.0, 4.0, 3.0, 11.4);
 * 
 * Creating a matrix (line breaks for readability, but are optional):
 * ila::Matrix<2, 2> myMatrix
 * (
 *     1.0, 2.0,
 *     3.0, 4.0
 * ); 
 * 
 * Vectors and matrices are 1-indexed (the first row is at index 1 rather than 0):
 *   - Accessing the 2nd component of a vector: 'myVector[2]'
 *   - Accessing the (2, 1)-entry of a matrix: 'myMatrix[2][1]'
 * 
 * Creating a user-defined set of vectors (line breaks for readability, but are optional):
 * std::array<ila::Vector<3>, 2> mySet = // A set of 2 vectors in R3
 * {
 *     Vector<3>(1.0, 2.0, 3.0),
 *     Vector<3>(4.0, 5.0, 6.0)
 * };
 * 
 * NOTE FOR SETS OF VECTORS:
 * 
 *  - This library handles sets of vectors using std::arrays, so use std::arrays for user-defined vector sets.
 * 
 *  - However, dimensions of bases for row, column, and null spaces cannot be calculated at compile time and so are
 *    returned via an std::vector. Manually hardcode the std::vector contents** into new std::arrays and recompile
 *    to be compatible with this library's functions.
 *      
 *    ** Use the print overload for std::vector of Vectors to view their contents as an augmented matrix easily.
 */

namespace ila { // Intro Linear Algebra

//------------------------------------------------------------------------------------------//
//NON-LINEAR ALGEBRA FUNCTIONS:

/* Setting that prints the first 'x' digits of every vector component and matrix entry. Default set to 5.
 */
std::streamsize precision = 4;

/* Sets the precision of vector components and matrix entries when printed.
 * Maximum precision allowed is 6.
 * @param _precision Sets printing all components and entries to the first '_precision' digits.
 */
void set_precision(std::streamsize _precision) {
    assert (_precision >= 0 && _precision <= 6);
    precision = _precision;
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

//------------------------------------------------------------------------------------------//
//VECTOR STRUCT AND METHODS:

/* A vector is an array of real numbers.
 * Vectors are represented as column vectors as row vectors are rarely used.
 * @param D The dimension of the vector.
 */
template <size_t D>
class Vector {
public:
    std::array<double, D> components;

    template <typename... Ts>
    Vector(Ts... _components);
    Vector(std::array<double, D> _components);
    Vector();

    double& operator[](size_t index);

    Vector<D>& operator+=(const Vector<D>& v);
    Vector<D>& operator-=(const Vector<D>& v);
    Vector<D>& operator*=(double scalar);
    Vector<D>& operator/=(double scalar);
};

/* Constructs a vector with a list of arguments. Example Vector object creation: 
 * Vector<5> vec(1.0, 2.0, 3.0, 4.0, 5.0);
 * 
 * @param _components The list of scalar arguments containing the components of the vector.
 */
template <size_t D>
template <typename... Ts>
Vector<D>::Vector(Ts... _components) 
    : components{ static_cast<double>(std::forward<Ts>(_components))... } 
{
    assert (D != 0);
    assert (sizeof...(_components) == D);
}

/* Allows a vector to be constructed directly from an array.
 * @param _components The array holding the vector's components.
 */
template <size_t D>
Vector<D>::Vector(std::array<double, D> _components)
    : components(_components) { }

/* Constructs a vector with zero-initialized components.
 */
template <size_t D>
Vector<D>::Vector() : components{} { }

/* Vector indexing is 1-based.
 * @param index The index of the vector wanted to be changed/accessed.
 * @returns The component at the argument index.
 */
template <size_t D>
double& Vector<D>::operator[](size_t index) {
    assert (index >= 1 && index <= D);
    return components[index - 1];
}

/* Vector<D> + Vector<D> operator overload is in Chapter 1.
 */
template <size_t D>
Vector<D>& Vector<D>::operator+=(const Vector<D>& v) {
    *this = *this + v;
    return *this;
}

/* Vector<D> - Vector<D> operator overload is in Chapter 1.
 */
template <size_t D>
Vector<D>& Vector<D>::operator-=(const Vector<D>& v) {
    *this = *this - v;
    return *this;
}

/* double * Vector<D> operator overload is in Chapter 1.
 */
template <size_t D>
Vector<D>& Vector<D>::operator*=(double scalar) {
    *this = scalar * *this;
    return *this;
}

/* Vector<D> / double operator overload is in Chapter 1.
 */
template <size_t D>
Vector<D>& Vector<D>::operator/=(double scalar) {
    *this = *this / scalar;
    return *this;
}

//------------------------------------------------------------------------------------------//
//MATRIX STRUCT AND METHODS:

/* A matrix is an array of arrays of real numbers.
 * @param R The number of rows in the matrix.
 * @param C The number of columns in the matrix.
 */
template <size_t R, size_t C>
class Matrix {
public:
    std::array<std::array<double, C>, R> entries;

    template <typename... Ts>
    Matrix(Ts... _entries);
    Matrix(std::array<std::array<double, C>, R> _entries);
    Matrix();

    Vector<C> row_vector(size_t row) const;
    Vector<R> column_vector(size_t col) const;

    class Proxy {
    public:
        Proxy(std::array<double, C> _row) : row(_row) { };
        double& operator[](size_t col);

    private:
        std::array<double, C> row;
    };

    Proxy operator[](size_t row);

    Matrix<R, C>& operator+=(const Matrix<R, C>& m);
    Matrix<R, C>& operator-=(const Matrix<R, C>& m);
    Matrix<R, C>& operator*=(double scalar);
    Matrix<R, C>& operator/=(double scalar);
};

/* Constructs a matrix. Example Matrix object creation:
 * Matrix<2, 3> mat
 * (
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * );
 *
 * @param _entries The list of scalar arguments containing the entries of the matrix.
 */
template <size_t R, size_t C>
template <typename... Ts>
Matrix<R, C>::Matrix(Ts... _entries) 
    : entries{}
{
    assert (R != 0 && C != 0);
    assert (sizeof...(_entries) == R * C);

    std::array<double, R * C> temp{ static_cast<double>(std::forward<Ts>(_entries))... };
    for (size_t row = 0; row < R; row++) {
        size_t beg = row * C;
        size_t end = row * C + C;
        std::copy(temp.begin() + beg, temp.begin() + end, entries[row].begin());
    }
}

/* Allows a matrix to be constructed directly from a 2D array.
 * @param _entries The 2D array holding the matrix's entries.
 */
template <size_t R, size_t C>
Matrix<R, C>::Matrix(std::array<std::array<double, C>, R> _entries)
    : entries(_entries) { }

/* Constructs a Matrix object with zero-intialized entries.
 */
template <size_t R, size_t C>
Matrix<R, C>::Matrix() : entries{} { }

/* A row vector of a matrix is a vector that contains the entries of a row.
 * @param row The argument row index (1-based).
 * @returns The row vector at the argument index.
 */
template <size_t R, size_t C>
Vector<C> Matrix<R, C>::row_vector(size_t row) const {
    assert (row >= 1 && row <= R);
    return Vector<C>(entries[row - 1]);
}

/* A column vector of a matrix is a vector that contains the entries of a column.
 * @param col The argument column index (1-based).
 * @returns The column vector at the argument index.
 */
template <size_t R, size_t C>
Vector<R> Matrix<R, C>::column_vector(size_t col) const {
    assert (col >= 1 && col <= C);
    std::array<double, R> column{};
    for (size_t row = 0; row < R; row++)
        column[row] = entries[row][col - 1];
    return Vector<R>(column);
}

/* Accesses the argument column of the proxy's row.
 * This overload is intended for the left hand side of an assignment.
 * @param col The column index of the proxy row wanted to be changed/accessed.
 * @returns The entry at the argument column of the proxy row.
 */
template <size_t R, size_t C>
double& Matrix<R, C>::Proxy::operator[](size_t col) {
    assert (col >= 1 && col <= C);
    return row[col - 1];
}

/* Allows access to the entries of the matrix. Doubly overloaded subscript operators require use of an intermediate proxy class.
 * @param row The row of the entry wanted to be accessed.
 * @returns A proxy class containing a pointer to the first entry of the row.
 * The proxy's subscript overloaders can then be used to access specific entries in the row.
 */
template <size_t R, size_t C>
typename Matrix<R, C>::Proxy Matrix<R, C>::operator[](size_t row) {
    assert (row >= 1 && row <= R);
    return Proxy(entries[row - 1]);
}

/* Matrix<R, C> + Matrix<R, C> operator overload is in Chapter 1.
 */
template <size_t R, size_t C>
Matrix<R, C>& Matrix<R, C>::operator+=(const Matrix<R, C>& m) {
    *this = *this + m;
    return *this;
}

/* Matrix<R, C> - Matrix<R, C> operator overload is in Chapter 1.
 */
template <size_t R, size_t C>
Matrix<R, C>& Matrix<R, C>::operator-=(const Matrix<R, C>& m) {
    *this = *this - m;
    return *this;
}

/* double * Matrix<R, C> operator overload is in Chapter 1.
 */
template <size_t R, size_t C>
Matrix<R, C>& Matrix<R, C>::operator*=(double scalar) {
    *this = scalar * *this;
    return *this;
}

/* Matrix<R, C> / double operator overload is in Chapter 1.
 */
template <size_t R, size_t C>
Matrix<R, C>& Matrix<R, C>::operator/=(double scalar) {
    *this = *this / scalar;
    return *this;
}

//------------------------------------------------------------------------------------------//
//PRINT FUNCTION OVERLOADS

/* Alternative to the print() function, allows the vector to be printed directly to the standard output.
 * Given a vector "v": "std::cout << v;"
 * @param v The vector to be printed to the standard output.
 */
template <size_t X>
std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
    std::streamsize original = os.precision();
    os << std::setprecision(precision);
    for (double component : v.components)
        os << "{\t" << (abs(component) < epsilon() ? 0 : component) << "\t}\n";
    os << "\n";
    os << std::setprecision(original);
    return os;
}

/* Prints a vector to the terminal. Each component is rounded to the nearest thousandth.
 * @param v The argument vector to be printed.
 */
template <size_t D>
void print(const Vector<D>& v) {
    std::cout << v;
}

/* Alternative to the print() function, allows the matrix to be printed directly to the standard output.
 * Given a matrix "m": "std::cout << m;"
 * @param m The matrix to be printed to the standard output.
 */
template <size_t X, size_t Y>
std::ostream& operator<<(std::ostream& os, const Matrix<X, Y>& m) {
    std::streamsize original = os.precision();
    os << std::setprecision(precision);
    for (size_t row = 0; row < X; row++) {
        os << "{\t";
        for (double entry : m.entries[row])
            os << (abs(entry) < epsilon() ? 0 : entry) << "\t";
        os << "}\n";
    }
    os << std::setprecision(original);
    return os;
}

/* Prints the matrix to the terminal. Each entry is rouded to the nearest thousandths.
 * @param m The argument matrix to be printed.
 */
template <size_t R, size_t C>
void print(const Matrix<R, C>& m) {
    std::cout << m;
}

/* Prints a set of vectors (std::array of vectors) as its augmented matrix.
 * @param set The argument std::array of vectors to be printed.
 */
template <size_t D, size_t S>
void print(const std::array<Vector<D>, S>& set) {
    print(augment_vector_set(set));
}

/* Prints a set of vectors (std::vector of vectors) as its augmented matrix.
 * This function prints something that resembles the set's augmented matrix, but no augmenting actually occurs
 * (because vector sizes are not constant, so they cannot be used in matrix template arguments).
 * @param set The argument std::vector of vectors to be printed.
 */
template <size_t D>
void print(const std::vector<Vector<D>>& set) {
    std::streamsize original = std::cout.precision();
    std::cout << std::setprecision(precision);
    for (size_t row = 0; row < D; row++) {
        std::cout << "{\t";
        for (size_t col = 0; col < set.size(); col++) {
            double component = set.at(col).components[row];
            std::cout << (abs(component) < epsilon() ? 0 : component) << "\t";
        }
        std::cout << "}\n";
    }
    std::cout << std::setprecision(original);
}

//------------------------------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

/* Two vectors are equal if all corresponding components are equal.
 * @returns True if the vector arguments v1 and v2 are equal. False if otherwise.
 */
template <size_t D>
bool operator==(const Vector<D>& lhs, const Vector<D>& rhs) {
    return std::equal(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), is_equal);
}

/* Opposite of vector equality.
 */
template <size_t D>
bool operator!=(const Vector<D>& lhs, const Vector<D>& rhs) {
    return !(lhs == rhs);
}

/* Two matrices are equal if all corresponding entries are equal.
 * @returns True if the matrix arguments m1 and m2 are equal. False if otherwise.
 */
template <size_t R, size_t C>
bool operator==(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    for (size_t row = 0; row < R; row++)
        if (!std::equal(lhs.entries[row].begin(), lhs.entries[row].end(), rhs.entries[row].begin()))
            return false;
        
    return true;
}

/* Opposite of matrix equality.
 */
template <size_t R, size_t C>
bool operator!=(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    return !(lhs == rhs);
}

/* The sum of two vectors is a vector of the same size with corresponding components added.
 * @returns A vector that is the sum of two argument vectors.
 */
template <size_t D>
Vector<D> operator+(const Vector<D>& lhs, const Vector<D> rhs) {
    Vector<D> sum{};
    std::transform(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), sum.components.begin(), std::plus<double>());
    return sum;
}

/* The sum of two matrices is a matrix of the same size with corresponding entries added.
 * @returns A matrix that is the sum of two argument matrices.
 */
template <size_t R, size_t C>
Matrix<R, C> operator+(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    Matrix<R, C> sum{};
    for (size_t row = 0; row < R; row++)
        std::transform(lhs.entries[row].begin(), lhs.entries[row].end(), rhs.entries[row].begin(), sum.entries[row].begin(), std::plus<double>());
    return sum;
}

/* The difference of two vectors is a vector of the same size with corresponding components subtracted.
 * @returns A vector that is the difference of two argument vectors.
 */
template <size_t D>
Vector<D> operator-(const Vector<D>& lhs, const Vector<D>& rhs) {
    Vector<D> diff{};
    std::transform(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), diff.components.begin(), std::minus<double>());
    return diff;
}

//Vector operator-(const Vector& v1, const Vector& v2);
/* The difference of two matrices is a matrix of the same size with corresponding entries subtracted.
 * @returns A matrix that is the difference of two argument matrices.
 */
template <size_t R, size_t C>
Matrix<R, C> operator-(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    Matrix<R, C> diff{};
    for (size_t row = 0; row < R; row++)
        std::transform(lhs.entries[row].begin(), lhs.entries[row].end(), rhs.entries[row].begin(), diff.entries[row].begin(), std::minus<double>());
    return diff;
}

/* The product of a scalar and a vector is a vector of the same size with all its components multiplied by the scalar.
 * @returns A vector that is the product of a scalar and a vector.
 */
template <size_t D>
Vector<D> operator*(double scalar, const Vector<D>& v) {
    Vector<D> product{};
    std::transform(v.components.begin(), v.components.end(), product.components.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return product;
}   

/* Scalar-vector multiplication is commutative.
 * @returns A vector that is the product of a vector and a scalar.
 */
template <size_t D>
Vector<D> operator*(const Vector<D>& v, double scalar) {
    return scalar * v;
}

/* The quotient of a vector and a scalar is a vector of the same size with all its components divided by the scalar.
 * @returns A vector that is the quotient of a vector and a scalar.
 */
template <size_t D>
Vector<D> operator/(const Vector<D>& v, double scalar) {
    double invScalar = 1 / scalar;
    return invScalar * v;
}

/* The product of a scalar and a matrix is a matrix of the same size with all its entries multiplied by the scalar.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> operator*(double scalar, const Matrix<R, C>& m) {
    Matrix<R, C> product{};
    for (size_t row = 0; row < R; row++)
        std::transform(m.entries[row].begin(), m.entries[row].end(), product.entries[row].begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return product;
}

/* Scalar-matrix multiplication is commutative.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> operator*(const Matrix<R, C>& m, double scalar) {
    return scalar * m;
}

/* The quotient of a matrix and a scalar is a matrix of the same size with all its entries divided by the scalar.
 * @returns A matrix that is the quotient of a matrix and a scalar.
 */
template <size_t R, size_t C>
Matrix<R, C> operator/(const Matrix<R, C>& m, double scalar) {
    double invScalar = 1 / scalar;
    return invScalar * m;
}

/* A matrix-vector product is the linear combination of the vector's components and the matrix's column vectors.
 * @returns The matrix-vector product of the argument matrix and vector.
 */
//TODO: Replace with dot product
template <size_t R, size_t C>
Vector<R> operator*(const Matrix<R, C>& m, const Vector<C>& v) {
    Vector<R> product{};
    for (size_t i = 0; i < R; i++) {
        double entry = 0.0;
        for (size_t j = 0; j < C; j++)
            entry += v.components[j] * m.entries[i][j];
        product.components[i] = entry;
    }
    return product;
}

/* A zero vector is a vector where all components are zero.
 * @param D The dimension of the zero matrix.
 * @returns A zero vector of the argument size.
 */
template <size_t D>
Vector<D> zero_vector() {
    assert (D != 0);
    return Vector<D>();
}

/* A zero matrix is a matrix where all entries are zero.
 * @param R The number of rows in the zero matrix.
 * @param C The number of columns in the zero matrix.
 * @returns A zero matrix of the argument size.
 */
template <size_t R, size_t C>
Matrix<R, C> zero_matrix() {
    assert (R != 0 && C != 0);
    return Matrix<R, C>();
}

/* A standard vector is a zero vector with one component being a one instead of a zero.
 * @param D The dimension of the standard vector.
 * @param one_component The location of the "one" component.
 * @returns A standard vector of the argument dimension with the 1 in the argument location.
 */
template <size_t D>
Vector<D> standard_vector(size_t one_component) {
    assert (one_component >= 1 && one_component <= D);

    Vector<D> standardVector{};
    standardVector.components[one_component - 1] = 1.0;
    return standardVector;
}

/* An identity matrix is a square zero matrix with diagonal entries being a one instead of a zero.
 * @param S The size: the number of rows and columns of the identity matrix.
 * @returns An identity matrix of the argument size.
 */
template <size_t S>
Matrix<S, S> identity_matrix() {
    Matrix<S, S> identityMatrix{};
    for (size_t i = 0; i < S; i++)
        identityMatrix.entries[i][i] = 1.0;
    return identityMatrix;
}

/* A rotation matrix is a 2x2 matrix that rotates an R^2 vector by some amount of degrees counter-clockwise.
 * Given rotation matrix A and some R^2 vector x, Ax = x' where x' is the rotated R^2 vector.
 * @param degrees The angle (in degrees) of the rotation matrix.
 * @returns The 2x2 rotation matrix of the argument angle in degrees.
 */
Matrix<2, 2> rotation_matrix(double degrees) {
    return Matrix<2, 2>(
        cos(deg_to_rad(degrees)), -sin(deg_to_rad(degrees)),
        sin(deg_to_rad(degrees)),  cos(deg_to_rad(degrees))
    );
}

/* The transpose of an nxm matrix is an mxn matrix where (i,j)-entries are transformed to (j,i)-entries.
 * @param m The matrix whose transpose is to be returned.
 * @returns The transpose of the argument matrix.
 */
template <size_t R, size_t C>
Matrix<C, R> transpose(const Matrix<R, C>& m) {
    Matrix<C, R> transpose{};
    for (size_t i = 0; i < C; i++)
    for (size_t j = 0; j < R; j++)
        transpose.entries[i][j] = m.entries[j][i];
    return transpose;
}

/* An elementary row operation where two rows are exchanged in a matrix: row1 <--> row2
 * @param m The matrix to be modified.
 * @param row1 The first row to be exchanged.
 * @param row2 The second row to be exchanged.
 */
template <size_t R, size_t C>
void ERO_row_swap(Matrix<R, C>& m, size_t row1, size_t row2) {
    assert (row1 >= 1 && row1 <= R);
    assert (row2 >= 1 && row2 <= R);

    std::swap(m.entries[row1 - 1], m.entries[row2 - 1]);
}

/* An elementary row operation where a row is multiplied by a constant in a matrix: scalar * row --> row
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the row by.
 * @param row The row to be multiplied.
 */
template <size_t R, size_t C>
void ERO_scalar_multiplication(Matrix<R, C>& m, double scalar, size_t row) {
    assert (row >= 1 && row <= R);
    std::transform(m.entries[row - 1].begin(), m.entries[row - 1].end(), m.entries[row - 1].begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
}

/* An elementary row where a multiple of one row is added to another row in a matrix: scalar * scaledRow + outputRow --> outputRow
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the scaled row by.
 * @param scaledRow The row to be scaled by.
 * @param outputRow The output row to add and to copy the results into.
 */
template <size_t R, size_t C>
void ERO_row_sum(Matrix<R, C>& m, double scalar, size_t rowToScale, size_t outputRow) {
    assert (rowToScale >= 1 && rowToScale <= R);
    assert (outputRow >= 1 && outputRow <= R);

    std::array<double, C> scaledRow{};
    std::transform(m.entries[rowToScale - 1].begin(), m.entries[rowToScale - 1].end(), scaledRow.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    std::transform(m.entries[outputRow - 1].begin(), m.entries[outputRow - 1].end(), scaledRow.begin(), m.entries[outputRow - 1].begin(), std::plus<double>());
}

/* This function transforms the matrix argument itself to reduced-row echelon form, avoiding an unnecessary copy.
 * This function should only be read for implementation details and *SHOULD NOT* be used outside this file.
 * Instead, use the function "ref(Matrix)" to get the row-echelon form of a matrix.
 * @param m The argument matrix to be changed to its row-echelon form.
 * @returns The first value of the returned pair returns the rank of the argument matrix. This is used in rank calculation.
 * @returns The second value of the returned pair returns the number of row swaps. This is used for calculating the determinant of a matrix.
 */
template <size_t R, size_t C>
std::pair<unsigned int, unsigned int> ref_by_reference(Matrix<R, C>& m) {
    size_t pivotRow = 0;
    unsigned int rowSwaps = 0;

    for (size_t col = 0; col < C; col++) {
        if (pivotRow > R - 1)
            break;

        double nonzeroFound = 0.0; /* 0.0 means a nonzero entry has not been found, else nonzero entry is set to this variable */
        for (size_t row = pivotRow; row < R; row++) {
            if (is_equal(m.entries[row][col], 0.0))
                continue;

            if (nonzeroFound != 0.0) {
                double scalar = -m.entries[row][col] / nonzeroFound;
                ERO_row_sum(m, scalar, pivotRow, row + 1);
                continue;
            }

            nonzeroFound = m.entries[row][col]; /* First nonzero of the row found */
            if (pivotRow != row) {
                ERO_row_swap(m, pivotRow + 1, row + 1);
                rowSwaps++;
            }
            pivotRow++;
        }
    }

    unsigned int rank = pivotRow;
    return std::make_pair(rank, rowSwaps);
}

/* The row-echelon form (ref) of a matrix is a matrix with the same solution set that follows 2 restrictions:
 *  1. Every nonzero row lies above all zero rows
 *  2. The leading entry of a nonzero row is in a column to the right of every leading entry of a nonzero row above
 * 
 * Matrices may have an infinite amount of row-echelon form. This function returns the one calculated by the forward pass of Gaussian Elimination.
 * For implementation details, read the function "ref_by_reference(Matrix&)".
 * @param m The argument matrix.
 * @returns The row-echelon form of the argument matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> ref(Matrix<R, C> m) {
    ref_by_reference(m);
    return m;
}

/* The reduced row-echelon form (rref) of a matrix is the matrix in row-echelon form (ref) with 2 additional conditions:
 *  1. If a column contains the leading entry of some nonzero row, all other entries in that column are 0.
 *  2. The leading entry of all nonzero rows is 1.
 * 
 * Matrices have one unique reduced row-echelon form. This function finds the rref via. Gaussian Elimination.
 * @param m The argument matrix.
 * @returns The reduced row-echelon form of the argument matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> rref(Matrix<R, C> m) {
    ref_by_reference(m);
    size_t pivotRow = 0;
    for (size_t col = 0; col < C; col++) {
        if (pivotRow > R - 1)
            break;

        if (is_equal(m.entries[pivotRow][col], 0.0))
            continue;

        if (m.entries[pivotRow][col] != 1.0) {
            double scalar = 1.0 / m.entries[pivotRow][col];
            ERO_scalar_multiplication(m, scalar, pivotRow + 1);
        }

        for (int row = pivotRow - 1; row >= 0; row--)
            if (m.entries[row][col] != 0) {
                double scalar = -m.entries[row][col];
                ERO_row_sum(m, scalar, pivotRow + 1, row + 1);
            }

        pivotRow++;
    }
    return m; 
}

/* The rank of a matrix is the dimension of the matrix's column space.
 * This is equivalent to the number of pivot rows in the matrix's row-echelon form (which is how this method calculates the rank).
 * @param m The argument matrix.
 * @returns The rank of the argument matrix.
 */
template <size_t R, size_t C>
constexpr unsigned int rank(const Matrix<R, C>& m) {
    Matrix<R, C> copy = m;
    unsigned int matRank = ref_by_reference(copy).first;
    return matRank;
}

/* The nullity of a matrix is the dimension of the matrix's null space.
 * This is equivalent to the rank of the matrix subtracted from its total number of columns.
 * @param m The argument matrix.
 * @returns The nullity of the argument matrix.
 */
template <size_t R, size_t C>
unsigned int nullity(const Matrix<R, C>& m) {
    return C - rank(m);
}

/* A matrix augmented with a vector produces a matrix containing the matrix on the left side and the vector on the right side.
 * Given matrix A and vector v, A augmented with v is written as [ A | v ].
 * @param m The argument matrix.
 * @param v The argument vector.
 * @returns The argument matrix augmented with the argument vector.
 */
template <size_t R, size_t C>
Matrix<R, C + 1> augment(const Matrix<R, C>& m, const Vector<R>& v) {
    Matrix<R, C + 1> augmentedMatrix{};
    for (size_t row = 0; row < R; row++)
    for (size_t col = 0; col < C + 1; col++)
        augmentedMatrix.entries[row][col] = (col == C) ? v.components[row] : m.entries[row][col];
    return augmentedMatrix;
}

/* A matrix augmented with another matrix appends the matrices together to form a larger matrix.
 * Given matrix A and matrix B, A augmented with B is written as [ A | B ].
 * @param lhs The matrix to be augmented to the left hand side of the result matrix.
 * @param rhs The matrix to be augmented to the right hand side of the result matrix.
 * @returns The lhs argument matrix augmented with the rhs augmented matrix.
 */
template <size_t R, size_t C1, size_t C2>
Matrix<R, C1 + C2> augment(const Matrix<R, C1>& lhs, const Matrix<R, C2>& rhs) {
    Matrix<R, C1 + C2> augmentedMatrix{};
    for (size_t row = 0; row < R; row++)
    for (size_t col = 0; col < C1 + C2; col++)
        augmentedMatrix.entries[row][col] = (col < C1) ? lhs.entries[row][col] : rhs.entries[row][col - C1];
    return augmentedMatrix;
}

/* Given a coefficient matrix A and a constant vector b, this solves Ax = b (where x is the solution vector).
 * The reduced row-echelon form of A augmented to b allows the general solution of x to be easily found.
 * 
 *                       x1 x2 x3 x4  x5 b
 * Given an rref row: [  0  1  0  2  -3  3  ]
 * The leading entry of the rref can be solved for: x2 + 2x4 - 3x5 = 3 --> x2 = 3 - 2x4 + 3x5
 * A solution variable whose corresponding column has no pivot spot is a free variable.
 * 
 * @param coeffMat The argument coefficient matrix.
 * @param constantVec The argument constant vector.
 * @returns The reduced row-echelon form of the augmented matrix of the coefficient matrix and the constant vector.
 */
template <size_t R, size_t C>
Matrix<R, C + 1> solve(const Matrix<R, C>& coeffMat, const Vector<R>& constantVec) {
    return rref(augment(coeffMat, constantVec));
}

/* A system is consistent if given a coefficient matrix A and a constant vector b, a solution x exists for Ax = b.
 * Given the reduced row-echelon form of the augmented matrix of A and b, the system is not consistent if a row as follows exists:
 *
 *    x1 x2      xn b
 * [  0  0  ...  0  c ] where c is a nonzero scalar
 * This is synonymous to 0x1 + 0x2 + ... + 0xn = c --> 0 = c. 0 cannot equal a nonzero number, so no solution exists for x and the system is inconsistent.
 * 
 * @param coeffMat The argument coefficient matrix (A).
 * @param constantVec The argument constant vector (b).
 * @returns True if Ax = b is a consistent solution. False if otherwise.
 */
template <size_t R, size_t C>
bool is_consistent(const Matrix<R, C>& coeffMat, const Vector<R>& constantVec) {
    Matrix<R, C + 1> rrefAugment = rref(solve(coeffMat, constantVec));
    for (size_t row = 0; row < R; row++) {
        if (rrefAugment.entries[row][C] == 0.0)
            continue;

        if (std::all_of(rrefAugment.entries[row].begin(), rrefAugment.entries[row].end() - 1, [](double entry){ return entry == 0.0; }))
            return false;
    }
    return true;
}

/* Augmenting a set of vectors merges all of the vectors to produce a matrix.
 * The nth column of the result augmented matrix corresponds to the nth vector in the set.
 * @param set The argument vector set.
 * @returns An augmented matrix of the set of vectors.
 */
template <size_t D, size_t S>
Matrix<D, S> augment_vector_set(const std::array<Vector<D>, S>& set) {
    Matrix<D, S> augmentedVectorSet{};
    for (size_t col = 0; col < S; col++)
    for (size_t row = 0; row < D; row++)
        augmentedVectorSet.entries[row][col] = set.at(col).components[row];
    return augmentedVectorSet;
}

/* A vector is in the span of a set of vectors if the vector can be written as a linear combination of the other vectors.
 * @param vec The argument vector.
 * @param set The argument set of vectors.
 * @returns True if the argument vector is in the span of the argument set of vectors. False if otherwise.
 */
template <size_t D, size_t S>
bool is_in_span(const Vector<D>& vec, const std::array<Vector<D>, S>& set) {
    Matrix<D, S> augmentedSet = augment_vector_set(set);
    return is_consistent(augmentedSet, vec);
}

/* A set of vectors is linearly independent if it is NOT possible to get a vector in the set as a linear combination of the other vectors in the set.
 * @param set The argument set of vectors.
 * @returns True if the argument set of vectors is linearly independent. False if otherwise.
 */
template <size_t D, size_t S>
bool is_linearly_independent(const std::array<Vector<D>, S>& set) {
    return rank(augment_vector_set(set)) == S;
}

/* A set of vectors is linearly dependent if it is possible to get a vector in the set as a linear combination of the other vectors in the set.
 * @param set The argument set of vectors.
 * @returns True if the argument set of vectors is linearly dependent. False if otherwise.
 */
template <size_t D, size_t S>
bool is_linearly_dependent(const std::array<Vector<D>, S>& set) {
    return !is_linearly_independent(set);
}

/* A homogenous system (given a coefficient matrix "A") is Ax = 0.
 * Similar to the solve() function, this will output the reduced row-echelon form of the coefficient matrix augmented with a zero vector.
 * @param coeffMat The argument coefficient matrix.
 * @returns The reduced row-echelon form of the coefficient matrix augmented with a zero vector.
 */
template <size_t R, size_t C>
Matrix<R, C + 1> solve_homogenous_system(const Matrix<R, C>& coeffMat) {
    return solve(coeffMat, zero_vector<R>());
}

//------------------------------------------------------------------------------------------//
//CHAPTER 2 - MATRICES AND LINEAR TRANSFORMATIONS

/* The product of an AxB matrix and a BxC matrix produces an AxC matrix.
 * The (i,j)-entry of the product matrix is the dot product of the ith row of the AxB matrix and the jth column of the BxC matrix.
 * @param lhs The left hand side argument matrix (of dimension AxB).
 * @param rhs The right hand side argument matrix (of dimension BxC).
 * @returns The product of the two argument matrices (of dimension AxC).
 */
//TODO: Replace with dot product
template <size_t A, size_t B, size_t C>
Matrix<A, C> operator*(const Matrix<A, B>& lhs, const Matrix<B, C>& rhs) {
    Matrix<A, C> product{};
    for (size_t row = 0; row < A; row++)
    for (size_t col = 0; col < C; col++) {
        double entry = 0.0;
        for (size_t dot = 0; dot < B; dot++)
            entry += lhs.entries[row][dot] * rhs.entries[dot][col];
        product.entries[row][col] = entry;
    }
    return product;
}

/* A matrix is diagonal if it is a square matrix and all non-diagonal entries are zero.
 * A diagonal entry is an (i,j)-entry where i = j.
 * @param m The argument matrix.
 * @returns True if the argument matrix is a diagonal matrix. False if otherwise.
 */
template <size_t R, size_t C>
bool is_diagonal(const Matrix<R, C>& m) {
    if (R != C)
        return false;

    for (size_t row = 0; row < R; row++)
    for (size_t col = 0; col < C; col++) {
        if (row == col)
            continue;

        if (!is_equal(m.entries[row][col], 0.0))
            return false;
    }
    return true;
}

/* A matrix is symmetric if (i,j)-entries equal the matrix's (j,i)-entries for all applicable i and j.
 * @param m The argument matrix.
 * @returns True if the argument matrix is symmetric. False if otherwise.
 */
template <size_t R, size_t C>
bool is_symmetric(const Matrix<R, C>& m) {
    if (R != C)
        return false;

    for (size_t row = 0; row < R; row++)
    for (size_t col = row + 1; col < C; col++)
        if (m.entries[row][col] != m.entries[R - row - 1][C - col - 1])
            return false;

    return true;
}

/* An elementary matrix has the same effect as an elementary row operation.
 * An elementary row operation applied to a matrix "A" yields the same result as EA where "E" is the corresponding elementary matrix.
 * An elementary matrix is found by applying the elementary row operation on the identity matrix.
 * @param S The size (number of rows and columns) of the elementary matrix.
 * @param row1 The first row to be exchanged.
 * @param row2 The second row to be exchanged.
 * @returns The corresponding elementary matrix.
 */
template <size_t S>
Matrix<S, S> EM_row_swap(size_t row1, size_t row2) {
    Matrix<S, S> elementaryMatrix = identity_matrix<S>();
    ERO_row_swap(elementaryMatrix, row1, row2);
    return elementaryMatrix;
}

/* An elementary matrix has the same effect as an elementary row operation.
 * An elementary row operation applied to a matrix "A" yields the same result as EA where "E" is the corresponding elementary matrix.
 * An elementary matrix is found by applying the elementary row operation on the identity matrix.
 * @param S The size (number of rows and columns) of the elementary matrix.
 * @param scalar The scalar to multiply the row by.
 * @param row The row to be multiplied.
 * @returns The corresponding elementary matrix.
 */
template <size_t S>
Matrix<S, S> EM_scalar_multiplication(double scalar, size_t row) {
    Matrix<S, S> elementaryMatrix = identity_matrix<S>();
    ERO_scalar_multiplication(elementaryMatrix, scalar, row);
    return elementaryMatrix;
}

/* An elementary matrix has the same effect as an elementary row operation.
 * An elementary row operation applied to a matrix "A" yields the same result as EA where "E" is the corresponding elementary matrix.
 * An elementary matrix is found by applying the elementary row operation on the identity matrix.
 * @param S The size (number of rows and columns) of the elementary matrix.
 * @param scalar The scalar to multiply the scaled row by.
 * @param scaledRow The row to be scaled by.
 * @param outputRow The output row to add and to copy the results into.
 * @returns The corresponding elementary matrix.
 */
template <size_t S>
Matrix<S, S> EM_row_sum(double scalar, size_t rowToScale, size_t outputRow) {
    Matrix<S, S> elementaryMatrix = identity_matrix<S>();
    ERO_row_sum(elementaryMatrix, scalar, rowToScale, outputRow);
    return elementaryMatrix;
}

/* An nxn matrix "A" is invertible if there exists a matrix "P" such that PA = AP = I (identity matrix).
 * @param m The argument matrix.
 * @returns True if the argument matrix is invertible. False if otherwise.
 */
template <size_t R, size_t C>
bool is_invertible(const Matrix<R, C>& m) {
    if (R != C)
        return false;

    return rref(m) == identity_matrix<R>();
}

/* Given a matrix "A" its inverse is "A^-1" where A*A^-1 = A^-1*A = I (identity matrix).
 * @param m The argument matrix.
 * @returns The inverse of the argument matrix if it is invertible.
 */
template <size_t S>
Matrix<S, S> inverse(const Matrix<S, S>& m) {
    assert (is_invertible(m));

    Matrix<S, 2 * S> rrefAugmented = rref(augment(m, identity_matrix<S>()));
    Matrix<S, S> inverse{};
    for (size_t row = 0; row < S; row++)
        std::copy(rrefAugmented.entries[row].begin() + S, rrefAugmented.entries[row].begin() + 2 * S, inverse.entries[row].begin());
    return inverse;
}

//------------------------------------------------------------------------------------------//
//CHAPTER 3 - DETERMINANTS

/* The determinant of a square matrix is a number that characterizes some properties of the matrix.
 * The determinant of a square matrix is 0 if and only if it is NOT invertible.
 * @param m The argument square matrix.
 * @returns The determinant of the argument square matrix.
 */
template <size_t S>
double det(const Matrix<S, S>& m) {
    Matrix<S, S> copy = m;
    unsigned int rowSwaps = ref_by_reference(copy).second;

    double determinant = 1;
    for (int i = 0; i < S; i++)
        determinant *= copy.entries[i][i];
    if (rowSwaps % 2 == 1)
        determinant *= -1;

    return determinant;
}

//------------------------------------------------------------------------------------------//
//CHAPTER 4 - SUBSPACES AND THEIR PROPERTIES

/* The dimension of a subspace is the number of vectors in its basis.
 * This function is only compatible with basis(), row(), col(), and null(), which already outputs its basis.
 * @param set The argument std::vector set.
 * @returns The size of the vector / dimension of given basis.
 */
template <size_t D>
unsigned int dim(const std::vector<Vector<D>>& set) {
    return set.size();
}

/* The row space of a matrix A is the span of its row vectors and is denoted as Row A.
 * A basis is the minimal set that spans the same space.
 * @param m The argument matrix.
 * @returns The basis of the argument matrix's row space as an std::vector.
 */
template <size_t R, size_t C>
std::vector<Vector<C>> row(const Matrix<R, C>& m) {
    Matrix<R, C> rrefM = rref(m);
    Vector<C> zeroVec = zero_vector<C>();
    std::vector<Vector<C>> rowBasis{};

    for (size_t row = 1; row <= R; row++) {
        Vector<C> rowVec = rrefM.row_vector(row);
        if (rowVec != zeroVec)
            rowBasis.emplace_back(rowVec);
    }
    return rowBasis;
}

/* The column space of matrix A is the span of its column vectors and is denoted as Col A.
 * A basis is the minimal set that spans the same space.
 * @param m The argument matrix.
 * @returns The basis of the argument matrix's column space as an std::vector.
 */
template <size_t R, size_t C>
std::vector<Vector<R>> col(const Matrix<R, C>& m) {
    Matrix<R, C> refM = ref(m);
    std::vector<Vector<R>> columnBasis{};

    size_t pivotRow = 0;
    for (size_t col = 0; col < C; col++) {
        if (pivotRow >= R)
            break;

        if (is_equal(refM.entries[pivotRow][col], 0.0))
            continue;

        columnBasis.emplace_back(m.column_vector(col + 1));
        pivotRow++;
    }
    return columnBasis;
}

/* The null space of matrix A is the set of all possible x's that satisfy the equation Ax = 0.
 * A basis is the minimal set that spans the same space.
 * @param m The argument matrix.
 * @returns The basis of the argument matrix's null space as an std::vector.
 */
template <size_t R, size_t C>
std::vector<Vector<C>> null(const Matrix<R, C>& m) {
    Matrix<R, C> rrefM = rref(m);
    std::vector<Vector<C>> nullBasis{};

    size_t pivotRow = 0;
    for (size_t col = 0; col < C; col++) {
        if (is_equal(rrefM.entries[pivotRow][col], 1.0)) {
            pivotRow++;
            continue;
        }
        
        Vector<C> vectorInBasis{};
        for (int i = 0; i < C; i++)
            vectorInBasis.components[i] = -rrefM.entries[i][col];
        vectorInBasis.components[col] = 1.0;
        nullBasis.emplace_back(vectorInBasis);
    }

    if (nullBasis.empty())
        nullBasis.emplace_back(Vector<C>()); /* zero vector */
    return nullBasis;
}

/* A basis is the minimal set that spans the same space.
 * @param set A user-defined set given as an std::array.
 * @returns The basis of the argument set as an std::vector.
 */
template <size_t D, size_t S>
std::vector<Vector<D>> basis(const std::array<Vector<D>, S>& set) {
    return col(augment_vector_set(set));
}

//------------------------------------------------------------------------------------------//
//CHAPTER 6 - ORTHOGONALITY

/* Many functions from Chapter 5 require functions from Chapter 6.
 * My Intro to Linear Algebra course covered only up until section 6.2. This chapter will cover that and QR Decomposition.
 */

/* The dot product of two vectors are the sum of squares of corresponding entries.
 * @param v1 The first argument vector.
 * @param v2 The second argument vector.
 * @returns The dot product of the two argument vectors.
 */
template <size_t D>
double dot(const Vector<D>& v1, const Vector<D>& v2) {
    return std::inner_product(v1.components.begin(), v1.components.end(), v2.components.begin(), 0.0);
}

/* The norm of a vector is its length. A vector v's norm is denoted as || v ||.
 * @param v The argument vector.
 * @returns The norm of the vector.
 */
template <size_t D>
double norm(const Vector<D>& v) {
    return sqrt(dot(v, v));
}

/* Normalizing a vector changes its norm/length to 1 (unit vector) while preserving its direction.
 * @param v The argument vector to be normalized.
 */
template <size_t D>
void normalize(Vector<D>& v) {
    v /= norm(v);
}

/* The distance between two vectors is the norm of their difference.
 * @param v1 The first argument vector.
 * @param v2 The second argument vector.
 * @returns The distance between the two vectors.
 */
template <size_t D>
double distance(const Vector<D>& v1, const Vector<D>& v2) {
    return norm(v1 - v2);
}

/* Two vectors are orthogonal if they are perpendicular when geometrically represented.
 * Two vectors are orthogonal if their dot product is 0.
 * @param v1 The first argument vector.
 * @param v2 The second argument vector.
 * @returns True if the two argument vectors are perpendicular. False if otherwise.
 */
template <size_t D>
bool is_orthogonal(const Vector<D>& v1, const Vector<D>& v2) {
    return is_equal(dot(v1, v2), 0.0);
}

/* CHAPTER 5:
 * Eigenvalue given square matrix and vector
 * Eigenspace basis given square matrix and eigenvalue
 * All eigenvalues of a matrix.
 * Is matrix diagonalizeable
 * If diagonalizeable, find invertible matrix P and diagonal matrix D
 */

} // namespace ila