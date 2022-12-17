#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>

#include <cmath>
#include <string>

#include <algorithm>
#include <functional>
#include <array>

//------------------------------------------------------------------------------------------//
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
    Vector();

    double operator[](size_t index);

    void print() const;
};

/* Constructs a vector with a list of arguments. Example Vector object creation: 
 * Vector<5> vec(1.0, 2.0, 3.0, 4.0, 5.0);
 * 
 * @param _components The list of scalar arguments containing the components of the vector.
 */
template <size_t D>
template <typename... Ts>
Vector<D>::Vector(Ts... _components) 
    : components{ (double)std::forward<Ts>(_components)... } 
{
    assert (D != 0);
    assert (sizeof...(_components) == D);
}

/* Constructs a vector with zero-initialized components.
 */
template <size_t D>
Vector<D>::Vector() : components{} { }

/* Vector indexing is 1-based.
 * @param index The index of the vector wanted to be changed/accessed.
 * @returns The component at the argument index.
 */
template <size_t D>
double Vector<D>::operator[](size_t index) {
    assert (index >= 1 && index <= D);
    return components[index - 1];
}

/* Alternative to the print() function, allows the vector to be printed directly to the standard output.
 * Given a vector "v": "std::cout << v;"
 * @param v The vector to be printed to the standard output.
 */
template <size_t X>
std::ostream& operator<<(std::ostream& os, const Vector<X>& v) {
    for (double component : v.components)
        os << "{\t" << (abs(component) < epsilon() ? 0 : component) << "\t}\n";
    os << "\n";
    return os;
}

/* Prints the vector to the terminal. Each component is rounded to the nearest thousandth.
 */
template <size_t D>
void Vector<D>::print() const {
    std::cout << *this;
}

//------------------------------------------------------------------------------------------//
//VECTORSET STRUCT:

//TODO: Remove/rework this

/* A VectorSet holds a set of vectors.
 * @param D The dimension of the vectors in the set.
 * @param S The number of vectors or "size" in the set. 
 */
template <size_t D, size_t S>
class VectorSet {
public:
    std::array<Vector<D>, S> set;

    template <typename... Ts>
    VectorSet(Ts&&... _set);

    Vector<D>&       operator[](size_t index);
    Vector<D> const& operator[](size_t index) const;
};

/* Constructs a VectorSet object.
 * Example VectorSet object creation:
 * VectorSet<3, 2> vectorSet(
 *     Vector<3>(1.0, 2.0, 3.0),
 *     Vector<3>(4.0, 5.0, 6.0)
 * );
 * @param _set The list holding the vectors in the set.
 */
template <size_t D, size_t S>
template <typename... Ts>
VectorSet<D, S>::VectorSet(Ts&&... _set) : set{_set...} {
    assert (D != 0);
    assert (sizeof...(_set) == S);
}

/* VectorSets are accessed like normal arrays; indexing is 0-based. This overload is primarily intended for the left-hand side of an assignment.
 * @param index The index of the vector wanted to be accessed.
 * @returns The vector at the argument index.
 */
template <size_t D, size_t S>
Vector<D>& VectorSet<D, S>::operator[](size_t index) {
    assert (index >= 0 && index < S);
    return set[index];
}

/* VectorSets are accessed like normal arrays; indexing is 0-based. This overload is exclusively intended for the right-hand side of an assingment.
 * @param index The index of the vector wanted to be accessed.
 * @returns The vector at the argument index.
 */
template <size_t D, size_t S>
Vector<D> const& VectorSet<D, S>::operator[](size_t index) const {
    assert (index >= 0 && index < S);
    return set[index];
}

//------------------------------------------------------------------------------------------//
//MATRIX STRUCT AND METHODS:

//TODO: Change this to an array of arrays.

/* A matrix is an array of arrays of real numbers.
 * @param R The number of rows in the matrix.
 * @param C The number of columns in the matrix.
 */
template <size_t R, size_t C>
class Matrix {
public:
    std::array<double, R * C> entries;

    template <typename... Ts>
    Matrix(Ts... _entries);
    Matrix();

    class Proxy {
    public:
        Proxy(double* _rowPtr) : rowPtr(_rowPtr) { };
        double operator[](size_t col);

    private:
        double* rowPtr;
    };

    Proxy operator[](size_t row);

    void print() const;
};

/* Constructs a matrix. Example Matrix object creation:
 * Matrix<2, 3> mat(
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * );
 *
 * @param _entries The list of scalar arguments containing the entries of the matrix.
 */
template <size_t R, size_t C>
template <typename... Ts>
Matrix<R, C>::Matrix(Ts... _entries) 
    : entries{ (double)std::forward<Ts>(_entries)... } 
{
    assert (R != 0 && C != 0);
    assert (sizeof...(_entries) == R * C);
}

/* Constructs a Matrix object with zero-intialized entries.
 */
template <size_t R, size_t C>
Matrix<R, C>::Matrix() : entries{} { }

/* Accesses the argument column of the proxy's row.
 * This overload is intended for the left hand side of an assignment.
 * @param col The column index of the proxy row wanted to be changed/accessed.
 * @returns The entry at the argument column of the proxy row.
 */
template <size_t R, size_t C>
double Matrix<R, C>::Proxy::operator[](size_t col) {
    assert (col >= 1 && col <= C);
    return rowPtr[col - 1];
}

/* Allows access to the entries of the matrix. Doubly overloaded subscript operators require use of an intermediate proxy class.
 * @param row The row of the entry wanted to be accessed.
 * @returns A proxy class containing a pointer to the first entry of the row.
 * The proxy's subscript overloaders can then be used to access specific entries in the row.
 */
template <size_t R, size_t C>
typename Matrix<R, C>::Proxy Matrix<R, C>::operator[](size_t row) {
    assert (row >= 1 && row <= R);
    return Proxy(&entries[(row - 1) * C]);
}

/* Alternative to the print() function, allows the matrix to be printed directly to the standard output.
 * Given a matrix "m": "std::cout << m;"
 * @param m The matrix to be printed to the standard output.
 */
template <size_t X, size_t Y>
std::ostream& operator<<(std::ostream& os, const Matrix<X, Y>& m) {
    for (int i = 0; i < X; i++) {
        os << "{\t";
        for (int j = 0; j < Y; j++)
            os << (abs(m.entries[i * Y + j]) < epsilon() ? 0 : m.entries[i * Y + j]) << "\t";
        os << "}\n";
    }
    os << "\n";
    return os;
}

/* Prints the matrix to the terminal. Each entry is rouded to the nearest thousandths.
 */
template <size_t R, size_t C>
void Matrix<R, C>::print() const {
    std::cout << *this;
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

/* Two matrices are equal if all corresponding entries are equal.
 * @returns True if the matrix arguments m1 and m2 are equal. False if otherwise.
 */
template <size_t R, size_t C>
bool operator==(const Matrix<R, C>& lhs, const Matrix<R, C>& rhs) {
    return std::equal(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), is_equal);
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
    std::transform(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), sum.entries.begin(), std::plus<double>());
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
    std::transform(lhs.entries.begin(), lhs.entries.end(), rhs.entries.begin(), diff.entries.begin(), std::minus<double>());
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

/* The product of a scalar and a matrix is a matrix of the same size with all its entries multiplied by the scalar.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> operator*(double scalar, const Matrix<R, C>& m) {
    Matrix<R, C> product{};
    std::transform(m.entries.begin(), m.entries.end(), product.entries.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return product;
}

/* Scalar-matrix multiplication is commutative.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
template <size_t R, size_t C>
Matrix<R, C> operator*(const Matrix<R, C>& m, double scalar) {
    return scalar * m;
}

/* A matrix-vector product is the linear combination of the vector's components and the matrix's column vectors.
 * @returns The matrix-vector product of the argument matrix and vector.
 */
template <size_t R, size_t C>
Vector<R> operator*(const Matrix<R, C>& m, const Vector<C>& v) {
    Vector<R> product{};
    for (size_t i = 0; i < R; i++) {
        double entry = 0.0;
        for (size_t j = 0; j < C; j++)
            entry += v.components[j] * m.entries[i * C + j];
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
        identityMatrix.entries[i * S + i] = 1.0;
    return Matrix<S, S>(identityMatrix);
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
        transpose.entries[i * R + j] = m.entries[j * C + i];
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

    for (size_t col = 0; col < C; col++)
        std::swap(m.entries[(row1 - 1) * C + col], m.entries[(row2 - 1) * C + col]);
}

/* An elementary row operation where a row is multiplied by a constant in a matrix: scalar * row --> row
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the row by.
 * @param row The row to be multiplied.
 */
template <size_t R, size_t C>
void ERO_scalar_multiplication(Matrix<R, C>& m, double scalar, size_t row) {
    assert (row >= 1 && row <= R);

    size_t beg = (row - 1) * C;
    size_t end = (row - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, m.entries.begin() + beg, std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
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
    size_t beg = (rowToScale - 1) * C;
    size_t end = (rowToScale - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));

    beg = (outputRow - 1) * C;
    end = (outputRow - 1) * C + C;
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), m.entries.begin() + beg, std::plus<double>());
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
            if (is_equal(m.entries[row * C + col], 0.0))
                continue;

            if (nonzeroFound != 0.0) {
                double scalar = -m.entries[row * C + col] / nonzeroFound;
                ERO_row_sum(m, scalar, pivotRow, row + 1);
                continue;
            }

            nonzeroFound = m.entries[row * C + col]; /* First nonzero of the row found */
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
 * For implementation details, read the function "ref_by_reference(Matrix&, unsigned int*)".
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

        if (is_equal(m.entries[pivotRow * C + col], 0.0))
            continue;

        if (m.entries[pivotRow * C + col] != 1.0) {
            double scalar = 1.0 / m.entries[pivotRow * C + col];
            ERO_scalar_multiplication(m, scalar, pivotRow + 1);
        }

        for (int row = pivotRow - 1; row >= 0; row--)
            if (m.entries[row * C + col] != 0) {
                double scalar = -m.entries[row * C + col];
                ERO_row_sum(m, scalar, pivotRow + 1, row + 1);
            }

        pivotRow++;
    }
    return m; 
}

/* The rank of a matrix is the dimension of the vector space generated by its column.
 * This is equivalent to the number of pivot rows in the matrix's row-echelon form (which is how this method calculates the rank).
 * @param m The argument matrix.
 * @returns The rank of the argument matrix.
 */
template <size_t R, size_t C>
unsigned int rank(const Matrix<R, C>& m) {
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
        augmentedMatrix.entries[row * (C + 1) + col] = (col == C) ? v.components[row] : m.entries[row * C + col];
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
        augmentedMatrix.entries[row * (C1 + C2) + col] = (col < C1) ? lhs.entries[row * C1 + col] : rhs.entries[row * C2 + (col - C1)];
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
        if (rrefAugment.entries[row * (C + 1) + C] == 0.0)
            continue;

        size_t beg = row * (C + 1);
        size_t end = row * (C + 1) + C;
        if (std::all_of(rrefAugment.entries.begin() + beg, rrefAugment.entries.begin() + end, [](double entry){ return entry == 0.0; }))
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
Matrix<D, S> augment_vector_set(const VectorSet<D, S>& set) {
    Matrix<D, S> augmentedVectorSet{};
    for (size_t col = 0; col < S; col++)
    for (size_t row = 0; row < D; row++)
        augmentedVectorSet.entries[row * S + col] = set.set[col].components[row];
    return augmentedVectorSet;
}

/* A vector is in the span of a set of vectors if the vector can be written as a linear combination of the other vectors.
 * @param vec The argument vector.
 * @param set The argument set of vectors.
 * @returns True if the argument vector is in the span of the argument set of vectors. False if otherwise.
 */
template <size_t D, size_t S>
bool is_in_span(const Vector<D>& vec, const VectorSet<D, S>& set) {
    Matrix<D, S> augmentedSet = augment_vector_set(set);
    return is_consistent(augmentedSet, vec);
}

/* A set of vectors is linearly independent if it is NOT possible to get a vector in the set as a linear combination of the other vectors in the set.
 * @param set The argument set of vectors.
 * @returns True if the argument set of vectors is linearly independent. False if otherwise.
 */
template <size_t D, size_t S>
bool is_linearly_independent(const VectorSet<D, S>& set) {
    return rank(augment_vector_set(set)) == S;
}

/* A set of vectors is linearly dependent if it is possible to get a vector in the set as a linear combination of the other vectors in the set.
 * @param set The argument set of vectors.
 * @returns True if the argument set of vectors is linearly dependent. False if otherwise.
 */
template <size_t D, size_t S>
bool is_linearly_dependent(const VectorSet<D, S>& set) {
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
template <size_t A, size_t B, size_t C>
Matrix<A, C> operator*(const Matrix<A, B>& lhs, const Matrix<B, C>& rhs) {
    Matrix<A, C> product{};
    for (size_t row = 0; row < A; row++)
    for (size_t col = 0; col < C; col++) {
        double entry = 0.0;
        for (size_t dot = 0; dot < B; dot++)
            entry += lhs.entries[row * B + dot] * rhs.entries[dot * C + col];
        product.entries[row * C + col] = entry;
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

        if (!is_equal(m.entries[row * C + col], 0.0))
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
    for (size_t col = row + 1; col < C; col++) {
        size_t original = row * C + col;
        size_t symmetric = (R - row - 1) * C + (C - col - 1);
        if (m.entries[original] != m.entries[symmetric])
            return false;
    }
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
    for (size_t row = 0; row < S; row++) {
        size_t beg = row * (2 * S) + S;
        size_t end = row * (2 * S) + (2 * S);
        std::copy(rrefAugmented.entries.begin() + beg, rrefAugmented.entries.begin() + end, inverse.entries.begin() + row * S);
    }
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
        determinant *= copy.entries[i * S + i];
    if (rowSwaps % 2 == 1)
        determinant *= -1;

    return determinant;
}

//------------------------------------------------------------------------------------------//
//CHAPTER 4 - SUBSPACES AND THEIR PROPERTIES

//TODO (after refactoring):

/* CHAPTER 4:
 * Row Space
 * Column Space
 * Null Space
 * Dimension of Space??
 * Are two matrices similar
 */

/* CHAPTER 5:
 * Eigenvalue given square matrix and vector
 * Eigenspace basis given square matrix and eigenvalue
 * All eigenvalues of a matrix.
 * Is matrix diagonalizeable
 * If diagonalizeable, find invertible matrix P and diagonal matrix D
 */