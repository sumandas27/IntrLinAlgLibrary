#pragma once
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>

//--- Standard Template Library ---//
#include <algorithm>
#include <functional>
#include <initializer_list>
#include <array>

//TODO: Fix all methods lol

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
template <size_t D>
struct Vector {
    std::array<double, D> components;

    template <typename... Ts>
    Vector(Ts... _components);
    Vector();

    double operator[](size_t index);

    void print() const;

    template <size_t X>
    friend std::ostream& operator<<(std::ostream& os, const Vector<X>& v);
};

/* Constructs a vector with a list of arguments. Example Vector object creation: 
 * Vector<5> vec(1.0, 2.0, 3.0, 4.0, 5.0);
 * 
 * @param _components The list of scalar arguments containing the components of the vector.
 */
template <size_t D>
template <typename... Ts>
Vector<D>::Vector(Ts... _components) : components{ (double)_components... } {
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

/* Prints the vector to the terminal. Each component is rounded to the nearest thousandth.
 */
template <size_t D>
void Vector<D>::print() const {
    for (double component : components)
        std::cout << "{\t" << (abs(component) < epsilon() ? 0 : component) << "\t}\n";
    std::cout << "\n";
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

//----------------------------------------------------------------------//
//MATRIX STRUCT AND METHODS:

/* A matrix is an array of arrays of real numbers.
 * @param R The number of rows in the matrix.
 * @param C The number of columns in the matrix.
 */
template <size_t R, size_t C>
struct Matrix {
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

    template <size_t X, size_t Y>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<X, Y>& m);
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
Matrix<R, C>::Matrix(Ts... _entries) : entries{ (double)_entries... } {
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

/* Prints the matrix to the terminal. Each entry is rouded to the nearest thousandths.
 */
template <size_t R, size_t C>
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

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

/* Two vectors are equal if all corresponding components are equal.
 * @returns true if the vector arguments v1 and v2 are equal, false if otherwise.
 */
template <size_t D>
bool operator==(const Vector<D>& lhs, const Vector<D>& rhs) {
    return std::equal(lhs.components.begin(), lhs.components.end(), rhs.components.begin(), is_equal);
}

/* Two matrices are equal if all corresponding entries are equal.
 * @returns true if the matrix arguments m1 and m2 are equal, false if otherwise.
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
    return Matrix<2, 2>
    (
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

    if (row1 == row2)
        return;

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
 * @param rank This pointer will be populated with the rank of the matrix. Pass in 'nullptr' if the rank isn't wanted.
 */
template <size_t R, size_t C>
void ref_by_reference(Matrix<R, C>& m, unsigned int* rank) {
    size_t pivotRow = 0;
    for (size_t col = 0; col < C; col++) {
        if (pivotRow > R - 1)
            break;

        double nonzeroFound = 0.0; /* 0.0 means a nonzero entry has not been found, else nonzero entry is set to this variable */
        for (size_t row = pivotRow; row < R; row++) {
            if (is_equal(m.entries[row * C + col], 0.0))
                continue;

            if (nonzeroFound == 0.0) {
                nonzeroFound = m.entries[row * C + col];
                ERO_row_swap(m, pivotRow + 1, row + 1);
                pivotRow++;
            }
            else {
                double scalar = -m.entries[row * C + col] / nonzeroFound;
                ERO_row_sum(m, scalar, pivotRow, row + 1);
            }
        }
    }

    if (rank != nullptr)
        *rank = pivotRow;
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
    ref_by_reference(m, nullptr);
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
    ref_by_reference(m, nullptr);
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
    unsigned int rank = 0;
    Matrix<R, C> copy = m;
    ref_by_reference(copy, &rank);
    return rank;
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
 * This is synonymous to 0x1 + 0x2 + ... + 0xn = c --> 0 = c. 0 cannot equal a nonzero number. Therefore, no solution exists for x and the system is inconsistent.
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

/*template <size_t D, size_t S>
Matrix<D, S> augment_vector_set(const std::array<Vector<D>, S>& set) {
    assert (S != 0);

    std::array<double, D * S> augmentedVectorSet
    for (size_t col = 0; col < S; col++)
    for (size_t row = 0; row < D; row++) {

    }
}*/