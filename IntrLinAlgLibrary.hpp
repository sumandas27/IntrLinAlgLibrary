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
 * Vectors will be represented as column vectors as row vectors are rarely used in Linear Algebra.
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

    private:
        double* rowPtr;

    public:
        Proxy(double* _rowPtr) : rowPtr(_rowPtr) { };

        double&       operator[](unsigned int col);
        double const& operator[](unsigned int col) const;
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
            std::cout << (abs(entries.at(i * C + j)) < epsilon() ? 0 : entries.at(i * C + j)) << "\t";
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
//bool operator==(const Vector& v1, const Vector& v2);
/* Two matrices are equal if all corresponding entries are equal.
 * @returns true if the matrix arguments m1 and m2 are equal, false if otherwise.
 */
//bool operator==(const Matrix& m1, const Matrix& m2);
/* The sum of two vectors is a vector of the same size with corresponding components added.
 * @returns A vector that is the sum of two argument vectors.
 */
//Vector operator+(const Vector& v1, const Vector& v2);
/* The sum of two matrices is a matrix of the same size with corresponding entries added.
 * @returns A matrix that is the sum of two argument matrices.
 */
//Matrix operator+(const Matrix& m1, const Matrix& m2);
/* The difference of two vectors is a vector of the same size with corresponding components subtracted.
 * @returns A vector that is the difference of two argument vectors.
 */
//Vector operator-(const Vector& v1, const Vector& v2);
/* The difference of two matrices is a matrix of the same size with corresponding entries subtracted.
 * @returns A matrix that is the difference of two argument matrices.
 */
//Matrix operator-(const Matrix& m1, const Matrix& m2);
/* The product of a scalar and a vector is a vector of the same size with all its components multiplied by the scalar.
 * @returns A vector that is the product of a scalar and a vector.
 */
//Vector operator*(double scalar, const Vector& v);
/* Scalar-vector multiplication is commutative.
 * @returns A vector that is the product of a vector and a scalar.
 */
//Vector operator*(const Vector& v, double scalar);
/* The product of a scalar and a matrix is a matrix of the same size with all its entries multiplied by the scalar.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
//Matrix operator*(double scalar, const Matrix& m);
/* Scalar-matrix multiplication is commutative.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
//Matrix operator*(const Matrix& m, double scalar);
/* A matrix-vector product is the linear combination of the vector's components and the matrix's column vectors.
 * @returns The matrix-vector product of the argument matrix and vector.
 */
//Vector operator*(const Matrix& m, const Vector& v);

/* A zero vector is a vector where all components are zero.
 * @param dim The dimension of the zero matrix.
 * @returns A zero vector of the argument size.
 */
//Vector zero_vector(unsigned int dim);
/* A zero matrix is a matrix where all entries are zero.
 * @param rows The number of rows in the zero matrix.
 * @param cols The number of columns in the zero matrix.
 * @returns A zero matrix of the argument size.
 */
//Matrix zero_matrix(unsigned int rows, unsigned int cols);

/* A standard vector is a zero vector with one component being a one instead of a zero.
 * @param dim The dimension of the standard vector.
 * @param one_component The location of the "one" component.
 * @returns A standard vector of the argument dimension with the 1 in the argument location.
 */
//Vector standard_vector(unsigned int dim, unsigned int one_component);
/* An identity matrix is a square zero matrix with diagonal entries being a one instead of a zero.
 * @param size The number of rows and columns of the identity matrix.
 * @returns An identity matrix of the argument size.
 */
//Matrix identity_matrix(unsigned int size);

/* A rotation matrix is a 2x2 matrix that rotates an R^2 vector by some amount of degrees counter-clockwise.
 * Given rotation matrix A and some R^2 vector x, Ax = x' where x' is the rotated R^2 vector.
 * @param degrees The angle (in degrees) of the rotation matrix.
 * @returns The 2x2 rotation matrix of the argument angle in degrees.
 */
//Matrix rotation_matrix(double degrees);

/* A matrix is a square if it has the same number of rows and columns.
 * @param m The matrix argument.
 * @returns true if the matrix argument is a square, false if otherwise.
 */
//bool is_square(const Matrix& m);

/* The transpose of an nxm matrix is an mxn matrix where (i,j)-entries are transformed to (j,i)-entries.
 * @param m The matrix whose transpose is to be returned.
 * @returns The transpose of the argument matrix.
 */
//Matrix transpose(const Matrix& m);

/* An elementary row operation where two rows are exchanged in a matrix: row1 <--> row2
 * @param m The matrix to be modified.
 * @param row1 The first row to be exchanged.
 * @param row2 The second row to be exchanged.
 */
//void ERO_row_swap(Matrix& m, unsigned int row1, unsigned int row2);
/* An elementary row operation where a row is multiplied by a constant in a matrix: scalar * row --> row
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the row by.
 * @param row The row to be multiplied.
 */
//void ERO_scalar_multiplication(Matrix& m, double scalar, unsigned int row);
/* An elementary row where a multiple of one row is added to another row in a matrix: scalar * scaledRow + outputRow --> outputRow
 * @param m The matrix to be modified.
 * @param scalar The scalar to multiply the scaled row by.
 * @param scaledRow The row to be scaled by.
 * @param outputRow The output row to add and to copy the results into.
 */
//void ERO_row_sum(Matrix& m, double scalar, unsigned int rowToScale, unsigned int outputRow);

/*//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

bool operator==(const Vector& v1, const Vector& v2) {
    if (v1.get_dim() != v2.get_dim())
        return false;

    return std::equal(v1.components.begin(), v1.components.end(), v2.components.begin(), is_equal);
}

bool operator==(const Matrix& m1, const Matrix& m2) {
    if (m1.get_rows() != m2.get_rows() || m1.get_cols() != m2.get_cols())
        return false;

    return std::equal(m1.entries.begin(), m1.entries.end(), m2.entries.begin(), is_equal);
}

Vector operator+(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());

    std::vector<double> sum(v1.components.size());
    std::transform(v1.components.begin(), v1.components.end(), v2.components.begin(), sum.begin(), std::plus<double>());
    return Vector(v1.get_dim(), sum);
}

Matrix operator+(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());

    std::vector<double> sum(m1.entries.size());
    std::transform(m1.entries.begin(), m1.entries.end(), m2.entries.begin(), sum.begin(), std::plus<double>());
    return Matrix(m1.get_rows(), m1.get_cols(), sum);
}

Vector operator-(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());

    std::vector<double> diff(v1.components.size());
    std::transform(v1.components.begin(), v1.components.end(), v2.components.begin(), diff.begin(), std::minus<double>());
    return Vector(v1.get_dim(), diff);
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());

    std::vector<double> diff(m1.entries.size());
    std::transform(m1.entries.begin(), m1.entries.end(), m2.entries.begin(), diff.begin(), std::minus<double>());
    return Matrix(m1.get_rows(), m1.get_cols(), diff);
}

Vector operator*(double scalar, const Vector& v) {
    std::vector<double> product(v.components.size());
    std::transform(v.components.begin(), v.components.end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Vector(v.get_dim(), product);
}

Vector operator*(const Vector& v, double scalar) {
    return scalar * v;
}

Matrix operator*(double scalar, const Matrix& m) {
    std::vector<double> product(m.entries.size());
    std::transform(m.entries.begin(), m.entries.end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Matrix(m.get_rows(), m.get_cols(), product);
}

Matrix operator*(const Matrix& m, double scalar) {
    return scalar * m;
}

Vector operator*(const Matrix& m, const Vector& v) {
    assert (m.get_cols() == v.get_dim());

    std::vector<double> product(m.get_rows());
    for (unsigned int i = 0; i < m.get_rows(); i++) {
        double entry = 0.0;
        for (unsigned int j = 0; j < v.get_dim(); j++)
            entry += v.components.at(j) * m.entries.at(i * m.get_cols() + j);
        product.at(i) = entry;
    }
    return Vector(m.get_rows(), product);
}

Vector zero_vector(unsigned int dim) {
    std::vector<double> zeroVector(dim, 0.0);
    return Vector(dim, zeroVector);
}

Matrix zero_matrix(unsigned int rows, unsigned int cols) {
    std::vector<double> zeroMatrix(rows * cols, 0.0);
    return Matrix(rows, cols, zeroMatrix);
}

Vector standard_vector(unsigned int dim, unsigned int one_component) {
    assert (one_component >= 1 && one_component <= dim);
    
    std::vector<double> standardVector(dim, 0.0);
    standardVector.at(one_component - 1) = 1.0;
    return Vector(dim, standardVector);
}

Matrix identity_matrix(unsigned int size) {
    std::vector<double> identityMatrix(size * size, 0.0);
    for (unsigned int i = 0; i < size; i++)
        identityMatrix.at(i * size + i) = 1.0;
    return Matrix(size, size, identityMatrix);
}

Matrix transpose(const Matrix& m) {
    unsigned int tRows = m.get_cols();
    unsigned int tCols = m.get_rows();
    std::vector<double> transpose(tRows * tCols);
    for (unsigned int i = 0; i < tRows; i++)
    for (unsigned int j = 0; j < tCols; j++)
        transpose.at(i * tCols + j) = m.entries.at(j * tRows + i);
    return Matrix(tRows, tCols, transpose);
}

bool is_square(const Matrix& m) {
    return m.get_rows() == m.get_cols();
}

Matrix rotation_matrix(double degrees) {
    std::vector<double> rotationMatrix{
         cos(deg_to_rad(degrees)),
        -sin(deg_to_rad(degrees)),
         sin(deg_to_rad(degrees)),
         cos(deg_to_rad(degrees))
    };
    return Matrix(2, 2, rotationMatrix);
}

void ERO_row_swap(Matrix& m, unsigned int row1, unsigned int row2) {
    assert (row1 >= 1 && row1 <= m.get_rows());
    assert (row2 >= 1 && row2 <= m.get_rows());
    assert (row1 != row2);

    for (unsigned int col = 0; col < m.get_cols(); col++)
        std::swap(m.entries.at((row1 - 1) * m.get_cols() + col), m.entries.at((row2 - 1) * m.get_cols() + col));
}

void ERO_scalar_multiplication(Matrix& m, double scalar, unsigned int row) {
    assert (row >= 1 && row <= m.get_rows());
    
    unsigned int beg = (row - 1) * m.get_cols();
    unsigned int end = (row - 1) * m.get_cols() + m.get_cols();
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, m.entries.begin() + beg, std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
}

void ERO_row_sum(Matrix& m, double scalar, unsigned int rowToScale, unsigned int outputRow) {
    assert (rowToScale >= 1 && rowToScale <= m.get_rows());
    assert (outputRow >= 1 && outputRow <= m.get_rows());
    assert (rowToScale != outputRow);

    std::vector<double> scaledRow(m.get_cols());
    unsigned int beg = (rowToScale - 1) * m.get_cols();
    unsigned int end = (rowToScale - 1) * m.get_cols() + m.get_cols();
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));

    beg = (outputRow - 1) * m.get_cols();
    end = (outputRow - 1) * m.get_cols() + m.get_cols();
    std::transform(m.entries.begin() + beg, m.entries.begin() + end, scaledRow.begin(), m.entries.begin() + beg, std::plus<double>());
}*/