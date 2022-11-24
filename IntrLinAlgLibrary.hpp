#pragma once
#include <iostream>
#include <cassert>
#include <cmath>
#include <string>

//--- Standard Template Library ---//
#include <algorithm>
#include <functional>
#include <vector>

//TODO: Rotation Matrix
//TODO: Matrix-Vector Product

/* A vector is an array of real numbers.
 * Vectors will be represented as column vectors as row vectors are rarely used in Linear Algebra.
 */
class Vector {

private:
    // The dimension of the vector.
    const unsigned int dim;
    // The structure containing the components of the vector.
    std::vector<double> components;

public:
    /* Constructs a vector.
     * @param _dim The input vector dimension.
     * @param _components The input structure containing the components of the vector.
     */
    Vector(unsigned int _dim, std::vector<double>& _components);
    /* A vector object shouldn't be created with the default constructor.
     */
    Vector() = delete;

    /* Gets the dimension of the vector.
     */
    unsigned int get_dim() const;
    /* Gets the structure containing the components of the vector.
     */
    std::vector<double> get_components() const;
    
    /* Vectors in linear algebra are traditionally one-indexed. 
     * @param index The index of the desired scalar.
     * @returns The scalar at the index argument.
     */
    double get(unsigned int index) const;
    /* Sets the components at the index argument to the desired scalar.
     * @param index The (1-based) component of the vector to be set.
     * @param value The scalar to set the desired component to. 
     */
    void set(unsigned int index, double value);

    /* Prints the vector to the terminal.
     * Only the first 8 characters of every component are printed.
     */
    void print();
};

// A matrix is an array of arrays of real numbers.
class Matrix {

private:
    // The number of rows in the matrix.
    const unsigned int rows;
    // The number of columns in the matrix.
    const unsigned int cols;

    // The structure containing the entries of the matrix.
    std::vector<double> entries;

public:
    /* Constructs a matrix.
     * @param _rows The input number of rows in the matrix.
     * @param _cols The input number of columns in the matrix.
     * @param _entries The input structure containing the entries of the vector.
     */
    Matrix(unsigned int _rows, unsigned int _cols, std::vector<double>& _entries);
    /* A matrix object shouldn't be created with the default constructor.
     */
    Matrix() = delete;

    /* Gets the number of rows in the matrix.
     */
    unsigned int get_rows() const;
    /* Gets the number of columns in the matrix.
     */
    unsigned int get_cols() const;
    /* Gets the structure (vector of vectors) containing the entries of the matrix.
     */
    std::vector<double> get_entries() const;

    /* Matrices (like vectors) in linear algebra are traditionally one-indexed.
     * @param row The row of the desired scalar.
     * @param col The column of the desired scalar.
     */
    double get(unsigned int row, unsigned int col) const;
    /* Sets the matrix entry at the row and column argument to the desired scalar.
     * @param row The row of the matrix entry to be set.
     * @param col The column of the matrix entry to be set.
     * @param value The scalar to set the desired matrix entry to.
     */
    void set(unsigned int row, unsigned int col, double value);

    /* Prints the matrix to the terminal.
     * Only the first 8 characters of every entry are printed.
     */
    void print();
};

//----------------------------------------------------------------------//
//NON-LINEAR ALGEBRA FUNCTIONS:

/* Checks equality between two doubles if their difference lies within a tolerance range.
 */
bool is_equal(double val1, double val2);

/* Converts an angle from degrees to radians.
 */
constexpr double deg_to_rad(double degrees);

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

/* Two vectors are equal if all corresponding components are equal.
 * @returns true if the vector arguments v1 and v2 are equal, false if otherwise.
 */
bool operator==(const Vector& v1, const Vector& v2);
/* Two matrices are equal if all corresponding entries are equal.
 * @returns true if the matrix arguments m1 and m2 are equal, false if otherwise.
 */
bool operator==(const Matrix& m1, const Matrix& m2);
/* The sum of two vectors is a vector of the same size with corresponding components added.
 * @returns A vector that is the sum of two argument vectors.
 */
Vector operator+(const Vector& v1, const Vector& v2);
/* The sum of two matrices is a matrix of the same size with corresponding entries added.
 * @returns A matrix that is the sum of two argument matrices.
 */
Matrix operator+(const Matrix& m1, const Matrix& m2);
/* The difference of two vectors is a vector of the same size with corresponding components subtracted.
 * @returns A vector that is the difference of two argument vectors.
 */
Vector operator-(const Vector& v1, const Vector& v2);
/* The difference of two matrices is a matrix of the same size with corresponding entries subtracted.
 * @returns A matrix that is the difference of two argument matrices.
 */
Matrix operator-(const Matrix& m1, const Matrix& m2);
/* The product of a scalar and a vector is a vector of the same size with all its components multiplied by the scalar.
 * @returns A vector that is the product of a scalar and a vector.
 */
Vector operator*(double scalar, const Vector& v);
/* Scalar-vector multiplication is commutative.
 * @returns A vector that is the product of a vector and a scalar.
 */
Vector operator*(const Vector& v, double scalar);
/* The product of a scalar and a matrix is a matrix of the same size with all its entries multiplied by the scalar.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
Matrix operator*(double scalar, const Matrix& m);
/* Scalar-matrix multiplication is commutative.
 * @returns A matrix that is the product of a scalar and a matrix.
 */
Matrix operator*(const Matrix& m, double scalar);
/* A matrix-vector product is the linear combination of the vector's components and the matrix's column vectors.
 * @returns The matrix-vector product of the argument matrix and vector.
 */
Vector operator*(const Matrix& m, const Vector& v);

/* A zero vector is a vector where all components are zero.
 * @param dim The dimension of the zero matrix.
 * @returns A zero vector of the argument size.
 */
Vector zero_vector(unsigned int dim);
/* A zero matrix is a matrix where all entries are zero.
 * @param rows The number of rows in the zero matrix.
 * @param cols The number of columns in the zero matrix.
 * @returns A zero matrix of the argument size.
 */
Matrix zero_matrix(unsigned int rows, unsigned int cols);

/* A standard vector is a zero vector with one component being a one instead of a zero.
 * @param dim The dimension of the standard vector.
 * @param one_component The location of the "one" component.
 * @returns A standard vector of the argument dimension with the 1 in the argument location.
 */
Vector standard_vector(unsigned int dim, unsigned int one_component);

/* An identity matrix is a square zero matrix with diagonal entries being a one instead of a zero.
 * @param size The number of rows and columns of the identity matrix.
 * @returns An identity matrix of the argument size.
 */
Matrix identity_matrix(unsigned int size);

/* A rotation matrix is a 2x2 matrix that rotates an R^2 vector by some amount of degrees counter-clockwise.
 * Given rotation matrix A and some R^2 vector x, Ax = x' where x' is the rotated R^2 vector.
 * @param degrees The angle (in degrees) of the rotation matrix.
 * @returns The 2x2 rotation matrix of the argument angle in degrees.
 */
Matrix rotation_matrix(double degrees);

/* The transpose of an nxm matrix is an mxn matrix where (i,j)-entries are transformed to (j,i)-entries.
 * @param m The matrix whose transpose is to be returned.
 * @returns The transpose of the argument matrix.
 */
Matrix transpose(const Matrix& m);

/* A matrix is a square if it has the same number of rows and columns.
 * @param m The matrix argument.
 * @returns true if the matrix argument is a square, false if otherwise.
 */
bool is_square(const Matrix& m);