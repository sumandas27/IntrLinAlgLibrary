#pragma once
#include <iostream>
#include <cassert>

//--- Containers ---//
#include <vector>

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

    /* Returns the dimension of the vector.
     */
    unsigned int get_dim();
    
    /* Vectors in linear algebra are traditionally one-indexed. 
     * @param index The index of the desired scalar.
     * @returns The scalar at the index argument.
     */
    double get(unsigned int index);
    /* Sets the components at the index argument to the desired scalar.
     * @param index The (1-based) component of the vector to be set.
     * @param value The scalar to set the desired component to. 
     */
    void set(unsigned int index, double value);
};

// A matrix is an array of arrays of real numbers.
class Matrix {

private:
    // The number of rows in the matrix.
    const unsigned int rows;
    // The number of columns in the matrix.
    const unsigned int cols;

    // The structure containing the entries of the matrix.
    std::vector< std::vector<double> > entries;

public:
    /* Constructs a matrix.
     * @param _rows The input number of rows in the matrix.
     * @param _cols The input number of columns in the matrix.
     * @param _entries The input structure containing the entries of the vector.
     */
    Matrix(unsigned int _rows, unsigned int _cols, std::vector< std::vector<double> >& _entries);
    /* A matrix object shouldn't be created with the default constructor.
     */
    Matrix() = delete;

    /* Returns the number of rows in the matrix.
     */
    unsigned int get_rows();
    /* Returns the number of columns in the matrix.
     */
    unsigned int get_cols();

    /* Matrices (like vectors) in linear algebra are traditionally one-indexed.
     * @param row The row of the desired scalar.
     * @param col The column of the desired scalar.
     */
    double get(unsigned int row, unsigned int col);
    /* Sets the matrix entry at the row and column argument to the desired scalar.
     * @param row The row of the matrix entry to be set.
     * @param col The column of the matrix entry to be set.
     * @param value The scalar to set the desired matrix entry to.
     */
    void set(unsigned int row, unsigned int col, double value);
};