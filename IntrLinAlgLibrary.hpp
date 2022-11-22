#pragma once
#include <iostream>
#include <cassert>

//--- Containers ---//
#include <vector>

//TODO NEXT: implement setters

// A vector is an array of real numbers.
struct Vector {
    // The dimension of the vector.
    const unsigned int dim;
    // The structure containing the components of the vector.
    std::vector<double> components;

    /* Constructs a vector.
     * @param _dim The input vector dimension.
     * @param _components The input structure containing the components of the vector.
     */
    Vector(unsigned int _dim, std::vector<double>& _components);
    
    /* Vectors in linear algebra are traditionally one-indexed. 
     * @param index The index of the desired scalar.
     * @returns The scalar at the index argument.
     */
    double get(unsigned int index);
};

// A matrix is an array of arrays of real numbers.
struct Matrix {
    // The number of rows in the matrix.
    const unsigned int rows;
    // The number of columns in the matrix.
    const unsigned int cols;

    // The structure containing the entries of the matrix.
    std::vector< std::vector<double> > entries;

    /* Constructs a matrix.
     * @param _rows The input number of rows in the matrix.
     * @param _cols The input number of columns in the matrix.
     * @param _entries The input structure containing the entries of the vector.
     */
    Matrix(unsigned int _rows, unsigned int _cols, std::vector< std::vector<double> >& _entries);

    /* Matrices (like vectors) in linear algebra are traditionally one-indexed.
     * @param row The row of the desired scalar.
     * @param col The column of the desired scalar.
     */
    double get(unsigned int row, unsigned int col);
};