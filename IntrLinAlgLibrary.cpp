#include "IntrLinAlgLibrary.hpp"

//----------------------------------------------------------------------//
//VECTOR METHODS:

Vector::Vector(unsigned int _dim, std::vector<double>& _components) : dim(_dim) {
    components.reserve(dim);
    components.insert(components.end(), _components.begin(), _components.end());
}

double Vector::get(unsigned int index) {
    assert (index >= 1 && index <= dim);
    return components.at(index - 1);
}

//----------------------------------------------------------------------//
//MATRIX METHODS:

Matrix::Matrix(unsigned int _rows, unsigned int _cols, std::vector< std::vector<double> >& _entries) : rows(_rows), cols(_cols) {
    entries.reserve(rows);
    entries.resize(rows);
    for (int i = 0; i < rows; i++) {
        entries.at(i).reserve(cols);
        entries.at(i).insert(entries.at(i).end(), _entries.at(i).begin(), _entries.at(i).end());
    }
}

double Matrix::get(unsigned int row, unsigned int col) {
    assert (row >= 1 && row <= rows);
    assert (col >= 1 && col <= cols);
    return entries.at(row - 1).at(col - 1);
}