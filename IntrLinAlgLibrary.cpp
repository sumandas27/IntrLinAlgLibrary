#include "IntrLinAlgLibrary.hpp"

//----------------------------------------------------------------------//
//VECTOR METHODS:

Vector::Vector(unsigned int _dim, std::vector<double>& _components) : dim(_dim) {
    assert (_components.size() == _dim);

    components.reserve(dim);
    components.insert(components.end(), _components.begin(), _components.end());
}

unsigned int Vector::get_dim() const {
    return dim;
}

double Vector::get(unsigned int index) const {
    assert (index >= 1 && index <= dim);
    return components.at(index - 1);
}

void Vector::set(unsigned int index, double value) {
    assert (index >= 1 && index <= dim);
    components.at(index - 1) = value;
}

//----------------------------------------------------------------------//
//MATRIX METHODS:

Matrix::Matrix(unsigned int _rows, unsigned int _cols, std::vector< std::vector<double> >& _entries) : rows(_rows), cols(_cols) {
    assert (_entries.size() == _rows);
    for (const std::vector<double>& col : _entries)
        assert(col.size() == _cols);

    entries.reserve(rows);
    entries.resize(rows);
    for (int i = 0; i < rows; i++) {
        entries.at(i).reserve(cols);
        entries.at(i).insert(entries.at(i).end(), _entries.at(i).begin(), _entries.at(i).end());
    }
}

unsigned int Matrix::get_rows() const {
    return rows;
}

unsigned int Matrix::get_cols() const {
    return cols;
}

double Matrix::get(unsigned int row, unsigned int col) const {
    assert (row >= 1 && row <= rows);
    assert (col >= 1 && col <= cols);
    return entries.at(row - 1).at(col - 1);
}

void Matrix::set(unsigned int row, unsigned int col, double value) {
    assert (row >= 1 && row <= rows);
    assert (col >= 1 && col <= cols);
    entries.at(row - 1).at(col - 1) = value;
}

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

bool operator==(const Vector& v1, const Vector& v2) {
    if (v1.get_dim() != v2.get_dim())
        return false;

    for (unsigned int index = 1; index <= v1.get_dim(); index++)
        if (v1.get(index) != v2.get(index))
            return false;

    return true;
}

bool operator==(const Matrix& m1, const Matrix& m2) {
    if (m1.get_rows() != m2.get_rows() || m1.get_cols() != m2.get_cols())
        return false;

    for (unsigned int row = 1; row <= m1.get_rows(); row++)
    for (unsigned int col = 1; col <= m1.get_cols(); col++)
        if (m1.get(row, col) != m2.get(row, col))
            return false;

    return true;
}

bool is_square(const Matrix& m) {
    return m.get_rows() == m.get_cols();
}