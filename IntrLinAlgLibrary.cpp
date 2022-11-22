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

std::vector<double> Vector::get_vector() const {
    return components;
}

double Vector::get(unsigned int index) const {
    assert (index >= 1 && index <= dim);
    return components.at(index - 1);
}

void Vector::set(unsigned int index, double value) {
    assert (index >= 1 && index <= dim);
    components.at(index - 1) = value;
}

void Vector::print() {
    for (double& component : components)
        std::cout << "{  " << std::to_string(component).substr(0, 7) << "  }\n";
    std::cout << "\n";
}

//----------------------------------------------------------------------//
//MATRIX METHODS:

Matrix::Matrix(unsigned int _rows, unsigned int _cols, std::vector< std::vector<double> >& _entries) : rows(_rows), cols(_cols) {
    assert (_entries.size() == _rows);
    for (const std::vector<double>& col : _entries)
        assert (col.size() == _cols);

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

std::vector< std::vector<double> > Matrix::get_matrix() const {
    return entries;
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

void Matrix::print() {
    for (std::vector<double>& row : entries) {
        std::cout << "{  ";
        for (double& entry : row)
            std::cout << std::to_string(entry).substr(0, 7) << "  ";
        std::cout << "}\n";
    }
    std::cout << "\n";
}

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

bool operator==(const Vector& v1, const Vector& v2) {
    return v1.get_vector() == v2.get_vector();
}

bool operator==(const Matrix& m1, const Matrix& m2) {
    for (unsigned int row = 1; row <= m1.get_rows(); row++)
        if (m1.get_matrix().at(row - 1) != m2.get_matrix().at(row - 1))
            return false;

    return true;
}

bool is_square(const Matrix& m) {
    return m.get_rows() == m.get_cols();
}