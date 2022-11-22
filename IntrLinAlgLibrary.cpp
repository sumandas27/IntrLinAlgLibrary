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
    for (unsigned int row = 0; row < m1.get_rows(); row++)
        if (m1.get_matrix().at(row) != m2.get_matrix().at(row))
            return false;

    return true;
}

Vector operator+(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());

    std::vector<double> sum(v1.get_dim());
    std::transform(v1.get_vector().begin(), v1.get_vector().end(), v2.get_vector().begin(), sum.begin(), std::plus<double>());
    return Vector(v1.get_dim(), sum);
}


// COME BACK TO MATRIX ADD AND SUBTRACT LATER:
Matrix operator+(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());

    std::vector< std::vector<double> > sum(m1.get_rows(), std::vector<double>(m1.get_cols()));
    for (int row = 0; row < m1.get_rows(); row++) {
        std::vector<double> r1 = m1.get_matrix().at(row);
        std::vector<double> r2 = m2.get_matrix().at(row);
        std::transform(r1.begin(), r1.end(), r2.begin(), sum.at(row).begin(), std::plus<double>());
    }
    return Matrix(m1.get_rows(), m1.get_cols(), sum);
}

Vector operator-(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());

    std::vector<double> diff(v1.get_dim());
    std::transform(v1.get_vector().begin(), v1.get_vector().end(), v2.get_vector().begin(), diff.begin(), std::minus<double>());
    return Vector(v1.get_dim(), diff);
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());

   std::vector< std::vector<double> > sum(m1.get_rows(), std::vector<double>(m1.get_cols()));
    for (int row = 0; row < m1.get_rows(); row++) {
        std::vector<double> r1 = m1.get_matrix().at(row);
        std::vector<double> r2 = m2.get_matrix().at(row);
        std::transform(r1.begin(), r1.end(), r2.begin(), sum.at(row).begin(), std::minus<double>());
    }
    return Matrix(m1.get_rows(), m1.get_cols(), sum);
}

Vector operator*(double scalar, const Vector& v) {
    std::vector<double> product = v.get_vector();
    for (double& component : product)
        component *= scalar;
    return Vector(v.get_dim(), product);
}

Vector operator*(const Vector& v, double scalar) {
    return scalar * v;
}

Matrix operator*(double scalar, const Matrix& m) {
    std::vector< std::vector<double> > product = m.get_matrix();
    for (std::vector<double>& row : product)
        for (double& entry : row)
            entry *= scalar;
    return Matrix(m.get_rows(), m.get_cols(), product);
}

Matrix operator*(const Matrix& m, double scalar) {
    return scalar * m;
}

bool is_square(const Matrix& m) {
    return m.get_rows() == m.get_cols();
}