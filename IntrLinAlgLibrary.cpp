#include "IntrLinAlgLibrary.hpp"

//----------------------------------------------------------------------//
//VECTOR METHODS:

Vector::Vector(unsigned int _dim, std::vector<double>& _components) : dim(_dim) {
    assert (_dim != 0);
    assert (_components.size() == _dim);
    components.insert(components.end(), _components.begin(), _components.end());
}

unsigned int Vector::get_dim() const {
    return dim;
}

std::vector<double> Vector::get_components() const {
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
        std::cout << "{  " << std::to_string(component).substr(0, 8) << "  }\n";
    std::cout << "\n";
}

//----------------------------------------------------------------------//
//MATRIX METHODS:

Matrix::Matrix(unsigned int _rows, unsigned int _cols, std::vector<double>& _entries) : rows(_rows), cols(_cols) {
    assert (_rows != 0);
    assert (_cols != 0);
    assert (_entries.size() == _rows * _cols);
    entries.insert(entries.end(), _entries.begin(), _entries.end());
}

unsigned int Matrix::get_rows() const {
    return rows;
}

unsigned int Matrix::get_cols() const {
    return cols;
}

std::vector<double> Matrix::get_entries() const {
    return entries;
}

double Matrix::get(unsigned int row, unsigned int col) const {
    assert (row >= 1 && row <= rows);
    assert (col >= 1 && col <= cols);
    return entries.at((row - 1) * cols + (col - 1));
}

void Matrix::set(unsigned int row, unsigned int col, double value) {
    assert (row >= 1 && row <= rows);
    assert (col >= 1 && col <= cols);
    entries.at((row - 1) * cols + (col - 1)) = value;
}

void Matrix::print() {
    for (int i = 0; i < rows; i++) {
        std::cout << "{  ";
        for (int j = 0; j < cols; j++)
            std::cout << std::to_string(entries.at(i * cols + j)).substr(0, 8) << "  ";
        std::cout << "}\n";
    }
    std::cout << "\n";
}

//----------------------------------------------------------------------//
//NON-LINEAR ALGEBRA FUNCTIONS:

constexpr double deg_to_rad(double degrees) {
    return degrees * M_PI / 180;
}

//----------------------------------------------------------------------//
//CHAPTER 1 - MATRICES, VECTORS, AND SYSTEMS OF LINEAR EQUATIONS

bool operator==(const Vector& v1, const Vector& v2) {
    return v1.get_components() == v2.get_components();
}

bool operator==(const Matrix& m1, const Matrix& m2) {
    if (m1.get_rows() != m2.get_rows() || m1.get_cols() != m2.get_cols())
        return false;
    return m1.get_entries() == m2.get_entries();
}

Vector operator+(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());
    std::vector<double> sum(v1.get_components().size());
    std::transform(v1.get_components().begin(), v1.get_components().end(), v2.get_components().begin(), sum.begin(), std::plus<double>());
    return Vector(v1.get_dim(), sum);
}

Matrix operator+(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());
    std::vector<double> sum(m1.get_entries().size());
    std::transform(m1.get_entries().begin(), m1.get_entries().end(), m2.get_entries().begin(), sum.begin(), std::plus<double>());
    return Matrix(m1.get_rows(), m1.get_cols(), sum);
}

Vector operator-(const Vector& v1, const Vector& v2) {
    assert (v1.get_dim() == v2.get_dim());
    std::vector<double> diff(v1.get_components().size());
    std::transform(v1.get_components().begin(), v1.get_components().end(), v2.get_components().begin(), diff.begin(), std::minus<double>());
    return Vector(v1.get_dim(), diff);
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
    assert (m1.get_rows() == m2.get_rows());
    assert (m1.get_cols() == m2.get_cols());
    std::vector<double> diff(m1.get_entries().size());
    std::transform(m1.get_entries().begin(), m1.get_entries().end(), m2.get_entries().begin(), diff.begin(), std::plus<double>());
    return Matrix(m1.get_rows(), m1.get_cols(), diff);
}

Vector operator*(double scalar, const Vector& v) {
    std::vector<double> product(v.get_components().size());
    std::transform(v.get_components().begin(), v.get_components().end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Vector(v.get_dim(), product);
}

Vector operator*(const Vector& v, double scalar) {
    return scalar * v;
}

Matrix operator*(double scalar, const Matrix& m) {
    std::vector<double> product(m.get_entries().size());
    std::transform(m.get_entries().begin(), m.get_entries().end(), product.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
    return Matrix(m.get_rows(), m.get_cols(), product);
}

Matrix operator*(const Matrix& m, double scalar) {
    return scalar * m;
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
    for (int i = 0; i < size; i++)
        identityMatrix.at(i * size + i) = 1.0;
    return Matrix(size, size, identityMatrix);
}

Matrix transpose(const Matrix& m) {
    unsigned int rows = m.get_cols();
    unsigned int cols = m.get_rows();
    std::vector<double> transpose(rows * cols);
    for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
        transpose.at(i * cols + j) = m.get_entries().at(j * rows + i);
    return Matrix(rows, cols, transpose);
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

bool is_square(const Matrix& m) {
    return m.get_rows() == m.get_cols();
}