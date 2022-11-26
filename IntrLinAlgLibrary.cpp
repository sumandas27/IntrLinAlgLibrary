#include "IntrLinAlgLibrary.hpp"

void IntrLinAlgLibrary_init() {
    std::cout << std::setprecision(3);
}

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

double Vector::get(unsigned int index) const {
    assert (index >= 1 && index <= dim);
    return components.at(index - 1);
}

void Vector::set(unsigned int index, double value) {
    assert (index >= 1 && index <= dim);
    components.at(index - 1) = value;
}

void Vector::print() {
    for (const double& component : components)
        std::cout << "{\t" << (component < epsilon() ? 0 : component) << "\t}\n";
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
        std::cout << "{\t";
        for (int j = 0; j < cols; j++)
            std::cout << (abs(entries.at(i * cols + j)) < epsilon() ? 0 : entries.at(i * cols + j)) << "\t";
        std::cout << "}\n";
    }
    std::cout << "\n";
}

//----------------------------------------------------------------------//
//NON-LINEAR ALGEBRA FUNCTIONS:

constexpr double epsilon() {
    return 0.00001;
}

bool is_equal(double val1, double val2) {
    return abs(val1 - val2) < epsilon();
}

constexpr double deg_to_rad(double degrees) {
    return degrees * M_PI / 180;
}

//----------------------------------------------------------------------//
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
}