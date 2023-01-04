# IntrLinAlgLibrary

A single header library with functionality for all concepts covered in Intr Lin Alg 250 @ rutgers.

* ```IntrLinAlgLibrary.hpp```
* Author: Sumanta Das (2022)
* Licence: MIT
* Textbook used: **Elementary Linear Algebra: A Matrix Approach by Lawrence E. Spence, Arnold J. Insel, Stephen H. Friedberg**

# Setup and Intro

To use, simply include the file: ```#include "IntrLinAlgLibrary.hpp"```

This file uses ```namespace ila``` for "Intro Linear Algebra".

* **NOTE:** There may be naming conflicts with other libraries if you use ```using namespace ila;```. For example, ```IntrLinAlgLibrary``` contains a ```print()``` method which may appear in other files. Instead, consider limiting to certain data types or functions such as ```using namespace ila::Matrix, ila::Vector;``` and/or any other functions you would like.

Vectors are arrays and matrices are arrays of arrays. Both contain ```float``` as their data type.

* Vector elements are called **components**. A vector with dimension ```D``` has ```D``` components.
* Matrix elements are called **entries**. A ```R x C``` matrix has ```R``` rows and ```C``` columns, thus having ```R x C``` entries.

Intr Lin Alg 250 @ rutgers covers 6 chapters (some not in its entirety):

* Chapter 1 - Matrices, Vectors, and Systems of Linear Equations
* Chapter 2 - Matrices and Linear Transformations
* Chapter 3 - Determinants
* Chapter 4 - Subspaces and their Properties
* Chapter 5 - Eigenvalues, Eigenvectors, and Diagonalization
* Chapter 6 - Orthogonality

# Creating a vector and/or a matrix, and print overloads:

Creating a vector: ```ila::Vector<D> myVec(1, 2, 3, ...);```

* ```D``` - the dimension/size of the vector.
* Add ```D``` number of ```float``` parameters for the contents of the vector.

**Example:** Creating the vector called ```myVec``` with 4 components: ```[ 1, 2, 3, 4 ]```

```cpp
ila::Vector<4> myVec(1.0, 2.0, 3.0, 4.0);
```

Vectors are traditionally one-indexed in linear algebra, and **NOT** zero-indexed like traditional array.
Access the **first** component by ```myVec[1]```, and ***NOT*** ```myVec[0] //Error```.

* If you really prefer zero-indexing, type ```myVec.components[0]``` to access the first component. 

A print overload is provided to you to print your vector to the standard output. ```ila::print(myVec)``` outputs:

```
{    1    }
{    2    }
{    3    }
{    4    }
```

Alternatively, ```std::cout << myVec;``` can also be used.

### Special Vectors

If your vector has 2 components (of type ```ila::Vector<2>```), you can use:

* ```myVec.x()``` to access the first component and ```myVec.y()``` to access the second component (use as a 2D point).

If your vector has 3 components (of type ```ila::Vector<3>```), you can use:

* ```myVec.x()``` for the first component, ```myVec.y()``` for the second component, and ```myVec.z()``` for the third component (use as a 3D point).
* ```myVec.r()``` for the first component, ```myVec.g()``` for the second component, and ```myVec.b()``` for the third component (use as a color).

If your vector has 4 components (of type ```ila::Vector<4>```), you can use:

* ```myVec.r()``` for the first component, ```myVec.g()``` for the second component, ```myVec.b()``` for the third component, and ```myVec.a()``` for the fourth component (use as a color).

---

Creating a matrix: ```ila::Matrix<R, C> myMat(1, 2, 3, ...);```

* ```R``` - the number of rows of the matrix.
* ```C``` - the number of columns of the matrix.
* Add ```R``` times ```C``` number of ```float``` parameters for the contents of the matrix.

**Example:** Creating a 3x3 (3 rows, 3 columns) matrix called ```myMat``` with entries: ```[[ 1, 2, 3 ], [4, 5, 6], [7, 8, 9]]```

```cpp
/* Line breaks are optional... included for readability */
ila::Matrix<3, 3> myMat
(
  1.0, 2.0, 3.0,
  4.0, 5.0, 6.0,
  7.0, 8.0, 9.0
);
/* There must be 9 parameters (3 x 3 = 9) */
```

Similar to vectors, matrices are also traditionally one-indexed in linear algebra.
Access the (1,1)-entry (entry at the 1st row and 1st column of the matrix) by ```myMat[1][1]``` and ***NOT*** ```myMat[0][0] //Error```.

* If you really prefer zero-indexing, type ```myMat.entries[0][0]``` to access the (1,1)-entry.

There is also a print overload for matrices to the standard output. ```ila::print(myMat);``` outputs:

```
{    1    2    3    }
{    4    5    6    }
{    7    8    9    }
```

Alternatively, ```std::cout << myMat;``` can also be used.

## Sets of vectors:

* This library handles sets of vectors using ```std::arrays```, so use ```std::arrays``` for user-defined vector sets. Example of a user-defined set of vectors:

```cpp
/* A set of 2 vectors with dimension 3 */
std::array<ila::Vector<3>, 2> mySet = 
{
    Vector<3>(1.0, 2.0, 3.0),
    Vector<3>(4.0, 5.0, 6.0)
};
```

* However, sizes of bases for row, column, and null spaces cannot be calculated at compile time and so are returned via an ```std::vector```. Manually hardcode the ```std::vector``` contents into new ```std::arrays``` and recompile to be compatible with this library's functions *(Read the Chapter 4 subsection for more details)*.

* ```print()``` overloads for sets represented as ```std::array```s and ```std::vector```s are given.

```cpp
std::array<ila::Vector<3>, 4> mySetAsStdArray  = { /* ... */ };
std::vector<ila::Vector<3>>   mySetAsStdVector = { /* ... */ };

/* Both print statements work */
ila::print(mySetAsStdArray);
ila::print(mySetAsStdVector);
```

* Sets printed to the standard output are printed with their augmented matrix representation. For example:

```cpp
std::array<ila::Vector<3>, 2> mySet = 
{
    Vector<3>(1.0, 2.0, 3.0),
    Vector<3>(4.0, 5.0, 6.0)
};

ila::print(mySet);

/* Outputs:

{    1    4    }
{    2    5    }
{    3    6    }

*/
```

---

All ```ila::print``` overloads:

* Overload for ```ila::Vector<D>```
* Overload for ```ila::Matrix<R, C>```
* Overload for ```std::array<Vector<D>, S>``` and ```std::vector<Vector<D>>```
* Overload for ```std::vector<ila::Eigenvalue>``` *(Read the Chapter 5 subsection for more details)*

# Library function details:

For all function descriptions, function and template arguments can be replaced with whatever you like. I am just providing examples.

* ```myMat```, ```myMat1```, and ```myMat2``` are sample names for matrices. ```mySquareMat``` signifies that the matrix parameter must be a square matrix.
* ```myVec```, ```myVec1```, and ```myVec2``` are sample names for vectors.
* ```mySetAsStdArray``` is a sample name for a set of vectors as an ```std::array``` whose size is known at compile time.
* ```mySetAsStdVector``` is a sample name for a set of vectors as an ```std::vector``` whose size is *not* known at compile time.

### Operator overloads exist for:

* Adding and subtracting vectors and matrices: ```myVec1 + myVec2```
* Multiplying and dividing a vector/matrix by a scalar: ```4.3 * myMat``` and ```myVec / 3.3```
* Multiplying a matrix by a vector: A ```A x B``` matrix times a vector of dimension ```B``` results in a vector with dimension ```A```: ```myMat * myVec```. This is *not* commutative.
* Multiplying a matrix by matrix: A ```A x B``` matrix times a ```B x C``` matrix results in a ```A x C``` matrix: ```myMat1 * myMat2```. This is *not* commutative.

Matrices and vectors *must* be of valid dimensions when performing operations on them, or else the program will not compile.

### Special Vectors and Matrices:

* ```ila::zero_vector<6>()``` produces a zero vector of dimension ```6```.
* ```ila::zero_matrix<7, 5>()``` produces a ```7 x 5``` zero matrix.
* ```ila::standard_vector<4>(3)``` produces a standard vector of dimension ```4``` whose ```3```rd component contains the 1.0: ```[ 0, 0, 1, 0 ]```
* ```ila::identity_matrix<6>()``` produces a ```6 x 6``` identity matrix.
* ```ila::rotation_matrix(45.0)``` produces a ```2 x 2``` rotation matrix which rotates a vector of dimension 2 by ```45.0``` degrees counter-clockwise.

## Chapter 1 - Matrices, Vectors, and Systems of Linear Equations

* ```ila::transpose(myMat)``` returns the transpose of ```myMat```.
* ```ila::ERO_row_swap(myMat, 3, 4)``` swaps the ```3```rd and ```4```th row of ```myMat```.
* ```ila::ERO_scalar_multiplication(myMat, 3.44, 2)``` multiplies the ```2```nd row of ```myMat``` by ```3.44```.
* ```ila::ERO_row_sum(myMat, 2.21, 2, 3)``` adds the ```2```nd row multiplied by ```2.21``` to the ```3```rd row of ```myMat``` (the ```2```nd row is left unchanged).
* ```ila::ref(myMat)``` returns the row-echelon form of ```myMat```.
* ```ila::rref(myMat)``` returns the reduced row-echelon form of ```myMat```.
* ```ila::rank(myMat)``` returns the rank of ```myMat```.
* ```ila::nullity(myMat)``` returns the nullity of ```myMat```.
* ```ila::augment(myMat, myVec)``` returns a matrix that augments ```myVec``` to the end of ```myMat```.
* ```ila::augment``` is overloaded to augment two matrices together as well.
* ```ila::augment_vector_set(mySetAsStdArray)``` merges all the vectors in ```mySetAsStdArray``` into a matrix and returns it.
* ```ila::solve(myMat, myVec)``` returns the reduced row-echelon form of ```myVec``` augmented to ```myMat```. rref is designed so that its solution to the system ```myMat * x = myVec``` is easy to extract:

```
                      x1 x2 x3 x4  x5    b
Given an rref row: [  0  1  0  2  -3  |  3  ]
The leading entry of the rref can be solved for: x2 + 2x4 - 3x5 = 3 --> x2 = 3 - 2x4 + 3x5
A solution variable whose corresponding column has no pivot spot is a free variable.
```

* ```ila::is_consistent(myMat, myVec)``` returns whether or not a solution vector ```x``` exists for the system ```myMat * x = myVec```.
* ```ila::is_in_span(myVec, mySetAsStdArray)``` returns whether or not ```myVec``` is in the span of ```mySetAsStdArray```.
* ```ila::is_linearly_independent(mySetAsStdArray)``` returns whether or not ```mySetAsStdArray``` is linearly independent.
* ```ila::is_linearly_dependent(mySetAsStdArray)``` is also provided.
* ```ila::solve_homogenous_system(myMat)``` solves the *homogenous system* ```myMat * x = zeroVec```. Similar to ```ila::solve()```, it returns the reduced row-echelon form of the zero vector augmented to ```myMat```.

## Chapter 2 - Matrices and Linear Transformations

* ```ila::is_diagonal(myMat)``` returns whether or not ```myMat``` is diagonal.
* ```ila::is_symmetric(myMat)``` returns whether or not ```myMat``` is symmetric.
* ```ila::EM_row_swap<4>(2, 3)``` returns a ```4 x 4``` elementary matrix which swaps the ```2```nd and ```3```rd row of any matrix multiplied by this elementary matrix.
* ```ila::EM_scalar_multiplication<5>(3.23, 2)``` returns a ```5 x 5``` elementary matrix which multiplies the ```2```nd row of any matrix multiplied by this elementary matrix by ```3.23```.
* ```ila::EM_row_sum<3>(2.21, 1, 2)``` returns a ```3 x 3``` elementary matrix which adds the ```1```st row multiplied by ```2.21``` to the ```2```nd row of any matrix multiplied by this elementary matrix (the ```1```st row is left unchanged).
* ```ila::is_invertible(myMat)``` returns whether or not ```myMat``` is invertible.
* ```ila::inverse(mySquareMat)``` returns the inverse matrix of ```mySquareMat```.

## Chapter 3 - Determinants

* ```ila::det(mySquareMat)``` returns the determinant of ```mySquareMat```.

## Chapter 4 - Subspaces and their Properties

* ```ila::row(myMat)``` returns the basis of ```myMat```'s row space as an ```std::vector``` of vectors.
* ```ila::col(myMat)``` returns the basis of ```myMat```'s column space as an ```std::vector``` of vectors.
* ```ila::null(myMat)``` returns the basis of ```myMat```'s null space as an ```std::vector``` of vectors.

**IMPORTANT:** ```std::vector```s are not supported as sets in this library. Print the ```std::vector```'s contents and rewrite them as ```std::array```s to use them for this library.

### Getting the dimension of bases:

* ```ila::dim(ila::row(myMat))``` to get the dimension of ```myMat```'s row space.
* ```ila::dim(ila::col(myMat))``` to get the dimension of ```myMat```'s column space.
* ```ila::dim(ila::null(myMat))``` to get the dimension of ```myMat```'s null space.

**NOTE:** ```ila::dim``` simply returns the size of the ```std::vector```. This function *will not* produce the correct dimension for a user-defined set of vectors as an ```std::vector```.

* ```ila::basis(mySetAsStdArray)``` returns the basis of ```mySetAsStdArray``` as an ```std::vector```.

## Chapter 5 - Eigenvalues, Eigenvectors, and Diagonalization

* ```ila::is_lower_triangular(myMat)``` returns whether or not ```myMat``` is lower-triangular.
* ```ila::is_upper_triangular(myMat)``` returns whether or not ```myMat``` is upper-triangular.
* ```ila::generate_eigenvalues(mySquareMat)``` returns all of ```mySquareMat```'s eigenvalues as an ```std::vector``` of ```ila::Eigenvalue``` structs.

The ```ila::Eigenvalue``` struct contains two attributes: ```eigenvalue``` and ```multiplicity```.

An ```ila::print``` is overloaded for ```std::vector<ila::Eigenvalue>```. To print all the eigenvalues of ```mySquareMat```, write ```ila::print(ila::generate_eigenvalues(mySquareMat))```. Example:

```cpp
ila::Matrix<3, 3> mySquareMat
(
    5, -10, -5,
    2, 14, 2,
    -4, -8, 6
);
    
ila::print(ila::generate_eigenvalues(mySquareMat));
```

This outputs:

```
Eigenvalue: 10	Multiplicity: 2
Eigenvalue: 5   Multiplicity: 1
```

**IMPORTANT:** This library can only handle real eigenvalues. Complex eigenvalues are not supported.

* ```ila::get_eigenvalue(mySquareMat, eigenvector)``` returns ```mySquareMat```'s corresponding eigenvalue of its ```eigenvector```.
* ```ila::get_eigenspace_basis(mySquareMat, eigenvalue)``` returns ```mySquareMat```'s corresponsing eigenspace basis of its ```eigenvalue``` as an ```std::vector```.

Like all other bases, write ```ila::dim(ila::get_eigenspace_basis(mySquareMat, eigenvalue))``` to get the dimension of the eigenspace basis.

* ```ila::is_diagonalizable(myMat)``` returns whether or not ```myMat``` is diagonalizable.
* ```ila::diagonalize(mySquareMat)``` returns an invertible matrix ```P``` and a diagonal matrix ```D``` (such that ```mySquareMat = PDP^-1```) as an ```std::pair```.
  * To get matrix ```p```, write: ```ila::Matrix<S, S> p = ila::diagonalize(mySquareMat).first;``` *(replace ```S``` with the size of ```mySquareMat```)*.
  * To get matrix ```d```, write: ```ila::Matrix<S, S> d = ila::diagonalize(mySquareMat).second;``` *(replace ```S``` with the size of ```mySquareMat```)*.
  * To get both ```p``` and ```d``` at once, write: ```auto [p, d] = ila::diagonalize(mySquareMat)```.

## Chapter 6 - Orthogonality

* ```ila::dot(myVec1, myVec2)``` returns the dot product of ```myVec1``` and ```myVec2```.
  * ```myVec1 * myVec2``` is operator overloaded to return their dot product.
* ```ila::unit_vector(myVec)``` returns the unit vector of ```myVec```.
* ```ila::norm(myVec)``` returns the norm or length of ```myVec```.
* ```ila::length``` is the same function as ```ila::norm``` and can be used interchangeably.
* ```ila::normalize(myVec)``` normalizes ```myVec``` into a norm/length of 1.
* ```ila::distance(myVec1, myVec2)``` returns the distance between ```myVec1``` and ```myVec2```.
* ```ila::is_orthogonal(myVec1, myVec2)``` returns whether or not ```myVec1``` and ```myVec2``` are orthogonal to each other.
  * This function is overloaded for ```ila::is_orthogonal(mySetAsStdArray)``` which returns whether or not ```mySetAsStdArray``` is orthogonal.
* ```ila::orthogonal_projection(myVec1, myVec2)``` returns the orthogonal projection of ```myVec1``` onto ```myVec2```.
* ```ila::orthonormal_basis(myBasisAsStdArray)``` returns the orthonormal basis that generates the same space as ```myBasisAsStdArray```.
  * This function is overloaded for ```std::vector``` arguments: ```ila::orthonormal_basis(myBasisAsStdVector)```
* ```ila::qr_factorization(myMat)``` returns an orthogonal/semiorthogonal matrix ```Q``` and an upper-triangular matrix ```R``` (such that ```myMat = QR```) as an ```std::pair```.
  * To get matrix ```q```, write: ```ila::Matrix<R, C> q = ila::qr_factorization(myMat).first;``` *(replace ```R``` and ```C``` with the number of rows and columns of ```myMat```)*
  * To get matrix ```r```, write: ```ila::Matrix<C, C> r = ila::qr_factorization(myMat).second;``` *(replace ```C``` with the number of columns of ```myMat```)*
  * To get both ```q``` and ```r``` at once, write: ```auto [q, r] = ila::qr_factorization(myMat);```

## Additional Functions

* ```ila::cross(myVec1, myVec2)``` returns the cross product of ```myVec1``` and ```myVec2```. **NOTE:** This function only works for vectors in R2 and R3:
  * For two vectors in R2, this function returns the scalar of the area of its enclosed parallelogram as a ```float```.
  * For two vectors in R3, this function returns the vector (as an ```ila::Vector<3>``` of appropriate magnitude that is orthogonal to both the argument vectors.