# IntrLinAlgLibrary
A single header library with functionality for all concepts covered in Intr Lin Alg 250 @ rutgers.

* ```IntrLinAlgLibrary.hpp```
* Author: Sumanta Das (2022)

## Setup
Simply include the file: ```#include "IntrLinAlgLibrary.hpp```

This file uses ```namespace ila``` for "Intro Linear Algebra". To disable this prefix, type: ```using namespace ila;```
NOTE: There may be naming conflicts with other libraries. For example, ```IntrLinAlgLibrary``` contains a ```print()``` method which may appear in other files.

Vectors are arrays while matrices are an array of arrays. They both contain ```double``` as their data type.

* Vector elements are called **components**. A vector with dimension ```D``` has ```D``` components.
* Matrix elements are called **entries**. A ```R x C``` matrix has ```R``` rows and ```C``` columns, thus having ```R x C``` entries.

Intr Lin Alg 250 @ rutgers covers 6 chapters (some not in its entirety):

* Chapter 1 - Matrices, Vectors, and Systems of Linear Equations
* Chapter 2 - Matrices and Linear Transformations
* Chapter 3 - Determinants
* Chapter 4 - Subspaces and their Properties
* Chapter 5 - Eigenvalues, Eigenvectors, and Diagonalization
* Chapter 6 - Orthogonality

## Creating a vector and/or a matrix

Creating a vector: ```ila::Vector<D> myVec(1, 2, 3, ...);```

* ```D``` - the dimension/size of the vector.
* Add ```D``` number of ```double``` parameters for the contents of the vector.

Example: Creating the vector called ```myVec``` with 4 components: ```[ 1, 2, 3, 4 ]```

```cpp
ila::Vector<4> myVec(1.0, 2.0, 3.0, 4.0);
```

Vectors are traditionally one-indexed in linear algebra, and **NOT** zero-indexed like traditional . 
Access the **first** component by ```myVec[1]```, and ***NOT*** ```myVec[0] //Error```.

A print overload is provided to you to print your vector to the standard output. ```ila::print(myVec)``` outputs:

```
{    1    }
{    2    }
{    3    }
{    4    }
```

Alternatively, ```std::cout << myVec;``` can also be used.

---

Creating a matrix: ```ila::Matrix<R, C> myMat(1, 2, 3, ...);```

* ```R``` - the number of rows of the matrix.
* ```C``` - the number of columns of the matrix.
* Add ```R``` times ```C``` number of ```double``` parameters for the contents of the matrix.

Example: Creating a 3x3 (3 rows, 3 columns) matrix called ```myMat``` with entries: ```[[ 1, 2, 3 ], [4, 5, 6], [7, 8, 9]]```

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

There is also a print overload for matrices to the standard output. ```ila::print(myMat);``` outputs:

```
{    1    2    3    }
{    4    5    6    }
{    7    8    9    }
```

Alternatively, ```std::cout << myMat;``` can also be used.

## Library function details:

For all function descriptions, function and template arguments can be replaced with whatever you like. I am just providing examples.

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
* ```ila::ERO_row_swap(myMat, 3, 4)``` performs the *row swap* elementary row operation on ```myMat```. It swaps the ```3```rd and ```4```th row of ```myMat```.
* ```ila::ERO_scalar_multiplication(myMat, 3.44, 2)``` performs the *scalar multiplication* elementary row operation on ```myMat```. It multiplies the ```2```nd row of ```myMat``` by ```3.44```.
