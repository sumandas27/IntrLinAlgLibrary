# IntrLinAlgLibrary
A single header library with functionality for all concepts covered in Intr Lin Alg 250 @ rutgers.

* IntrLinAlgLibrary.hpp
* Author: Sumanta Das (2022)

## Setup
Simply include the file: ```#include "IntrLinAlgLibrary.hpp```

This file uses ```namespace ila``` for "Intro Linear Algebra". To disable this prefix, type: ```using namespace ila;```
NOTE: There may be naming conflicts with other libraries. For example, ```IntrLinAlgLibrary``` contains a ```print()``` method which may appear in other files.

## Creating a vector and/or a matrix

Vectors are arrays while matrices are an array of arrays. They both contain ```double``` as their data type.

Creating a vector: ```ila::Vector<D> myVec(1, 2, 3, ...);```

* ```D``` - the dimension/size of the vector.
* Add ```D``` number of ```double``` parameters for the contents of the vector.

Example: Creating the vector with 4 components: ```[ 1, 2, 3, 4 ]```

```cpp
ila::Vector<4> myVec(1, 2, 3, 4);
```

Vectors are traditionally one-indexed in linear algebra, and **NOT** zero-indexed like traditional . Access the **first** element (component) by ```myVec[1]```, and ***NOT*** ```myVec[0] // Error```.

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

Example: Creating a 3x3 (3 rows, 3 columns) matrix: ```[[ 1, 2, 3 ], [4, 5, 6], [7, 8, 9]]```

```cpp
/* Line breaks are optional... included for readability */
ila::Matrix<3, 3> myMat
(
  1, 2, 3,
  4, 5, 6,
  7, 8, 9
);
```

There is also a print overload for matrices to the standard output. ```ila::print(myMat);``` outputs:

```
{    1    2    3    }
{    4    5    6    }
{    7    8    9    }
```

Alternatively, ```std::cout << myMat;``` can also be used.
