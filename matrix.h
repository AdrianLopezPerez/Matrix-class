//  ********************************************
// *                DISCLAIMER                  *
// *                                            *
// * BASED IN THE WORK OF MICHAEL L HALLS-MOORE *
// *      AS FOUND IN QUANTSTART.COM AND        *
// *  HIS BOOK "C++ FOR QUANTITATIVE FINANCE"   *
// *  SLIGHTLY MODIFIED BY ME (ADRIAN LOPEZ)    *
// *          FOR MY OWN CONVENIENCE            *
//  ********************************************
// 
// MY OWN CONTRIBUTIONS
// Changes:
// - Rows and columns capped to 256
// - Changed the name of some original methods and vars
// New functions:
// - Determinant
// - Inversion
// - Trace
// - Print
// - Access rows and cols

#ifndef __MATRIX_H
#define __MATRIX_H 

#include <iostream>
#include <vector>

template <typename T> 
class Matrix {
private:
    std::vector<std::vector<T>> mat;
    unsigned short nrows;
    unsigned short ncols;

public:
    Matrix(int _rows, int _cols);
    Matrix(const Matrix<T>& rhs);
    virtual ~Matrix();

    // Operator overloading, for "standard" mathematical matrix operations                                                                                                                                                          
    Matrix<T>& operator=(const Matrix<T>& rhs);

    // MATRIX-MATRIX OPERATIONS                                                                                                                                                                                               
    Matrix<T> operator+(const Matrix<T>& rhs);          // ADDITION
    Matrix<T>& operator+=(const Matrix<T>& rhs);        // CUMULATIVE ADDITION
    Matrix<T> operator-(const Matrix<T>& rhs);          // SUBSTRACTION
    Matrix<T>& operator-=(const Matrix<T>& rhs);        // CUMULATIVE SUBSTRACTION
    Matrix<T> operator*(const Matrix<T>& rhs);          // MATRIX PRODUCT
    Matrix<T>& operator*=(const Matrix<T>& rhs);        // CUMULATIVE MATRIX PRODUCT
    Matrix<T> t();                                      // TRANSPOSITION
    Matrix<T> inv();                                    // INVERSION
    void identity();                                    // CONSTRUCTION OF N-DIMENSIONAL IDENTITY

    // MATRIX-SCALAR OPERATIONS                                                                                                                                                                                                    
    Matrix<T> operator+(const T& rhs);                  // ADDITION OF SCALAR TO EACH ELEMENT 
    Matrix<T> operator-(const T& rhs);                  // SUBSTRACTION OF SCALAR TO EACH ELEMENT
    Matrix<T> operator*(const T& rhs);                  // PRODUCT OF SCALAR TO EACH ELEMENT
    Matrix<T> operator/(const T& rhs);                  // DIVISION BY SCALAR TO EACH ELEMENT
    T det();                                            // DETERMINANT
    T trace();                                          // TRACE

    // Matrix/vector operations                                                                                                                                                                                                     
    std::vector<T> diagonal();                          // ELEMENTS IN DIAGONAL

    // Access the elements
    T& operator()(const int& row, const int& col);
    const T& operator()(const int& row, const int& col) const;

    // Access individual rows and columns
    std::vector<T> row(const int& row) const;
    std::vector<T> col(const int& col) const;

    // Access the row and column sizes                                                                                                                                                                                              
    unsigned short nrow() const;
    unsigned short ncol() const;

    // Print
    void print();
};

// Parameter Constructor                                                                                                                                                      
template<typename T>
Matrix<T>::Matrix(int _rows, int _cols) {
     if (_rows > 0 && _cols > 0 && _rows <= 256 && _cols <= 256) {
        mat.resize(_rows);
        for (unsigned short i = 0; i < mat.size(); ++i) {
            mat[i].resize(_cols, 0);
        }
        nrows = _rows;
        ncols = _cols;
    }
    else throw "A matrix with non-positive or more than 256 rows or columns is not allowed";
}

// Copy Constructor                                                                                                                                                           
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs) {
    mat = rhs.mat;
    nrows = rhs.nrow();
    ncols = rhs.ncol();
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Matrix<T>::~Matrix() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
    if (&rhs == this) return *this;
    unsigned short new_rows = rhs.nrow();
    unsigned short new_cols = rhs.ncol();
    mat.resize(new_rows);
    for (unsigned short i = 0; i < mat.size(); ++i) {
        mat[i].resize(new_cols);
    }
    for (unsigned short i = 0; i < new_rows; ++i) {
        mat[i] = rhs.row(i);
    }
    nrows = new_rows;
    ncols = new_cols;
    return *this;
}


// Addition of two matrices                                                                                                                                                   
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& rhs) {
    Matrix result(nrows, ncols);
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            result(i, j) = this->mat[i][j] + rhs(i, j);
        }
    }
    return result;
}

// Cumulative addition of this matrix and another                                                                                                                             
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.nrow(), cols = rhs.ncol();
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < cols; ++j) {
            this->mat[i][j] += rhs(i, j);
        }
    }
    return *this;
}

// Subtraction of this matrix and another                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.nrow(), cols = rhs.ncol();
    Matrix result(nrows, ncols);
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            result(i, j) = this->mat[i][j] - rhs(i, j);
        }
    }
    return result;
}

// Cumulative subtraction of this matrix and another                                                                                                                          
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.nrow(), cols = rhs.ncol();
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            this->mat[i][j] -= rhs(i, j);
        }
    }
    return *this;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) {
    if (ncols == rhs.nrow()) {
        Matrix result(nrows, rhs.ncol());
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                for (unsigned short k = 0; k < nrows; k++) {
                    result(i, j) += this->mat[i][k] * rhs(k, j);
                }
            }
        }
        return result;
    }
    else throw "Cannot be multiplied.";
}

// Cumulative left multiplication of this matrix and another                                                                                                                  
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& rhs) {
    Matrix result = (*this) * rhs;
    (*this) = result;
    return *this;
}

// TRANSPOSE                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::t() {
    Matrix result(nrows, ncols);
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            result(i, j) = this->mat[j][i];
        }
    }
    return result;
}

// INVERSION
template<typename T>
Matrix<T> Matrix<T>::inv() {
    if (nrows != ncols) throw "Non-invertible";
    else {
        Matrix<T> result(nrows, nrows);
        result.identity();
        Matrix<T> up_triang(nrows, nrows);
        up_triang.mat = mat;
        std::vector<T> v;
        T coeff;
        for (int i = 0; i < nrows - 1; ++i) {
            if (up_triang.mat[i][i] != 0) {
                for (unsigned short j = 1 + i; j < nrows; ++j) {
                    coeff = up_triang.mat[j][i] / up_triang.mat[i][i];
                    for (unsigned short k = 0; k < ncols; ++k) {
                        up_triang.mat[j][k] = up_triang.mat[j][k]
                            - coeff * up_triang.mat[i][k];
                        result.mat[j][k] = result.mat[j][k]
                            - coeff * result.mat[i][k];
                    }
                }
            }
            else if (i < nrows - 1) {
                std::vector<T> v;
                for (unsigned short j = i; j < nrows - 1; ++j) {
                    if (up_triang.mat[j + 1][i] != 0) {
                        v = up_triang.mat[i];
                        up_triang.mat[i] = up_triang.mat[j + 1];
                        up_triang.mat[j + 1] = v;
                        v = result.mat[i];
                        result.mat[i] = result.mat[j + 1];
                        result.mat[j + 1] = v;
                        break;
                    }
                    else throw "Non-invertible";
                }
                --i;
            }
            else if (i = nrows - 1) throw "Non-invertible";
        }
        // Get all pivots equal to one on lhs
        for (unsigned short i = 0; i < nrows; ++i) {
            coeff = up_triang.mat[i][i];
            for (unsigned short j = 0; j < nrows; ++j) {
                up_triang.mat[i][j] = up_triang.mat[i][j] / coeff;
                result.mat[i][j] = result.mat[i][j] / coeff;
            }
        }
        // Convert lhs to upper triangular
        for (unsigned short i = nrows - 1; i > 0 && i < nrows; --i) {
            for (unsigned short j = i - 1; j >= 0 && j < nrows - 1; --j) {
                coeff = up_triang.mat[j][i];
                for (unsigned short k = nrows - 1; k >= 0 && k < nrows; --k) {
                    up_triang.mat[j][k] = up_triang.mat[j][k]
                        - coeff * up_triang.mat[i][k];
                    result.mat[j][k] = result.mat[j][k]
                        - coeff * result.mat[i][k];
                }
            }
        }

        return result;
    }
}

// CONSTRUCTION OF N-DIMENSIONAL IDENTITY
template<typename T>
void Matrix<T>::identity() {
    if (nrows == ncols) {
        for (unsigned short i = 0; i < nrows; ++i) {
            mat[i][i] = 1;
        }
    }
    else throw "Number of rows and columns must coincide.";
}

// Matrix/scalar addition                                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::operator+(const T& rhs) {
    Matrix result(nrows, ncols);
    if (rhs == 0) {
        result.mat = mat;
        return result;
    }
    else {
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                result(i, j) = this->mat[i][j] + rhs;
            }
        }
        return result;
    }
}

// Matrix/scalar subtraction                                                                                                                                                  
template<typename T>
Matrix<T> Matrix<T>::operator-(const T& rhs) {
    Matrix result(nrows, ncols);
    if (rhs == 0) {
        result.mat = mat;
        return result;
    }
    else {
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                result(i, j) = this->mat[i][j] - rhs;
            }
        }
        return result;
    }
}

// Matrix/scalar multiplication                                                                                                                                               
template<typename T>
Matrix<T> Matrix<T>::operator*(const T& rhs) {
    Matrix result(nrows, ncols);
    if (rhs == 0) return result;
    else if (rhs == 1) {
        result.mat = mat;
        return result;
    }
    else {
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                result(i, j) = this->mat[i][j] * rhs;
            }
        }
        return result;
    }
}

// Matrix/scalar division                                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::operator/(const T& rhs) {
    if (rhs == 0) throw "Cannot divide by zero";
    else {
        Matrix result(nrows, ncols);
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                result(i, j) = this->mat[i][j] / rhs;
            }
        }
        return result;
    }
}

// DETERMINANT
template<typename T>
T Matrix<T>::det() {
    if (nrows != ncols) throw "Matrix is non-square"; 
    else {                                            
        switch (nrows) {
        // Simplified method for matrices of dimension 1, 2, and 3
        case 1:
            return mat[0][0];
        case 2:
            return mat[0][0] * mat[1][1]
                 - mat[0][1] * mat[1][0];
        case 3:
            return mat[0][0] * mat[1][1] * mat[2][2]
                 + mat[0][1] * mat[1][2] * mat[2][0]
                 + mat[0][2] * mat[1][0] * mat[2][1]
                 - mat[0][2] * mat[1][1] * mat[2][0]
                 - mat[0][1] * mat[1][0] * mat[2][2]
                 - mat[0][0] * mat[1][2] * mat[2][1];
        // Complex method for n-dimensional matrices
        // Convert matrix to upper triangular, then multiply elements
        // in the diagonal to obtain the determinant.
        default:
            // Declare working matrix up_triang, initialize as the original matrix
            Matrix<T> up_triang(nrows,ncols);
            up_triang.mat = mat;
            // Variable to change the sign of the determinant when rows are swapped.
            short sign_by_row_swap = 1; 
            T coeff;               
            for (int i = 0; i < nrows - 1; ++i) { 
                if (up_triang.mat[i][i] != 0) {
                    for (unsigned short j = 1 + i; j < nrows; ++j) {
                        coeff = up_triang.mat[j][i] / up_triang.mat[i][i];
                        for (unsigned short k = i; k < ncols; ++k) {
                            up_triang.mat[j][k] = up_triang.mat[j][k]
                                - coeff * up_triang.mat[i][k];
                        }
                    }
                }
                else if (i < nrows - 1) {
                    std::vector<T> v;
                    for (unsigned short j = i; j < nrows - 1; ++j) {
                        if (up_triang.mat[j + 1][i] != 0) {
                            v = up_triang.mat[i];
                            up_triang.mat[i] = up_triang.mat[j + 1];
                            up_triang.mat[j + 1] = v;
                            sign_by_row_swap *= -1;
                            break;
                        }
                        else return 0;
                    }
                    --i;
                }
                else if (i = nrows - 1) return 0;
            }
            std::vector<T> v = up_triang.diagonal();
            T result = 1;
            for (unsigned short i = 0; i < v.size(); ++i) {
                result *= v[i];
            }
            result = sign_by_row_swap * result;
            return result;
        }
    }
}

// TRACE
template<typename T>
T Matrix<T>::trace() {
    if (nrows != ncols) throw "Matrix is non-square";
    else {
        T result = 0;
        std::vector<T> v = diagonal();
        for (unsigned short i = 0; i < v.size(); ++i) {
            result += v[i];
        }
        return result;
    }
}

// Obtain a vector of the diagonal elements                                                                                                                                   
template<typename T>
std::vector<T> Matrix<T>::diagonal() {
    std::vector<T> result(nrows, 0);
    for (unsigned short i = 0; i < nrows; ++i) {
        result[i] = this->mat[i][i];
    }
    return result;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& Matrix<T>::operator()(const int& row, const int& col) {
    if (row >= 0 && col >= 0 && row < nrows && col < ncols) {
        return this->mat[row][col];
    }
    else throw "Out of limits";
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& Matrix<T>::operator()(const int& row, const int& col) const {
    if (row >= 0 && col >= 0 && row < nrows && col < ncols) {
        return this->mat[row][col];
    }
    else throw "Out of limits";
}

// Access row
template<typename T>
std::vector<T> Matrix<T>::row(const int& row) const{
    if (row >= 0 && row < nrows) {
        return this->mat[row];
    }
    else throw "Out of limits";
}

template<typename T>
std::vector<T> Matrix<T>::col(const int& col) const{
    if (col >= 0 && col < ncols) {
        std::vector<T> v;
        for (unsigned short i = 0; i < nrows; ++i) {
            v.push_back(mat[i][col]);
        }
        return v;
    }
    else throw "Out of limits";
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
inline unsigned short Matrix<T>::nrow() const {
    return this->nrows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
inline unsigned short Matrix<T>::ncol() const {
    return this->ncols;
}

// PRINT MATRIX
template<typename T>
void Matrix<T>::print() {
    std::cout << std::endl;
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            if (mat[i][j] >= 0) std::cout << ' ';
            std::cout << mat[i][j] << '\t';
            if (j == ncols - 1) std::cout << '\n';
        }
    }
}

#endif
