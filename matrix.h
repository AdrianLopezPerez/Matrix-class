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
    Matrix(signed int _rows, signed int _cols);
    Matrix(const Matrix<T>& rhs);
    virtual ~Matrix();

    // Operator overloading, for "standard" mathematical matrix operations                                                                                                                                                          
    Matrix<T>& operator=(const Matrix<T>& rhs);

    // Matrix-matrix operations                                                                                                                                                                                             
    Matrix<T> operator+(const Matrix<T>& rhs);   // Addition
    Matrix<T>& operator+=(const Matrix<T>& rhs); // Cumulative addition
    Matrix<T> operator-(const Matrix<T>& rhs);   // Substraction
    Matrix<T>& operator-=(const Matrix<T>& rhs); // Cumulative substraction
    Matrix<T> operator*(const Matrix<T>& rhs);   // Matrix product
    Matrix<T>& operator*=(const Matrix<T>& rhs); // Cumulative product
    Matrix<T> t();                               // Transposition
    Matrix<T> inv();                             // Inversion
    void identity();                             // Construction of n-dimensional matrix

    // Matrix-scalar operations                                                                                                                                                                                                    
    Matrix<T> operator+(const T& rhs);           // Addition of scalar to each element
    Matrix<T> operator-(const T& rhs);           // Substraction of scalar to each element
    Matrix<T> operator*(const T& rhs);           // Product by scalar
    Matrix<T> operator/(const T& rhs);           // Division by scalar
    T det();                                     // Determinant
    T trace();                                   // Trace

    // Matrix/vector operations                                                                                                                                                                                                     
    std::vector<T> diagonal();                   // Elements in diagonal

    // Access the elements by rows
    std::vector<T>& operator[](const int& row);
    const std::vector<T>& operator[](const int& row) const;

    // Access the row and column sizes                                                                                                                                                                                              
    unsigned short rows() const;
    unsigned short cols() const;

    // Print the matrix in console
    void print();
};

// Parameter Constructor                                                                                                                                                      
template<typename T>
Matrix<T>::Matrix(signed int _rows, signed int _cols) {
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
    nrows = rhs.rows();
    ncols = rhs.cols();
}

// (Virtual) Destructor                                                                                                                                                       
template<typename T>
Matrix<T>::~Matrix() {}

// Assignment Operator                                                                                                                                                        
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
    if (&rhs == this) return *this;
    unsigned short new_rows = rhs.rows();
    unsigned short new_cols = rhs.cols();
    mat.resize(new_rows);
    for (unsigned short i = 0; i < mat.size(); ++i) {
        mat[i].resize(new_cols);
    }
    for (unsigned short i = 0; i < new_rows; ++i) {
        mat[i] = rhs[i];
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
            result[i][j] = this->mat[i][j] + rhs[i][j];
        }
    }
    return result;
}

// Cumulative addition of this matrix and another                                                                                                                             
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.rows(), ncols = rhs.cols();
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < cols; ++j) {
            this->mat[i][j] += rhs[i][j];
        }
    }
    return *this;
}

// Subtraction of this matrix and another                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.rows(), ncols = rhs.cols();
    Matrix result(nrows, ncols);
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            result[i][j] = this->mat[i][j] - rhs[i][j];
        }
    }
    return result;
}

// Cumulative subtraction of this matrix and another                                                                                                                          
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& rhs) {
    unsigned short nrows = rhs.rows(), cols = rhs.cols();
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            this->mat[i][j] -= rhs[i][j];
        }
    }
    return *this;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& rhs) {
    if (ncols == rhs.rows()) {
        Matrix result(nrows, rhs.cols());
        for (unsigned short i = 0; i < nrows; ++i) {
            for (unsigned short j = 0; j < ncols; ++j) {
                for (unsigned short k = 0; k < nrows; k++) {
                    result[i][j] += this->mat[i][k] * rhs[k][j];
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

// Transpose                                                                                                                                     
template<typename T>
Matrix<T> Matrix<T>::t() {
    Matrix result(nrows, ncols);
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            result[i][j] = this->mat[j][i];
        }
    }
    return result;
}

// Inversion
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
        for (signed int i = 0; i < nrows - 1; ++i) {
            if (up_triang[i][i] != 0) {
                for (unsigned short j = i + 1; j < nrows; ++j) {
                    coeff = up_triang[j][i] / up_triang[i][i];
                    for (unsigned short k = 0; k < ncols; ++k) {
                        up_triang[j][k] = up_triang[j][k]
                            - coeff * up_triang[i][k];
                        result[j][k] = result[j][k]
                            - coeff * result[i][k];
                    }
                }
            }
            else {
                std::vector<T> v;
                for (unsigned short j = i + 1; j < nrows; ++j) {
                    if (up_triang[j][i] != 0) {
                        v = up_triang[i];
                        up_triang[i] = up_triang[j];
                        up_triang[j] = v;
                        v = result[i];
                        result[i] = result[j];
                        result[j] = v;
                        --i;
                        break;
                    }
                    else if (j = nrows - 1) throw "Non-invertible";
                }
            }
        }
        // Get all pivots equal to one on lhs
        for (unsigned short i = 0; i < nrows; ++i) {
            coeff = up_triang[i][i];
            for (unsigned short j = 0; j < nrows; ++j) {
                up_triang[i][j] = up_triang[i][j] / coeff;
                result[i][j] = result[i][j] / coeff;
            }
        }
        // Convert lhs to upper triangular
        for (unsigned short i = nrows - 1; i > 0 && i < nrows; --i) {
            for (unsigned short j = i - 1; j >= 0 && j < nrows - 1; --j) {
                coeff = up_triang[j][i];
                for (unsigned short k = nrows - 1; k >= 0 && k < nrows; --k) {
                    up_triang[j][k] = up_triang[j][k]
                        - coeff * up_triang[i][k];
                    result[j][k] = result[j][k]
                        - coeff * result[i][k];
                }
            }
        }

        return result;
    }
}

// Construction of n-dimensional identity
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
                result[i][j] = this->mat[i][j] + rhs;
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
                result[i][j] = this->mat[i][j] - rhs;
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
                result[i][j] = this->mat[i][j] * rhs;
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
                result[i][j] = this->mat[i][j] / rhs;
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
            Matrix<T> up_triang(nrows, ncols);
            up_triang.mat = mat;
            // Variable to change the sign of the determinant when rows are swapped.
            short sign_by_row_swap = 1;
            T coeff;
            // DO NOT change 'signed int' to 'unsigned short' or smth else below
            // 'signed' bc when up_triang[0][0]=0, then i = -1.
            // 'int' bc the max positive value of signed int is >= that of unsigned short.
            for (signed int i = 0; i < nrows - 1; ++i) { 
                if (up_triang[i][i] != 0) {
                    for (unsigned short j = i + 1; j < nrows; ++j) {
                        coeff = up_triang[j][i] / up_triang[i][i];
                        for (unsigned short k = 0; k < ncols; ++k) {
                            up_triang[j][k] = up_triang[j][k]
                                - coeff * up_triang[i][k];
                        }
                    }
                }
                else {
                    std::vector<T> v;
                    for (unsigned short j = i + 1; j < nrows; ++j) {
                        if (up_triang[j][i] != 0) {
                            v = up_triang[i];
                            up_triang[i] = up_triang[j];
                            up_triang[j] = v;
                            sign_by_row_swap *= -1;
                            --i;
                            break;
                        }
                        else if (j = nrows - 1) return 0;
                    }
                }
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

// Trace
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
    std::vector<T> result;
    for (unsigned short i = 0; i < nrows; ++i) {
        result.push_back(this->mat[i][i]);
    }
    return result;
}

// Access the elements or rows
template<typename T>
std::vector<T>& Matrix<T>::operator[](const int& row) {
    if (row >= 0 && row < nrows) {
        return this->mat[row];
    }
    else throw "Out of limits";
}

template<typename T>
const std::vector<T>& Matrix<T>::operator[](const int& row) const {
    if (row >= 0 && row < nrows) {
        return this->mat[row];
    }
    else throw "Out of limits";
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
inline unsigned short Matrix<T>::rows() const {
    return this->nrows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
inline unsigned short Matrix<T>::cols() const {
    return this->ncols;
}

// Print matrix to console
template<typename T>
void Matrix<T>::print() {
    std::cout << std::endl;
    for (unsigned short i = 0; i < nrows; ++i) {
        for (unsigned short j = 0; j < ncols; ++j) {
            // Round very small doubles to zero
            // Avoids spurious values which should be zero
            if (std::abs(mat[i][j]) < 1.0e-15) mat[i][j] = 0;
            // Align positives with negatives (aesthetic)
            if (mat[i][j] >= 0) std::cout << ' ';
            // Print element
            std::cout << mat[i][j];
            // Jump to next column
            std::cout << '\t';
            // Jump to next row when columns end
            if (j == ncols - 1) std::cout << '\n';
        }
    }
}

#endif
