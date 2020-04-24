# Matrix
## Warning
- All operations and functions are intended for numeric types only.
- Be warned that, generally, determinant and inversion do not work properly for integers as a consequence of division taking place in these calculations. Use float or an equivalent in these cases.
- Matrix row and column indices ```x[row][col]``` start at 0 and end at 255.
- Row and column size limited to 256 for perfomance. If planning to increase them, have in mind that inversion and determinant algorithms execute in cubic time.
## Declaration and construction
```
Matrix<type> x(rows, columns);
```
## Copy constructor
```
Matrix<type> y(rows, columns);
Matrix<type> x = y;
```
## Assignment operator
```
Matrix<type> x(rows,columns), y(rows_y,columns_y);
x = y;
```
Now the ```x``` matrix has ```rows_y``` rows and ```columns_y``` columns.
## Access/assign
Access whole rows.
```
// Access
x[row];
// Assign (vector variable or value)
std::vector<type> y;
x[row] = y;
x[row] = {1, 2, 3, 4}
```
Access/assign elements.
```
// Access
x[row][column];
// Assign (scalar variable or value)
type y;
x[row][column] = y;
x[row][column] = 1;
```
A clever way to access columns is to store the transpose matrix, then access its rows.
```
Matrix y = x.t();
y[column];
```
The row of the ```y``` matrix is the column of the x matrix.
DO NOT access with ```x.t()[column]```. While it is possible, it will transpose the whole matrix each time, which is inefficient. Just store the transpose and access the rows as done above.
## Get the number of rows or columns
```
x.rows();
x.cols();
```
## Operations with matrices or scalars
### Addition, substraction, product, and cumulative variants
Let ```y``` be a matrix or a scalar,
```
Matrix<type> x(rows,columns);
x + y, x - y, x * y;
x += y, x -= y, x *= y;
```
There is also the division by scalar,
```
x / y;
```
All scalars must be on the right hand side, never on the left hand side.
## Transposition
Returns an object of class ```Matrix```.
```
x.t();
```
## Inversion
Returns an object of class ```Matrix```.
```
x.inv();
```
## Identity
Does not return anything. Transforms a square matrix into a n-dimensional identity matrix.
```
x.identity();
```
## Determinant
Returns a value of the same type of the elements in the matrix.
```
x.det();
```
## Trace
Returns a value of the same type of the elements in the matrix.
```
x.trace();
```
## Diagonal
Returns a vector whose elements are of the same type of the elements in the matrix.
```
x.diagonal();
```
## Print the matrix
```
x.print();
```
###### Disclaimer
This is based on the original work of Michael L. Halls-Moore as found in Quantstart.com and his book "C++ for Quantitative Finance. I have modified it thoroughly to fit my needs and tastes for academic purposes.
