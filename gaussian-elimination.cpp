/* The workhorse of linear algebra and computational math, Gaussian Elimination is actually
quite a simple algorithm. The following implementation is used on augmented matrices that DO NOT
require 'pivoting'. The pseudocode describes the rough implementation to solve a system of linear
equations using Gaussian-Elimination.

Calculations done on an augmented N x (N+1) matrix W (from the matrix equation Ax = b, aug -> Ab)

// elimination first, then substitution

for i=1 to N-1: // iterate over rows
    // pivot here if needed, but not covered in this algorithm
    for j = i+1 to N: // position of row directly under i
        if W_ii too small: // what is 'too small'?
            throw error -- divisor too small, cannot continue
        m = W_ji / W_ii // create the multiplier which we will use for row subtraction
        for k=(i+1) to (N): // iterating over columns
            W_jk = W_jk - m*W_ik // reassign the value to zero it out

        // here, we have completed ONE full row swap. We can use the total number of row swaps to
        // easily calculate the determinant of this matrix!

if W_NN is too small: // again-- what is too small?
    throw error -- singluar matrix

// back-substitution
Given an N-vector x for storing solution values, we perform the following

x[N] = W_N,N+1 / W_NN
for i=N-1 to 1 (backwards):

    // accumulate the sum of SUM_j=i+1^n (W_ij * x[j])
    s = 0;
    for j=i+1 to N:
        s+=W_ij * x[j] 
    x[i] = (W_i,N+1 - s) / W_ii

*/

#include <iostream>
#include <vector>

void printVector(std::vector<std::vector<double>> &vec){
    /* Print the contents of a given 2-dimensional vector, vec
    in the way it would be viewed on paper.

    Important to note: this function only works using a 2-D
    vector--that is, a vector of vectors

        [ [ ... ],
            ...
          [ ... ] ]

    resembling an N-dimensional matrix. This function does not
    handle any more than 1 nesting of vectors-inside-vectors.
    This can be accomplished using a recursive function to
    operate on each vector within each vector, but is far more
    complex than this short demo entails.

    @param std::vector<double> reference vec
    @rtype void */

    std::size_t vecSize = vec.size();

    for (std::size_t i=0; i<vecSize; i++){
        std::printf("Row %ld  ", i);
        for (std::size_t j=0; j<vec[i].size(); j++){
              std::printf("%f ", vec[i][j]);
        }
        std::printf("%s", "\n");
    }
    std::printf("%s", "\n");
}

void gauss_elim(std::vector<std::vector<double>> &mtx){
    /* Perform Gaussian-Elimination without pivoting on an augmented
    N x N+1 matrix. Modify the matrix in place. The final form of
    this matrix is in echelon form, but not reduced echelon form.
    @param std::vector<std::vector<double>> N x N+1 vector object
    @rtype NULL */

    std::size_t ROWS = mtx.size();    // N
    std::size_t COLS = mtx[0].size(); // N+1

    // N-2 because 0-based and we reduce by going to the next row
    for (std::size_t i=0; i<ROWS-1; i++){ // iterate over rows (i = starting row)
        // pivot here if needed, but not covered in this algorithm
        for (std::size_t j = i+1; j<ROWS; j++){ // position of row directly under i
            //if mtx_ii too small: // what is 'too small'?
            //    throw error -- divisor too small, cannot continue
            double m = mtx[j][i] / mtx[i][i]; // create the multiplier which we will use for row subtraction

            // iterate over each element in the row to complete the row swap
            for (std::size_t k = 0; k<COLS; k++){ // iterating over columns-- COLS-1 = N
                mtx[j][k] = mtx[j][k] - m*(mtx[i][k]); // reassign the value to zero it out

                /* An interesting note is that this procedure 'zeros-out' the entire row, from
                column position k=0 to COLS-1, but it does not need to. Once the matrix has
                been reduced, other operations that rely on this reduced form
                (like computing the determinant, or solving this simple system of linear equations)
                don't actually need those other values, so we can optimize the calculation
                speed by just leaving those values starting at k=i as is. */

            // here, we have completed ONE full row swap. We can use the total number of row swaps to
            // easily calculate the determinant of this matrix!

            }
        }
    }

    //if mtx_NN is too small: // again-- what is too small?
    //    throw error -- singluar matrix
}

int main(){

    // first define the matrix we seek to reduce
    std::vector< std::vector<double> > A = {
        {8, 1, 6, 1},
        {3, 5, 7, 4},
        {4, 9, 2, 2}
    };

    // print out original
    printVector(A);

    // perform elimination
    gauss_elim(A);

    // print out reduced matrix
    std::printf("%s", "Echelon-form by Gaussian elimination:\n");
    printVector(A);

    // another matrix, this time 2 x 3
    std::vector< std::vector<double> > B = {
        {1.0, 2.5,     -1.18},
        {0.0, -5.4837, 3.69 },
    };

    std::printf("%s", "Another matrix:\n");
    printVector(B);
    gauss_elim(B);
    std::printf("%s", "Echelon-form by Gaussian elimination:\n");
    printVector(B);

    return 0;
}
