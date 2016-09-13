#include <Rcpp.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

/*
*
* The inverse of merge sort a.k.a magic
* 
* 3  1
* 2  2
* -  -
* 1  4
* 4  3
*
*    2
*    3
*    -
*    1
*    4
*
*/

void Reorder(int split_col, int ncol,
             int split_n,
             int start, int end,
             IntegerMatrix &ordering) {
  const int n = end - start; // the total number of observations

  // column j of ordering
  for (int j=0; j < ncol; j++) {
    // Rprintf("Reorder: column %d with size %d\n", j, ordering[j].size());
    // we need to take start away from everything
    int leftpos = 0;
    int rightpos = split_n - start;
    int tempvec[n];
    // making sure we're not in the column that's already sorted
    if (j == split_col) continue;

    // data row i in column j
    for (int i=start; i < end; i++) {
      // Rprintf("Reorder: data row %d, we are searching for %d \n", i, ordering[j][i]);
      bool wentright = true;
      // we only need to check from start to split_n (the first group)
      for (int k=start; k < split_n; k++) {
        // do we have a match
        // Rprintf("Reorder: match row %d with value %d\n", k, ordering[split_col][k]);
        if (ordering(i,j) == ordering(k,split_col)) {
          // Rprintf("Reorder: a match! %d\n", ordering[split_col][k]);
          // if so, then that means this guy is in the left group
          tempvec[leftpos] = ordering(i,j);
          leftpos++;
          wentright = false;
          break;
        }
      }
      // if we didn't put them in the left group
      if (wentright) {
        tempvec[rightpos] = ordering(i,j);
        rightpos++;
      }
    }
    // Rprintf("start %d leftpos %d rightpos %d split_n %d end %d\n",
    //         start, leftpos, rightpos, split_n, end);

    // replace column j with tempvec
    for (int i = start; i < end; i++) {
      ordering(i,j) = tempvec[i-start]; // originally was doing [i]
    }
  }
}