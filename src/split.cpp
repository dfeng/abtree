#include <Rcpp.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

/*
*
* @param splitnode pointer to the node that we wish to split
*
* @return ?
*
* @pseudo
*   for columns of X
*     run function eval() on column i restricted to n_1, n_2 (this returns the maximum)
*   pick the column with the best max
*   run function split(column i, split_n) which will do the actual splitting
*     this involves reordering the ordering() matrix to keep groupings intact
*   recursively run Partition on left and right branch
*
*
*/
void Partition(Node *splitnode,
               const NumericVector &y, const NumericMatrix &x,
               const IntegerVector &trt, IntegerMatrix &ordering,
               int ntrt, const IntegerVector &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int level) {
  Rprintf("Entering Partition\n");
  Rprintf("start/end %d %d\n", start, end);
  int split_col, split_n;
  int split_first, split_last; // first and last for categorical
  int split_type; // optimal split type (ncat)

  double opt_Q = splitnode->blok.opt_Q; // Do better than parents!! aka Chinese Motto
  Block opt_left, opt_right;
  double split_tau;

  if (level == max_depth)
    return;

  // looping over the columns
  for (int i = 0; i < ncol; i++) {
    double current_tau;
    int current_split_n = -1;
    int current_split_first = -1, current_split_last = -1;
    Block current_left, current_right;
    Rprintf("split on col %d\n", i);
    if (ncat[i] == 0) {
      if (!BestSplitNum(y, x(_,i), trt, ordering(_,i), ntrt,
                        start, end, min_bucket, min_split,
                        current_left, current_right,
                        current_tau, current_split_n))
        continue;
    } else { // ncat[i] > 0 <==> categorical
      if (!BestSplitCat(y, x(_,i), trt, ordering(_,i), ntrt, 
                        ncat[i],
                        start, end, min_bucket, min_split,
                        current_left, current_right,
                        current_tau,
                        current_split_first, current_split_last))
        continue;
    }

    double current_Q = current_left.opt_Q + current_right.opt_Q;
    Rprintf("current_Q:%0.2f opt_Q:%0.2f split_tau:%0.2f col: %d\n",
            current_Q, opt_Q, current_tau, i);
    if (current_Q > opt_Q) {
      split_type = ncat[i]; // categorical?
      opt_Q = current_Q;
      opt_left = current_left;
      opt_right = current_right;

      split_tau = current_tau;
      split_col = i;
      if (split_type == 0) {
        split_n = current_split_n;
      } else {
        split_first = current_split_first;
        split_last = current_split_last;
      }
    }
  }
  Rprintf("finished loop for possible splits\n");
  // if every attempted split returned false
  // which should happen when min_bucket/min_split fails every time
  if (opt_Q == splitnode->blok.opt_Q) {
    // splitnode->is_leaf = true;
    Rprintf("but we didn't find anything good\n");
    return;
  }
  Rprintf("And we found good stuff\n");
  // Rprintf("split_col %d split_tau %0.2f start %d end %d\n",
  //         split_col, split_tau, start, end);

  if (split_type > 0) {
    int tempvec[end - start];
    for (int i = 0; i < split_last - split_first; i++) {
      tempvec[i] = ordering(split_first + i, split_col);
    }
    for (int i = 0; i < split_first - start; i++) {
      tempvec[split_last - split_first + i] = ordering(start + i, split_col);
    }
    for (int i = 0; i < end - split_last; i++) {
      tempvec[split_last - start + i] = ordering(split_last + i, split_col);
    }
    for (int i = 0; i < end - start; i++) {
      ordering(start + i,split_col) = tempvec[i];
    }
    split_n = start + split_last - split_first;
  }

  // Reorder the ordering matrix
  // Rprintf("ORDERING START");
  Reorder(split_col, ncol, split_n, start, end, ordering);
  // Rprintf("ORDERING STOP");

  // Rprintf("opt_Q %0.2f \n", opt_Q);

  // populating current node
  splitnode->opt_Q = opt_Q;
  splitnode->split_tau = split_tau;
  // Rprintf("split: %d\n", split_col);
  splitnode->split_col = split_col;

  // splitnode->split_n = split_n;

  // initializing left, right branches
  splitnode->left = new Node();
  splitnode->left->blok = opt_left;
  // Rprintf("left: opt_Q %0.2f opt_trt %d\n", opt_left.opt_Q, opt_left.opt_trt);
  splitnode->right = new Node();
  splitnode->right->blok = opt_right;
  // Rprintf("right: opt_Q %0.2f opt_trt %d\n", opt_right.opt_Q, opt_right.opt_trt);

  Rprintf("recurse!\n");
  // now recurse!
  // if (opt_left.total_n > min_split)
  Partition(splitnode->left, y, x, trt, ordering, ntrt, ncat, ncol, start, split_n, min_bucket, min_split, max_depth, level+1);
  // if (opt_right.total_n > min_split)
  Partition(splitnode->right, y, x, trt, ordering, ntrt, ncat, ncol, split_n, end, min_bucket, min_split, max_depth, level+1);
}

bool BestSplitNum(const NumericVector &y, const NumericVector &x,
                  const IntegerVector &trt, const IntegerVector &ordering,
                  int ntrt, int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_n) {
  Rprintf("Start Numeric Best Split\n");
  Rprintf("start/end %d %d\n", start, end);
  const int len = end - start;

  // we will use ncum to keep track of counts of each treatment
  // and response for each treatment as we proceed through the data
  IntegerMatrix ncum(len, ntrt);
  NumericMatrix ycum(len, ntrt);

  Rprintf("Calculating Cumulative Sums\n");
  // calculate cumulative sums for use later
  for (int i=start; i < end; i++) {
    int o = ordering[i];
    if (i == start) {
      ncum(0,trt[o])++;
      ycum(0,trt[o]) = y[o];
    } else {
      for (int j = 0; j < ntrt; j++) {
        ncum(i-start,j) = ncum(i-start-1,j);
        ycum(i-start,j) = ycum(i-start-1,j);
      }
      ncum(i-start,trt[o])++;
      ycum(i-start,trt[o]) += y[o];
    }
  }

  // start looping through the data
  Rprintf("Start Loop\n");
  double opt_Q = -DBL_MAX;
  int n = 0; // convenient to keep track of n thus far, same as i-start+1
  double prevX = x[ordering[start]];

  for (int i = start; i < end; i++) {
    // Rprintf("i:%d o:%d prevX:%0.2f  currX:%0.2f, ncum:%d %d\n", i, o, prevX, x[o], ncum(n,0), ncum(n,1));

    // check if current x is different, in which case we have passed enough
    // data to consider a split
    if (x[ordering[i]] > prevX) {
      bool flags_pass = true;
      // bucket test
      for (int j = 0; j < ntrt; j++) {
        if ((ncum(n-1,j) < min_bucket) ||
            (ncum(len-1,j)-ncum(n-1,j) < min_bucket)) {
          flags_pass = false;
          continue;
        }
      }
      // split test
      if ((sum(ncum(n-1,_))) < min_split ||
          (sum(ncum(len-1,_)-ncum(n-1,_))) < min_split) {
        flags_pass = false;
      }
      
      if (flags_pass) {
        Rprintf("bucket pass\n");
        // now evaluate the split at prevX
        Block b_left((NumericVector) ycum(n-1,_),
                     (IntegerVector) ncum(n-1,_));
        Block b_right((NumericVector) ycum(len-1,_)-ycum(n-1,_),
                      (IntegerVector) ncum(len-1,_)-ncum(n-1,_));

        if (b_left.opt_Q + b_right.opt_Q > opt_Q) {
          opt_Q = b_left.opt_Q + b_right.opt_Q;
          opt_left = b_left;
          opt_right = b_right;
          split_tau = prevX;
          split_n = i;
          // Rprintf("tau: %0.2f split_n %d i %d n:%d start:%d, left: n(%d,%d) p(%0.2f, %0.2f), right:n (%d, %d) p(%0.2f, %0.2f) opt: %d %d\n", split_tau, split_n, i, n, start,
          //         b_left.n[0], b_left.n[1],
          //         b_left.y[0], b_left.y[1],
          //         b_right.n[0], b_right.n[1],
          //         b_right.y[0], b_right.y[1],
          //         b_left.opt_trt, b_right.opt_trt);
        }
      }
      // now we have a new prevX;
      prevX = x[ordering[i]];
    }
    n++;
  }

  // if we didn't make any splits (because no buckets), return false
  return (opt_Q != -DBL_MAX);
}


bool BestSplitCat(const NumericVector &y, const NumericVector &x,
                  const IntegerVector &trt, const IntegerVector &ordering, 
                  int ntrt, int K,
                  int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau,
                  int &split_left, int &split_right) {
  Rprintf("Start Categorical Best Split\n");
  Rprintf("start/end %d %d\n", start, end);
  IntegerMatrix ncat(K, ntrt); // defaults to 0
  NumericMatrix ycat(K, ntrt); // defaults to 0.0

  IntegerVector ntot(ntrt); // defaults to 0
  NumericVector ytot(ntrt); // defaults to 0.0

  Rprintf("Calculating Sums\n");
  int current_cat = -INT_MAX;
  double prev_cat = -DBL_MAX;
  for (int i = start; i < end; i++) {
    int o = ordering[i];
    if (x[o] != prev_cat) { // new category
      prev_cat = x[o];
      current_cat = (int) prev_cat;
    }
    ncat(current_cat,trt[o])++;
    ycat(current_cat,trt[o]) += y[o];
    ntot[trt[o]]++;
    ytot[trt[o]] += y[o];
  }

  for (int i = 0; i < K; i++) {
    for (int j = 0; j < ntrt; j++) {
      Rprintf("cat,%d trt,%d (n %d, y %0.2f)\n",
              i, j, ncat(i,j), ycat(i,j));
    }
  }

  double opt_Q = -DBL_MAX;

  Rprintf("Start Loop\n");
  // possible splits are 0, 1, ..., K-1
  int left = start;
  for (int i = 0; i < K; i++) {
    // Rprintf("i:%d o:%d prevX:%0.2f  currX:%0.2f, ncum:%d %d\n", i, o, prevX, x[o], ncum[n][0], ncum[n][1]);
    // check that we have at least the right minimum number of observations
    bool bucket_pass = true;
    for (int j = 0; j < ntrt; j++) {
      if ((ncat(i,j) < min_bucket) || (ntot(j)-ncat(i,j) < min_bucket))
        bucket_pass = false;
    }
    if (bucket_pass) {

      // for convenience, left is the category being split on and
      // throw everyone else into right branch)
      Block b_left((NumericVector) ycat(i,_),
                   (IntegerVector) ncat(i,_));
      Block b_right((NumericVector) ytot-ycat(i,_),
                    (IntegerVector) ntot-ncat(i,_));

      Rprintf("cat:%d left: %0.2f right: %0.2f \n", i, b_left.opt_Q, b_right.opt_Q);

      if (b_left.opt_Q + b_right.opt_Q > opt_Q) {
        opt_Q = b_left.opt_Q + b_right.opt_Q;
        opt_left = b_left;
        opt_right = b_right;
        split_tau = (double) i;
        split_left = left;
        split_right = left + sum(ncat(i,_));
        // Rprintf("tau: %0.2f split_n %d i %d start:%d\n", split_tau, split_n, i, start);
      }
    }
    left += sum(ncat(i,_));
  }

  // if we didn't make any splits, return false
  return (opt_Q != -DBL_MAX);
}