#include <Rcpp.h>
// #include <Rcpp/Benchmark/Timer.h>
#include <float.h> // DBL_MAX

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
  int split_col, split_n;
  int split_first, split_last; // first and last for categorical
  int split_type; // optimal split type (ncat)

  double opt_Q = splitnode->blok.opt_Q; // Do better than parents!! aka Chinese Motto
  Block opt_left, opt_right;
  double split_tau;

  // Timer timer;

  if (level == max_depth)
    return;

  // looping over the columns
  for (int i = 0; i < ncol; i++) {
    double current_tau;
    int current_split_n = -1;
    int current_split_first = -1, current_split_last = -1;
    Block current_left, current_right;
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
    // Rprintf("current_Q:%0.2f opt_Q:%0.2f split_tau:%0.2f col: %d\n",
    //         current_Q, opt_Q, current_tau, i);
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

  // timer.step("loop");

  // if every attempted split returned false
  // which should happen when min_bucket/min_split fails every time
  if (opt_Q == splitnode->blok.opt_Q) {
    // splitnode->is_leaf = true;
    // Rprintf("but we didn't find anything good\n");
    return;
  }

  if (split_type > 0) { // categorical
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
      ordering(start + i, split_col) = tempvec[i];
    }
    split_n = start + split_last - split_first;
  }

  // timer.step("reorder1");

  // construct a 'lookup' table to determine which rows are going left or right
  bool which[y.size()];
  std::fill_n(which, y.size(), 0);
  for (int i = start; i < split_n; i++) {
    which[ordering(i, split_col)] = true;
  }

  Reorder(split_col, ncol, split_n, start, end, ordering, which);
  // Rcpp::Rcout << "ordering" << std::endl << ordering << std::endl;
  // Rprintf("ORDERING STOP");

  // Rprintf("opt_Q %0.2f \n", opt_Q);

  // populating current node
  splitnode->opt_Q = opt_Q;
  splitnode->split_tau = split_tau;
  // Rprintf("split: %d\n", split_col);
  splitnode->split_col = split_col;

  // splitnode->split_n = split_n;

  // initializing left, right branches
  splitnode->left = new Node(ntrt);
  splitnode->left->blok = opt_left;
  // Rprintf("left: opt_Q %0.2f opt_trt %d\n", opt_left.opt_Q, opt_left.opt_trt);
  splitnode->right = new Node(ntrt);
  splitnode->right->blok = opt_right;
  // Rprintf("right: opt_Q %0.2f opt_trt %d\n", opt_right.opt_Q, opt_right.opt_trt);

  // timer.step("reorder2");

  // NumericVector res(timer);
  // for (int i=0; i<res.size(); i++) {
  //   res[i] = res[i] / 1000000;
  // }
  // Rcpp::Rcout << res << std::endl;
  // Rprintf("recurse!\n");
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
  // Rprintf("Start Numeric Best Split\n");
  // Rcpp::Rcout << "start num" << std::endl;
  // Rprintf("start/end %d %d\n", start, end);
  const int len = end - start;
  // we will use ncum to keep track of counts of each treatment
  // and response for each treatment as we proceed through the data
  IntegerMatrix ncum(len, ntrt);
  NumericMatrix ycum(len, ntrt);

  // Rprintf("Calculating Cumulative Sums\n");
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
  // Rprintf("Start Loop\n");
  double opt_Q = -DBL_MAX;
  int n = 0; // convenient to keep track of n thus far, same as i-start+1
  double prev_x = x[ordering[start]];
  for (int i = start; i < end; i++) {
    // Rprintf("i:%d o:%d prev_x:%0.2f  currX:%0.2f, ncum:%d %d\n", i, ordering[i], prev_x, x[ordering[i]], ncum(n,0), ncum(n,1));

    // check if current x is different, in which case we have passed enough
    // data to consider a split
    if (x[ordering[i]] > prev_x) {
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
        // now evaluate the split at prev_x
        Block b_left((NumericVector) ycum(n-1,_),
                     (IntegerVector) ncum(n-1,_));
        Block b_right((NumericVector) ycum(len-1,_)-ycum(n-1,_),
                      (IntegerVector) ncum(len-1,_)-ncum(n-1,_));

        if (b_left.opt_Q + b_right.opt_Q > opt_Q) {
          opt_Q = b_left.opt_Q + b_right.opt_Q;
          opt_left = b_left;
          opt_right = b_right;
          split_tau = prev_x;
          split_n = i;
        }
      }
      // now we have a new prev_x;
      prev_x = x[ordering[i]];
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
  // Rprintf("Start Categorical Best Split\n");
  // Rprintf("start/end %d %d\n", start, end);
  IntegerMatrix ncat(K,ntrt); // defaults to 0
  NumericMatrix ycat(K,ntrt); // defaults to 0.0

  IntegerVector ntot(ntrt); // defaults to 0
  NumericVector ytot(ntrt); // defaults to 0.0

  // Rprintf("Calculating Sums\n");
  int current_cat = -1;
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
  // Rcpp::Rcout << "ncat " << ncat << std::endl;
  // Rcpp::Rcout << "ycat " << ycat << std::endl;

  // Rprintf("Start Loop\n");
  double opt_Q = -DBL_MAX;
  int left = start;
  for (int i = 0; i < K; i++) {
    // Rprintf("i: %d\n", i);
    bool flags_pass = true;
    // bucket test
    for (int j = 0; j < ntrt; j++) {
      if ((ncat(i,j) < min_bucket) ||
          (ntot(j)-ncat(i,j) < min_bucket)) {
        flags_pass = false;
        continue;
      }
    }
    // split test
    if ((sum(ncat(i,_)) < min_split) ||
        (sum(ntot-ncat(i,_)) < min_split)) {
      flags_pass = false;
    }
    if (flags_pass) {
      // Rprintf("Passed\n");
      // for convenience, left is the category being split on and
      // throw everyone else into right branch)
      Block b_left((NumericVector) ycat(i,_),
                   (IntegerVector) ncat(i,_));
      Block b_right((NumericVector) (ytot - (NumericVector) ycat(i,_)),
                    (IntegerVector) (ntot - (IntegerVector) ncat(i,_)));

      // Rcpp::Rcout << "ly " << (NumericVector) ycat(i,_) << std::endl;
      // Rcpp::Rcout << "ln " << (IntegerVector) ncat(i,_) << std::endl;
      // Rcpp::Rcout << "ry " << (NumericVector) (ytot - (NumericVector) ycat(i,_)) << std::endl;
      // Rcpp::Rcout << "rn " << (IntegerVector) (ntot - (IntegerVector) ncat(i,_)) << std::endl;

      // Rprintf("cat:%d left: %0.2f right: %0.2f \n", i, b_left.opt_Q, b_right.opt_Q);

      if (b_left.opt_Q + b_right.opt_Q > opt_Q) {
        opt_Q = b_left.opt_Q + b_right.opt_Q;
        opt_left = b_left;
        opt_right = b_right;
        split_tau = (double) i;
        split_left = left;
        split_right = left + sum(ncat(i,_));
      }
    }
    left += sum(ncat(i,_));
  }
  // if we didn't make any splits, return false
  return (opt_Q != -DBL_MAX);
}