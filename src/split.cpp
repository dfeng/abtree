#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
#include <float.h> // DBL_MAX
#include <math.h> // floor

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
*     this involves reordering the ordering matrix to keep groupings intact
*   recursively run Partition on left and right branch
*
*
*/
void Partition(Node *splitnode,
               NumericVector y, NumericMatrix x,
               IntegerVector trt, IntegerMatrix ordering,
               int ntrt, const IntegerVector &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int mtry,
               int level) {
  int split_col, split_n;
  int split_first, split_last; // first and last for categorical
  int split_type; // optimal split type (ncat)

  // double opt_Q = splitnode->blok.opt_Q; // Do better than parents!! aka Chinese Motto
  double opt_Q = -DBL_MAX;
  Block opt_left, opt_right;
  double split_tau;

  // Timer timer;

  if (level == max_depth)
    return;

  // mtry bit (http://gallery.rcpp.org/articles/stl-random-shuffle/)
  IntegerVector col_order = seq_len(ncol) - 1;
  // Rcout << col_order << std::endl;
  std::random_shuffle(col_order.begin(), col_order.end(), randWrapper);
  // Rcout << col_order << std::endl;
  
  // looping over the columns
  for (int t = 0; t < mtry; t++) {
    int i = col_order(t);
    // Rcout << i << t << std::endl;
    double current_tau;
    double current_Q;
    int current_split_n = -1;
    int current_split_first = -1, current_split_last = -1;
    Block current_left, current_right;
    if (ncat[i] == 0) {
      current_Q = BestSplitNum(y, x(_,i), trt, ordering(_,i), ntrt,
                        start, end, min_bucket, min_split,
                        current_left, current_right,
                        current_tau, current_split_n);
      if (current_Q == -DBL_MAX)
        continue;
    } else { // ncat[i] > 0 <==> categorical
      current_Q = BestSplitCat(y, x(_,i), trt, ordering(_,i), ntrt, 
                        ncat[i],
                        start, end, min_bucket, min_split,
                        current_left, current_right,
                        current_tau,
                        current_split_first, current_split_last);
      if (current_Q == -DBL_MAX)
        continue;
    }

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
    // timer.step("loop");
  }

  // timer.step("loop");

  // if every attempted split returned false
  // which should happen when min_bucket/min_split fails every time
  if (opt_Q == -DBL_MAX) {
  // if (opt_Q == splitnode->blok.opt_Q) {
    // splitnode->is_leaf = true;
    // Rprintf("but we didn't find anything good\n");
    return;
  }

  // if things are categorical, then we have to manually move the selected category
  if (split_type > 0) { // categorical
    int tempvec[end - start]; // WARNING
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

  // construct a 'lookup' table to determine which rows are going left or right
  bool which[y.size()]; // WARNING
  std::fill_n(which, y.size(), 0); // false
  for (int i = start; i < split_n; i++) {
    which[ordering(i, split_col)] = true;
  }

  Reorder(split_col, ncol, split_n, start, end, ordering, which);
  // Rcout << "ordering" << std::endl << ordering << std::endl;

  // populating current node
  splitnode->opt_Q = opt_Q;
  splitnode->split_tau = split_tau;
  splitnode->split_col = split_col;

  // splitnode->split_n = split_n;

  // initializing left, right branches
  splitnode->left = new Node(ntrt);
  splitnode->left->blok = opt_left;
  // Rprintf("left: opt_Q %0.2f opt_trt %d\n", opt_left.opt_Q, opt_left.opt_trt);
  splitnode->right = new Node(ntrt);
  splitnode->right->blok = opt_right;
  // Rprintf("right: opt_Q %0.2f opt_trt %d\n", opt_right.opt_Q, opt_right.opt_trt);

  // timer.step("reorder");

  // NumericVector res(timer);
  // res = res / 1000000;
  // Rcout << "split: " << res << std::endl << std::endl;
  
  if (isinf(opt_Q)) return;
  
  // now recurse!
  Partition(splitnode->left, y, x, trt, ordering, ntrt, ncat, ncol, start, split_n, min_bucket, min_split, max_depth, mtry, level+1);
  Partition(splitnode->right, y, x, trt, ordering, ntrt, ncat, ncol, split_n, end, min_bucket, min_split, max_depth, mtry, level+1);
}

int randWrapper(const int n) { return floor(unif_rand() * n); }

double splitCriteria(Block &left, Block &right) {
  // [ (yAL - YBL)^2 + (yAR - YBR)^2 ] / [1/n + ...] / (v + ...)
  // Rcout << "left " << left.mean[1] << left.var[1] << std::endl;
  // Rcout << "left.p: " << left.p << std::endl;
  
  // double bss = sum(pow(left.total_p - left.p, 2)) + sum(pow(right.total_p - right.p, 2));
  double bss = left.ntot * pow(left.mean[1] - left.mean[0], 2) + right.ntot * pow(right.mean[1] - right.mean[0], 2);
  // double nf = 1.0 / left.n[0] + 1.0 / left.n[1] + 1.0 / right.n[0] + 1.0 / right.n[1];
  double wss = left.ntot * (left.var[0] + left.var[1]) + right.ntot * (right.var[0] + right.var[1]);
  // Rprintf("n: %d, %d, %d, %d | var: %0.2f, %0.2f, %0.2f, %0.2f\n", left.n[0], left.n[1], right.n[0],
          // right.n[1], left.var[0], left.var[1], right.var[0], right.var[1]);
  // double wss = left.n[0]*left.var[0] + left.n[1]*left.var[1] + 
     // right.n[0] * right.var[0] + right.n[1]*right.var[1];
  // double wss = sum(((NumericVector) left.n) * left.p * (1-left.p)) + sum(((NumericVector) right.n) * right.p * (1-right.p));
  // Rcout << "bss: " << bss << std::endl;
  // Rcout << "wss: " << wss << std::endl;
  // Rcout << "nf: " << nf << std::endl;
  return bss / wss;
}

double BestSplitNum(NumericVector y, NumericMatrix::Column x,
                  IntegerVector trt, IntegerMatrix::Column ordering,
                  int ntrt, int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_n) {
  const int len = end - start;
  // we will use ncum to keep track of counts of each treatment
  // and response for each treatment as we proceed through the data
  IntegerMatrix ncum(len, ntrt);
  NumericMatrix ycum(len, ntrt);
  NumericMatrix y2cum(len, ntrt);

  // Timer timer;

  // Rprintf("Calculating Cumulative Sums\n");
  // calculate cumulative sums for use later
  for (int i = start; i < end; i++) {
    int o = ordering[i];
    if (i == start) {
      ncum(0,trt[o])++;
      ycum(0,trt[o]) = y[o];
      y2cum(0,trt[o]) = pow(y[o],2);
    } else {
      for (int j = 0; j < ntrt; j++) {
        ncum(i-start,j) = ncum(i-start-1,j);
        ycum(i-start,j) = ycum(i-start-1,j);
        y2cum(i-start,j) = y2cum(i-start-1,j);
      }
      ncum(i-start,trt[o])++;
      ycum(i-start,trt[o]) += y[o];
      y2cum(i-start,trt[o]) += pow(y[o],2);
    }
  }

  
  // timer.step("sums");

  // start looping through the data
  // Rprintf("Start Loop\n");
  double opt_Q = -DBL_MAX;
  int n = 0; // convenient to keep track of n thus far, same as i-start+1
  double prev_x = x[ordering[start]];
  for (int i = start; i < end; i++) {
    // check if current x is different, in which case we have passed enough
    // data to consider a split
    if (x[ordering[i]] > prev_x) {
      bool flags_pass = true;
      // alpha test
      double alpha = n / (double) len;
      if (alpha < 0.1 || alpha > 0.9) flags_pass = false;
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
        // Rcout << "ncum" << (IntegerVector) ncum(len-1,_) << std::endl;
        // Rcout << "ycum" << (NumericVector) ycum(len-1,_) << std::endl;
        // Rcout << "y2cum" << (NumericVector) y2cum(len-1,_) << std::endl;
        Block b_left((NumericVector) ycum(n-1,_),
                     (NumericVector) y2cum(n-1,_),
                     (IntegerVector) ncum(n-1,_));
        Block b_right((NumericVector) ycum(len-1,_)-ycum(n-1,_),
                      (NumericVector) y2cum(len-1,_)-y2cum(n-1,_),
                      (IntegerVector) ncum(len-1,_)-ncum(n-1,_));

        // Rcout << "left" << b_left.mean[0] << b_left.mean[1] << std::endl;
        double Q = splitCriteria(b_left, b_right);
        // Rcout << "Q: " << Q << std::endl;
        if (Q > opt_Q) {
          opt_Q     = Q;
          opt_left  = b_left;
          opt_right = b_right;
          split_tau = (x[ordering[i]] + prev_x)/2; // picking the midpoint
          split_n   = i;
        }
      }
      // now we have a new prev_x;
      prev_x = x[ordering[i]];
    }
    n++;
  }

  // timer.step("stop");

  // NumericVector res(timer);
  // res = res / 1000000;
  // Rcout << "num: " << res << std::endl;


  // if we didn't make any splits (because no buckets), return false
  // return (opt_Q != -DBL_MAX);
  return opt_Q;
}

double BestSplitCat(NumericVector y, NumericMatrix::Column x,
                  IntegerVector trt, IntegerMatrix::Column ordering, 
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
  NumericMatrix y2cat(K,ntrt); // defaults to 0.0
  IntegerVector ntot(ntrt);
  NumericVector ytot(ntrt);
  NumericVector y2tot(ntrt);

  // Timer timer;

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
    y2cat(current_cat,trt[o]) += pow(y[o],2);
    ntot[trt[o]]++;
    ytot[trt[o]] += y[o];
    y2tot[trt[o]] += pow(y[o],2);
  }
  // Rcout << "ncat " << ncat << std::endl;
  // Rcout << "ycat " << ycat << std::endl;

  // timer.step("sums");

  // Rprintf("Start Loop\n");
  double opt_Q = -DBL_MAX;
  int left = start;
  for (int i = 0; i < K; i++) {
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
                   (NumericVector) y2cat(i,_),
                   (IntegerVector) ncat(i,_));
      Block b_right((NumericVector) (ytot - (NumericVector) ycat(i,_)),
                    (NumericVector) (y2tot - (NumericVector) y2cat(i,_)),
                    (IntegerVector) (ntot - (IntegerVector) ncat(i,_)));

      // Rcout << "ly " << (NumericVector) ycat(i,_) << std::endl;
      // Rcout << "ln " << (IntegerVector) ncat(i,_) << std::endl;
      // Rcout << "ry " << (NumericVector) (ytot - (NumericVector) ycat(i,_)) << std::endl;
      // Rcout << "rn " << (IntegerVector) (ntot - (IntegerVector) ncat(i,_)) << std::endl;

      // Rprintf("cat:%d left: %0.2f right: %0.2f \n", i, b_left.opt_Q, b_right.opt_Q);
      double Q = splitCriteria(b_left, b_right);
      if (Q > opt_Q) {
        opt_Q = Q;
        opt_left = b_left;
        opt_right = b_right;
        split_tau = (double) i;
        split_left = left;
        split_right = left + sum(ncat(i,_));
      }
    }
    left += sum(ncat(i,_));
  }

  // timer.step("stop");

  // NumericVector res(timer);
  // res = res / 1000000;
  // Rcout << "cat: " << res << std::endl;

  // if we didn't make any splits, return false
  // return (opt_Q != -DBL_MAX);
  return opt_Q;
}