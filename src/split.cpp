#include <Rcpp.h>
#include <vector>

#include "node.h"
#include "abtree.h"

/*
*
* @param splitnode pointer to the node that we wish to split
*
* @return ?
*
* @comment should follow the same form as our BuildTree() function
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
               const DoubleVec &y, const DoubleMat &x, const IntVec &trt,
               IntMat &ordering, const IntVec &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int level) {
  int split_col, split_n;
  int split_first, split_last; // first and last for categorical
  int split_type; // optimal split type (ncat)
  // double opt_Q = -DBL_MAX;
  double opt_Q = splitnode->blok.opt_Q; // Do better than parents!! Chinese Motto
  Block opt_left, opt_right;
  double split_tau;

  // id for node
//   splitnode->id = id;
//   id++;

  // calculate current profit and best treatment (actually, we did that already)
  // Rprintf("start/end %d %d\n", start, end);

  // TODO: Levels?
  if (level == max_depth)
    return;

  // looping over the columns
  for (int i = 0; i < ncol; i++) {
    double current_tau;
    int current_split_n = -1;
    int current_split_first = -1, current_split_last = -1;
    Block current_left, current_right;
    // Rprintf("split on col %d\n", i);
    if (ncat[i] == 0) {
      if (!BestSplitNum(y, x[i], trt, ordering[i],
                        start, end, min_bucket,
                        current_left, current_right,
                        current_tau, current_split_n))
        continue;
    } else { // ncat[i] > 0 <==> categorical
      if (!BestSplitCat(y, x[i], trt, ordering[i],
                        ncat[i], start, end, min_bucket,
                        current_left, current_right,
                        current_tau, current_split_first, current_split_last))
        continue;
    }

    double current_Q = current_left.opt_Q + current_right.opt_Q;

    if (current_Q > opt_Q) {
      split_type = ncat[i];
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
      // Rprintf("current_Q:%0.2f split_tau:%0.2f col: %d\n",
      //         current_Q, current_tau, i);
    }
  }
  // Rprintf("finished loop for possible splits\n");
  if (opt_Q == -DBL_MAX || opt_left.total_n < min_split || opt_right.total_n < min_split) {
    // splitnode->is_leaf = true;
    // Rprintf("but we didn't find anything good\n");
    return;
  }
  // Rprintf("And we found good stuff\n");
  // Rprintf("split_col %d split_tau %0.2f start %d end %d\n",
  //         split_col, split_tau, start, end);

  if (split_type > 0) {
    int tempvec[end - start];
    for (int i=0; i < split_last - split_first; i++) {
      tempvec[i] = ordering[split_col][split_first + i];
    }
    for (int i=0; i < split_first - start; i++) {
      tempvec[split_last - split_first + i] = ordering[split_col][start + i];
    }
    for (int i=0; i < end - split_last; i++) {
      tempvec[split_last - start + i] = ordering[split_col][split_last + i];
    }
    for (int i=0; i < end - start; i++) {
      ordering[split_col][start + i] = tempvec[i];
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
  splitnode->left = new Node(opt_left);
  // Rprintf("left: opt_Q %0.2f opt_trt %d\n", opt_left.opt_Q, opt_left.opt_trt);
  splitnode->right = new Node(opt_right);
  // Rprintf("right: opt_Q %0.2f opt_trt %d\n", opt_right.opt_Q, opt_right.opt_trt);

  // Rprintf("recurse!\n");
  // now recurse!
  // if (opt_left.total_n > min_split)
  Partition(splitnode->left, y, x, trt, ordering, ncat, ncol, start, split_n, min_bucket, min_split, max_depth, level+1);
  // if (opt_right.total_n > min_split)
  Partition(splitnode->right, y, x, trt, ordering, ncat, ncol, split_n, end, min_bucket, min_split, max_depth, level+1);
}

bool BestSplitCat(const DoubleVec &y, const DoubleVec &x, const IntVec &trt,
                  const IntVec &ordering, int K, // ncat[i]
                  int start, int end, int min_bucket,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_left, int &split_right) {
  // if I recall correctly, categorical variables are simply converted
  // into integers... which are probably stored as doubles here.
  // we will keep track of counts of each treatment
  // and response for each treatment per category
  int ncat[K][2]; // K by 2
  double ycat[K][2]; // K by 2

  // calculate category values for use later
  int ntot[2] = {0, 0};
  double ytot[2] = {0.0, 0.0};

  for (int k = 0; k < K; k++) {
    ncat[k][0] = 0;
    ncat[k][1] = 0;
    ycat[k][0] = 0.0;
    ycat[k][1] = 0.0;
  }
  int cur_cat = -1;
  double cat_dbl = -1.0;
  for (int i=start; i < end; i++) {
    int o = ordering[i];
    if (x[o] != cat_dbl) { // new category
      cat_dbl = x[o];
      cur_cat = (int) cat_dbl;
    }
    ncat[cur_cat][trt[o]]++;
    ycat[cur_cat][trt[o]] += y[o];
    ntot[trt[o]]++;
    ytot[trt[o]] += y[o];
  }

  // for (int i = 0; i < K; i++) {
  //   Rprintf("currCat:%d, catA: (n %d, y %0.2f), catB: (n %d, y %0.2f)\n", i,
  //           ncat[i][0], ycat[i][0], ncat[i][1], ycat[i][1]);
  // }
  // start looping through the data (ahem, the categories)

  double opt_Q = -DBL_MAX;
  // possible splits are 0, 1, ..., K-1
  int left = start;
  for (int i = 0; i < K; i++) {
    // Rprintf("i:%d o:%d prevX:%0.2f  currX:%0.2f, ncum:%d %d\n", i, o, prevX, x[o], ncum[n][0], ncum[n][1]);
    // check that we have at least the right minimum number of observations
    if (ncat[i][0] > min_bucket & ncat[i][1] > min_bucket &
        ((ntot[0]-ncat[i][0]) > min_bucket) &
        ((ntot[1]-ncat[i][1]) >  min_bucket)) {
      // Rprintf("n:%d left: (A %d,B %d), right: (A %d,B %d)\n", i,
      //         ncat[i][0],
      //         ncat[i][1],
      //         (ntot[0]-ncat[i][0]),
      //         (ntot[1]-ncat[i][1])
      //       );
      // Rprintf("y:%d left: (A %0.2f,B %0.2f), right: (A %0.2f,B %0.2f)\n", i,
      //         ycat[i][0],
      //         ycat[i][1],
      //         (ytot[0]-ycat[i][0]),
      //         (ytot[1]-ycat[i][1])
      //       );

      // for convenience, left is the category being split on and
      // throw everyone else into right branch)
      Block b_left(ycat[i][0], ycat[i][1],
                              ncat[i][0], ncat[i][1]);
      Block b_right(ytot[0]-ycat[i][0],
                               ytot[1]-ycat[i][1],
                               ntot[0]-ncat[i][0],
                               ntot[1]-ncat[i][1]);

      // Rprintf("cat:%d left: %0.2f right: %0.2f \n", i, b_left.opt_Q, b_right.opt_Q);

      if (b_left.opt_Q + b_right.opt_Q > opt_Q) {
        // Rprintf("good?\n");
        opt_Q = b_left.opt_Q + b_right.opt_Q;
        opt_left = b_left;
        opt_right = b_right;
        split_tau = (double) i;
        split_left = left;
        split_right = left+ncat[i][0]+ncat[i][1];
        // Rprintf("tau: %0.2f split_n %d i %d start:%d\n", split_tau, split_n, i, start);
      }
    }
    left += ncat[i][0]+ncat[i][1];
  }

  // if we didn't make any splits, return false
  return (opt_Q != -DBL_MAX);
}

bool BestSplitNum(const DoubleVec &y, const DoubleVec &x,
               const IntVec &trt, const IntVec &ordering,
               int start, int end, int min_bucket,
               Block &opt_left, Block &opt_right,
               double &split_tau, int &split_n) {
  const int len = end - start;

  // we will use ncum to keep track of counts of each treatment
  // and response for each treatment as we proceed through the data
  int ncum[len][2];
  double ycum[len][2];

  // initialize (like defaultdict in python)
  for (int j = 0; j < 2; j++) {
    ncum[0][j] = 0;
    ycum[0][j] = 0.0;
  }

  // calculate cumulative sums for use later
  for (int i=start; i < end; i++) {
    const int o = ordering[i];
    if (i == start) {
      ncum[0][trt[o]]++;
      ycum[0][trt[o]] = y[o];
    } else {
      for (int j = 0; j < 2; j++) {
        ncum[i-start][j] = ncum[i-start-1][j];
        ycum[i-start][j] = ycum[i-start-1][j];
      }
      ncum[i-start][trt[o]]++;
      ycum[i-start][trt[o]] += y[o];
    }
  }

  // start looping through the data

  double opt_Q = -DBL_MAX;
  int n = 0; // convenient to keep track of n thus far, same as i-start+1
  double prevX = x[ordering[start]];

  for (int i = start; i < end; i++) {
    const int o = ordering[i];
    // Rprintf("i:%d o:%d prevX:%0.2f  currX:%0.2f, ncum:%d %d\n", i, o, prevX, x[o], ncum[n][0], ncum[n][1]);

    // check if current x is different, in which case we have passed enough
    // data to consider a split
    if (x[o] > prevX) {
      if (ncum[n-1][0] > min_bucket & ncum[n-1][1] > min_bucket &
            (ncum[len-1][0]-ncum[n-1][0] > min_bucket) &
            (ncum[len-1][1]-ncum[n-1][1]) > min_bucket) {
        // now evaluate the split at prevX

        Block b_left(ycum[n-1][0], ycum[n-1][1],
                                ncum[n-1][0], ncum[n-1][1]);
        Block b_right(ycum[len-1][0]-ycum[n-1][0],
                                 ycum[len-1][1]-ycum[n-1][1],
                                 ncum[len-1][0]-ncum[n-1][0],
                                 ncum[len-1][1]-ncum[n-1][1]);

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
      prevX = x[o];
    }
    n++;
  }

  // if we didn't make any splits, return false
  return (opt_Q != -DBL_MAX);
}
