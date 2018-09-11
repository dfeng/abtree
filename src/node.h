#ifndef NODE_H
#define NODE_H

#include <Rcpp.h>
using namespace Rcpp;

/*
  Data Structure defining blocks and nodes of the tree
 */

/*
  Blocks of data that correspond to splits in the decision tree
 */

struct Block {
  NumericVector mean; // mean / trt
  NumericVector var; // var / trt
  // NumericVector y; // response / trt
  // NumericVector yy; // response^2 / trt
  IntegerVector n; // count / trt

  // inferred values
  // NumericVector mean; // average response / trt
  // NumericVector var; // variance / trt
  int ntot;
  // double meantot;

  // optimal
  int opt_trt;
  // double opt_Q;
  double opt_mean;

  // Constructors
  Block() {}; // default constructor, never used
  Block(NumericVector y0, NumericVector yy0, IntegerVector n0);
};

/*
  Node of a decision tree
 */
struct Node {
  // Split Information
  int split_col; // column index of X matrix to split on
  // int split_n; // the index to split on
  double split_tau; // X_col <= tau
  double opt_Q; // the profit corresponding to the best split

  // PRUNING
  double total_Q; // the total Q from here on down (assuming this node is the root)
  int num_leaves;

  double complexity; // complexity value
  bool pruned; // whether or not this non-leaf node is a pruned beginning of branch
  double branch; // which branch this node is a part of (id'd by complexity)

  Block blok;
  Block test_blok;

  NumericVector prune_y;
  NumericVector prune_y2;
  IntegerVector prune_n;
  NumericVector predict_y;
  NumericVector predict_y2;
  IntegerVector predict_n;

  Node *left; // pointer to left branch
  Node *right; // right branch

  // Constructors
  Node() {};
  Node(int ntrt);

  // Helper Functions
  void print();
};

#endif