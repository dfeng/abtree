#ifndef NODE_H
#define NODE_H

/*
  Data Structure defining blocks and nodes of the tree
 */

#include <Rcpp.h>
using namespace Rcpp;

/*
  Blocks of data that correspond to splits in the decision tree
 */

struct Block {
  NumericVector y; // response / trt
  IntegerVector n; // count / trt

  // inferred values
  NumericVector p; // average response / trt
  int total_n;

  // optimal
  int opt_trt;
  double opt_Q;
  double opt_prob;

  // Constructors
  Block() {}; // default constructor, never used
  Block(NumericVector y0, IntegerVector n0);
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
  IntegerVector prune_n;
  NumericVector predict_y;
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