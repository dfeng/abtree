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
  Block(NumericVector y0, IntegerVector n0)  {
  // assigning to block struct
  y = y0; n = n0;
  int ntrt = y0.size();
  NumericVector p(ntrt);
  total_n = 0;
  
  // calculating p
  opt_prob = -1.0;
  double tot_prob = 0.0;
  for (int i = 0; i < ntrt; i++) {
    p[i] = y[i] / n[i];
    tot_prob += p[i];
    total_n += n[i];
    if (p[i] > opt_prob) {
      opt_prob = p[i];
      opt_trt = i;
    }
  }
  // TODO: loss types
  // calculate optimal Q
  switch(1) {
    // regret/L1
    case 0: opt_Q = total_n * sum(opt_prob - p);
            break;
    // LS/L2
    case 1: opt_Q = total_n * sum(pow(opt_prob - p, 2));
            break;
  }
};
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

  // POST-PRUNING

  double prune_y[2];
  int prune_n[2];

  double predict_y[2];
  int predict_n[2];

  Block blok;
  Block test_blok;

  // bool isTerminal; // is this a leaf?
  Node *left; // pointer to left branch
  Node *right; // right branch

  // Constructors
  Node() {
    this->left = nullptr;
    this->right = nullptr;
    this->pruned = false;
    this->branch = -1;
    this->complexity = -1;
    this->split_col = -1;
    this->opt_Q = -1.0;
  };
  // Node(Block b);

  void print();
};

#endif