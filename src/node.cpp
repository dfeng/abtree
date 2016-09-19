#include "vector.h"
#include "node.h"
#include "abtree.h"

extern int loss_type; // variable to determine which loss function to use

Block::Block(NumericVector y0, IntegerVector n0) {
  // assigning to block struct
  y = y0; n = n0;
  int ntrt = y0.size();
  NumericVector p(ntrt);
  total_n = 0;
  
  // calculating p
  opt_prob = -1.0;
  double tot_prob = 0;
  for (int i = 0; i < ntrt; i++) {
    p[i] = y[i] / n[i];
    tot_prob += p[i];
    total_n += n[i];
    if (p[i] > opt_prob) {
      opt_prob = p[i];
      opt_trt = i;
    }
  }
  int loss_type = 1;
  // TODO: loss types
  // calculate optimal Q
  switch(loss_type) {
    // regret/L1
    case 0: opt_Q = total_n * sum(opt_prob - p);
            break;
    // LS/L2
    case 1: opt_Q = total_n * sum(pow(opt_prob - p, 2));
            break;
  }
};

Node::Node() {
  this->left = nullptr;
  this->right = nullptr;
  this->pruned = false;
  this->branch = -1;
  this->complexity = -1;
  this->split_col = -1;
  this->opt_Q = -1.0;
  Block b;
  this->blok = b;
  Block t_b;
  this->test_blok = t_b;
};