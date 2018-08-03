#include "vector.h"
#include "node.h"
#include "abtree.h"

extern int loss_type; // variable to determine which loss function to use

Block::Block(NumericVector y0, NumericVector yy0, IntegerVector n0) {
  // assigning to block struct
  y = y0; yy = yy0; n = n0;
  int ntrt = y0.size();
  NumericVector _mean(ntrt);
  NumericVector _var(ntrt);
  ntot = 0;
  
  // calculating prob
  // opt_prob = -1.0;
  opt_mean = -DBL_MAX;
  double ytot = 0;
  for (int i = 0; i < ntrt; i++) {
    // Rcpp::Rcout << "y" << y[i] << "n" << n[i] << std::endl;
    _mean[i] = y[i] / n[i];
    _var[i] = yy[i] - pow(_mean[i],2) * n[i];
    ytot += y[i];
    ntot += n[i];
    if (_mean[i] > opt_mean) {
      opt_mean = _mean[i];
      opt_trt = i;
    }
  }
  meantot = ytot / ntot;
  int loss_type = 2;
  // TODO: loss types
  // calculate optimal Q
  switch(loss_type) {
    // regret/L1
    case 0: opt_Q = ntot * sum(opt_mean - _mean);
            break;
    // LS/L2
    case 1: opt_Q = ntot * sum(pow(opt_mean - _mean, 2));
            break;
    // WSS/BSS
    case 2: opt_Q = ( sum( pow(meantot - _mean, 2)) ) / sum(((NumericVector) n) * _mean * (1-_mean));
    // case 2: opt_Q = sum(pow(opt_prob - _mean, 2)) / ( (n[0] * n[1] / pow(ntot, 2)) * sum(pow(opt_prob - _mean, 2)) + sum(((NumericVector) n) * _mean * (1-_mean))/pow(ntot, 2) + 1);
    // case 2: opt_Q = ( sum(((NumericVector) n) * _mean * (1-_mean)) / ntot ) / ( sum( ((NumericVector) n) * pow(meantot - _mean, 2)) );
    // case 2: opt_Q = ( sum( ((NumericVector) n) * pow(opt_prob - _mean, 2)) ) / ( sum(y)  - sum(((NumericVector) n) * _mean * _mean) );
    // case 2: opt_Q = -sum(((NumericVector) n) * _mean * (1-_mean));
            break;
  }
  mean = _mean;
  var = _var;
};

Node::Node(int ntrt) {
  this->left = nullptr;
  this->right = nullptr;
  this->pruned = false;
  this->branch = -1;
  this->complexity = -1;
  this->split_col = -1;
  this->opt_Q = -1.0;
  // Block b;
  // Block t_b;
  // this->blok = b;
  // this->test_blok = t_b;
  NumericVector p_y(ntrt);
  this->prune_y = p_y;
  IntegerVector p_n(ntrt);
  this->prune_n = p_n;
  NumericVector r_y(ntrt);
  this->predict_y = r_y;
  IntegerVector r_n(ntrt);
  this->predict_n = r_n;
};

void Node::print() {
  // Rprintf("block n_A: %d, p_B: %0.2f, tau: %0.2f, col: %d, Q: %0.2f\n",
  //         blok.n[0], blok.p[1], split_tau, split_col, total_Q);
  // Rprintf("complexity: %0.2f, branch: %0.2f, pruned: %d\n",
  //         complexity, branch, (int) pruned);
  // Rprintf("prune: (%0.2f, %0.2f) (%d, %d)\n\n", prune_y[0], prune_y[1], prune_n[0], prune_n[1]);
}