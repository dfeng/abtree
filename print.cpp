#include <Rcpp.h>
#include <vector>

#include "node.h"
#include "abtree.h"

using namespace Rcpp;

void Node::print() {
  Rprintf("block n_A: %d, p_B: %0.2f, tau: %0.2f, col: %d, Q: %0.2f\n",
          blok.n[0], blok.p[1], split_tau, split_col, total_Q);
  Rprintf("complexity: %0.2f, branch: %0.2f, pruned: %d\n",
          complexity, branch, (int) pruned);
  Rprintf("prune: (%0.2f, %0.2f) (%d, %d)\n\n", prune_y[0], prune_y[1], prune_n[0], prune_n[1]);

}
