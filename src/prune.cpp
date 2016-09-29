#include <Rcpp.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

// x: test data
void PredictPrune(Node *root,
                  const NumericVector &y, const NumericMatrix &x,
                  const IntegerVector &trt, IntegerVector &ncat,
                  NumericMatrix &cp_table) {
  int n = cp_table.nrow();
  // std::map<int, double> cp_values;

  // looping through every data point
  for (int j = 0; j < x.nrow(); j++) {
    // node travels down where the data point goes
    Node *node = root;
    // while not at leaf
    while (node->left) {
      if (node->pruned) {
        for (int i = 0; i < n; i++) {
          // if the complexity value is higher
          if (cp_table(i,0) >= node->complexity &&
              (cp_table(i,0) < node->branch || node->branch == -1.0)) {
            node->prune_y[trt[j]] += y[j];
            node->prune_n[trt[j]]++;
            // Rcout << "prune_y " << node->prune_y << std::endl;
          }
        }
        // cp_values[node->id] += node->blok.opt_prob;
      }
      if (ncat[node->split_col] == 0) {
        if (x(j,node->split_col) <= node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      } else {
        if (x(j,node->split_col) == node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      }
    }
    // we've reached a leaf
    for (int i = 0; i < cp_table.nrow(); i++) {
      if (cp_table(i,0) < node->branch) {
        node->prune_y[trt[j]] += y[j];
        node->prune_n[trt[j]]++;
        // Rprintf("node->blok.opt_prob %0.2f\n", node->blok.opt_prob);
      }
    }
    // cp_values[]
  }
  // loop through again, this time by node
  FillComplexity(root, cp_table);
}

void FillComplexity(Node *node, NumericMatrix &cp_table) {
  // node->print();
  if (node->pruned) {
    for (int i = 0; i < cp_table.nrow(); i++) {
      if (cp_table(i,0) >= node->complexity &&
          (cp_table(i,0) < node->branch || node->branch == -1.0)) {
        Block b = Block(node->prune_y, node->prune_n);
        cp_table(i,1) += b.opt_Q;
      }
    }
  }
  if (node->left) {
    FillComplexity(node->left, cp_table);
    FillComplexity(node->right, cp_table);    
  } else {
    for (int i = 0; i < cp_table.nrow(); i++) {
      if (cp_table(i,0) < node->branch) {
        Block b = Block(node->prune_y, node->prune_n);
        cp_table(i,1) += b.opt_Q;
      }
    }
  }
}

// given a complexity parameter, return pruned tree
void PruneTree(Node *node, double complexity) {
  if (node->pruned) {
    // Rprintf("node %d comp %0.2f\n", node->id+1, node->complexity);
    if (node->complexity <= complexity) {
      Rcout << "truncating" << std::endl;
      TruncateNode(node);
      return;
    }
  }
  // Rprintf("if\n");
  if (node->left) {
    // Rprintf("going down %d\n", node->id+1);
    PruneTree(node->left, complexity);
    PruneTree(node->right, complexity);
  }
}

void TruncateNode(Node *node) {
  // Rprintf("truncating %d\n", node->id+1);
  // Rprintf("  left %d\n", node->left->id+1);
  // Rprintf("  right %d\n", node->right->id+1);
  node->left = nullptr;
  node->right = nullptr;
  // Rprintf("done\n");
}
