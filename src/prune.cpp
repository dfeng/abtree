#include <Rcpp.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

// x: test data
void PredictPrune(Node *root, const DoubleVec &y, const DoubleMat &x,
                  const IntVec &trt, IntVec &ncat,
                  DoubleMat &cp_table) {
  int n = cp_table.size();
  // Rprintf("size %d", x.size());
  // std::map<int, double> cp_values;

  // looping through every data point
  for (int j = 0; j < x[0].size(); j++) {
    // node travels down where the data point goes
    Node *node = root;
    // while not at leaf
    while (node->left) {
      if (node->pruned) {
        for (int i = 0; i < n; i++) {
          // Rprintf("cp %d %0.2f. node complexity %0.2f \n", i, cp_table[i][0], node->complexity);
          // if the complexity value is higher
          if (cp_table[i][0] >= node->complexity &&
              (cp_table[i][0] < node->branch || node->branch == -1.0)) {
            // cp_table[i][1] += node->blok.opt_prob;
            if (trt[j] == 0) {
              node->prune_y[0] += y[j];
              node->prune_n[0]++;
            } else {
              node->prune_y[1] += y[j];
              node->prune_n[1]++;
            }
            // Rprintf("node->blok.opt_prob %0.2f\n", node->blok.opt_prob);
          }
        }
        // cp_values[node->id] += node->blok.opt_prob;
      }
      if (ncat[node->split_col] == 0) {
        if (x[node->split_col][j] <= node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      } else {
        if (x[node->split_col][j] == node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      }
    }
    // we've reached a leaf
    for (int i = 0; i < cp_table.size(); i++) {
      if (cp_table[i][0] < node->branch) {
        if (trt[j] == 0) {
          node->prune_y[0] += y[j];
          node->prune_n[0]++;
        } else {
          node->prune_y[1] += y[j];
          node->prune_n[1]++;
        }
        // Rprintf("node->blok.opt_prob %0.2f\n", node->blok.opt_prob);
      }
    }
    // cp_values[]
  }
  // loop through again, this time by node
  FillComplexity(root, cp_table);
}

void FillComplexity(Node *node, DoubleMat &cp_table) {
  // node->print();
  if (node->pruned) {
    for (int i = 0; i < cp_table.size(); i++) {
      if (cp_table[i][0] >= node->complexity &&
          (cp_table[i][0] < node->branch || node->branch == -1.0)) {
        Block b = Block(node->prune_y[0], node->prune_y[1], node->prune_n[0], node->prune_n[1]);
        cp_table[i][1] += b.opt_Q;
      }
    }
  }
  if (node->left) {
    FillComplexity(node->left, cp_table);
    FillComplexity(node->right, cp_table);    
  } else {
    for (int i = 0; i < cp_table.size(); i++) {
      if (cp_table[i][0] < node->branch) {
        Block b = Block(node->prune_y[0], node->prune_y[1], node->prune_n[0], node->prune_n[1]);
        cp_table[i][1] += b.opt_Q;
      }
    }
  }
}

void PredictPrune2(Node *root, const DoubleVec &y, const DoubleMat &x, const IntVec &trt, IntVec &ncat, DoubleMat &cp_table) {

  // we need to keep track of the numbers of treatment match vs non. matches
  // so that we can take an average later
  int n = cp_table.size();

  IntVec match_count(n);
  IntVec no_match_count(n);
  DoubleVec cp_match(n);
  DoubleVec cp_no_match(n);

  // Rprintf("size %d", x.size());
  // std::map<int, double> cp_values;

  // looping through every data point
  for (int j = 0; j < x[0].size(); j++) {
    // node travels down where the data point goes
    Node *node = root;
    // while not at leaf
    while (node->left) {
      if (node->pruned) {
        for (int i = 0; i < n; i++) {
          // Rprintf("cp %d %0.2f. node complexity %0.2f \n", i, cp_table[i][0], node->complexity);
          // if the complexity value is higher
          // branch == -1.0 & pruned is the root
          if (cp_table[i][0] >= node->complexity &&
              (cp_table[i][0] < node->branch || node->branch == -1.0)) {
            if (trt[j] == node->blok.opt_trt) {
              // cp_table[i][1] += y[j];
              cp_match[i] += y[j];
              match_count[i]++;
            } else {
              // cp_table[i][1] -= y[j];
              cp_no_match[i] += y[j];
              no_match_count[i]++;
            }
            // Rprintf("node->blok.opt_prob %0.2f\n", node->blok.opt_prob);
          }
        }
        // cp_values[node->id] += node->blok.opt_prob;
      }
      if (ncat[node->split_col] == 0) {
        if (x[node->split_col][j] <= node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      } else {
        if (x[node->split_col][j] == node->split_tau) {
          node = node->left;
        } else {
          node = node->right;
        }
      }
    }
    for (int i = 0; i < n; i++) {
      if (cp_table[i][0] < node->branch) {
        if (trt[j] == node->blok.opt_trt) {
          // cp_table[i][1] += y[j];
          cp_match[i] += y[j];
          match_count[i]++;
        } else {
          // cp_table[i][1] -= y[j];
          cp_no_match[i] += y[j];
          no_match_count[i]++;
        }
        // Rprintf("node->blok.opt_prob %0.2f\n", node->blok.opt_prob);
      }
    }
    // we've reached the end (leaf)
    // cp_values[]
  }
  for (int i = 0; i < n; i++) {
    cp_table[i][1] = cp_match[i]/match_count[i] - cp_no_match[i]/no_match_count[i];
    // Rprintf("%d: match cp %0.3f %0.3f\n", i, cp_match[i], cp_no_match[i]);
    // Rprintf("%d: match counts %d %d\n", i, match_count[i], no_match_count[i]);
  }
}

// given a complexity parameter, return pruned tree
void PruneTree(Node *node, double complexity) {
  if (node->pruned) {
    // Rprintf("node %d comp %0.2f\n", node->id+1, node->complexity);
    if (node->complexity <= complexity) {
      TruncateNode(node);
      return;
      // Rprintf("shouldn't happen\n");
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
  // TODO: maybe we should delete the other information
}
