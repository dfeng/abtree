#include <Rcpp.h>
#include <float.h> // DBL_MAX
#include <assert.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

// determine the ordered list of complexity parameters/pruned trees
DoubleMat TreeComplexity(Node *root) {
  // init complexity table
  DoubleMat cp_table;

  // fill out the full tree
  DoubleVec cp_row(2);
  cp_row[0] = 0.0; // the complexity parameter for the full tree is 0!
  cp_row[1] = 0.0;
  cp_table.push_back(cp_row);

  while (!root->pruned) {
    double max_complexity = DBL_MAX;
    Node *max_node = nullptr;
    Rcout << "Loop: " << std::endl;

    SetComplexity(root, max_complexity, &max_node);
    Rcout << "Set Branch" << max_node->complexity << std::endl;
    SetBranch(max_node, max_complexity); // we need complexity instead of branch id, because then we have an ordering
    // SetBranch(max_node, max_node->id);

    Rcout << "Done Setting" << std::endl;
    max_node->pruned = true;
    max_node->complexity = max_complexity;

    DoubleVec cp_row(2);
    // cp_row[0] = (double) 1.0; // delete plz
    cp_row[0] = (double) max_complexity;
    cp_row[1] = 0.0; // placeholder where the predicted Q will go

    // push into complexity table
    cp_table.push_back(cp_row);

    Rprintf("max_complexity %0.2f\n", max_complexity);
    Rprintf("max_node %0.2f\n", max_node->complexity);
  }
  return cp_table;
}

// Populating nodes with complexity information
void SetComplexity(Node *node, double &max_complexity, Node **max_node) {
  double complexity;
  Rcout << "max_complexity" << max_complexity << std::endl;
  // if we are not a leaf
  if (node->left && !node->pruned) {
    SetComplexity(node->left, max_complexity, max_node);
    SetComplexity(node->right, max_complexity, max_node);
    node->total_Q = node->left->total_Q + node->right->total_Q;
    Rcout << "total_Q" << node->total_Q << std::endl;
    node->num_leaves = node->left->num_leaves + node->right->num_leaves;
    Rcout << "num_leaves" << node->num_leaves << std::endl;
    complexity = (node->total_Q - node->blok.opt_Q)/(node->num_leaves - 1);
    Rcout << "complexity" << complexity << std::endl;
    assert(complexity >= 0); // with the new measure, not sure if this is true
    // if (complexity < 0) {
    //   Rprintf("tQ %0.2f oQ %0.2f nL %d \n", node->total_Q, node->blok.opt_Q, node->num_leaves);
    // }
    if (complexity < max_complexity) {
      max_complexity = complexity;
      *max_node = node;
    }
  // if this node has already been pruned in a previous iteration,
  // treat it as a leaf
  } else if (node->pruned) {
    node->num_leaves = 1;
    node->total_Q = node->blok.opt_Q;
  // else, it's a leaf
  } else {
    // node->complexity = -1;
    node->num_leaves = 1;
    node->total_Q = node->blok.opt_Q;
  }
}

// find all leaves in this branch that are unique to it (identified by the complexity value)
void SetBranch(Node *node, double branch) {
  // if we are a leaf
  if (!node->left || node->pruned) {
    assert(node->branch != -1); // sanity check
    // if this node hasn't been visited
    // if (node->branch == -1) {
      // that means this
    node->branch = branch;
    // }
  } else {
    if (!node->pruned) {
      SetBranch(node->left, branch);
      SetBranch(node->right, branch);
    }
  }
}