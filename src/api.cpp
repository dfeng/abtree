#include <Rcpp.h>

#include "node.h"
#include "abtree.h"

using namespace Rcpp;

void DeleteTree(Node *node) {
  if (node->left) {
    DeleteTree(node->left);
    DeleteTree(node->right);
  }
  delete node;
}

// ========================  Tree -> DF  ========================  // 

void ExportTree(Node *node, DoubleMat &tree_df) {
  NodeToRow(node, tree_df, 0);
}

/* Converting a BST into an Array
 *  - id is the first empty slot (avaiable)
 *  @return the new first empty slot
 */
int NodeToRow(Node *node, DoubleMat &tree_df, int id) {
  // if (!node)
  //   return id+1;
  int ntrt = node->blok.n.size();
  DoubleVec node_row(10 + 2*ntrt);
  node_row[0] = id;
  node_row[1] = (double) node->split_col; // var
  node_row[2] = (double) node->split_tau; //tau
  node_row[3] = (double) node->blok.opt_trt; //opt_trt
  node_row[4] = (double) node->opt_Q;
  node_row[5] = (double) node->complexity;
  node_row[6] = (double) node->branch;
  node_row[7] = (double) node->level; // (before) pruned
  for (int i = 0; i < ntrt; i++) {
    node_row[10+i] = (double) node->blok.mean[i]; // mean first, then n!
    node_row[10+ntrt+i] = (double) node->blok.n[i];
  }
  tree_df.push_back(node_row);
  // int rownum = tree_df.size();
  if (!node->left)
    return id+1;
  int left_id = (double) NodeToRow(node->left, tree_df, id+1);
  tree_df[id][8] = id+1;
  int right_id = (double) NodeToRow(node->right, tree_df, left_id);
  tree_df[id][9] = left_id;
  return right_id;
}