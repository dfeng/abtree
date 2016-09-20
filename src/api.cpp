#include <Rcpp.h>

#include "node.h"
#include "abtree.h"

using namespace Rcpp;

void ExportTree(Node *node, DoubleMat &tree_df) {
  NodeToRow(node, tree_df, 1);
}

/* Converting a BST into an Array */
int NodeToRow(Node *node, DoubleMat &tree_df, int id) {
  if (!node)
    return -1;

  int ntrt = node->blok.y.size();
  DoubleVec node_row(10 + 2*ntrt);

  node_row[0] = id;
  node_row[1] = (double) node->split_col; // var
  node_row[2] = (double) node->split_tau; //tau
  node_row[3] = (double) node->blok.opt_trt; //opt_trt
  node_row[4] = (double) node->total_Q;
  node_row[5] = (double) node->complexity;
  node_row[6] = (double) node->branch;
  node_row[7] = (double) node->pruned;
  for (int i = 0; i < ntrt; i++) {
    node_row[10+i] = (double) node->blok.y[i]; // y first, then n!
    node_row[10+ntrt+i] = (double) node->blok.n[i];
  }
  tree_df.push_back(node_row);
  int rownum = tree_df.size();
  tree_df[rownum-1][8] = (double) NodeToRow(node->left, tree_df, id*2);
  tree_df[rownum-1][9] = (double) NodeToRow(node->right, tree_df, id*2+1);
  return rownum;
}

// ========================  DF -> Tree  ========================  // 


Node ImportTree(NumericMatrix &tree_df, int ntrt) {
  Node root = Node();
  RowToNode(&root, tree_df, 0, ntrt);
  return root;
}

void RowToNode(Node *parent, NumericMatrix &tree_df, int row) {
  Block b((NumericVector) tree_df(row, Range(10,10+)),
          (NumericVector) tree_df(row, Range()));
  // Rprintf("optprob %0.2f \n", b.opt_prob);
  // parent->id = id;
  parent->blok = b;
  parent->split_tau = tree_df(row,2);
  parent->split_col = tree_df(row,1);

  parent->total_Q = tree_df(row,8);
  parent->complexity = tree_df(row,9);
  parent->branch = tree_df(row,10);
  parent->pruned = (bool) tree_df(row,11);

  int childleft_row = (int) tree_df(row,12) - 1;
  if (childleft_row == -2)
    return;
  int childright_row = (int) tree_df(row,13) - 1;

  parent->left = new Node();
  parent->right = new Node();
//   parent->print();
//   Rprintf("left:%d, right:%d\n", childleft_row, childright_row);
  RowToNode(parent->left, tree_df, childleft_row);
  RowToNode(parent->right, tree_df, childright_row);
}