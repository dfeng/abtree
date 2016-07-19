#include <Rcpp.h>
#include <vector>

#include "node.h"
#include "abtree.h"


using namespace Rcpp;

/* Converting a BST into an Array */
int NodeToRow(Node *node, DoubleMat &tree_df, double id, bool test) {
  if (node) {
    DoubleVec node_row;
    if (test) {
      node_row = DoubleVec(18);
    } else {
      node_row = DoubleVec(14);
    }

    node_row[0] = id;
    node_row[1] = (double) node->split_col; // var
    node_row[2] = (double) node->split_tau; //tau
    node_row[3] = (double) node->blok.opt_trt; //opt_trt
    node_row[4] = (double) node->blok.n[0];
    node_row[5] = (double) node->blok.y[0];
    node_row[6] = (double) node->blok.n[1];
    node_row[7] = (double) node->blok.y[1];
    node_row[8] = (double) node->total_Q;
    node_row[9] = (double) node->complexity;
    node_row[10] = (double) node->branch;
    node_row[11] = (double) node->pruned;
    tree_df.push_back(node_row);
    int rownum = tree_df.size();
    tree_df[rownum-1][12] = (double) NodeToRow(node->left, tree_df, id*2, test);
    tree_df[rownum-1][13] = (double) NodeToRow(node->right, tree_df, id*2+1, test);
    if (test) {
      tree_df[rownum-1][14] = (double) node->test_blok.n[0];
      tree_df[rownum-1][15] = node->test_blok.p[0];
      tree_df[rownum-1][16] = (double) node->test_blok.n[1];
      tree_df[rownum-1][17] = node->test_blok.p[1];
    }
    return rownum;
  } else {
    return -1;
  }
}

void ExportTree(Node *node, DoubleMat &tree_df, bool test) {
  NodeToRow(node, tree_df, 1.0, test);
}



// void ExportTree(Node *root, DoubleMat &tree_df) {
//   std::queue<Node*> q;
//   int count = 0;
//   for (q.push(root); !q.empty(); q.pop()) {
//     Node* node = q.front();
//     NodeToRow(node, tree_df, count);
//     count++;

//     if (node->left) {
//       q.push(temp_node->left);
//       q.push(temp_node->right);
//     }
//   }
// }

void RowToNode(Node *parent, NumericMatrix &tree_df, int row) {
  Block b(tree_df(row,5), tree_df(row,7),
          tree_df(row,4), tree_df(row,6));
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

Node * ImportTree(NumericMatrix &tree_df) {
  Node *root = new Node();
  RowToNode(root, tree_df, 0);
  return root;
}
