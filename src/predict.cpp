#include <Rcpp.h>
#include <vector>

#include "node.h"
#include "abtree.h"

// given a tree (root), return the leaf that this data ends up at
// xrow is a row of the X matrix, not a column
Node * PredictNode(Node *node, NumericMatrix::Row xrow, IntVec &ncat) {
  if (!node->left) {
    return node;
  } else {
    // if numeric
    if (ncat[node->split_col] == 0) {
      if (xrow[node->split_col] <= node->split_tau) {
        return PredictNode(node->left, xrow, ncat);
      } else {
        return PredictNode(node->right, xrow, ncat);
      }
    // if categorical
    } else {
      if (xrow[node->split_col] == node->split_tau) {
        return PredictNode(node->left, xrow, ncat);
      } else {
        return PredictNode(node->right, xrow, ncat);
      }
    }
  }
}

IntegerVector Predict(Node *root, const DoubleVec &y, NumericMatrix x,
             const IntVec &trt, IntVec &ncat) {
  // IntegerVector trt(x.nrow());
  // for each data point
  IntegerVector pred_trt(x.nrow());
  for (int i = 0; i < x.nrow(); i++) {
    Node *leaf = PredictNode(root, x(i,_), ncat);
    pred_trt[i] = leaf->blok.opt_trt;
    if (trt[i] == 0) {
      leaf->predict_y[0] += y[i];
      leaf->predict_n[0]++;
    } else {
      leaf->predict_y[1] += y[i];
      leaf->predict_n[1]++;
    }
  }
  // for each leaf
  FillTest(root);
  return pred_trt;
}

void FillTest(Node *node) {
  if (!node->left) {
    Block b = Block(node->predict_y[0], node->predict_y[1], node->predict_n[0], node->predict_n[1]);
    node->test_blok = b;
  } else {
    FillTest(node->left);
    FillTest(node->right);
  }
}

// int PredictTreatment(Node *root, NumericMatrix::Row xrow, IntVec &ncat)  {
//   Node *leaf = PredictNode(root, xrow, ncat);
//   return leaf->blok.opt_trt;
// }
