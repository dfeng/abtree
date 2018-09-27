#include <Rcpp.h>

#include "vector.h"
#include "node.h"
#include "abtree.h"

// IntegerVector Predict(Node *root,
//                       NumericVector y, NumericMatrix x,
//                       IntegerVector trt, IntegerVector ncat,
//                       int ntrt) {
//   IntegerVector pred_trt(x.nrow());
//   NumericMatrix pred_prob(x.nrow(), ntrt);
//   // for each data point
//   for (int i = 0; i < x.nrow(); i++) {
//     Node *leaf = PredictNode(root, x(i,_), ncat);
//     pred_prob(i,_) = leaf->blok.p;
//     pred_trt[i] = leaf->blok.opt_trt;
//     leaf->predict_y[trt[i]] += y[i];
//     leaf->predict_n[trt[i]]++;
//   }
//   // for each leaf
//   // TODO: remind myself what this is for? don't think it does anything useful right now...
//   // FillTest(root);
//   return pred_trt;
// }

// given a tree (root), return the leaf that this data ends up at
// xrow is a row of the X matrix, not a column
Node * PredictNode(Node *node,
                   NumericMatrix::Row xrow,
                   IntegerVector ncat, int pred_max_depth) {
  if (!node->left || node->level >= pred_max_depth) {
    return node;
  } else {
    // if numeric
    if (ncat[node->split_col] == 0) {
      if (xrow[node->split_col] <= node->split_tau) {
        return PredictNode(node->left, xrow, ncat, pred_max_depth);
      } else {
        return PredictNode(node->right, xrow, ncat, pred_max_depth);
      }
    // if categorical
    } else {
      if (xrow[node->split_col] == node->split_tau) {
        return PredictNode(node->left, xrow, ncat, pred_max_depth);
      } else {
        return PredictNode(node->right, xrow, ncat, pred_max_depth);
      }
    }
  }
}

// void FillTest(Node *node) {
//   if (!node->left) {
//     Block b = Block(node->predict_y, node->predict_n);
//     node->test_blok = b;
//   } else {
//     FillTest(node->left);
//     FillTest(node->right);
//   }
// }