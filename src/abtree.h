#ifndef ABTREE_H
#define ABTREE_H

#include "node.h"
#include "vector.h"

using namespace Rcpp;

// ========================  Split  ========================  // 

void Partition(Node *splitnode,
               const NumericVector &y, const NumericMatrix &x,
               const IntegerVector &trt, IntegerMatrix &ordering,
               int ntrt, const IntegerVector &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int level);

bool BestSplitNum(const NumericVector &y, const NumericVector &x,
                  const IntegerVector &trt, const IntegerVector &ordering,
                  int ntrt, int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_n);

bool BestSplitCat(const NumericVector &y, const NumericVector &x,
                  const IntegerVector &trt, const IntegerVector &ordering,
                  int ntrt, int K,
                  int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau,
                  int &split_left, int &split_right);

// ========================  Reorder  ========================  // 

void Reorder(int split_col, int ncol,
             int split_n,
             int start, int end,
             IntegerMatrix &ordering,
             bool *which);

// ========================  Complexity  ========================  // 


DoubleMat TreeComplexity(Node *root);
void SetComplexity(Node *node, double &max_complexity, Node **max_node);
void SetBranch(Node *node, double branch);

// ========================  Prune  ========================  // 

void PredictPrune(Node *root,
                  const NumericVector &y, const NumericMatrix &x,
                  const IntegerVector &trt, IntegerVector &ncat,
                  NumericMatrix &cp_table);

void FillComplexity(Node *node, NumericMatrix &cp_table);
void PruneTree(Node *root, double complexity);
void TruncateNode(Node *node);

// ========================  API  ========================  // 

void DeleteTree(Node *node);
void ExportTree(Node *node, DoubleMat &tree_df);
int NodeToRow(Node *node, DoubleMat &tree_df, int id);

// ========================  Predict  ========================  // 

IntegerVector Predict(Node *root,
                      const NumericVector &y, const NumericMatrix x,
                      const IntegerVector &trt, const IntegerVector &ncat);
Node * PredictNode(Node *node,
                   const NumericVector &xrow,
                   const IntegerVector &ncat);
void FillTest(Node *node);

#endif
