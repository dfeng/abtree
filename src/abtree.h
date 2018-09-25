#ifndef ABTREE_H
#define ABTREE_H

#include "node.h"
#include "vector.h"

using namespace Rcpp;

// ========================  Split  ========================  // 

void Partition(Node *splitnode,
               NumericVector y, NumericMatrix x,
               IntegerVector trt, IntegerMatrix ordering,
               int ntrt, const IntegerVector &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int mtry, int split_cond,
               int level);

int randWrapper(const int n);
double splitCriteria(Block &left, Block &right, int &split_cond);

double BestSplitNum(NumericVector y, NumericMatrix::Column x,
                  IntegerVector trt, IntegerMatrix::Column ordering,
                  int ntrt, int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_n, int &split_cond);

double BestSplitCat(NumericVector y, NumericMatrix::Column x,
                  IntegerVector trt, IntegerMatrix::Column ordering,
                  int ntrt, int K,
                  int start, int end,
                  int min_bucket, int min_split,
                  Block &opt_left, Block &opt_right,
                  double &split_tau,
                  int &split_left, int &split_right, int &split_cond);

// ========================  Reorder  ========================  // 

void Reorder(int split_col, int ncol,
             int split_n, int start, int end,
             IntegerMatrix ordering, bool *which);

// ========================  Complexity  ========================  // 


DoubleMat TreeComplexity(Node *root);
void SetComplexity(Node *node, double &max_complexity, Node **max_node);
void SetBranch(Node *node, double branch);

// ========================  Prune  ========================  // 

void PredictPrune(Node *root,
                  NumericVector y, NumericMatrix x,
                  IntegerVector trt, IntegerVector ncat,
                  NumericMatrix cp_table);

void FillComplexity(Node *node, NumericMatrix cp_table);
void PruneTree(Node *root, double complexity);
void TruncateNode(Node *node);

// ========================  API  ========================  // 

void DeleteTree(Node *node);
void ExportTree(Node *node, DoubleMat &tree_df);
int NodeToRow(Node *node, DoubleMat &tree_df, int id);

// ========================  Predict  ========================  // 

// IntegerVector Predict(Node *root,
//                       NumericVector y, NumericMatrix x,
//                       IntegerVector trt, IntegerVector ncat,
//                       int ntrt);
Node * PredictNode(Node *node,
                   NumericMatrix::Row xrow,
                   IntegerVector ncat);
// void FillTest(Node *node);

#endif
