#ifndef ABTREE_H
#define ABTREE_H

#include "node.h"
#include "vector.h"

using namespace Rcpp;

// =======================================  //
// ===  Declaring Auxiliary Functions  ===  //
// =======================================  //

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
             IntegerMatrix &ordering);

// ========================  Complexity  ========================  // 


DoubleMat TreeComplexity(Node *root);
void SetComplexity(Node *node, double &max_complexity, Node **max_node);
void SetBranch(Node *node, double branch);

// ========================  Prune  ========================  // 


void PredictPrune(Node *root, const DoubleVec &y, const DoubleMat &x,
                  const IntVec &trt, IntVec &ncat,
                  DoubleMat &cp_table);
void FillComplexity(Node *node, DoubleMat &cp_table);
// void PredictPrune2(Node *root, const DoubleVec &y, const DoubleMat &x, const IntVec &trt, IntVec &ncat, DoubleMat &cp_table);
void PruneTree(Node *root, double complexity);
void TruncateNode(Node *node);

// ========================  API  ========================  // 

void DeleteTree(Node *node);
void ExportTree(Node *node, DoubleMat &tree_df);
int NodeToRow(Node *node, DoubleMat &tree_df, int id);
Node * ImportTree(NumericMatrix &tree_df, int ntrt);
void RowToNode(Node *parent, NumericMatrix &tree_df, int row, int ntrt);

DoubleMat NumToDoubleMat(NumericMatrix m);

// predict
// Node * PredictNode(Node *node, NumericMatrix::Row xrow, IntVec &ncat);
// IntegerVector Predict(Node *root, const DoubleVec &y, NumericMatrix x,
//              const IntVec &trt, IntVec &ncat);
// void FillTest(Node *root);

#endif
