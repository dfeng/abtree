#ifndef ABTREE_H
#define ABTREE_H

typedef std::vector<double> DoubleVec;
typedef std::vector<int> IntVec;
typedef std::vector<DoubleVec> DoubleMat;
typedef std::vector<IntVec> IntMat;

#include "node.h"


using namespace Rcpp;

DoubleMat NumToDoubleMat(NumericMatrix m);

// =======================================  //
// ===  Declaring Auxiliary Functions  ===  //
// =======================================  //

// split

void Partition(Node *splitnode,
               const DoubleVec &y, const DoubleMat &x, const IntVec &trt,
               IntMat &ordering, const IntVec &ncat,
               int ncol, int start, int end,
               int min_bucket, int min_split, int max_depth,
               int level);

bool BestSplitCat(const DoubleVec &y, const DoubleVec &x, const IntVec &trt,
                  const IntVec &ordering, int K,
                  int start, int end, int min_bucket,
                  Block &opt_left, Block &opt_right,
                  double &split_tau, int &split_left, int &split_right);

bool BestSplitNum(const DoubleVec &y, const DoubleVec &x,
               const IntVec &trt, const IntVec &ordering,
               int start, int end, int min_bucket,
               Block &opt_left, Block &opt_right,
               double &split_tau, int &split_n);

// reorder

void Reorder(int split_col, int ncol,
             int split_n,
             int start, int end,
             IntMat &ordering);

// prune

DoubleMat TreeComplexity(Node *root);
void SetComplexity(Node *node, double &max_complexity, Node **max_node);
void SetBranch(Node *node, double branch);

void PredictPrune(Node *root, const DoubleVec &y, const DoubleMat &x,
                  const IntVec &trt, IntVec &ncat,
                  DoubleMat &cp_table);
void FillComplexity(Node *node, DoubleMat &cp_table);
void PredictPrune2(Node *root, const DoubleVec &y, const DoubleMat &x, const IntVec &trt, IntVec &ncat, DoubleMat &cp_table);
void PruneTree(Node *root, double complexity);
void TruncateNode(Node *node);

// api

int NodeToRow(Node *node, DoubleMat &tree_df, double id, bool test);
void ExportTree(Node *node, DoubleMat &tree_df, bool test = FALSE);
void RowToNode(Node *parent, NumericMatrix &tree_df, int id);
Node * ImportTree(NumericMatrix &tree_df);

// predict
Node * PredictNode(Node *node, NumericMatrix::Row xrow, IntVec &ncat);
IntegerVector Predict(Node *root, const DoubleVec &y, NumericMatrix x,
             const IntVec &trt, IntVec &ncat);
void FillTest(Node *root);

#endif
