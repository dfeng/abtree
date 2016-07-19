#ifndef TREE_H
#define TREE_H

/*
 * A block contains y0, y1, n0, n1
 *
 * There should be a way to compare blocks based on the opt_Q
 *
 */
struct Block {
  double y[2];
  int n[2];
  double p[2]; // probability of success
  int total_n;
  // optimal
  int opt_trt;
  double opt_Q;
  double opt_prob;

  // Constructors
  Block() {
    y[0] = -1.0;
    y[1] = -1.0;
    n[0] = -1;
    n[1] = -1;
    p[0] = -1.0;
    p[1] = -1.0;
  }
  Block(double y0, double y1, int n0, int n1, bool max=FALSE) {
    y[0] = y0;
    y[1] = y1;
    n[0] = n0;
    n[1] = n1;
    p[0] = y0 / n0;
    p[1] = y1 / n1;
    total_n = n0 + n1;
    // if (y1 > y0) {
    if (p[1] > p[0]) {
      opt_trt = 1;
      if (max)
        opt_Q = p[1] * total_n;
      opt_prob = p[1];
    } else {
      opt_trt = 0;
      if (max)
        opt_Q = p[0] * total_n;
      opt_prob = p[0];
    }
    if (!max)
      opt_Q = (p[0] - p[1])*(p[0] - p[1])*total_n;
  }
};


// we use pointers here because pointers can point to NULL, which means we don't actually need
// a flag of whether or not it's a leaf, though for pruning purposes it might be helpful
struct Node {
  //int id; // I succumbed – we need an id for each node, to make it easy to reference it later
            // SW: No we don't.
  // Split Information
  int split_col; // column index of X matrix to split on
  // int split_n; // the index to split on
  double split_tau; // X_col <= tau
  double opt_Q; // the profit corresponding to the best split

  // PRUNING
  double total_Q; // the total Q from here on down (assuming this node is the root)
  int num_leaves;

  double complexity; // complexity value
  bool pruned; // whether or not this non-leaf node is a pruned beginning of branch
  double branch; // which branch this node is a part of (id'd by complexity)

  // POST-PRUNING

  double prune_y[2];
  int prune_n[2];

  double predict_y[2];
  int predict_n[2];

  Block blok;
  Block test_blok;

  // bool isTerminal; // is this a leaf?
  Node *left; // pointer to left branch
  Node *right; // right branch

  // Constructors
  Node() {
    this->left = nullptr;
    this->right = nullptr;
    this->pruned = false;
    this->branch = -1;
    this->complexity = -1;
    this->split_col = -1;
    this->opt_Q = -1.0;
    this->prune_y[0] = 0.0;
    this->prune_y[1] = 0.0;
    this->prune_n[0] = 0;
    this->prune_n[1] = 0;
    this->predict_y[0] = 0.0;
    this->predict_y[1] = 0.0;
    this->predict_n[0] = 0;
    this->predict_n[1] = 0;
    Block b;
    this->blok = b;
    Block t_b;
    this->test_blok = t_b;
  }
  // there was a * here
  // but now it's gone
  // and yet
  // nothing has changed
  // why, c++
  // why?
  Node(Block b) {
    this->left = nullptr;
    this->right = nullptr;
    this->pruned = false;
    this->branch = -1;
    this->complexity = -1;
    this->blok = b; // adding the block
    this->split_col = -1;
    this->opt_Q = -1.0;
    this->prune_y[0] = 0.0;
    this->prune_y[1] = 0.0;
    this->prune_n[0] = 0;
    this->prune_n[1] = 0;
    this->predict_y[0] = 0.0;
    this->predict_y[1] = 0.0;
    this->predict_n[0] = 0;
    this->predict_n[1] = 0;
    Block t_b;
    this->test_blok = t_b;
  }

public:
  void print();
};

#endif
