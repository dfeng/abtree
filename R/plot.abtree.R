#' plot.abtree
#'
#' @param x an object of class 'abtree'
#' @param margin some margin settings?
#' @param digits the number of digits to show for split points of quantitative variables
#' @param scientific FALSE if tree is to be plotted without scientific notation for split points
#' @param binaryY whether we have binary outcome
#' @param ... other optional arguments
#' @return a plot of the tree
#' @import graphics
#' @export
#'
#' @examples
#' \dontrun{
#' tree <- MakeTree(y ~ grp | hour + browser, data=x)
#' plot(tree)
#' }
plot.abtree <- function (x, margin = 0, digits=3, scientific=FALSE, binaryY = FALSE,
                         ...)
{
  tree <- FormatTree(x)#, child_id=TRUE, digits=digits, scientific=scientific)
  if (nrow(tree) <= 1L)
    stop("A tree with just a root; nothing to plot.")

  temp <- abtreeco(tree)
  xx <- temp$x
  yy <- temp$y
  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- c(0, max(yy)) + max(yy) * c(-margin, margin)
  plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "",
       ...)
  node <- as.numeric(rownames(tree))
  temp <- abtree.branch(xx, yy, node, 1)
  text(xx[1L], yy[1L], "|")
  lines(c(temp$x), c(temp$y))
  
  # write text
  n_leaves <- sum(tree$split_Var == "<leaf>")
  n_parents <- nrow(tree)-n_leaves
  label <- xvals <- yvals <- rep(NA, n_leaves+n_parents*3)
  pos <- rep(1L, n_leaves+n_parents*3)
  i <- 1
  for (j in 1:nrow(tree)) {
    if (tree$split_var[j] == "<leaf>") {
      label[i] <- formatLeaf(tree$optimal_trt[j], 
                             tree[j, 10:(10+2*length(x$trt.levels)-1)], 
                             binaryY, digits)
      xvals[i] <- xx[j]
      yvals[i] <- yy[j]
    } else {
      #split_val <- format(tree$split_value[j], digits=digits)
      label[i] <- tree$split_var[j]
      xvals[i] <- xx[j]
      thisid <- node[j]
      parentid <- match(floor(thisid/2), node)
      yvals[i] <- ifelse(j==1, yy[j], (yy[j]+yy[parentid])*0.5)
      pos[i] <- 3
      i <- i + 1
      label[i] <- ifelse(tree$var_type[j] == "factor",
                      paste0(tree$split_factor[j]),
                      paste0("<=", round(tree$split_value[j], digits=digits)))
      xvals[i] <- xx[tree$childleft_id[j]] + 0.5*(xx[j]-xx[tree$childleft_id[j]])
      yvals[i] <- yy[j]
      pos[i] <- 3
      i <- i + 1
      label[i] <- ifelse(tree$var_type[j] == "factor",
                         paste0("!=", tree$split_factor[j]),
                         paste0(">", round(tree$split_value[j], digits=digits)))
      xvals[i] <- xx[j] + 0.5*(xx[j]-xx[tree$childleft_id[j]])
      yvals[i] <- yy[j]
      pos[i] <- 3
    }
    i <- i + 1
  }
  text(xvals, yvals, pos=pos, label, ...)
  invisible(list(x = xx, y = yy))
}
# 
# plot.legacy <- function(obj) {
#   # TODO: put require as part of package
#   require(igraph)
#   tree <- FormatTree(obj, child_id=TRUE)
#   isParent <- tree$split_var != "<leaf>"
#   parentRows <- tree[isParent,]
#   childRows <- tree[!isParent,]
#   parents <- which(isParent)
#   vertices <- data.frame(name=1:nrow(tree),
#                          label=tree$split_var,
#                          size2=15,
#                          label.cex = 1,
#                          stringsAsFactors=FALSE)
#   vertices$label[!isParent] <- paste0("A: ",
#                                       round(100*childRows$p_A,2),
#                                       "%\n (", childRows$n_A, ")\n B: ",
#                                       round(100*childRows$p_B, 2),
#                                       "%\n (", childRows$n_B, ")\n ",
#                                       childRows$optimal_trt)
#   vertices$label.cex[!isParent] <- 0.8
#   vertices$size2[!isParent] <- 30
# 
#   edge_labels_left <- rep(NA, length(parents))
#   edge_labels_right <- rep(NA, length(parents))
#   edge_labels_left[parentRows$var_type %in% c(0,2)] <- paste0("<=", parentRows$split_value[parentRows$var_type %in% c(0,2)])
#   edge_labels_left[parentRows$var_type == 1] <- parentRows$split_value[parentRows$var_type == 1]
#   edge_labels_right[parentRows$var_type==1] <- paste0("!=", parentRows$split_value[parentRows$var_type == 1])
#   edge_labels_right[parentRows$var_type %in% c(0,2)] <- paste0(">", parentRows$split_value[parentRows$var_type %in% c(0,2)])
#   edges <- data.frame(from=c(parents, parents),
#                       to = c(parentRows$childleft_id,
#                              parentRows$childright_id),
#                       label=c(edge_labels_left, edge_labels_right))
#   g <- graph_from_data_frame(edges, directed=FALSE,
#                              vertices=vertices)
#   par(mar=c(0,0,0,0))
#   #co <- layout.reingold.tilford(g, root=1)
#   co <- layout_(g, as_tree(root=1, rootlevel=1))
#   plot(g, layout=co, vertex.size=30, vertex.shape="rectangle")
#   invisible(g)
# }

# =================  #
# ===  Helpers  ===  #
# =================  #
#' @import utils
abtree_segments <- function (x, ...)
{
  dat <- data.frame(stack(as.data.frame(x$x)), stack(as.data.frame(x$y)))[,
                                                                          c("ind", "values", "values.1")]
  dat <- cbind(head(dat, -1), tail(dat, -1))
  dat <- dat[complete.cases(dat), -4]
  names(dat) <- c("n", "x", "y", "xend", "yend")
  dat
}
abtree.branch <- function (x, y, node, branch)
{
  is.left <- (node%%2L == 0L)
  node.left <- node[is.left]
  parent <- match(node.left/2L, node)
  sibling <- match(node.left + 1L, node)
  temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
  xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp,
              x[sibling], NA)
  yy <- rbind(y[is.left], y[parent], y[parent], y[sibling],
              NA)
  list(x = xx, y = yy)
}

abtreeco <- function (frame)
{
  node <- as.numeric(rownames(frame))
  depth <- tree.depth(node)
  maxDepth <- max(depth)
  is.leaf <- frame$split_var == "<leaf>"
  y <- (1 +maxDepth - depth)/max(depth, 4L)
  ymax <- y[depth==maxDepth][1]
  ymax.1 <- y[depth==maxDepth-1][1]
  ymax <- ymax.1*(2/3)+ymax/3
  y[depth==maxDepth] <- ymax
  x <- double(length(node))
  x[is.leaf] <- seq(sum(is.leaf))
  left.child <- match(node * 2L, node)
  right.child <- match(node * 2L + 1L, node)
  temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
  for (i in rev(temp))
    x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  return(list(x = x, y = y))
}


tree.depth <- function (nodes)
{
  depth <- floor(log(nodes, base = 2) + 1e-07)
  depth - min(depth)
}

formatLeaf <- function(opt_trt, x, binary=TRUE, ndigits=2) {
  ret <- opt_trt
  ntrt <- length(x)/2
  
  for (i in 1:ntrt) {
    pcount <- i+ntrt
    if(binary) {
      ret <- c(ret, paste0(LETTERS[i],": ",
             round(100*x[i]/x[pcount], ndigits),
             "%\n (", x[pcount], ")"))
    } else {
      ret <- c(ret, paste0(LETTERS[i],": ",
                           round(x[i]/x[pcount], ndigits),
                           "\n (", x[pcount], ")"))
    }
  }
  return(paste(ret, collapse="\n"))
  
}
