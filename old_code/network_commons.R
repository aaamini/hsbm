# n <- 100
# y <- sample(1:4,n,replace=T)
# z <- sample(1:4,n,replace=T)
# 
# # label_vec2mat(z,4)
# # label_vec2mat(y,4)
# # t(label_vec2mat(z,4)) %*% label_vec2mat(y,4)
# compute_confusion_matrix(z,y)
# compute_mutual_info(z,y)

compute_confusion_matrix <- function (z, y, K=NULL) {
  # Compute the confusion matrix between labels "y" and "z"
  # z,y Two sets of labels
  # K   number of labels in both "c" and "e"
  
  if (is.null(K)) K = max(c(z,y))
  t(label_vec2mat(z,K)) %*% label_vec2mat(y,K)
}


# M = label_vec2mat(c,K)'*label_vec2mat(e,K);

label_vec2mat <- function(z, K=NULL) {
  if (is.null(K)) K <- max(z)
  temp <- diag(K)
  temp[z,]
}


compute_mutual_info  <- function(z,y) {
  # normMUI Computes the normalized mutual information between two clusters
  #  Labels should be either vectors or n x k matrices

  # c = turn_into_column_ifvec(c);
  # e = turn_into_column_ifvec(e);

  # if ( !is.null(dim(z)) ) z = label_mat2vec(z)
  # if ( !is.null(dim(y)) ) z = label_mat2vec(y)

  CM = compute_confusion_matrix(z,y)
  normCM = CM / sum(CM); # normalized confusion matrix
  IDX = CM == 0 # index of zero elements of CM so that we can avoid them

  jointEnt = - sum( (normCM[!IDX])*log(normCM[!IDX]) )
  indpt = matrix(rowSums(normCM),ncol=1) %*% matrix(colSums(normCM),nrow=1)
  muInfo = sum(normCM[!IDX] * log(normCM[!IDX] / indpt[!IDX]) )

  muInfo / jointEnt
}
