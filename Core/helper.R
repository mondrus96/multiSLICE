# Extended Bayes Information Criterion
ebic = function(likl, p, n, k, gamma = 0.5){
  return(-2 * likl + k * log(n) + 2 * gamma * k * log(p))
}

# Bayes Information Criterion
bic = function(likl, n, k){
  return(-2*likl + k*log(n))
}

# Log likelihood function
logL = function(Sigma, invSigmahat){
  if(!isPD(invSigmahat)){
    invSigmahat <- makePD(invSigmahat)
  }
  return(log(det(invSigmahat)) - sum(diag(Sigma %*% (invSigmahat))))
}

# Make a matrix positive definite by adding a small value to diagonal
makePD = function(mat){
  p = ncol(mat)
  eigvals = suppressWarnings(eigs(mat, ncol(mat), opts = list(retvec = FALSE))$values)
  perturb = max(max(eigvals) - p*min(eigvals), 0)/(p-1)
  mat = mat+diag(p)*perturb
  return(mat)
}

# Check if a matrix is PD through Cholesky decomposition (faster than full eigendecomp)
isPD = function(mat){
  tryCatch({
    chol(mat)
    return(TRUE)
  }, error = function(e){
    return(FALSE)
  })
}

# Get sin angle between two vectors
sintheta <- function(v_hat, v){
  return(sqrt(1 - sum(v_hat * v)^2))
}

# Calculates the mutual information score
nmi <- function(labels_true, labels_pred) {
  # Adapted from sklearn \url{https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html}
  
  contingency <- table(labels_true, labels_pred) # Create contingency table
  contingency_sum <- sum(contingency) # Total sum of the contingency table
  row_sums <- rowSums(contingency) # Row and column sums
  col_sums <- colSums(contingency)
  
  mi <- 0 # Initialize mutual information
  for (i in seq_len(nrow(contingency))) {
    for (j in seq_len(ncol(contingency))) {
      if (contingency[i, j] > 0) {
        log_contingency_nm <- log(contingency[i, j])
        contingency_nm <- contingency[i, j] / contingency_sum
        outer <- row_sums[i] * col_sums[j]
        log_outer <- -log(outer) + log(sum(row_sums)) + log(sum(col_sums))
        mi_contrib <- contingency_nm * (log_contingency_nm - log(contingency_sum)) + contingency_nm * log_outer
        mi <- mi + mi_contrib
      }
    }
  }
  
  if (mi == 0) { # No need to proceed if mi == 0
    return(0)
  }
  
  h_true <- entropy(labels_true); h_pred <- entropy(labels_pred) # Calculate entropies
  norm <- mean(h_true, h_pred) # Normalization
  return(as.numeric(mi/norm))
}

# Calculates the adjusted rand index
ari <- function(labels_true, labels_pred){
  # Adapted from sklearn \url{https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html}
  
  contingency <- pair_confusion_matrix(labels_true, labels_pred)
  num <- sum(diag(contingency))
  den <- sum(contingency)
  
  if (num == den || den == 0){ 
    return(1)
  } else {
    return(num/den)
  }
}

# Calculate cluster entropy
entropy <- function(labels) {
  # Adapted from sklearn \url{https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/metrics/cluster/_supervised.py#L1261}
  
  if (length(labels) == 0) {
    return(1)
  }
  
  label_counts <- table(labels) # Count the occurrence of each label
  if (length(label_counts) == 1) { # If there's only one label, entropy is 0
    return(0)
  }
  
  total_count <- sum(label_counts) # Calculate probabilities
  probs <- label_counts / total_count
  entropy <- -sum(probs * log(probs)) # Calculate entropy
  return(entropy)
}

# Matching labels through contingency matrix
pair_confusion_matrix <- function(labels_true, labels_pred) {
  # Adapted from sklearn \url{https://github.com/scikit-learn/scikit-learn/blob/main/sklearn/metrics/cluster/_supervised.py#L190}
  
  labels_true <- factor(labels_true) # Ensure labels are factors and have the same levels
  labels_pred <- factor(labels_pred, levels = levels(labels_true))
  
  contingency <- table(labels_true, labels_pred) # Compute contingency matrix
  
  n_samples <- sum(contingency)
  n_c <- rowSums(contingency)
  n_k <- colSums(contingency)
  contingency_values <- as.vector(contingency) # Convert contingency matrix to vector of its values
  sum_squares <- sum(contingency_values^2) # Compute sum of squares of contingency matrix
  
  cmat <- matrix(numeric(4), nrow = 2) # Initialize and compute the pair confusion matrix
  cmat[2, 2] <- sum_squares - n_samples
  cmat[1, 2] <- sum(contingency %*% n_k) - sum_squares
  cmat[2, 1] <- sum(t(contingency) %*% n_c) - sum_squares
  cmat[1, 1] <- n_samples^2 - cmat[1, 2] - cmat[2, 1] - sum_squares
  
  return(cmat) # Return contingency table
}

F1score <- function(S_star, S_hat){
  labels_true <- S_star[upper.tri(S_star)]; labels_pred <- S_hat[upper.tri(S_hat)]
  labels_true <- labels_true != 0; labels_pred <- labels_pred != 0
  
  TP <- sum(labels_true == 1 & labels_pred == 1)
  FP <- sum(labels_true == 0 & labels_pred == 1)
  FN <- sum(labels_true == 1 & labels_pred == 0)
  
  prec <- TP/(TP + FP)
  rec <- TP/(TP + FN)
  
  return(2 * (prec * rec) / (prec + rec))
}

TPrate <- function(S_star, S_hat){
  labels_true <- S_star[upper.tri(S_star)]; labels_pred <- S_hat[upper.tri(S_hat)]
  labels_true <- labels_true != 0; labels_pred <- labels_pred != 0
  
  TP <- sum(labels_true == 1 & labels_pred == 1)
  return(TP/sum(labels_true == 1))
}

TNrate <- function(S_star, S_hat){
  labels_true <- S_star[upper.tri(S_star)]; labels_pred <- S_hat[upper.tri(S_hat)]
  labels_true <- labels_true != 0; labels_pred <- labels_pred != 0
  
  TN <- sum(labels_true == 0 & labels_pred == 0)
  return(TN/sum(labels_true == 0))
}

# For making simulations
makeseq <- function(p, layers){
  seqs <- vector("list", layers)  # Initialize the list to store sequences
  start <- 1  # Starting index for the first sequence
  
  for (i in 1:layers) {
    end <- start + p - 1  # Calculate the ending index
    seqs[[i]] <- start:end  # Assign the sequence to the list
    start <- end + 1  # Update start for the next sequence
  }
  
  return(seqs)
}

# Makes list from matrix data
mat2list <- function(mat, labs) {
  uniqlabs <- unique(unlist(labs)) # Get unique labels
  # List to store decomposed matrices
  matlist <- vector("list", length(labs))
  
  # Function to extract relevant matrix part
  extractMat <- function(lab) {
    inds <- match(lab, uniqlabs) # Find indices of labels in uniqlabs
    submat <- mat[inds, inds] # Extract submatrix for given labels
    submat
  }
  
  # Apply extractMat to each set of labels in labs
  matlist <- lapply(labs, extractMat)
  
  return(matlist)
}

# Converts list to matrix
list2mat <- function(matlist, labs){
  uniqlabs <- unique(unlist(labs)) # Get unique labels
  
  mat <- countmat <- matrix(0, nrow = length(uniqlabs), ncol = length(uniqlabs)) # Create matrices
  labmats <- vector("list", length(matlist)) # Create an empty list to collect index matrices
  
  # Fill labmats with the corresponding index matrices
  for(i in seq_along(matlist)){
    indmat <- matrix(FALSE, nrow = length(uniqlabs), ncol = length(uniqlabs))
    inds <- match(labs[[i]], uniqlabs)
    indmat[inds, inds] <- TRUE
    labmats[[i]] <- indmat
  }
  for(i in seq_along(matlist)) {
    mat[labmats[[i]]] <- mat[labmats[[i]]] + matlist[[i]]
    countmat[labmats[[i]]] <- countmat[labmats[[i]]] + 1
  }
  
  mat <- mat / countmat # Average overlapping values
  mat[is.nan(mat)] <- NA # Replace NaN w/ NA
  colnames(mat) <- rownames(mat) <- uniqlabs
  
  return(mat)
}

# Function to insert a separator row and column
insert_separator <- function(mat, pos) {
  cbind(
    rbind(mat[1:pos,], -1000, mat[(pos+1):nrow(mat),])[, 1:pos],
    -1000,
    rbind(mat[1:pos,], -1000, mat[(pos+1):nrow(mat),])[, (pos+1):ncol(mat)]
  )
}

make_palette <- function(n, s_range = c(0.5, 0.9), v_range = c(0.6, 0.9)) {
  # Generate evenly spaced hues
  hues <- seq(0, 360 * (n - 1) / n, length.out = n)
  
  # Create colors in HSV space with varying saturation and value
  H <- hues
  S <- seq(s_range[1], s_range[2], length.out = n)
  V <- seq(v_range[2], v_range[1], length.out = n)  # Reverse V to make first colors brighter
  
  # Shuffle S and V to avoid predictable patterns
  set.seed(42)  # For reproducibility
  S <- sample(S)
  V <- sample(V)
  
  # Convert to RGB
  rgb_colors <- HSV(H, S, V)
  
  # Ensure colors are distinct enough
  for (i in 2:n) {
    while (color_distance(rgb_colors[i], rgb_colors[i-1]) < 20) {
      S[i] <- runif(1, s_range[1], s_range[2])
      V[i] <- runif(1, v_range[1], v_range[2])
      rgb_colors[i] <- HSV(H[i], S[i], V[i])
    }
  }
  
  return(hex(rgb_colors))
}

# Function to calculate Euclidean distance
eucl_dist <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}

# Calculates Euclidean distance between X1 and X2
dist.mat <- function(X1, X2) {
  X <- rbind(X1, X2) # Combine coordinates
  dist_mat <- as.matrix(dist(X, method = "euclidean")) # Distance matrix
  n1 <- nrow(X1); n2 <- nrow(X2) # Narrow to just interactions b/w 1 and 2
  dist_mat <- dist_mat[1:n1, (n1+1):(n1+n2)] 
  return(dist_mat)
}
# Function to create unique names for coordinates
make.names <- function(coords, prefix, start_index) {
  coords$names <- paste0(prefix, start_index:(start_index + nrow(coords) - 1))
  return(coords)
}
# Function to update names based on matches
update.names <- function(matches, source_coords, target_coords) {
  for (i in 1:nrow(matches)) {
    source_idx <- matches[i, 1]
    target_idx <- matches[i, 2]
    common_name <- source_coords$names[source_idx]
    target_coords$names[target_idx] <- common_name
  }
  return(target_coords)
}

color_distance <- function(color1, color2) {
  rgb1 <- as(color1, "RGB")@coords
  rgb2 <- as(color2, "RGB")@coords
  sqrt(sum((rgb1 - rgb2)^2)) * 255
}