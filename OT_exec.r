






##################################latestfunctions###
run_nmf <- function(list_mat, k,nrun,seed) {
  list_mat <- lapply(list_mat, function(x) as.data.frame(x))
  common_feat <- Reduce(intersect, lapply(list_mat, rownames))
  subset_list_mat <- lapply(list_mat, function(x) x[common_feat, , drop = FALSE])
  norm_list <- lapply(subset_list_mat, function(x) x / rowSums(x))
  norm_list <- lapply(norm_list, function(x) x / norm(as.matrix(x), type = "F"))
  combined_matrix <- do.call(rbind, lapply(norm_list, t))
  nmf_object <- nmf(combined_matrix, k,nrun=nrun,seed=seed, method = "Frobenius")
  return(nmf_object)
}

selecting_nmf_model <- function(nmf_fits, k) {
  if (!is.numeric(k) || k <= 0 || k > length(nmf_fits$fit)) {
    stop("Invalid value for k. Please enter an integer to select the k")
  }
  
  selected_nmf <- nmf_fits$fit[[as.character(k)]]
  
  H_matrix <- basis(selected_nmf)
  W_matrix <- coef(selected_nmf)
  
  nmf.clusters <- max.col(H_matrix)
  names(nmf.clusters) <- rownames(H_matrix)
  
  return(list(H_matrix = H_matrix, 
              W_matrix = W_matrix, 
              nmf.clusters = nmf.clusters))
}


max_clusters=50

elbow_kmeans<-function(H_matrix,max_clusters){
max_clusters <- max_clusters  
wcss <- numeric(max_clusters)  

# Calculate WCSS for each number of clusters
for (k in 1:max_clusters) {
  kmeans_result <- kmeans(H_matrix, centers = k, nstart = 25,iter.max = 100)
  wcss[k] <- kmeans_result$tot.withinss  # Store the total within-cluster sum of squares
}

# Plot the elbow plot
elbow_plot<-plot(1:max_clusters, wcss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "Total Within-Cluster Sum of Squares (WCSS)",
     main = "Elbow Plot for Optimal Number of Clusters")
 return(list(
    wcss = wcss,
    elbow_plot=elbow_plot
  ))
}


Cluster_Representations_bulk <- function(list_SE, kmeans_result, regulizer = 0.2, normalization = TRUE) {
  
  all_clusters <- unlist(lapply(list_SE, function(SE) colData(SE)$clusters))
  
  # Unique labels for patients across all SE objects
  unique_labels <- unique(all_clusters)
  
  # Initialize a matrix to store proportions of each label in each kmeans cluster
  kmeans_clusters <- as.factor(kmeans_result$cluster)  # Get kmeans cluster assignments
  label_proportions <- matrix(0, nrow = length(unique(kmeans_clusters)), ncol = length(unique_labels))
  colnames(label_proportions) <- unique_labels
  rownames(label_proportions) <- paste0("kmeans_cluster_", 1:length(unique(kmeans_clusters)))
  
    for (k in levels(kmeans_clusters)) {
      indices <- which(kmeans_clusters == k)
    
      label_counts <- table(all_clusters[indices])
    
    for (label in unique_labels) {
      label_proportions[as.numeric(k), label] <- ifelse(label %in% names(label_counts), 
                                                         label_counts[label], 0)
    }
    
    # Normalize proportions if desired
    if (normalization) {
      label_proportions[as.numeric(k), ] <- (label_proportions[as.numeric(k), ] + regulizer) /
                                            (sum(label_proportions[as.numeric(k), ]) + length(unique_labels) * regulizer)
    }
  }
  
  return(label_proportions)
}







cost_matrix_centroid_new <- function(list_SE, kmeans_dat, nmf_H_matrix, metric = "cosine") {
  
  # Extract all cluster labels from colData and identify unique labels
  all_clusters <- unlist(lapply(list_SE, function(SE) colData(SE)$clusters))
  unique_labels <- unique(all_clusters)
  kmeans_clusters <- kmeans_dat$cluster
  num_kmeans_clusters <- length(unique(kmeans_clusters))
  
  # Prepare a list to store centroids for each label and each KMeans cluster
  label_centroids <- list()
  
  # Loop over each unique label
  for (label in unique_labels) {
    # Find indices for this label
    label_indices <- which(all_clusters == label)
    
    # Initialize a matrix to store centroids for each KMeans cluster under this label
    label_centroids[[label]] <- matrix(0, nrow = num_kmeans_clusters, ncol = ncol(nmf_H_matrix))
    
    # Loop over each KMeans cluster to calculate its centroid within this label
    for (k in 1:num_kmeans_clusters) {
      # Find indices of samples that belong to both the label and KMeans cluster
      cluster_indices <- label_indices[which(kmeans_clusters[label_indices] == k)]
      
      # Calculate centroid for this label and KMeans cluster by taking the median across rows of the H matrix
      if (length(cluster_indices) > 0) {
        label_centroids[[label]][k, ] <- apply(nmf_H_matrix[cluster_indices, ], 2, median)
      }
    }
    
    # Dynamically assign each label's centroid matrix to a variable in the list
    label_centroids[[label]] <- label_centroids[[label]]
  }
  
  # Initialize the cost matrices between each pair of unique labels
  cost_matrices <- list()
  
  # Generate pairwise cost matrices between each unique label
  for (i in seq_along(unique_labels)) {
    for (j in seq_along(unique_labels)) {
      if (i != j) {
        label_i <- unique_labels[i]
        label_j <- unique_labels[j]
        
        # Retrieve centroids for the two labels
        centroids_i <- label_centroids[[label_i]]
        centroids_j <- label_centroids[[label_j]]
        
        # Initialize cost matrix for this label pair
        cost_matrix <- matrix(NA, nrow = nrow(centroids_i), ncol = nrow(centroids_j))
        
        # Calculate pairwise distances between centroids
        for (m in 1:nrow(centroids_i)) {
          for (n in 1:nrow(centroids_j)) {
            if (metric == "cosine") {
              # Cosine distance calculation
              cost_matrix[m, n] <- 1 - sum(centroids_i[m, ] * centroids_j[n, ]) / 
                                        (sqrt(sum(centroids_i[m, ]^2)) * sqrt(sum(centroids_j[n, ]^2)))
            } else {
              # Euclidean distance calculation
              cost_matrix[m, n] <- sqrt(sum((centroids_i[m, ] - centroids_j[n, ])^2))
            }
          }
        }
        
        # Convert cost matrix to data frame for interpretation
        cost_df <- as.data.frame(cost_matrix)
        
        # Set row and column names dynamically
        rownames(cost_df) <- paste(label_i, "KMeans_Cluster", 1:nrow(centroids_i), sep = "_")
        colnames(cost_df) <- paste(label_j, "KMeans_Cluster", 1:nrow(centroids_j), sep = "_")
        
        # Store cost matrix in the list with a descriptive name
        cost_matrices[[paste(label_i, label_j, sep = "_to_")]] <- list(dis = cost_matrix, cost = cost_df)
        
      }
      
    }
  }
  cost_matrices<-cost_matrices[[1]] #prevents dup
  return(cost_matrices)
}



get_transport_matrix <- function(mass_x, mass_y, cost, p = 2, method = "exact", cost_a = NULL, cost_b = NULL, ...) {
  # Run transport_plan_given_C to get the transport plan
  tplan <- transport_plan_given_C(mass_x = mass_x, 
                                  mass_y = mass_y, 
                                  p = p, 
                                  cost = cost, 
                                  method = method, 
                                  cost_a = cost_a, 
                                  cost_b = cost_b, 
                                  ...)
  
  # Determine dimensions for the transport matrix
  n1 <- length(mass_x)
  n2 <- length(mass_y)
  
  # Initialize an n1 x n2 matrix to store the transported masses
  transport_matrix <- matrix(0, nrow = n1, ncol = n2)
  
  # Populate the matrix with transported masses
  for (i in seq_along(tplan$mass)) {
    row_index <- tplan$from[i]
    col_index <- tplan$to[i]
    transport_matrix[row_index, col_index] <- tplan$mass[i]
  }
  
  # Assign row and column names to indicate which values are being compared
  rownames(transport_matrix) <- paste0("mass_x_", seq_len(n1))
  colnames(transport_matrix) <- paste0("mass_y_", seq_len(n2))
  
  return(transport_matrix)
}









norm_tp<-tp/max(tp)
