#' Cluster Bright Spots and Plot Results
#'
#' Performs clustering on embedded bright spots using k-means,
#' selects the optimal number of clusters using silhouette scores,
#' and generates visualizations for both embedding space and
#' original image space.
#'
#' This function is the final step in the starFindr pipeline:
#' read_fits() -> plot_fits() -> find_bright_spots() ->
#' embed_bright_spots() -> cluster_and_plot()
#'
#' @param img_matrix Numeric matrix representing the image
#' @param embed_df Data frame containing embedded coordinates with columns dim1 and dim2
#' @param max_k Maximum number of clusters to try (default: 8)
#'
#' @return A list containing:
#' \itemize{
#'   \item embedding_plot: ggplot object of clustered embedding
#'   \item image_plot: ggplot object with clusters overlaid on image
#'   \item clusters: Data frame with cluster assignments
#'   \item optimal_k: Selected number of clusters
#' }
#'
#' @import ggplot2
#' @importFrom cluster silhouette
#' @importFrom stats kmeans dist aggregate
#' @importFrom utils head
#' @importFrom grDevices hcl
#' @import viridis
#' @importFrom grid unit
#' @export
#' 
#' @examples
#' # result <- cluster_and_plot(img, embed_df)
#' # result$embedding_plot
#' # result$image_plot
#' 
cluster_and_plot <- function(img_matrix, embed_df, max_k = 8) {
  
  data <- embed_df[, c("dim1", "dim2")]
  
  # Silhouette score for k = 2 to max_k
  sil_scores <- sapply(2:max_k, function(k) {
    km  <- kmeans(data, centers = k, nstart = 25, iter.max = 50)
    ss  <- silhouette(km$cluster, dist(data))
    mean(ss[, 3])
  })
  
  best_k <- which.max(sil_scores) + 1  # +1 because we started at k=2
  cat("=== Silhouette Scores ===\n")
  for (i in seq_along(sil_scores)) {
    cat("  k =", i + 1, ":", round(sil_scores[i], 3),
        if (i + 1 == best_k) " <- BEST" else "", "\n")
  }
  cat("Optimal k selected:", best_k, "\n\n")
  
  # Final clustering with best k
  set.seed(42)
  km <- kmeans(data, centers = best_k, nstart = 50, iter.max = 100)
  embed_df$cluster <- as.factor(km$cluster)
  
  # Distinct colors for any k  never repeats
  colors <- hcl(h   = seq(0, 300, length.out = best_k),
                c   = 80,
                l   = 60)
  
  # Plot 1: UMAP embedding with clusters
  p_embed <- ggplot(embed_df, aes(x = dim1, y = dim2, color = cluster)) +
    geom_point(size = 3, alpha = 0.9) +
    stat_ellipse(aes(group = cluster), linetype = "dashed", level = 0.8) +
    scale_color_manual(values = colors, name = "Cluster") +
    theme_dark() +
    labs(title = paste("UMAP ", best_k, "Clusters"),
         x = "UMAP 1", y = "UMAP 2")
  
  # Plot 2: Clusters on original image
  img_log <- log1p(img_matrix - min(img_matrix))
  img_df  <- expand.grid(x = 1:nrow(img_log), y = 1:ncol(img_log))
  img_df$value <- as.vector(img_log)
  
  # Centroid of each cluster for labels
  centroids <- aggregate(cbind(x, y) ~ cluster, data = embed_df, FUN = mean)
  
  p_image <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_viridis(option = "inferno", guide = "none") +
    geom_point(data = embed_df,
               aes(x = x, y = y, color = cluster),
               inherit.aes = FALSE,
               shape = 1, size = 5, stroke = 1.5) +
    geom_label(data = centroids,
               aes(x = x, y = y,
                   label = paste("C", cluster),
                   color = cluster),
               inherit.aes   = FALSE,
               fill          = "black",
               size          = 3,
               fontface      = "bold",
               label.padding = grid::unit(0.15, "lines")) +
    scale_color_manual(values = colors, name = "Cluster") +
    coord_fixed() + theme_void() +
    theme(plot.background   = element_rect(fill = "black"),
          plot.title        = element_text(color = "white",
                                           hjust = 0.5, size = 14),
          legend.text       = element_text(color = "white"),
          legend.title      = element_text(color = "white"),
          legend.background = element_rect(fill = "black")) +
    labs(title = paste("Bright Object Clusters ", best_k, "Groups"))
  
  list(embedding_plot = p_embed,
       image_plot     = p_image,
       clusters       = embed_df,
       optimal_k      = best_k)
}
utils::globalVariables(c("x", "y", "value", "dim1", "dim2", "cluster"))
