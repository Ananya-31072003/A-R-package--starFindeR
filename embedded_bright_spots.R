#' Embed Bright Sources in 2D Space using UMAP or t-SNE
#'
#' Extracts a small square pixel patch around each detected bright source,
#' normalises each patch independently, and projects all patches into 2D
#' using either UMAP or t-SNE. The resulting 2D coordinates reveal
#' structure among the sources based on the shape of their local pixel
#' neighbourhood -- for example separating point sources (stars) from
#' extended sources (galaxies, nebulae) or artefacts.
#'
#' @param img_matrix numeric matrix. The image, as returned by
#'   \code{\link{read_fits}}.
#' @param spots_df data frame. Detected bright sources, as returned by
#'   \code{\link{find_bright_spots}}. Must contain columns \code{x},
#'   \code{y}, and \code{value}.
#' @param method character. Dimensionality reduction algorithm to use:
#'   \itemize{
#'     \item \code{"umap"} (default) -- Uniform Manifold Approximation
#'       and Projection. Faster than t-SNE, better at preserving global
#'       structure. Good first choice for most images.
#'     \item \code{"tsne"} -- t-distributed Stochastic Neighbour
#'       Embedding. Slower but often produces tighter, more visually
#'       separated clusters. Use when UMAP clusters overlap too much.
#'   }
#' @param patch_size integer. Side length in pixels of the square patch
#'   extracted around each source. Must be an odd number.
#'   Default \code{7} gives a 7x7 = 49-pixel feature vector per source.
#'   \itemize{
#'     \item Larger patch (e.g. \code{11}) -- captures more context
#'       around each source. Better for separating extended objects
#'       from point sources.
#'     \item Smaller patch (e.g. \code{5}) -- faster, focuses only on
#'       the immediate peak of the source.
#'   }
#'   Sources whose patch would fall outside the image boundary are
#'   silently dropped.
#'
#' @return A data frame with one row per successfully embedded source
#'   and the following columns:
#'   \itemize{
#'     \item \code{dim1}  -- first embedding coordinate (UMAP 1 or tSNE 1).
#'     \item \code{dim2}  -- second embedding coordinate (UMAP 2 or tSNE 2).
#'     \item \code{value} -- pixel intensity of the source from
#'       \code{spots_df}.
#'     \item \code{x}     -- row position of the source in the image.
#'     \item \code{y}     -- column position of the source in the image.
#'   }
#'   Sources at the image edge (within \code{floor(patch_size / 2)}
#'   pixels of any border) are excluded and will not appear in the output.
#'
#' @details
#' \strong{Pipeline position:}
#' This function sits between \code{\link{find_bright_spots}} and
#' \code{cluster_and_plot} in the shiny_point workflow:
#' \preformatted{
#'   img      <- read_fits("image.fits")
#'   spots    <- find_bright_spots(img)
#'   embed_df <- embed_bright_spots(img, spots)
#'   result   <- cluster_and_plot(img, embed_df)
#' }
#'
#' \strong{How the embedding works:}
#' \enumerate{
#'   \item A \code{patch_size x patch_size} pixel patch is cut from
#'     \code{img_matrix} centred on each source in \code{spots_df}.
#'   \item Each patch is flattened to a vector of length
#'     \code{patch_size^2}.
#'   \item Each vector is min-max normalised to \code{[0, 1]} so that
#'     absolute brightness does not dominate the feature space -- only
#'     the local shape of the source matters.
#'   \item All normalised vectors are stacked into a matrix and passed
#'     to UMAP or t-SNE to produce 2D coordinates.
#' }
#'
#' \strong{Required packages:}
#' \itemize{
#'   \item UMAP requires the \code{uwot} package.
#'   \item t-SNE requires the \code{Rtsne} package.
#' }
#' Install both with:
#' \preformatted{install.packages(c("uwot", "Rtsne"))}
#'
#' @examples
#' \dontrun{
#' img   <- read_fits("/Users/ananyaghosh/Downloads/crab.fits")
#' spots <- find_bright_spots(img, sigma = 3, top_n = 50)
#'
#' # Default: UMAP with 7x7 patches
#' embed_df <- embed_bright_spots(img, spots)
#'
#' # Use t-SNE instead
#' embed_df <- embed_bright_spots(img, spots, method = "tsne")
#'
#' # Larger patch to capture more context around each source
#' embed_df <- embed_bright_spots(img, spots, patch_size = 11)
#'
#' # Quick plot coloured by source intensity
#' library(ggplot2)
#' library(viridis)
#' ggplot(embed_df, aes(x = dim1, y = dim2, colour = value)) +
#'   geom_point(size = 2) +
#'   scale_colour_viridis(option = "inferno") +
#'   theme_dark() +
#'   labs(title = "UMAP of bright sources",
#'        x = "UMAP 1", y = "UMAP 2",
#'        colour = "Intensity")
#'
#' # Pass to cluster_spots for automatic clustering
#' result <- cluster_spots(img, embed_df)
#' result$image_plot
#' }
#'
#' @seealso
#' \code{\link{find_bright_spots}} to produce the \code{spots_df} input.
#' \code{\link{cluster_and_plot}} to cluster the embedding and overlay
#' results on the original image.
#'
#' @importFrom uwot umap
#' @importFrom Rtsne Rtsne
#' @export
embed_bright_spots <- function(img_matrix,
                               spots_df,
                               method     = c("umap", "tsne"),
                               patch_size = 7L) {
  
  # ---- 1. Validate inputs --------------------------------------------------
  if (!is.matrix(img_matrix)) {
    stop(
      "`img_matrix` must be a numeric matrix from read_fits().\n",
      "  Run read_fits() first:\n",
      "    img      <- read_fits('your_file.fits')\n",
      "    embed_df <- embed_bright_spots(img, spots)"
    )
  }
  
  if (!is.data.frame(spots_df)) {
    stop(
      "`spots_df` must be a data frame from find_bright_spots().\n",
      "  Run find_bright_spots() first:\n",
      "    spots    <- find_bright_spots(img)\n",
      "    embed_df <- embed_bright_spots(img, spots)"
    )
  }
  
  missing_cols <- setdiff(c("x", "y", "value"), names(spots_df))
  if (length(missing_cols) > 0L) {
    stop(
      "`spots_df` is missing required columns: ",
      paste(missing_cols, collapse = ", "), "\n",
      "  Make sure spots_df comes from find_bright_spots()."
    )
  }
  
  method <- match.arg(method)
  
  patch_size <- as.integer(patch_size)
  if (patch_size %% 2L == 0L) {
    patch_size <- patch_size + 1L
    message("`patch_size` must be odd -- adjusted to ", patch_size)
  }
  if (patch_size < 3L) {
    stop("`patch_size` must be at least 3.")
  }
  
  # ---- 2. Extract patches around each source -------------------------------
  half      <- floor(patch_size / 2L)
  nr        <- nrow(img_matrix)
  nc        <- ncol(img_matrix)
  patches   <- list()
  valid_idx <- integer(0L)
  
  for (i in seq_len(nrow(spots_df))) {
    px <- spots_df$x[i]
    py <- spots_df$y[i]
    
    # Skip sources whose patch would fall outside the image boundary
    in_bounds <- px > half &&
      px <= (nr - half) &&
      py > half &&
      py <= (nc - half)
    
    if (!in_bounds) next
    
    patch <- img_matrix[(px - half):(px + half),
                        (py - half):(py + half)]
    
    patches[[length(patches) + 1L]] <- as.vector(patch)
    valid_idx <- c(valid_idx, i)
  }
  
  n_valid <- length(valid_idx)
  
  cat("=== 2D Embedding ===\n")
  cat(sprintf("  Method           : %s\n",   toupper(method)))
  cat(sprintf("  Patch size       : %d x %d pixels\n",
              patch_size, patch_size))
  cat(sprintf("  Sources input    : %d\n",   nrow(spots_df)))
  cat(sprintf("  Valid patches    : %d\n",   n_valid))
  cat(sprintf("  Dropped (edge)   : %d\n",   nrow(spots_df) - n_valid))
  
  if (n_valid < 4L) {
    stop(
      "Too few valid patches (", n_valid, ") to embed.\n\n",
      "To fix this, try one or more of:\n",
      "  1. Increase top_n in find_bright_spots() to detect more sources.\n",
      "  2. Reduce patch_size so fewer sources are dropped at the edge.\n",
      "  3. Lower sigma in find_bright_spots() to detect fainter sources."
    )
  }
  
  # ---- 3. Stack patches into a matrix --------------------------------------
  patch_matrix <- do.call(rbind, patches)
  
  # ---- 4. Normalise each patch to [0, 1] -----------------------------------
  # Normalise by local range so brightness does not dominate the embedding --
  # only the shape of the source neighbourhood matters.
  patch_matrix <- t(apply(patch_matrix, 1L, function(r) {
    rng <- max(r) - min(r)
    if (rng == 0) r else (r - min(r)) / rng
  }))
  
  # ---- 5. Run embedding ----------------------------------------------------
  if (method == "umap") {
    
    nn  <- min(5L, n_valid - 1L)   # n_neighbors cannot exceed n_sources - 1
    cat(sprintf("  UMAP n_neighbors : %d\n", nn))
    
    embedding <- uwot::umap(patch_matrix,
                            n_neighbors  = nn,
                            min_dist     = 0.1,
                            n_components = 2L,
                            verbose      = FALSE)
    
  } else {
    
    perp <- min(10L, n_valid - 1L) # perplexity cannot exceed n_sources - 1
    cat(sprintf("  tSNE perplexity  : %d\n", perp))
    
    tsne_out  <- Rtsne::Rtsne(patch_matrix,
                              dims             = 2L,
                              perplexity       = perp,
                              verbose          = FALSE,
                              check_duplicates = FALSE)
    embedding <- tsne_out$Y
  }
  
  # ---- 6. Build output data frame ------------------------------------------
  result <- data.frame(
    dim1  = embedding[, 1L],
    dim2  = embedding[, 2L],
    value = spots_df$value[valid_idx],
    x     = spots_df$x[valid_idx],
    y     = spots_df$y[valid_idx]
  )
  
  cat(sprintf("  Embedding done.\n"))
  cat(sprintf("  dim1 range       : [%.3f, %.3f]\n",
              min(result$dim1), max(result$dim1)))
  cat(sprintf("  dim2 range       : [%.3f, %.3f]\n",
              min(result$dim2), max(result$dim2)))
  
  result
}