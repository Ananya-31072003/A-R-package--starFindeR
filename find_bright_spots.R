#' Detect Bright Sources in an Astronomical Image
#'
#' Identifies bright sources (stars, galaxies, nebulae, cosmic rays) in
#' an image matrix using robust sigma-clipping. The detection threshold
#' is computed from the median and median absolute deviation (MAD) of all
#' pixel values, making it resistant to extreme outliers. Returns the
#' top brightest spatially-detected pixels as a data frame.
#'
#' @param img_matrix numeric matrix. The image to analyse, as returned
#'   by \code{\link{read_fits}}.
#' @param sigma numeric. Detection threshold in units of MAD above the
#'   median background. Higher values find fewer but more certain sources:
#'   \itemize{
#'     \item \code{3} (default) -- standard detection threshold. Finds
#'       all sources more than 3 MAD above the background. Good starting
#'       point for most images.
#'     \item \code{5} -- stricter. Only the brightest, most obvious sources.
#'       Use when \code{sigma = 3} returns too many false detections.
#'     \item \code{2} -- looser. Picks up faint sources too.
#'       Use when \code{sigma = 3} finds no sources (\code{top_n = 0}).
#'   }
#' @param top_n integer. Maximum number of sources to return, taken as
#'   the \code{top_n} brightest pixels above the threshold.
#'   Default \code{50}.
#'
#' @return A data frame with \code{top_n} rows (or fewer if less than
#'   \code{top_n} pixels exceed the threshold), sorted by brightness
#'   descending (brightest source first). Columns:
#'   \itemize{
#'     \item \code{x}     -- row index (pixel x-coordinate in the matrix).
#'     \item \code{y}     -- column index (pixel y-coordinate in the matrix).
#'     \item \code{value} -- pixel intensity at that position.
#'   }
#'   Also prints a summary to the console showing the threshold, total
#'   bright pixels found, and the position and value of the brightest source.
#'
#' @details
#' \strong{Why median and MAD instead of mean and SD?}
#' Astronomical images typically have a small number of extremely bright
#' pixels (bright stars, cosmic rays, galactic cores). The standard
#' deviation is inflated by these outliers, raising the threshold so
#' high that real faint sources are missed. The median and MAD are
#' robust to outliers and give a stable estimate of the background level
#' regardless of how bright the brightest source is.
#'
#' \strong{Threshold formula:}
#' \preformatted{threshold = median(img) + sigma * mad(img)}
#'
#' \strong{Tip -- if too many sources are found near one bright region:}
#' Use \code{\link{find_dense_regions}} after this function to understand
#' where the sources are spatially clustered, or increase \code{sigma}.
#'
#' @examples
#' \dontrun{
#' img <- read_fits("/Users/ananyaghosh/Downloads/crab.fits")
#'
#' # Default: 3-sigma threshold, return top 50 brightest pixels
#' spots <- find_bright_spots(img)
#' head(spots, 10)   # top 10 brightest sources
#' nrow(spots)       # how many were returned
#'
#' # Stricter threshold -- only the very brightest sources
#' spots <- find_bright_spots(img, sigma = 5, top_n = 20)
#'
#' # More sources from a faint image
#' spots <- find_bright_spots(img, sigma = 2, top_n = 200)
#'
#' # Check the brightest source position
#' spots$x[1]      # row in the matrix
#' spots$y[1]      # column in the matrix
#' spots$value[1]  # intensity value
#'
#' # Pass to find_dense_regions to compute spatial density
#' density_result <- find_dense_regions(img, spots)
#' }
#'
#' @seealso
#' \code{\link{read_fits}} to load the image matrix.
#' \code{\link{find_dense_regions}} to compute the spatial density of
#' the detected sources and find the peak density region.
#' \code{\link{plot_fits}} to visualise the image.
#'
#' @importFrom stats median mad
#' @export
find_bright_spots <- function(img_matrix,
                              sigma = 3,
                              top_n = 50L) {
  
  # ---- 1. Validate inputs --------------------------------------------------
  if (!is.matrix(img_matrix)) {
    stop(
      "`img_matrix` must be a numeric matrix.\n",
      "  Load your FITS file first:\n",
      "    img   <- read_fits('your_file.fits')\n",
      "    spots <- find_bright_spots(img)"
    )
  }
  
  if (!is.numeric(sigma) || length(sigma) != 1L || sigma <= 0) {
    stop("`sigma` must be a single positive number. Example: sigma = 3")
  }
  
  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1L) {
    stop("`top_n` must be a single positive integer. Example: top_n = 50")
  }
  
  top_n <- as.integer(top_n)
  
  # ---- 2. Compute robust background and threshold --------------------------
  bg_median <- stats::median(img_matrix)
  bg_mad    <- stats::mad(img_matrix)
  threshold <- bg_median + sigma * bg_mad
  
  # ---- 3. Find all pixels above threshold ----------------------------------
  bright_idx <- which(img_matrix >= threshold, arr.ind = TRUE)
  
  if (nrow(bright_idx) == 0L) {
    stop(
      "No bright sources found at sigma = ", sigma, ".\n\n",
      "To fix this:\n",
      "  Try a lower sigma value, for example:\n",
      "    find_bright_spots(img, sigma = 2)"
    )
  }
  
  # ---- 4. Build data frame and sort by brightness --------------------------
  bright_df <- data.frame(
    x     = bright_idx[, 1L],
    y     = bright_idx[, 2L],
    value = img_matrix[bright_idx]
  )
  
  bright_df <- bright_df[order(bright_df$value, decreasing = TRUE), ]
  rownames(bright_df) <- NULL
  
  # ---- 5. Return only the top N brightest ----------------------------------
  top_spots <- head(bright_df, top_n)
  
  # ---- 6. Print summary to console -----------------------------------------
  cat("=== Bright Source Detection ===\n")
  cat(sprintf("  Background median   : %.4g\n",   bg_median))
  cat(sprintf("  Background MAD      : %.4g\n",   bg_mad))
  cat(sprintf("  Threshold (%g-sigma) : %.4g\n",  sigma, threshold))
  cat(sprintf("  Total pixels found  : %d\n",     nrow(bright_df)))
  cat(sprintf("  Returning top       : %d spots\n", nrow(top_spots)))
  cat(sprintf("  Brightest value     : %.4g\n",   top_spots$value[1L]))
  cat(sprintf("  Brightest position  : (%d, %d)\n",
              top_spots$x[1L], top_spots$y[1L]))
  cat(sprintf("  Faintest returned   : %.4g\n",
              top_spots$value[nrow(top_spots)]))
  
  top_spots
}


#' Compute Spatial Density of Bright Sources
#'
#' Estimates the 2D spatial density of detected bright sources using
#' kernel density estimation (KDE) and identifies the image region where
#' sources are most concentrated. Useful for finding the brightest
#' structural region of a galaxy, nebula, or star cluster.
#'
#' @param img_matrix numeric matrix. The image, as returned by
#'   \code{\link{read_fits}}. Used only for its dimensions to set the
#'   KDE boundary.
#' @param spots_df data frame. Detected bright sources, as returned by
#'   \code{\link{find_bright_spots}}. Must contain columns \code{x}
#'   and \code{y}.
#' @param bandwidth numeric. KDE smoothing bandwidth in pixels. Controls
#'   how broadly each source contributes to the density estimate:
#'   \itemize{
#'     \item \code{50} (default) -- good for images where sources span
#'       a large area of the frame.
#'     \item Smaller values (e.g. \code{20}) -- tighter, more localised
#'       density peaks. Use when sources are well separated.
#'     \item Larger values (e.g. \code{100}) -- broader, smoother
#'       density. Use when sources are sparse.
#'   }
#'
#' @return A named list with three elements:
#'   \itemize{
#'     \item \code{$kde}    -- the full KDE object from
#'       \code{\link[MASS]{kde2d}}, containing \code{$x}, \code{$y},
#'       and \code{$z} (the density surface). Can be used for custom
#'       plotting with \code{image()} or \code{contour()}.
#'     \item \code{$peak_x} -- x coordinate (row) of the peak density
#'       location in the image.
#'     \item \code{$peak_y} -- y coordinate (column) of the peak density
#'       location in the image.
#'   }
#'   Also prints the peak density coordinates to the console.
#'
#' @details
#' Kernel density estimation places a smooth Gaussian kernel over each
#' detected source position and sums them. The result is a continuous
#' density surface showing where sources are most concentrated.
#' The peak of this surface is the single pixel location that has the
#' highest local concentration of bright sources.
#'
#' The KDE is computed on a 200 x 200 grid spanning the full image
#' extent. The peak coordinates are therefore approximate to within
#' \code{max(nrow, ncol) / 200} pixels.
#'
#' @examples
#' \dontrun{
#' img    <- read_fits("/Users/ananyaghosh/Downloads/crab.fits")
#' spots  <- find_bright_spots(img, sigma = 3, top_n = 50)
#'
#' # Default bandwidth
#' density_result <- find_dense_regions(img, spots)
#'
#' # Access results
#' density_result$peak_x   # x coordinate of densest region
#' density_result$peak_y   # y coordinate of densest region
#' density_result$kde      # full KDE object for custom plotting
#'
#' # Tighter bandwidth for a compact cluster
#' density_result <- find_dense_regions(img, spots, bandwidth = 20)
#'
#' # The KDE z matrix can be plotted directly
#' image(density_result$kde, main = "Source Density")
#' }
#'
#' @seealso
#' \code{\link{find_bright_spots}} to detect the sources passed to
#' this function.
#' \code{embed_spots} to embed sources in 2D for clustering.
#'
#' @importFrom MASS kde2d
#' @export
find_dense_regions <- function(img_matrix,
                               spots_df,
                               bandwidth = 50) {
  
  # ---- 1. Validate inputs --------------------------------------------------
  if (!is.matrix(img_matrix)) {
    stop("`img_matrix` must be a numeric matrix from read_fits().")
  }
  
  if (!is.data.frame(spots_df)) {
    stop(
      "`spots_df` must be a data frame from find_bright_spots().\n",
      "  Run find_bright_spots() first:\n",
      "    spots          <- find_bright_spots(img)\n",
      "    density_result <- find_dense_regions(img, spots)"
    )
  }
  
  missing_cols <- setdiff(c("x", "y"), names(spots_df))
  if (length(missing_cols) > 0L) {
    stop(
      "`spots_df` is missing required columns: ",
      paste(missing_cols, collapse = ", "), "\n",
      "  Make sure spots_df comes from find_bright_spots()."
    )
  }
  
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || bandwidth <= 0) {
    stop("`bandwidth` must be a single positive number. Example: bandwidth = 50")
  }
  
  if (nrow(spots_df) < 2L) {
    stop(
      "Need at least 2 sources to compute density. ",
      "Only ", nrow(spots_df), " source(s) found.\n",
      "  Lower sigma in find_bright_spots() to detect more sources."
    )
  }
  
  # ---- 2. Compute 2D kernel density estimate -------------------------------
  kde <- MASS::kde2d(
    x    = spots_df$x,
    y    = spots_df$y,
    n    = 200L,
    h    = bandwidth,
    lims = c(0, nrow(img_matrix), 0, ncol(img_matrix))
  )
  
  # ---- 3. Find the peak density location -----------------------------------
  peak_idx <- which(kde$z == max(kde$z), arr.ind = TRUE)
  peak_x   <- kde$x[peak_idx[1L, 1L]]
  peak_y   <- kde$y[peak_idx[1L, 2L]]
  
  # ---- 4. Print summary to console -----------------------------------------
  cat("=== Spatial Density Summary ===\n")
  cat(sprintf("  Sources used        : %d\n",         nrow(spots_df)))
  cat(sprintf("  Bandwidth           : %g pixels\n",  bandwidth))
  cat(sprintf("  KDE grid size       : 200 x 200\n"))
  cat(sprintf("  Peak density at     : (%.0f, %.0f)\n", peak_x, peak_y))
  cat(sprintf("  Peak as fraction    : (%.2f, %.2f) of image size\n",
              peak_x / nrow(img_matrix),
              peak_y / ncol(img_matrix)))
  
  list(kde    = kde,
       peak_x = peak_x,
       peak_y = peak_y)
}