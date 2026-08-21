#' Plot a FITS Image in Colour and Greyscale
#'
#' Takes the numeric matrix returned by \code{\link{read_fits}} and
#' produces two plots -- one false-colour plot using a viridis palette
#' and one greyscale plot. Both are returned as ggplot objects so you
#' can display, save, or customise them further.
#'
#' @param img_matrix numeric matrix. The image to plot, as returned by
#'   \code{\link{read_fits}}. Each cell is one pixel intensity value.
#' @param scale character. Intensity scaling applied before plotting.
#'   \itemize{
#'     \item \code{"log"} (default) -- \code{log1p} stretch. Best for
#'       most astronomical images where a bright core or star greatly
#'       outshines the surrounding emission.
#'     \item \code{"sqrt"} -- square-root stretch. Softer than log,
#'       good when bright and faint features are closer in brightness.
#'     \item \code{"linear"} -- no stretch. Use for calibrated flux
#'       maps where the absolute scale matters.
#'   }
#' @param colormap character. Any viridis palette name. Options:
#'   \itemize{
#'     \item \code{"inferno"} (default) -- black to red to yellow.
#'       Classic look for solar and high-energy images.
#'     \item \code{"magma"}   -- black to purple to white.
#'     \item \code{"plasma"}  -- purple to orange to yellow.
#'     \item \code{"viridis"} -- purple to green to yellow.
#'     \item \code{"cividis"} -- blue to yellow, colour-blind safe.
#'   }
#' @param title character. Base title shown in both plots. The suffix
#'   \code{"-- Colour"} or \code{"-- Grayscale"} is appended
#'   automatically. Default \code{"Astronomical Image"}.
#'
#' @return A named list with two ggplot objects:
#'   \itemize{
#'     \item \code{$color} -- false-colour plot using the chosen
#'       \code{colormap} on a black background.
#'     \item \code{$gray}  -- greyscale plot on a black background.
#'   }
#'   Display by typing \code{plots$color} or \code{plots$gray}.
#'   Save with \code{ggplot2::ggsave()}.
#'
#' @details
#' Before plotting, the image is shifted so its minimum value is zero
#' (\code{img - min(img)}), then the chosen \code{scale} is applied.
#' This ensures \code{log} and \code{sqrt} never receive negative input.
#'
#' Both plots use \code{theme_void()} with a black background so the
#' image fills the panel with no axis clutter, and \code{coord_fixed()}
#' keeps the pixel aspect ratio 1:1 so the image is not distorted.
#'
#' @examples
#' \dontrun{
#' # Step 1: read the image
#' img <- read_fits("/Users/ananyaghosh/Downloads/crab.fits")
#'
#' # Step 2: plot with defaults (log scale, inferno colormap)
#' plots <- plot_fits(img)
#' plots$color    # show colour plot
#' plots$gray     # show greyscale plot
#'
#' # Change scale and colormap
#' plots <- plot_fits(img,
#'                    scale    = "sqrt",
#'                    colormap = "magma",
#'                    title    = "Crab Nebula")
#' plots$color
#'
#' # Linear scaling for calibrated flux images
#' plots <- plot_fits(img, scale = "linear", title = "M101 UV")
#' plots$gray
#'
#' # Save to PNG
#' ggplot2::ggsave("crab_colour.png", plots$color,
#'                 width = 6, height = 6, dpi = 300)
#'
#' # Solar Orbiter EUI with plasma colormap
#' img   <- read_fits("/Users/ananyaghosh/Downloads/solo_L2_eui-hrieuv174-image_20231013T090050209_V01.fits")
#' plots <- plot_fits(img,
#'                    colormap = "plasma",
#'                    title    = "Solar Orbiter EUI 174 Ang")
#' plots$color
#' }
#'
#' @seealso
#' \code{\link{read_fits}} to read the FITS image matrix before plotting.
#' \code{\link{find_bright_spots}} to detect bright sources in the image.
#'
#' @import ggplot2
#' @import viridis
#' @importFrom grDevices gray.colors
#' @export

plot_fits <- function(img_matrix,
                      scale    = c("log", "sqrt", "linear"),
                      colormap = "inferno",
                      title    = "Astronomical Image") {
  
  # ---- 1. Validate inputs --------------------------------------------------
  if (!is.matrix(img_matrix)) {
    stop(
      "`img_matrix` must be a numeric matrix.\n",
      "  Load your FITS file first with read_fits():\n",
      "    img   <- read_fits('your_file.fits')\n",
      "    plots <- plot_fits(img)"
    )
  }
  
  scale <- match.arg(scale)
  
  valid_colormaps <- c("inferno", "magma", "plasma", "viridis", "cividis")
  if (!colormap %in% valid_colormaps) {
    stop(
      "`colormap` must be one of: ",
      paste(valid_colormaps, collapse = ", "), "\n",
      "  You supplied: '", colormap, "'"
    )
  }
  
  # ---- 2. Shift to zero then apply scaling ---------------------------------
  img_pos <- img_matrix - min(img_matrix)
  
  img_scaled <- switch(scale,
                       "log"    = log1p(img_pos),
                       "sqrt"   = sqrt(img_pos),
                       "linear" = img_pos
  )
  
  # ---- 3. Build long-format data frame for ggplot --------------------------
  img_df       <- expand.grid(x = seq_len(nrow(img_scaled)),
                              y = seq_len(ncol(img_scaled)))
  img_df$value <- as.vector(img_scaled)
  
  # ---- 4. Shared black background theme ------------------------------------
  dark_theme <- ggplot2::theme(
    plot.background   = ggplot2::element_rect(fill = "black", colour = NA),
    panel.background  = ggplot2::element_rect(fill = "black", colour = NA),
    plot.title        = ggplot2::element_text(colour = "white",
                                              hjust  = 0.5,
                                              size   = 13),
    legend.background = ggplot2::element_rect(fill = "black"),
    legend.text       = ggplot2::element_text(colour = "white"),
    legend.title      = ggplot2::element_text(colour = "white")
  )
  
  # ---- 5. Colour plot -------------------------------------------------------
  p_color <- ggplot2::ggplot(img_df,
                             ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    viridis::scale_fill_viridis(option = colormap,
                                name   = "Intensity") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    dark_theme +
    ggplot2::labs(title = paste(title, "\u2014 Colour"))
  
  # ---- 6. Greyscale plot ---------------------------------------------------
  p_gray <- ggplot2::ggplot(img_df,
                            ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(colours = grDevices::gray.colors(256),
                                  name    = "Intensity") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    dark_theme +
    ggplot2::labs(title = paste(title, "\u2014 Grayscale"))
  
  # ---- 7. Return both as named list ----------------------------------------
  list(color = p_color,
       gray  = p_gray)
}