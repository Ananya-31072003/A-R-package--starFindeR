#' Read a FITS Image File
#'
#' Reads any astronomical FITS file and returns a 2D numeric matrix.
#' Uses Python's astropy library via reticulate to handle all FITS formats
#' including standard images, compressed images, binary tables, and 3D cubes.
#'
#' @param fits_path character. Full path to the FITS file on your computer.
#'   Example: \code{"/Users/ananyaghosh/Downloads/crab.fits"}
#' @param hdu integer or NULL. Which HDU (Header Data Unit) index to read.
#'   Uses 0-based indexing (Python convention):
#'   \itemize{
#'     \item \code{NULL} (default) -- automatically scans all HDUs and returns
#'       the first one that contains a 2D image.
#'     \item \code{0L} -- reads the primary HDU.
#'     \item \code{1L} -- reads the first extension HDU (e.g. compressed images
#'       like Solar Orbiter EUI files where the image is in COMPRESSED_IMAGE).
#'   }
#'
#' @return A numeric matrix of class \code{"double"} with dimensions
#'   \code{nrow x ncol} where each cell is the pixel intensity value
#'   (raw instrument counts or calibrated flux from the file).
#'   \itemize{
#'     \item \code{NaN} pixels are replaced with \code{0}.
#'     \item \code{Inf} pixels are replaced with the maximum finite value.
#'   }
#'
#' @details
#' \strong{Python setup required.} Run once in Terminal before using this
#' function:
#' \preformatted{pip3 install astropy numpy}
#'
#' Then set your Python path permanently by adding this line to
#' \code{~/.Rprofile}:
#' \preformatted{
#' Sys.setenv(ASTROPY_PYTHON =
#'   "/Library/Frameworks/Python.framework/Versions/3.14/bin/python3")
#' }
#'
#' \strong{Why astropy and not FITSio?}
#' The base R package \code{FITSio} cannot read tile-compressed FITS files
#' (\code{CompImageHDU}), variable-length binary table columns
#' (\code{P}/\code{Q} format), or non-zero \code{pcount} extensions --
#' all of which are common in modern telescope pipelines such as
#' Solar Orbiter and LCO. Astropy handles every one of these transparently.
#'
#' \strong{Supported FITS formats:}
#' \itemize{
#'   \item Standard \code{PrimaryHDU} images (e.g. \code{crab.fits})
#'   \item Tile-compressed \code{CompImageHDU}
#'     (e.g. Solar Orbiter EUI \code{solo_L2_eui-hrieuv174-*.fits})
#'   \item Float32, Float64, Int16, Int32, UInt8 pixel data types
#'   \item 3D data cubes -- first frame is returned automatically
#'   \item Files with \code{NaN} or \code{Inf} values -- auto-cleaned
#' }
#'
#' \strong{Python is initialised once per R session.} The first call to
#' \code{read_fits()} starts Python and imports astropy. Subsequent calls
#' in the same session reuse the existing Python process.
#'
#' @examples
#' \dontrun{
#' # Set Python path once (add to ~/.Rprofile to make it permanent)
#' Sys.setenv(ASTROPY_PYTHON =
#'   "/Library/Frameworks/Python.framework/Versions/3.14/bin/python3")
#'
#' # 1. Standard image -- auto-detects the correct HDU
#' img <- read_fits("/Users/ananyaghosh/Downloads/crab.fits")
#' dim(img)   # 1024 1024
#'
#' # 2. UV image
#' img <- read_fits("/Users/ananyaghosh/Downloads/m101_uv.fits")
#' dim(img)   # 4096 4096  (or similar)
#'
#' # 3. LCO telescope observation
#' img <- read_fits("/Users/ananyaghosh/Downloads/tfn0m419-sq32-20250128-0353-e91.fits")
#' dim(img)   # 2400 2400
#'
#' # 4. Solar Orbiter EUI compressed image (CompImageHDU in HDU 1)
#' img <- read_fits("/Users/ananyaghosh/Downloads/solo_L2_eui-hrieuv174-image_20231013T090050209_V01.fits")
#' dim(img)   # 2048 2048
#'
#' # 5. Force a specific HDU when auto-detect picks the wrong one
#' img <- read_fits("/Users/ananyaghosh/Downloads/image.fits", hdu = 1L)
#'
#' # 6. Check what is inside before reading
#' read_fits("/Users/ananyaghosh/Downloads/image.fits")
#' # prints the full HDU table so you can choose the right hdu value
#' }
#'
#' @seealso
#' \code{\link{read_fits_header}} to extract metadata (telescope name,
#' filter, exposure time, observation date) from the FITS header.
#' \code{\link{plot_fits}} to display the image in colour and greyscale.
#' \code{\link{find_bright_spots}} to detect bright sources in the image.
#'
#' @importFrom reticulate import use_python py_available
#' @export
read_fits <- function(fits_path,
                      hdu = NULL) {
  
  # ---- 1. Validate inputs --------------------------------------------------
  if (!is.character(fits_path) || length(fits_path) != 1L) {
    stop(
      "`fits_path` must be a single character string.\n",
      "  Example: read_fits('/Users/ananyaghosh/Downloads/crab.fits')"
    )
  }
  
  if (!file.exists(fits_path)) {
    stop(
      "File not found: ", fits_path, "\n",
      "  Check the path is correct and the file exists."
    )
  }
  
  if (!is.null(hdu) && !is.numeric(hdu)) {
    stop("`hdu` must be an integer (e.g. 0L or 1L) or NULL.")
  }
  
  # ---- 2. Resolve Python path ----------------------------------------------
  python_path <- Sys.getenv("ASTROPY_PYTHON", unset = "python3")
  
  # ---- 3. Initialise Python (once per session) -----------------------------
  if (!reticulate::py_available()) {
    tryCatch(
      reticulate::use_python(python_path, required = TRUE),
      error = function(e) {
        stop(
          "Cannot start Python at: ", python_path, "\n\n",
          "To fix this:\n",
          "  1. Open Terminal and run:  which python3\n",
          "  2. Copy the path it prints.\n",
          "  3. Add this to ~/.Rprofile:\n",
          "       Sys.setenv(ASTROPY_PYTHON = '/paste/your/path/here')\n",
          "  4. Restart R (Session > Restart R) and try again."
        )
      }
    )
  }
  
  # ---- 4. Import astropy ---------------------------------------------------
  fits_module <- tryCatch(
    reticulate::import("astropy.io.fits"),
    error = function(e) {
      stop(
        "astropy is not installed.\n\n",
        "To fix this:\n",
        "  Open Terminal and run:  pip3 install astropy numpy\n",
        "  Then restart R and try again."
      )
    }
  )
  
  # ---- 5. Import numpy -----------------------------------------------------
  np <- tryCatch(
    reticulate::import("numpy"),
    error = function(e) {
      stop("numpy not found. In Terminal run:  pip3 install numpy")
    }
  )
  
  # ---- 6. Open FITS file ---------------------------------------------------
  hdul <- tryCatch(
    fits_module$open(fits_path, memmap = FALSE),
    error = function(e) {
      stop("Could not open FITS file: ", fits_path, "\n", e$message)
    }
  )
  
  # Always close file on exit, even if an error occurs below
  on.exit(try(hdul$close(), silent = TRUE))
  
  # Print full HDU table so user can see what is inside the file
  cat("=== FITS File Structure ===\n")
  hdul$info()
  cat("\n")
  
  # ---- 7. Locate image data ------------------------------------------------
  img_data  <- NULL
  found_hdu <- NA_integer_
  
  if (!is.null(hdu)) {
    
    # User told us exactly which HDU to use
    img_data <- tryCatch(
      hdul[[as.integer(hdu)]]$data,
      error = function(e) {
        stop(
          "Cannot read HDU ", hdu, ": ", e$message, "\n",
          "Check the structure printed above and try a different hdu value."
        )
      }
    )
    
    if (is.null(img_data)) {
      stop(
        "HDU ", hdu, " exists but contains no data.\n",
        "Check the structure printed above and try a different hdu value."
      )
    }
    found_hdu <- as.integer(hdu)
    
  } else {
    
    # Auto-scan: loop through every HDU and take the first 2D image
    n_hdu <- length(hdul)
    
    for (i in seq.int(0L, n_hdu - 1L)) {
      tryCatch({
        
        d    <- hdul[[as.integer(i)]]$data
        ndim <- length(dim(d))
        
        if (!is.null(d) && ndim == 2L) {
          cat(sprintf("  HDU %d : 2D image  (%d x %d)  -- selected\n",
                      i, dim(d)[1L], dim(d)[2L]))
          img_data  <- d
          found_hdu <- i
          break
          
        } else if (!is.null(d) && ndim == 3L) {
          cat(sprintf(
            "  HDU %d : 3D cube  (%d x %d x %d)  -- first frame selected\n",
            i, dim(d)[1L], dim(d)[2L], dim(d)[3L]))
          img_data  <- d[1L, , ]
          found_hdu <- i
          break
          
        } else {
          cat(sprintf("  HDU %d : no 2D image -- skipping\n", i))
        }
        
      }, error = function(e) {
        cat(sprintf("  HDU %d : could not read -- skipping\n", i))
      })
    }
  }
  
  # ---- 8. Error if nothing found -------------------------------------------
  if (is.null(img_data)) {
    stop(
      "No readable 2D image found in any HDU.\n\n",
      "To fix this:\n",
      "  Look at the structure printed above.\n",
      "  Find an HDU labelled IMAGE or CompImageHDU.\n",
      "  Pass its index manually, for example:\n",
      "    read_fits('your_file.fits', hdu = 1L)"
    )
  }
  
  # ---- 9. Convert to R matrix via numpy ------------------------------------
  # numpy handles all pixel types: int16, float32, float64, uint8 etc.
  mat <- tryCatch({
    arr <- np$array(img_data, dtype = np$float64)
    matrix(
      data = as.numeric(arr),
      nrow = dim(arr)[1L],
      ncol = dim(arr)[2L]
    )
  }, error = function(e) {
    # Direct fallback if numpy conversion fails
    matrix(
      data = as.numeric(img_data),
      nrow = nrow(img_data),
      ncol = ncol(img_data)
    )
  })
  
  # ---- 10. Clean non-finite values -----------------------------------------
  n_nan <- sum(is.nan(mat))
  n_inf <- sum(is.infinite(mat))
  
  if (n_nan > 0L) {
    message("  Note: ", n_nan, " NaN pixel(s) replaced with 0.")
    mat[is.nan(mat)] <- 0
  }
  if (n_inf > 0L) {
    message("  Note: ", n_inf, " Inf pixel(s) set to max finite value.")
    mat[is.infinite(mat)] <- max(mat[is.finite(mat)])
  }
  
  # ---- 11. Force double storage mode ---------------------------------------
  storage.mode(mat) <- "double"
  
  # ---- 12. Print summary ---------------------------------------------------
  cat("=== Image Loaded Successfully ===\n")
  cat(sprintf("  HDU used    : %d\n",             found_hdu))
  cat(sprintf("  Dimensions  : %d x %d pixels\n", nrow(mat), ncol(mat)))
  cat(sprintf("  Value range : %.4g  to  %.4g\n", min(mat),  max(mat)))
  cat(sprintf("  Mean value  : %.4g\n",            mean(mat)))
  cat(sprintf("  File        : %s\n",              basename(fits_path)))
  
  mat
}


#' Read FITS Header Metadata
#'
#' Returns all header keywords from a FITS HDU as a named R list and
#' prints the most common astronomical keywords to the console.
#'
#' @param fits_path character. Full path to the FITS file.
#'   Example: \code{"/Users/ananyaghosh/Downloads/crab.fits"}
#' @param hdu integer. Which HDU header to read (0-based index).
#'   Default \code{0L} -- the primary HDU which usually holds the
#'   main metadata (telescope, filter, exposure time, coordinates).
#'
#' @return A named list where every element is one FITS header keyword.
#'   The list is returned \strong{invisibly} to avoid flooding the
#'   console -- capture it with \code{hdr <- read_fits_header(...)}.
#'
#'   Commonly useful keywords:
#'   \itemize{
#'     \item \code{hdr[["OBJECT"]]}   -- target name
#'     \item \code{hdr[["TELESCOP"]]} -- telescope name
#'     \item \code{hdr[["INSTRUME"]]} -- instrument name
#'     \item \code{hdr[["DATE-OBS"]]} -- observation date and time (UTC)
#'     \item \code{hdr[["EXPTIME"]]}  -- exposure time in seconds
#'     \item \code{hdr[["FILTER"]]}   -- filter or bandpass used
#'     \item \code{hdr[["NAXIS1"]]}   -- image width in pixels
#'     \item \code{hdr[["NAXIS2"]]}   -- image height in pixels
#'     \item \code{hdr[["RA"]]}       -- right ascension of pointing
#'     \item \code{hdr[["DEC"]]}      -- declination of pointing
#'   }
#'
#' @examples
#' \dontrun{
#' # Print the key fields and capture the full list
#' hdr <- read_fits_header("/Users/ananyaghosh/Downloads/crab.fits")
#'
#' # Access individual keywords
#' hdr[["OBJECT"]]    # e.g. "Crab Nebula"
#' hdr[["EXPTIME"]]   # e.g. 120   (seconds)
#' hdr[["TELESCOP"]]  # e.g. "HST"
#' hdr[["DATE-OBS"]]  # e.g. "2023-10-13T09:00:50.209"
#' hdr[["FILTER"]]    # e.g. "V"
#'
#' # Read header from extension HDU (e.g. compressed image metadata)
#' hdr <- read_fits_header(
#'   "/Users/ananyaghosh/Downloads/solo_L2_eui-hrieuv174-image_20231013T090050209_V01.fits",
#'   hdu = 1L
#' )
#' }
#'
#' @seealso \code{\link{read_fits}} to read the image pixel data.
#'
#' @importFrom reticulate import use_python py_available py_to_r
#' @export
read_fits_header <- function(fits_path,
                             hdu = 0L) {
  
  # ---- Validate ------------------------------------------------------------
  if (!file.exists(fits_path)) {
    stop("File not found: ", fits_path)
  }
  
  # ---- Python setup --------------------------------------------------------
  python_path <- Sys.getenv("ASTROPY_PYTHON", unset = "python3")
  
  if (!reticulate::py_available()) {
    reticulate::use_python(python_path, required = TRUE)
  }
  
  fits_module <- reticulate::import("astropy.io.fits")
  
  hdul <- fits_module$open(fits_path, memmap = FALSE)
  on.exit(try(hdul$close(), silent = TRUE))
  
  header <- hdul[[as.integer(hdu)]]$header
  
  # ---- Extract every keyword -----------------------------------------------
  keys <- tryCatch(
    as.character(reticulate::py_to_r(header$keys())),
    error = function(e) character(0L)
  )
  keys <- keys[nchar(trimws(keys)) > 0L]   # drop blank keyword names
  
  header_list <- lapply(keys, function(k) {
    tryCatch(reticulate::py_to_r(header$get(k)), error = function(e) NA)
  })
  names(header_list) <- keys
  
  # ---- Print common astronomical keywords ----------------------------------
  important_keys <- c(
    "OBJECT",   "TELESCOP", "INSTRUME", "DATE-OBS",
    "EXPTIME",  "FILTER",   "WAVELNTH", "NAXIS1",
    "NAXIS2",   "BUNIT",    "OBSERVER", "RA",
    "DEC",      "EQUINOX",  "ORIGIN"
  )
  
  cat("=== FITS Header Keywords ===\n")
  found_any <- FALSE
  for (k in important_keys) {
    v <- header_list[[k]]
    if (!is.null(v) && !all(is.na(v))) {
      cat(sprintf("  %-12s : %s\n", k, as.character(v)))
      found_any <- TRUE
    }
  }
  if (!found_any) {
    cat("  (none of the standard keywords found in this HDU)\n")
  }
  cat(sprintf("\n  Total header cards : %d\n", length(keys)))
  cat(sprintf("  HDU index read     : %d\n",  as.integer(hdu)))
  
  invisible(header_list)
}
