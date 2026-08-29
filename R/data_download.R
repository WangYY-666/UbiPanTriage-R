# ============================================================================
# data_download.R - fetch & install the bundled pan-cancer data
# (下载并安装内置的 33 癌种预计算数据)
#
# The source package on GitHub stays small; the pre-computed 33-cancer caches
# (~4.2 GB unpacked) are shipped separately as split archives in the `data`
# branch of the GitHub repository. These helpers download the parts, verify
# checksums, reassemble the tarball and unpack it into the installed package.
# ============================================================================

#' Check whether the bundled pan-cancer data are installed (检查内置数据是否就绪)
#'
#' Verifies that the pre-computed cancer caches are present under
#' \code{data_dir}. The caches are required by \code{\link{load_cancer_data}},
#' most scoring functions and the Shiny app; install them once with
#' \code{\link{download_extdata}}.
#' 检查预计算癌种缓存是否已就位；缺失时用 \code{download_extdata()} 安装。
#'
#' @param data_dir directory containing the \code{<CANCER>_cache.rds} files
#'   (缓存目录，默认包内置 extdata)。
#' @param min_cancers minimum number of cancer caches expected (最小预期癌种数)。
#' @return \code{TRUE} if enough caches are found, otherwise \code{FALSE}
#'   (invisibly); a status message is printed.
#' @examples
#' check_extdata()
#' @export
check_extdata <- function(data_dir = system.file("extdata", package = "UbiPanTriage"),
                          min_cancers = 33) {
  n <- length(available_cancers(data_dir))
  ok  <- n >= min_cancers
  if (ok) {
    message("UbiPanTriage data OK: ", n, " cancer caches found under ", data_dir)
  } else {
    message("Bundled pan-cancer data are missing or incomplete: found ", n,
            " of ", min_cancers, " cancer caches.",
            " Run UbiPanTriage::download_extdata() once to install them.")
  }
  invisible(ok)
}

#' Default GitHub repository holding the split data parts (数据分片所在仓库)
#' @noRd
.data_repo <- function() "WangYY-666/UbiPanTriage-R"

#' Default git ref (branch/tag) holding the data parts (数据分片所在分支)
#' @noRd
.data_ref <- function() "data"

#' Read the parts manifest shipped with the package (读取随包的分片清单)
#' @noRd
.read_manifest <- function() {
  f <- system.file("extdata_manifest.csv", package = "UbiPanTriage")
  if (!file.exists(f))
    stop("parts manifest not found in the installed package; ",
         "please reinstall UbiPanTriage (remotes::install_github(\"",
         .data_repo(), "\"))")
  utils::read.csv(f, stringsAsFactors = FALSE)
}

#' Persistent cache directory for downloaded parts (分片本地缓存目录)
#' @noRd
.parts_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("UbiPanTriage", "cache"), "parts")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Download and install the bundled pan-cancer data (下载并安装内置泛癌数据)
#'
#' The GitHub repository only ships the package source. The pre-computed
#' 33-cancer caches (expression / clinical / differential / immune /
#' metabolism / per-gene statistics, about 4.2 GB unpacked) are distributed as
#' split archives in the \code{data} branch of the repository. This function
#' downloads every part once (with checksum verification and resumable parts),
#' reassembles the tarball and unpacks it into \code{dest_dir}, after which
#' \code{\link{available_cancers}()}, \code{\link{load_cancer_data}()} and
#' \code{\link{run_shiny_app}()} work out of the box.
#'
#' GitHub 仓库只保留包源码，33 癌种预计算缓存（解压后约 4.2 GB）以分片形式
#' 存放在仓库的 \code{data} 分支。运行本函数一次即可完成下载、校验、拼接与
#' 安装；之后 \code{available_cancers()} / \code{load_cancer_data()} /
#' \code{run_shiny_app()} 均可直接使用。
#'
#' @param repo GitHub repository in \code{"owner/name"} form (GitHub 仓库名)。
#' @param ref git ref (branch/tag) holding the data parts (数据所在分支/标签)。
#' @param dest_dir directory into which the \code{extdata} folder will be
#'   unpacked; defaults to the installed package directory (数据安装目录，
#'   默认已安装包目录)。
#' @param quiet suppress progress messages (是否静默)。
#' @param timeout_sec per-part download timeout in seconds (单分片下载超时秒数)。
#' @return invisible path to the installed \code{extdata} directory.
#' @examples
#' \dontrun{
#' download_extdata()          # default repo + installed package
#' download_extdata(dest_dir = tempdir())   # install to a custom location
#' }
#' @export
download_extdata <- function(repo = .data_repo(),
                             ref = .data_ref(),
                             dest_dir = system.file(package = "UbiPanTriage"),
                             quiet = FALSE,
                             timeout_sec = 3600) {
  ext_dir <- file.path(dest_dir, "extdata")
  if (length(available_cancers(ext_dir)) >= 33) {
    if (!quiet) message("Bundled pan-cancer data already installed under ", ext_dir)
    return(invisible(ext_dir))
  }
  if (!dir.exists(dest_dir))
    stop("dest_dir does not exist: ", dest_dir)
  if (file.access(dest_dir, 2) != 0)
    stop("dest_dir is not writable: ", dest_dir, "\n",
         "If the package was installed into a system library, reinstall it into ",
         "your personal library (e.g. remotes::install_github(...)) or pass ",
         "dest_dir = \"", normalizePath(tempdir(), winslash = "/"), "\".")

  man   <- .read_manifest()
  parts <- man$part
  tmp   <- .parts_cache_dir()
  dl    <- file.path(tmp, paste0(".dl_", Sys.getenv("USERNAME", "user")))
  dir.create(dl, showWarnings = FALSE, recursive = TRUE)
  tgz   <- file.path(tmp, paste0("UbiPanTriage_extdata_", ref, ".tar.gz"))

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout_sec)

  if (!quiet) message("Downloading ", length(parts), " data parts (",
                      sprintf("%.1f", sum(man$bytes) / 1e9), " GB) from ",
                      repo, "@", ref, " ...")
  for (i in seq_along(parts)) {
    pf <- file.path(dl, parts[i])
    need <- !file.exists(pf) ||
      !identical(as.character(tools::md5sum(pf)), man$md5[i])
    if (need) {
      url <- file.path("https://raw.githubusercontent.com", repo, ref,
                       "extdata_parts", parts[i])
      if (!quiet) message(sprintf("  [%d/%d] %s (%.0f MB)", i, length(parts),
                                  parts[i], man$bytes[i] / 1e6))
      utils::download.file(url, pf, mode = "wb", quiet = TRUE)
    }
    got <- as.character(tools::md5sum(pf))
    if (!identical(got, man$md5[i]))
      stop("checksum mismatch for ", parts[i], " (got ", got, ")")
  }
  if (!quiet) message("Reassembling archive ...")
  con <- file(tgz, "wb")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(parts)) {
    con2 <- file(file.path(dl, parts[i]), "rb")
    repeat {
      b <- readBin(con2, raw(), n = 8 * 1024 * 1024)
      if (length(b) == 0) break
      writeBin(b, con)
    }
    close(con2)
  }
  close(con)
  on.exit(NULL, add = FALSE)

  utils::untar(tgz, exdir = dest_dir, tar = "internal")
  if (!dir.exists(ext_dir))
    stop("unpacking failed: ", ext_dir, " not created")
  if (length(available_cancers(ext_dir)) < 33)
    stop("unpacked data are incomplete; please re-run download_extdata()")
  if (!quiet) message("Done. ", length(available_cancers(ext_dir)),
                      " cancer caches installed under ", ext_dir)
  invisible(ext_dir)
}

