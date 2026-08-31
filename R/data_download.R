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

#' Read the release parts manifest shipped with the package (读取随包的发布分卷清单)
#' @noRd
.read_release_manifest <- function() {
  f <- system.file("extdata_release_manifest.csv", package = "UbiPanTriage")
  if (!file.exists(f))
    stop("release manifest not found in the installed package; ",
         "please reinstall UbiPanTriage (remotes::install_github(\"",
         .data_repo(), "\"))")
  utils::read.csv(f, stringsAsFactors = FALSE)
}

#' Candidate download URLs for one data part (单个分片的候选下载地址)
#'
#' Primary source is the \code{data} branch served by raw.githubusercontent.com;
#' when it fails, user-supplied \code{mirror} bases and the built-in
#' mainland-China proxies are tried in order (主源为 raw.githubusercontent.com
#' 的 data 分支；失败后依次尝试用户镜像与国内常用代理)。
#' @noRd
.part_urls <- function(repo, ref, release, part, mirror) {
  if (!is.null(release))
    return(paste0("https://github.com/", repo, "/releases/download/",
                  release, "/", part))
  primary <- paste0("https://raw.githubusercontent.com/", repo, "/", ref)
  builtin <- c(
    paste0("https://gh-proxy.com/https://raw.githubusercontent.com/", repo, "/", ref),
    paste0("https://ghproxy.net/https://raw.githubusercontent.com/", repo, "/", ref),
    paste0("https://mirror.ghproxy.com/https://raw.githubusercontent.com/", repo, "/", ref)
  )
  bases <- sub("/+$", "", unique(c(primary, mirror, builtin)))
  paste0(bases[nzchar(bases)], "/extdata_parts/", part)
}

#' Persistent cache directory for downloaded parts (分片本地缓存目录)
#'
#' Falls back to the session temp directory when the user cache dir cannot be
#' created or written (e.g. a locked-down \code{R_USER_CACHE_DIR}), so a cache
#' location problem never blocks the download.
#' @noRd
.parts_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("UbiPanTriage", "cache"), "parts")
  ok <- tryCatch({
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.access(d, 2L) == 0L
  }, error = function(e) FALSE)
  if (isTRUE(ok)) return(d)
  alt <- file.path(tempdir(), "UbiPanTriage_parts")
  dir.create(alt, recursive = TRUE, showWarnings = FALSE)
  warning("R user cache dir is not writable; using temporary cache: ", alt,
          call. = FALSE)
  alt
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
#' @param local_dir optional path to a complete \code{extdata} folder;
#'   when given, the data are copied from this folder instead of being
#'   downloaded (offline / local install, 本地 extdata 目录，提供时直接复制安装)。
#' @param release optional GitHub Release tag; when given, the data parts
#'   are downloaded from \code{https://github.com/<repo>/releases/download/<release>/}
#'   using \code{extdata_release_manifest.csv} instead of the \code{data} branch
#'   (备选发布源：GitHub Release 网页上传的分卷，data 分支不可用时使用)。
#' @param mirror optional character vector of base URLs that serve the same
#'   \code{extdata_parts/} directory; tried automatically after the primary
#'   source fails (镜像源；主源失败后自动依次尝试，每个元素应能访问
#'   \code{<mirror>/extdata_parts/<part>})。
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
                             local_dir = NULL,
                             release = NULL,
                             mirror = NULL,
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

  # ---- offline / local install: copy from an existing extdata folder ----
  if (!is.null(local_dir)) {
    if (!dir.exists(local_dir))
      stop("local_dir does not exist: ", local_dir)
    n_src <- length(available_cancers(local_dir))
    if (n_src < 33)
      stop("local_dir does not contain all 33 cancer caches (found ", n_src,
           "); expected <CANCER>_cache.rds files under ", local_dir)
    if (!quiet) message("Installing pre-computed data from local_dir: ", local_dir)
    dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
    items <- list.files(local_dir, full.names = TRUE, all.files = FALSE,
                        recursive = FALSE)
    ok <- file.copy(items, ext_dir, recursive = TRUE, overwrite = TRUE)
    if (!all(ok))
      stop("failed to copy some files from local_dir: ", local_dir)
    if (length(available_cancers(ext_dir)) < 33)
      stop("local copy is incomplete; please re-run download_extdata()")
    if (!quiet) message("Done. ", length(available_cancers(ext_dir)),
                        " cancer caches installed under ", ext_dir)
    return(invisible(ext_dir))
  }

  man   <- if (is.null(release)) .read_manifest() else .read_release_manifest()
  parts <- man$part
  tmp   <- .parts_cache_dir()
  dl    <- file.path(tmp, paste0(".dl_", Sys.getenv("USERNAME", "user")))
  dir.create(dl, showWarnings = FALSE, recursive = TRUE)
  tgz   <- file.path(tmp, paste0("UbiPanTriage_extdata_", ref, ".tar.gz"))

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout_sec)

  src_label <- if (is.null(release)) paste0(repo, "@", ref) else paste0(repo, " release:", release)
  if (!quiet) message("Downloading ", length(parts), " data parts (",
                      sprintf("%.1f", sum(man$bytes) / 1e9), " GB) from ",
                      src_label, " ...")
  ok_base <- NULL
  for (i in seq_along(parts)) {
    pf <- file.path(dl, parts[i])
    need <- !file.exists(pf) ||
      !identical(as.character(tools::md5sum(pf)), man$md5[i])
    if (need) {
      urls <- .part_urls(repo, ref, release, parts[i], mirror)
      if (!is.null(ok_base))
        urls <- unique(c(paste0(ok_base, "/extdata_parts/", parts[i]), urls))
      if (!quiet) message(sprintf("  [%d/%d] %s (%.0f MB)", i, length(parts),
                                  parts[i], man$bytes[i] / 1e6))
      ok <- FALSE
      errs <- character(0)
      for (u in urls) {
        if (!quiet && length(urls) > 1) message("    trying: ", u)
        ok <- tryCatch({
          utils::download.file(u, pf, mode = "wb", quiet = TRUE)
          TRUE
        }, error = function(e) {
          errs <<- c(errs, conditionMessage(e))
          FALSE
        })
        if (isTRUE(ok)) {
          ok_base <- sub("/extdata_parts/.*$", "", u)
          if (!quiet && length(urls) > 1)
            message("    ok (", sub("/extdata_parts/.*$", "", u), ")")
          break
        }
      }
      if (!isTRUE(ok))
        stop("failed to download data part '", parts[i],
             "' from any source.\n",
             "Tried ", length(urls), " URL(s):\n  ",
             paste(urls, collapse = "\n  "), "\n",
             "Last error: ", tail(errs, 1), "\n",
             "Fix: (1) use a VPN/proxy or retry later;\n",
             "(2) pass a mirror that hosts extdata_parts/, e.g.\n",
             "    download_extdata(mirror = \"https://gh-proxy.com/https://raw.githubusercontent.com/",
             repo, "/", ref, "\")\n",
             "(3) or install offline from a complete extdata folder with\n",
             "    download_extdata(local_dir = \"path/to/extdata\").")
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
