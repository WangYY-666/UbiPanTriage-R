# ============================================================================
# shiny.R - launch the interactive Shiny app (启动交互界面)
# ============================================================================

#' Probe an available local port (探测一个可用的本地端口)
#' @param start first port to try (起始端口)
#' @param n how many consecutive ports to probe (连续探测个数)
#' @return an integer port that is currently free, or NULL (返回可用端口号)
#' @noRd
.find_free_port <- function(start = 3839L, n = 100L) {
  for (p in seq.int(start, start + n - 1L)) {
    free <- tryCatch({
      con <- socketConnection("127.0.0.1", p, open = "r", timeout = 0.2)
      close(con)
      FALSE
    }, error = function(e) TRUE)
    if (free) return(p)
  }
  NULL
}

# Package-level flag: at most one keep-alive timer per R session (包级标记，每会话最多一个保活定时器)
.ubi_keepalive <- new.env(parent = emptyenv())

#' Refresh the modification time of every file under the session temp dir
#' (刷新会话临时目录下所有文件的修改时间)
#'
#' The modification time of an entry is what age-based \code{/tmp} cleaners
#' (\code{systemd-tmpfiles}, \code{tmpwatch}, cloud security agents) look at, so
#' touching the files keeps the whole tree out of the next cleanup pass.
#' @param d directory to refresh (要刷新的目录，默认会话临时目录)
#' @return invisible number of entries touched (被刷新的条目数)
#' @noRd
.touch_tempdir_files <- function(d = tempdir()) {
  if (!dir.exists(d)) return(invisible(0L))  # already gone: nothing to keep (目录已丢失)
  fs <- c(d, list.files(d, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, no.. = TRUE))
  now <- Sys.time()
  n <- 0L
  for (f in fs) {
    if (isTRUE(tryCatch({ Sys.setFileTime(f, now); TRUE },
                        error = function(e) FALSE, warning = function(w) FALSE)))
      n <- n + 1L
  }
  invisible(n)
}

#' Keep the R session temp dir alive while the app runs (运行期间保持临时目录存活)
#'
#' \code{bslib} compiles the app theme into a sub-directory of
#' \code{tempdir()} and serves it from there. On long-running servers the OS
#' cleans \code{/tmp} (systemd-tmpfiles / tmpwatch / cloud security agents drop
#' entries whose access, modification and status-change times are all older than
#' ~10 days); once that happens every page request fails with \dQuote{The output
#' directory '.../bslib-...' does not exist}, because the already-registered CSS
#' cannot be re-created in the same R process. Refresh the file times
#' periodically and the tree never becomes old enough to be collected. The
#' durable fix on a server is to point \code{TMPDIR} at a directory outside
#' \code{/tmp} (see \code{deploy/Dockerfile}); this timer is the in-process
#' safety net for hosts where that is not possible.
#'
#' The timer is armed once per R session and re-arms itself, so it keeps
#' touching the temp dir even when no browser is connected.
#'
#' Do NOT call this from the app's global code: \code{shiny::testServer()}
#' drains the \code{later} loop while a session is simulated, and a pending
#' timer that is hours away makes it wait for the whole delay (i.e. the test
#' suite would hang). It is called from \code{run_shiny_app()} only, which the
#' tests never run.
#' @param every_hours how often to refresh the file times (刷新间隔，小时)
#' @return invisible TRUE if a timer was started, FALSE otherwise (是否启动)
#' @noRd
.start_tempdir_keepalive <- function(every_hours = 6) {
  if (isTRUE(.ubi_keepalive$running)) return(invisible(FALSE))
  if (!requireNamespace("later", quietly = TRUE)) return(invisible(FALSE))
  touch_once <- function() {
    .touch_tempdir_files()
    later::later(touch_once, delay = every_hours * 3600)
    invisible(NULL)
  }
  .ubi_keepalive$running <- TRUE
  later::later(touch_once, delay = every_hours * 3600)
  invisible(TRUE)
}

#' Does the directory hold the web lookup library? (是否存在网页用的 per-gene 文库)
#'
#' The Shiny app reads the pre-computed per-gene library \code{lookup/}
#' (\code{pcs_meta.rds}, \code{ubi_univ.rds}, \code{genes_index.rds}, ...),
#' NOT the legacy per-cancer \code{<CANCER>_cache.rds} files. Installation
#' deliberately drops those caches to keep the image small, so the launch check
#' must test the lookup library (otherwise the app refuses to start on a
#' correctly built image).
#' 网页读取的是 \code{lookup/} 文库而非旧的按癌种缓存；镜像会删掉旧缓存以缩小体积，
#' 因此启动自检必须检查 \code{lookup/}。
#' @param data_dir candidate data directory (数据目录)
#' @return TRUE when the lookup library is present (文库是否齐全)
#' @noRd
.has_web_data <- function(data_dir) {
  ld <- file.path(data_dir, "lookup")
  dir.exists(ld) && all(file.exists(file.path(ld, c("pcs_meta.rds", "ubi_univ.rds",
                                                   "genes_index.rds"))))
}

#' Launch the Shiny web interface (启动 Shiny 交互界面)
#'
#' Starts the interactive web app with three pages: (1) ubiquitin gene-set
#' pan-cancer scoring & ranking, (2) single-gene pan-cancer query with
#' radar / correlation / boxplot / ROC / KM plots and a bilingual inference,
#' (3) custom gene-set scoring with optional user feature scores.
#' 启动交互式网页界面：泛素基因集泛癌评分排序、单基因泛癌查询（雷达图、
#' 相关性、箱线图、ROC、生存曲线与自动推论）以及自定义基因集评分。
#'
#' @param data_dir directory containing \code{pancancer_ubi_scores.rds} and
#'   the \code{lookup/} per-gene library (web data directory, 默认包内置数据)
#' @param port port for the shiny server (端口号；被占用时自动改用随机端口)
#' @param ... passed to \code{shiny::runApp}
#' @return invisible NULL; the app blocks until closed (阻塞直至关闭)
#' @examples
#' if (interactive()) run_shiny_app()
#' @export
run_shiny_app <- function(data_dir = system.file("extdata", package = "UbiPanTriage"),
                          port = 3838, ...) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("package 'shiny' is required to run the app")
  app_dir <- system.file("shiny", package = "UbiPanTriage")
  if (!nzchar(app_dir) || !dir.exists(app_dir))
    stop("shiny app not found in the installed package")
  if (!.has_web_data(data_dir) && length(available_cancers(data_dir)) == 0)
    stop("Web data not found under ", data_dir, ". ",
         "Run UbiPanTriage::download_extdata() once to install the pan-cancer ",
         "data (~4.2 GB), then launch the app again.")
  options(ubi.data_dir = data_dir)
  options(sass.cache = FALSE)          # avoid sass cache permission issues (避免 sass 缓存权限问题)
  options(sass.cache_dir = tempdir())  # compile bootstrap CSS under the temp dir (编译目录放到临时目录)
  # Keep bslib's temp CSS alive on long-running servers. Arm the timer here and
  # nowhere else (仅在启动真实服务时启用：放在 app.R 全局代码里会让 testServer() 等待定时器)
  .start_tempdir_keepalive()
  dots <- list(...)
  if (is.null(dots$launch.browser)) {
    # open the browser when possible, but never crash the app if it fails (浏览器打开失败不崩溃)
    dots$launch.browser <- function(url) {
      ok <- tryCatch({ utils::browseURL(url); TRUE }, error = function(e) FALSE)
      if (!ok) message("Could not open the browser automatically; open this URL manually: ", url)
    }
  }
  args <- c(list(app_dir, port = port), dots)
  started <- tryCatch({
    do.call(shiny::runApp, args)
    TRUE
  }, error = function(e) {
    if (grepl("address already in use|failed to create server|in use|bind|port",
              conditionMessage(e), ignore.case = TRUE)) {
      FALSE
    } else {
      stop(e)
    }
  })
  if (!started) {
    # Port occupied -> pick a free port and restart (端口被占用时自动改用空闲端口)
    port2 <- .find_free_port()
    if (is.null(port2))
      stop("no free port found; please close other Shiny sessions and retry ",
           "(未找到空闲端口，请先关闭其他 Shiny 会话后重试)")
    message("Port ", port, " is in use; switching to port ", port2,
            " (端口 ", port, " 被占用，已自动改用端口 ", port2, ")")
    args$port <- port2
    do.call(shiny::runApp, args)
  }
  invisible(NULL)
}

