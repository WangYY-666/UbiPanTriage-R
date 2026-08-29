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
  if (length(available_cancers(data_dir)) == 0)
    stop("No pre-computed cancer caches found under ", data_dir, ". ",
         "Run UbiPanTriage::download_extdata() once to install the pan-cancer ",
         "data (~4.2 GB), then launch the app again.")
  options(ubi.data_dir = data_dir)
  options(sass.cache = FALSE)          # avoid sass cache permission issues (避免 sass 缓存权限问题)
  options(sass.cache_dir = tempdir())  # compile bootstrap CSS under the temp dir (编译目录放到临时目录)
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

