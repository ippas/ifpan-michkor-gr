#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ieugwasr)
  library(dplyr)
})

# =========================================================
# KONFIGURACJA
# =========================================================

base_output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/gwas_sampleSize10000_withoutMissingSubcategory"

download_report <- TRUE
overwrite <- FALSE
min_sample_size <- 10000
target_population <- "European"
require_nonmissing_subcategory <- TRUE
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJtYXRldXN6emllYmE5N0BnbWFpbC5jb20iLCJpYXQiOjE3NzQwMzgxMjYsImV4cCI6MTc3NTI0NzcyNn0.Z-7jUbDu8epZTn7drhoPHJtw9W9t2aQFg5SurDxJo-wL910H-hyohAofqkIShihx3eVSb_6QbJxYhB_UvEVx7SUe0dSYQIsR98lqxRTZXlsQ3S2k7mNu-s8M_aaTjgMfAQ1lgL-xBg7lv5lLxB9BO1WhynQMiTpSDkHDyUFEoYxC72T3FmOMsU6cjmHIDU7etzcRaA3d8_0d6We3WD_xHJ69KWoZwkeDGe3_0cqszz59o3DFF4rV9nDcsBkjua8VNccMCu2ykZqf9IitcmOeRtoz73OiObUEEPQNILoqFzwekdXflHo63f60Eg4j7VJZnZJjI5wL4JGcXUHXz438qA")


print("hello")

# =========================================================
# FUNKCJE POMOCNICZE
# =========================================================

msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

clean_for_path <- function(x) {
  x <- as.character(x)
  x <- gsub("[^a-zA-Z0-9]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == "" | is.na(x)] <- "UNKNOWN"
  x
}

format_seconds_pretty <- function(seconds) {
  if (length(seconds) == 0 || is.na(seconds) || is.infinite(seconds)) {
    return(NA_character_)
  }
  
  seconds <- round(seconds)
  h <- seconds %/% 3600
  m <- (seconds %% 3600) %/% 60
  s <- seconds %% 60
  
  sprintf("%02d:%02d:%02d", h, m, s)
}

print("hello2")

save_runtime_logs <- function(log_tsv_file,
                              log_rdata_file,
                              results_list,
                              iteration_times_sec,
                              download_df,
                              global_start_time,
                              global_end_time = Sys.time(),
                              current_index = NA_integer_) {
  results_df_tmp <- dplyr::bind_rows(results_list)
  
  utils::write.table(
    results_df_tmp,
    file = log_tsv_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  save(
    results_list,
    results_df_tmp,
    iteration_times_sec,
    download_df,
    global_start_time,
    global_end_time,
    current_index,
    file = log_rdata_file
  )
}

# =========================================================
# FUNKCJA POBIERAJĄCA PLIKI DLA JEDNEGO GWAS
# =========================================================

download_ieugwasr_gwas_files <- function(id,
                                         output_dir = ".",
                                         force_dir = FALSE,
                                         download_report = FALSE,
                                         overwrite = FALSE,
                                         method = c("wget", "download.file"),
                                         verbose = TRUE) {
  method <- match.arg(method)
  
  .msg <- function(...) {
    if (isTRUE(verbose)) {
      msg(...)
    }
  }
  
  if (!dir.exists(output_dir)) {
    if (isTRUE(force_dir)) {
      .msg("Creating output directory: ", output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
      stop(
        "Output directory does not exist: ", output_dir, "\n",
        "Set force_dir = TRUE to create it."
      )
    }
  }
  
  .msg("Querying OpenGWAS files for ID: ", id)
  files <- ieugwasr::gwasinfo_files(id)
  
  if (nrow(files) == 0) {
    stop("No files available for GWAS ID: ", id)
  }
  
  urls <- files[[1]]
  n_found <- length(urls)
  
  .msg("Found ", n_found, " file(s) for ID ", id, ":")
  if (isTRUE(verbose)) {
    for (i in seq_along(urls)) {
      message("    [", i, "/", n_found, "] ", basename(urls[i]))
    }
  }
  
  vcf_url <- urls[grepl("\\.vcf\\.gz$", urls)]
  tbi_url <- urls[grepl("\\.vcf\\.gz\\.tbi$", urls)]
  report_url <- urls[grepl("_report\\.html$", urls)]
  
  if (length(vcf_url) == 0) {
    stop("No .vcf.gz file found for GWAS ID: ", id)
  }
  
  vcf_url <- vcf_url[1]
  if (length(tbi_url) > 0) tbi_url <- tbi_url[1]
  if (length(report_url) > 0) report_url <- report_url[1]
  
  vcf_dest <- file.path(output_dir, paste0(id, ".vcf.gz"))
  tbi_dest <- file.path(output_dir, paste0(id, ".vcf.gz.tbi"))
  report_dest <- file.path(output_dir, paste0(id, "_report.html"))
  
  download_jobs <- list(
    list(type = "vcf", url = vcf_url, dest = vcf_dest)
  )
  
  if (length(tbi_url) > 0) {
    download_jobs[[length(download_jobs) + 1]] <- list(
      type = "tbi",
      url = tbi_url,
      dest = tbi_dest
    )
  }
  
  if (isTRUE(download_report) && length(report_url) > 0) {
    download_jobs[[length(download_jobs) + 1]] <- list(
      type = "report",
      url = report_url,
      dest = report_dest
    )
  }
  
  n_to_download <- length(download_jobs)
  .msg("Prepared ", n_to_download, " file(s) for download for ID: ", id)
  
  download_one <- function(url, dest, idx, total, type) {
    if (file.exists(dest) && !isTRUE(overwrite)) {
      .msg("[", idx, "/", total, "] Skipping existing file: ", basename(dest))
      return(invisible(list(
        file = dest,
        downloaded = FALSE,
        exists = TRUE,
        size = file.info(dest)$size
      )))
    }
    
    .msg("[", idx, "/", total, "] Downloading ", type, ": ", basename(dest))
    
    if (identical(method, "wget")) {
      cmd <- sprintf(
        "wget --progress=bar:force -c -O %s %s",
        shQuote(dest),
        shQuote(url)
      )
      status <- system(cmd)
      if (!identical(status, 0L)) {
        stop("Download failed for file: ", basename(dest))
      }
    } else {
      utils::download.file(
        url = url,
        destfile = dest,
        mode = "wb",
        quiet = !verbose
      )
    }
    
    if (!file.exists(dest)) {
      stop("File was not created: ", dest)
    }
    
    size <- file.info(dest)$size
    .msg("Finished: ", basename(dest), " (", format(size, big.mark = " "), " bytes)")
    
    invisible(list(
      file = dest,
      downloaded = TRUE,
      exists = TRUE,
      size = size
    ))
  }
  
  results <- vector("list", length(download_jobs))
  
  for (i in seq_along(download_jobs)) {
    job <- download_jobs[[i]]
    results[[i]] <- download_one(
      url = job$url,
      dest = job$dest,
      idx = i,
      total = n_to_download,
      type = job$type
    )
  }
  
  .msg("All downloads completed for GWAS ID: ", id)
  
  downloaded_files <- vapply(results, function(x) x$file, character(1))
  file_sizes <- vapply(results, function(x) x$size, numeric(1))
  
  list(
    id = id,
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    n_files_found = n_found,
    n_files_prepared = n_to_download,
    found_files = urls,
    downloaded_files = downloaded_files,
    file_sizes = file_sizes,
    vcf = if (file.exists(vcf_dest)) vcf_dest else NA_character_,
    tbi = if (file.exists(tbi_dest)) tbi_dest else NA_character_,
    report = if (file.exists(report_dest)) report_dest else NA_character_
  )
}

# =========================================================
# ŚCIEŻKI LOGÓW
# =========================================================

safe_dir_create(base_output_dir)

log_tsv_file <- file.path(base_output_dir, "download_log.tsv")
log_rdata_file <- file.path(base_output_dir, "download_log.RData")
final_rdata_file <- file.path(base_output_dir, "final_download_state.RData")
metadata_tsv_file <- file.path(base_output_dir, "metadata_filtered.tsv")
metadata_rdata_file <- file.path(base_output_dir, "metadata_filtered.RData")
console_log_file <- file.path(base_output_dir, "console.log")

# =========================================================
# LOGOWANIE KONSOLI DO PLIKU
# =========================================================

log_con <- file(console_log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

# =========================================================
# START
# =========================================================

msg("Starting OpenGWAS download pipeline")
msg("Base output directory: ", base_output_dir)
msg("download_report = ", download_report)
msg("overwrite = ", overwrite)
msg("min_sample_size = ", min_sample_size)
msg("target_population = ", target_population)
msg("require_nonmissing_subcategory = ", require_nonmissing_subcategory)

token <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (!nzchar(token)) {
  stop("Environment variable OPENGWAS_JWT is not set.")
}

jwt_value <- ieugwasr::get_opengwas_jwt()
if (!nzchar(jwt_value)) {
  stop("OpenGWAS JWT is not visible to ieugwasr.")
}

msg("Token detected successfully")

wget_path <- Sys.which("wget")
if (!nzchar(wget_path)) {
  stop("wget is not available in PATH.")
}
msg("wget detected: ", wget_path)

global_start_time <- Sys.time()

# =========================================================
# POBIERANIE METADANYCH
# =========================================================

msg("Downloading metadata with gwasinfo() ... this may take some time")
info <- ieugwasr::gwasinfo()
info_df <- info %>% as.data.frame()

msg("Initial number of datasets: ", nrow(info_df))

# =========================================================
# FILTROWANIE
# =========================================================

download_df <- info_df %>%
  filter(population == target_population) %>%
  filter(!is.na(sample_size)) %>%
  filter(sample_size > min_sample_size)

if (isTRUE(require_nonmissing_subcategory)) {
  download_df <- download_df %>%
    filter(!is.na(subcategory)) %>%
    filter(subcategory != "NA") %>%
    filter(subcategory != "")
}

download_df <- download_df %>%
  mutate(
    clean_trait = clean_for_path(trait),
    folder_name = paste0(id, "_", clean_trait),
    output_dir = file.path(base_output_dir, folder_name)
  )

n_total <- nrow(download_df)

if (n_total == 0) {
  stop("No rows to process after filtering.")
}

msg("Datasets after filtering: ", n_total)

# =========================================================
# ZAPIS METADATA
# =========================================================

save(
  info_df,
  download_df,
  file = metadata_rdata_file
)

utils::write.table(
  download_df,
  file = metadata_tsv_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

msg("Saved metadata:")
msg("  ", metadata_tsv_file)
msg("  ", metadata_rdata_file)

# =========================================================
# PĘTLA POBIERANIA
# =========================================================

results_list <- vector("list", length = n_total)
iteration_times_sec <- numeric(0)

for (i in seq_len(n_total)) {
  row_i <- download_df[i, ]
  
  id_i <- row_i$id[[1]]
  trait_i <- row_i$trait[[1]]
  category_i <- row_i$category[[1]]
  subcategory_i <- row_i$subcategory[[1]]
  output_dir_i <- row_i$output_dir[[1]]
  
  n_done_before <- i - 1
  n_left_before <- n_total - n_done_before
  
  avg_sec_before <- if (length(iteration_times_sec) > 0) mean(iteration_times_sec) else NA_real_
  eta_sec_before <- if (!is.na(avg_sec_before)) avg_sec_before * (n_total - i + 1) else NA_real_
  
  msg("====================================================")
  msg("Downloading GWAS ", i, "/", n_total)
  msg("ID: ", id_i)
  msg("Trait: ", trait_i)
  msg("Category: ", category_i)
  msg("Subcategory: ", subcategory_i)
  msg("Done before iteration: ", n_done_before, "/", n_total)
  msg("Remaining before iteration: ", n_left_before)
  msg("Average iteration time so far: ", format_seconds_pretty(avg_sec_before))
  msg("Estimated remaining time: ", format_seconds_pretty(eta_sec_before))
  msg("Output dir: ", output_dir_i)
  msg("====================================================")
  
  iter_start_time <- Sys.time()
  
  status <- "OK"
  error_message <- NA_character_
  res <- NULL
  
  tryCatch({
    res <- download_ieugwasr_gwas_files(
      id = id_i,
      output_dir = output_dir_i,
      force_dir = TRUE,
      download_report = download_report,
      overwrite = overwrite,
      method = "wget",
      verbose = TRUE
    )
  }, error = function(e) {
    status <<- "ERROR"
    error_message <<- conditionMessage(e)
    msg("ERROR for ", id_i, ": ", error_message)
  })
  
  iter_end_time <- Sys.time()
  iter_time_sec <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
  iteration_times_sec <- c(iteration_times_sec, iter_time_sec)
  
  avg_sec_after <- mean(iteration_times_sec)
  n_done_after <- i
  n_left_after <- n_total - i
  eta_sec_after <- avg_sec_after * n_left_after
  
  msg("---- Iteration finished ----")
  msg("Status: ", status)
  msg("Iteration time: ", format_seconds_pretty(iter_time_sec))
  msg("Average time per GWAS: ", format_seconds_pretty(avg_sec_after))
  msg("Processed: ", n_done_after, "/", n_total)
  msg("Remaining: ", n_left_after)
  msg("Estimated time remaining: ", format_seconds_pretty(eta_sec_after))
  
  results_list[[i]] <- data.frame(
    index = i,
    id = id_i,
    trait = trait_i,
    category = category_i,
    subcategory = subcategory_i,
    clean_trait = row_i$clean_trait[[1]],
    folder_name = row_i$folder_name[[1]],
    output_dir = output_dir_i,
    status = status,
    error_message = error_message,
    iteration_time_sec = iter_time_sec,
    avg_iteration_time_sec = avg_sec_after,
    eta_remaining_sec = eta_sec_after,
    timestamp_start = as.character(iter_start_time),
    timestamp_end = as.character(iter_end_time),
    stringsAsFactors = FALSE
  )
  
  save_runtime_logs(
    log_tsv_file = log_tsv_file,
    log_rdata_file = log_rdata_file,
    results_list = results_list,
    iteration_times_sec = iteration_times_sec,
    download_df = download_df,
    global_start_time = global_start_time,
    global_end_time = Sys.time(),
    current_index = i
  )
}

# =========================================================
# KOŃCOWE PODSUMOWANIE
# =========================================================

results_df <- dplyr::bind_rows(results_list)
global_end_time <- Sys.time()
total_time_sec <- as.numeric(difftime(global_end_time, global_start_time, units = "secs"))

save(
  info_df,
  download_df,
  results_list,
  results_df,
  iteration_times_sec,
  global_start_time,
  global_end_time,
  total_time_sec,
  file = final_rdata_file
)

msg("############################################")
msg("ALL DONE")
msg("Total datasets processed: ", n_total)
msg("Successful: ", sum(results_df$status == "OK", na.rm = TRUE))
msg("Failed: ", sum(results_df$status == "ERROR", na.rm = TRUE))
msg("Total elapsed time: ", format_seconds_pretty(total_time_sec))
msg("Average per GWAS: ", format_seconds_pretty(mean(results_df$iteration_time_sec, na.rm = TRUE)))
msg("Logs saved to:")
msg("  ", log_tsv_file)
msg("  ", log_rdata_file)
msg("  ", final_rdata_file)
msg("  ", metadata_tsv_file)
msg("  ", metadata_rdata_file)
msg("  ", console_log_file)
msg("############################################")