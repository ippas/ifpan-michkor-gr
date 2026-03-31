# library(ieugwasr)
# library(dplyr)

Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJtYXRldXN6emllYmE5N0BnbWFpbC5jb20iLCJpYXQiOjE3NzQwMzgxMjYsImV4cCI6MTc3NTI0NzcyNn0.Z-7jUbDu8epZTn7drhoPHJtw9W9t2aQFg5SurDxJo-wL910H-hyohAofqkIShihx3eVSb_6QbJxYhB_UvEVx7SUe0dSYQIsR98lqxRTZXlsQ3S2k7mNu-s8M_aaTjgMfAQ1lgL-xBg7lv5lLxB9BO1WhynQMiTpSDkHDyUFEoYxC72T3FmOMsU6cjmHIDU7etzcRaA3d8_0d6We3WD_xHJ69KWoZwkeDGe3_0cqszz59o3DFF4rV9nDcsBkjua8VNccMCu2ykZqf9IitcmOeRtoz73OiObUEEPQNILoqFzwekdXflHo63f60Eg4j7VJZnZJjI5wL4JGcXUHXz438qA")

ieugwasr::get_opengwas_jwt()

info %>%
  as.data.frame() %>% 
  filter(population == "European") -> ieu_open_gwas_project_EUR

ieu_open_gwas_project_EUR %>% 
  filter(sample_size > 10000) %>%
  filter(subcategory != "NA") -> ieu_open_gwas_project_EUR_filtered


ieu_open_gwas_project_EUR_filtered %>% 
  .$subcategory %>% table

ieu_open_gwas_project_EUR_filtered %>% 
  filter(subcategory == "Psychiatric / neurological") %>% 
  filter(trait == "Major Depressive Disorder")

ieu_open_gwas_project_EUR_filtered %>% 
  filter(subcategory == "Psychiatric / neurological") %>% 
  filter(grepl("bipolar disorder", trait, ignore.case = T))


files <- gwasinfo_files("ieu-b-5110")
files


download_ieugwasr_gwas_files <- function(id,
                                         output_dir = ".",
                                         force_dir = FALSE,
                                         use_wget = TRUE,
                                         download_report = FALSE,
                                         verbose = TRUE) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("Package 'ieugwasr' is required but not installed.")
  }
  
  # ---- helper do komunikatów ----
  .msg <- function(...) {
    if (isTRUE(verbose)) {
      message(...)
    }
  }
  
  # ---- sprawdzenie output_dir ----
  if (!dir.exists(output_dir)) {
    if (isTRUE(force_dir)) {
      .msg("Output directory does not exist. Creating: ", output_dir)
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    } else {
      stop("Output directory does not exist: ", output_dir,
           "\nSet force_dir = TRUE to create it.")
    }
  }
  
  # ---- pobranie listy plików ----
  .msg("Querying OpenGWAS file list for ID: ", id)
  files <- ieugwasr::gwasinfo_files(id)
  
  if (nrow(files) == 0) {
    stop("No files available for GWAS ID: ", id)
  }
  
  urls <- files[[1]]
  n_files <- length(urls)
  
  .msg("Found ", n_files, " file(s) for GWAS ID: ", id)
  if (isTRUE(verbose)) {
    for (i in seq_along(urls)) {
      message("  [", i, "/", n_files, "] ", basename(urls[i]))
    }
  }
  
  # ---- wybór plików ----
  vcf_url <- urls[grepl("\\.vcf\\.gz$", urls)]
  tbi_url <- urls[grepl("\\.vcf\\.gz\\.tbi$", urls)]
  report_url <- urls[grepl("_report\\.html$", urls)]
  
  if (length(vcf_url) == 0) {
    stop("No .vcf.gz file found for GWAS ID: ", id)
  }
  
  # zwykle jest po jednym, ale bierzemy pierwszy jeśli byłoby więcej
  vcf_url <- vcf_url[1]
  if (length(tbi_url) > 0) tbi_url <- tbi_url[1]
  if (length(report_url) > 0) report_url <- report_url[1]
  
  # ---- ścieżki docelowe ----
  vcf_dest <- file.path(output_dir, paste0(id, ".vcf.gz"))
  tbi_dest <- file.path(output_dir, paste0(id, ".vcf.gz.tbi"))
  report_dest <- file.path(output_dir, paste0(id, "_report.html"))
  
  # ---- funkcja pobierająca ----
  download_one_file <- function(url, dest, file_index, total_files) {
    .msg("Downloading file ", file_index, "/", total_files, ": ", basename(dest))
    .msg("  URL: ", url)
    .msg("  DEST: ", dest)
    
    if (isTRUE(use_wget)) {
      cmd <- sprintf("wget -O %s %s", shQuote(dest), shQuote(url))
      status <- system(cmd)
      if (!identical(status, 0L)) {
        stop("Download failed for file: ", basename(dest))
      }
    } else {
      utils::download.file(url = url, destfile = dest, mode = "wb", quiet = !verbose)
    }
    
    file_size <- if (file.exists(dest)) file.info(dest)$size else NA
    .msg("Finished: ", basename(dest),
         if (!is.na(file_size)) paste0(" (", format(file_size, big.mark = " "), " bytes)") else "")
  }
  
  # ---- lista plików do pobrania ----
  download_jobs <- list(
    list(url = vcf_url, dest = vcf_dest, type = "vcf")
  )
  
  if (length(tbi_url) > 0) {
    download_jobs[[length(download_jobs) + 1]] <- list(
      url = tbi_url,
      dest = tbi_dest,
      type = "tbi"
    )
  }
  
  if (isTRUE(download_report) && length(report_url) > 0) {
    download_jobs[[length(download_jobs) + 1]] <- list(
      url = report_url,
      dest = report_dest,
      type = "report"
    )
  }
  
  n_to_download <- length(download_jobs)
  .msg("Prepared ", n_to_download, " file(s) for download.")
  
  # ---- pobieranie ----
  for (i in seq_along(download_jobs)) {
    download_one_file(
      url = download_jobs[[i]]$url,
      dest = download_jobs[[i]]$dest,
      file_index = i,
      total_files = n_to_download
    )
  }
  
  .msg("All downloads completed for GWAS ID: ", id)
  
  return(list(
    id = id,
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    n_files_found = n_files,
    n_files_downloaded = n_to_download,
    found_files = urls,
    downloaded_files = vapply(download_jobs, function(x) x$dest, character(1)),
    vcf = if (file.exists(vcf_dest)) vcf_dest else NA_character_,
    tbi = if (file.exists(tbi_dest)) tbi_dest else NA_character_,
    report = if (file.exists(report_dest)) report_dest else NA_character_
  ))
}


res <- download_ieugwasr_gwas_files(
  id = "ieu-b-5110",
  output_dir = "/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/sandbox/ieu-b-5110",
  force_dir = TRUE,
  download_report = TRUE
)

download_ieugwasr_gwas_files <- function(id,
                                         output_dir = ".",
                                         force_dir = FALSE,
                                         download_report = FALSE,
                                         overwrite = FALSE,
                                         method = c("wget", "download.file"),
                                         verbose = TRUE) {
  method <- match.arg(method)
  
  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("Package 'ieugwasr' is required but not installed.")
  }
  
  .msg <- function(...) {
    if (isTRUE(verbose)) {
      message(...)
    }
  }
  
  # ---- sprawdzenie katalogu wyjściowego ----
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
  
  # ---- pobranie listy plików ----
  .msg("Querying OpenGWAS files for ID: ", id)
  files <- ieugwasr::gwasinfo_files(id)
  
  if (nrow(files) == 0) {
    stop("No files available for GWAS ID: ", id)
  }
  
  urls <- files[[1]]
  n_found <- length(urls)
  
  .msg("Found ", n_found, " file(s):")
  if (isTRUE(verbose)) {
    for (i in seq_along(urls)) {
      message("  [", i, "/", n_found, "] ", basename(urls[i]))
    }
  }
  
  # ---- wybór interesujących plików ----
  vcf_url <- urls[grepl("\\.vcf\\.gz$", urls)]
  tbi_url <- urls[grepl("\\.vcf\\.gz\\.tbi$", urls)]
  report_url <- urls[grepl("_report\\.html$", urls)]
  
  if (length(vcf_url) == 0) {
    stop("No .vcf.gz file found for GWAS ID: ", id)
  }
  
  vcf_url <- vcf_url[1]
  if (length(tbi_url) > 0) tbi_url <- tbi_url[1]
  if (length(report_url) > 0) report_url <- report_url[1]
  
  # ---- ścieżki docelowe ----
  vcf_dest <- file.path(output_dir, paste0(id, ".vcf.gz"))
  tbi_dest <- file.path(output_dir, paste0(id, ".vcf.gz.tbi"))
  report_dest <- file.path(output_dir, paste0(id, "_report.html"))
  
  # ---- lista zadań pobierania ----
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
  .msg("Prepared ", n_to_download, " file(s) for download.")
  
  # ---- funkcja pomocnicza do pobrania jednego pliku ----
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
      # resume + ładny progress bar
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
  
  # ---- pobieranie ----
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
  
  return(list(
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
  ))
}

res <- download_ieugwasr_gwas_files(
  id = "ieu-b-5110",
  output_dir = "/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/sandbox/ieu-b-5110_v2",
  force_dir = TRUE,
  overwrite = TRUE
)


ieu_open_gwas_project_EUR_filtered %>% 
  head %>% 
  mutate(
    clean_trait = trait %>%
      gsub("[^a-zA-Z0-9]", "_", .) %>%
      gsub("_+", "_", .) %>%
      gsub("^_|_$", "", .)
  ) %>% 
  mutate(folder_name = paste0(id, "_", clean_trait))


# ##############################################################################
library(dplyr)

# =========================
# helper: czyszczenie nazw
# =========================
clean_for_path <- function(x) {
  x %>%
    gsub("[^a-zA-Z0-9]", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}

# =========================
# helper: format czasu
# =========================
format_seconds_pretty <- function(seconds) {
  if (is.na(seconds) || is.infinite(seconds)) return(NA_character_)
  
  seconds <- round(seconds)
  h <- seconds %/% 3600
  m <- (seconds %% 3600) %/% 60
  s <- seconds %% 60
  
  sprintf("%02d:%02d:%02d", h, m, s)
}

# =========================
# katalog bazowy
# =========================
base_output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/gwas_sampleSize10000_withoutMissingSubcategory"

# =========================
# przygotowanie tabeli
# =========================
download_df <- ieu_open_gwas_project_EUR_filtered %>%
  mutate(
    clean_trait = clean_for_path(trait),
    clean_category = clean_for_path(category),
    folder_name = paste0(id, "_", clean_trait),
    output_dir = file.path(base_output_dir, clean_category, folder_name)
  )

# opcjonalnie:
# download_df <- download_df %>% slice(1:10)

n_total <- nrow(download_df)

if (n_total == 0) {
  stop("No rows to process.")
}

message("Total GWAS datasets to download: ", n_total)

# =========================
# obiekty do logowania
# =========================
results_list <- vector("list", length = n_total)
iteration_times_sec <- numeric(0)

global_start_time <- Sys.time()

# =========================
# pętla
# =========================
for (i in seq_len(n_total)) {
  row_i <- download_df[i, ]
  
  id_i <- row_i$id[[1]]
  trait_i <- row_i$trait[[1]]
  category_i <- row_i$category[[1]]
  output_dir_i <- row_i$output_dir[[1]]
  
  n_done_before <- i - 1
  n_left_before <- n_total - n_done_before
  
  avg_sec <- if (length(iteration_times_sec) > 0) mean(iteration_times_sec) else NA_real_
  eta_sec <- if (!is.na(avg_sec)) avg_sec * (n_total - i + 1) else NA_real_
  
  message("\n====================================================")
  message("Downloading GWAS ", i, "/", n_total)
  message("ID: ", id_i)
  message("Trait: ", trait_i)
  message("Category: ", category_i)
  message("Done: ", n_done_before, "/", n_total)
  message("Remaining: ", n_left_before)
  message("Average iteration time: ", format_seconds_pretty(avg_sec))
  message("Estimated time remaining: ", format_seconds_pretty(eta_sec))
  message("Output dir: ", output_dir_i)
  message("====================================================")
  
  iter_start_time <- Sys.time()
  
  status <- "OK"
  error_message <- NA_character_
  res <- NULL
  
  tryCatch({
    res <- download_ieugwasr_gwas_files(
      id = id_i,
      output_dir = output_dir_i,
      force_dir = TRUE,
      overwrite = FALSE,
      download_report = FALSE,
      method = "wget",
      verbose = TRUE
    )
  }, error = function(e) {
    status <<- "ERROR"
    error_message <<- conditionMessage(e)
  })
  
  iter_end_time <- Sys.time()
  iter_time_sec <- as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
  iteration_times_sec <- c(iteration_times_sec, iter_time_sec)
  
  avg_sec_after <- mean(iteration_times_sec)
  n_done_after <- i
  n_left_after <- n_total - i
  eta_sec_after <- avg_sec_after * n_left_after
  
  message("---- Iteration finished ----")
  message("Status: ", status)
  message("Iteration time: ", format_seconds_pretty(iter_time_sec))
  message("Average time per GWAS: ", format_seconds_pretty(avg_sec_after))
  message("Downloaded: ", n_done_after, "/", n_total)
  message("Remaining: ", n_left_after)
  message("Estimated time remaining: ", format_seconds_pretty(eta_sec_after))
  
  results_list[[i]] <- data.frame(
    index = i,
    id = id_i,
    trait = trait_i,
    category = category_i,
    clean_trait = row_i$clean_trait[[1]],
    clean_category = row_i$clean_category[[1]],
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
  
  # zapis loga po każdej iteracji
  results_df_tmp <- bind_rows(results_list)
  write.table(
    results_df_tmp,
    file = file.path(base_output_dir, "download_log.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# =========================
# końcowe podsumowanie
# =========================
results_df <- bind_rows(results_list)

total_time_sec <- as.numeric(difftime(Sys.time(), global_start_time, units = "secs"))

message("\n############################################")
message("ALL DONE")
message("Total datasets processed: ", n_total)
message("Successful: ", sum(results_df$status == "OK", na.rm = TRUE))
message("Failed: ", sum(results_df$status == "ERROR", na.rm = TRUE))
message("Total elapsed time: ", format_seconds_pretty(total_time_sec))
message("Average per GWAS: ", format_seconds_pretty(mean(results_df$iteration_time_sec, na.rm = TRUE)))
message("Log saved to: ", file.path(base_output_dir, "download_log.tsv"))
message("############################################")
  
