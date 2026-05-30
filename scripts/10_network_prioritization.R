#!/usr/bin/env Rscript

# Neutral entry point retained for manuscript provenance.
# The original implementation file is kept unchanged to preserve run history.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

source(file.path(script_dir, "10_zhao_method_network_prioritization.R"), local = FALSE)
