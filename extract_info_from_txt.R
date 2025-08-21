# PDF_to_excel_Zuzka/extract_info_from_txt.R
# extract data from all all txt files in selected folder
library("writexl")

files <- list.files()

results <- do.call(rbind, lapply(files, get_data_from_txt))

get_data_from_txt <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  lines <- lines[lines != ""]
  
  values <- list(
    "Total TMB" = NA,
    "Coding Region Size Mb" = NA,
    "Passing Eligible Variants" = NA,
    "Usable MSI Sites" = NA,
    "Unstable Microsatellite Sites" = NA,
    "Percent Unstable Sites" = NA
  )
  
  for (line in lines) {
    if (grepl("Total TMB", line)) {
      values[["Total TMB"]] <- sub(".*:\t", "", line)
    } else if (grepl("Coding Region Size", line)) {
      values[["Coding Region Size Mb"]] <- sub(".*:\t", "", line)
    } else if (grepl("Number of Passing", line)) {
      values[["Passing Eligible Variants"]] <- sub(".*:\t", "", line)
    } else if (grepl("Usable MSI Sites", line)) {
      values[["Usable MSI Sites"]] <- sub(".*:\t", "", line)
    } else if (grepl("Total Microsatellite Sites Unstable", line)) {
      values[["Unstable Microsatellite Sites"]] <- sub(".*:\t", "", line)
    } else if (grepl("Percent Unstable Sites", line)) {
      values[["Percent Unstable Sites"]] <- sub(".*:\t", "", line)
    }
  }

  result <- as.data.frame(values, stringsAsFactors = FALSE, check.names = FALSE)
  result[["Patient ID"]] <- sub("_Biomarker.*", "", basename(file_path))
  result <- result[, c("Patient ID", setdiff(names(result), "Patient ID"))]
  return(result)
}

write_xlsx(results, "biomarker_summary.xlsx")
