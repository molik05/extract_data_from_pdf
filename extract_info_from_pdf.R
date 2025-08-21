# PDF_to_excel_Zuzka/extract_info_from_pdf.R
# extract data from specific section in all pdf reports in selected folder

library(pdftools)
library(dplyr)
library(tidyr)
library(stringr)
library(tools)
library(purrr)

# Extract text from PDF file
pdf_file <- pdf_data("2_24-0151-DG-16977.pdf")[1:2]
pdf_file <- pdf_data("1_24-0151-DG-16976.pdf")[1:2]
pdf_file <- pdf_data("data/245_37_24-0151-DG-17383.pdf")[1:2]

vysledek <- extract_data_from_pdf(pdf_file)

pdf_folder <- getwd()

pdf_files <- list.files(path= pdf_folder, full.names = TRUE)

#extract all results
all_results <- map_dfr(pdf_files, function(pdf_path) {
  patient_id <- file_path_sans_ext(basename(pdf_path))
  message("Processing",patient_id)
  result <- extract_data_from_pdf(pdf_data(pdf_path)[1:2])
  
  if (nrow(result) == 0) {
    message("  -> No data extracted, skipping.")
    return(NULL)
  }
  result$Patient_ID <- patient_id
  result <- result %>% select(Patient_ID, everything())
  
  return(result)
})

#extract data from pdf
extract_data_from_pdf <- function(pdf_file){
  columns <- setNames(vector("list", 4), c("IA", "IB", "IIC", "IID"))
  position_data <- setNames(vector("list", 4), c("IA", "IB", "IIC", "IID"))
  result <- data.frame(Category = character(),
                       Gene = character(),
                       Variant = character(),
                       stringsAsFactors = FALSE)
  
  for (page in seq_along(pdf_file)){
    page_data <- pdf_file[[page]]
    if (page == 1){
      from_y <- page_data %>%
        filter(str_detect(page_data$text, regex("IID", ignore_case = TRUE))) %>%
        pull(y) %>% unique()
      to_y_page <- page_data %>%
        filter(str_detect(page_data$text, regex("Page", ignore_case = TRUE))) %>%
        pull(y) %>% unique()
      to_y_clin <- page_data %>%
        filter(str_detect(page_data$text, regex("CLINICALLY", ignore_case = TRUE))) %>%
        pull(y) %>% unique()
      to_y <- min(to_y_clin,to_y_page)
      filtered_data <- page_data %>%
        filter(y > from_y & y < to_y)
    }
    else{
      from_y <- page_data %>%
        filter(str_detect(page_data$text, regex("Final", ignore_case = TRUE))) %>%
        pull(y) %>% unique()
      to_y_clin <- page_data %>%
        filter(str_detect(page_data$text, regex("CLINICALLY", ignore_case = TRUE))) %>%
        pull(y) %>% unique()
      filtered_data <- page_data %>%
        filter(y > from_y & y < to_y)
    } 
    
    columns$IA  <- filtered_data %>% filter(x >= 22  & x < 163) %>% arrange(y) %>% pull(text)
    position_data$IA <- filtered_data %>% filter(x >= 22 & x < 163) %>% arrange(y)
    columns$IB  <- filtered_data %>% filter(x >= 163 & x < 304) %>% arrange(y) %>% pull(text)
    position_data$IB <- filtered_data %>% filter(x >= 163 & x < 304) %>% arrange(y)
    columns$IIC <- filtered_data %>% filter(x >= 304 & x < 446) %>% arrange(y) %>% pull(text)
    position_data$IIC <- filtered_data %>% filter(x >= 304 & x < 446) %>% arrange(y)
    columns$IID <- filtered_data %>% filter(x >= 446) %>% arrange(y) %>% pull(text)
    position_data$IID <- filtered_data %>% filter(x >= 446) %>% arrange(y)
    
    to_remove <- c("No", "variants", "reported.")
    for (col in names(columns)) {
      keep_idx <- !columns[[col]] %in% to_remove
      columns[[col]] <- columns[[col]][keep_idx]
      position_data[[col]] <- position_data[[col]][keep_idx, ]
      vec <- columns[[col]]
      n <- length(vec)
      if (n >= 3){
        columns[[col]] <- vec[1:(n - 3)]
        position_data[[col]] <- position_data[[col]][1:(n - 3), ]  
      }
    }
    
    
    for (col in seq_along(columns)) {
      for (item in seq_along(columns[[col]])){
        x_value <- position_data[[col]]$x[item]
        y_value <- position_data[[col]]$y[item]
        if (x_value %in% c(22, 163, 304, 446)){
          gene <- columns[[col]][item]
          result <- rbind(result, data.frame(
            Category = names(columns)[col],
            Gene = gene,
            Variant = y_value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    for (col in names(columns)) {
      keep_idx <- !position_data[[col]]$x %in% c(22, 163, 304, 446)
      columns[[col]] <- columns[[col]][keep_idx]
      position_data[[col]] <- position_data[[col]][keep_idx, ]
    }
  
    for (col in seq_along(columns)) {
      gene_list_by_category <- split(result$Gene, result$Category)
      genes <- gene_list_by_category[[names(columns)[col]]]
      gene_pos_list_by_category <- split(result$Variant, result$Category)
      genes_positions <- gene_pos_list_by_category[[names(columns)[col]]]
      variants <- c()
      id <- 1
      for (item in seq_along(columns[[col]])){
        y_value <- position_data[[col]]$y[item]
        if (id == length(genes)){
          if (item == length(columns[[col]])){
            variants <- paste0(variants,columns[[col]][item]," ")
            result$Variant[result$Category == names(columns)[col] & result$Gene == genes[id]] <- variants 
          }
          else{
            variants <- paste0(variants,columns[[col]][item]," ")
          }
        }
        else if (y_value < (as.numeric(genes_positions[id+1])-3)){
          variants <- paste0(variants,columns[[col]][item]," ")
        }
        else {
          result$Variant[result$Category == names(columns)[col] & result$Gene == genes[id]] <- variants
          id <- id+1
          variants <- paste0(columns[[col]][item]," ")
        }
      }
    }
  result$Variant <- gsub("c\\. ", "c.", result$Variant)
  result <- result[order(result$Category), ]
  return(result)
  }
}


#extract the most common variants
all_results_2 <- all_results %>% mutate(simplified_variants = str_extract(
  Variant,"p\\.[A-Za-z0-9*_/-]+|Copy number gain in [A-Za-z0-9]+"))

top_20_df <- as.data.frame(sort(table(simplified_variants), decreasing = TRUE)[1:10])

variant_genes <- all_results_2 %>%
  filter(simplified_variants %in% top_20_df$simplified_variants) %>%
  group_by(simplified_variants) %>%
  summarise(Genes = paste(sort(unique(Gene)), collapse = ", "),
            Patients = paste(sort(unique(Patient_ID)), collapse = ", "),
            Original_Variants = paste(sort(unique(Variant)), collapse = "; "), 
            Categories = paste(sort(unique(Category)), collapse = ", "))


top_20_with_genes <- left_join(top_20_df, variant_genes, by = "simplified_variants")

top_20_full <- top_20_with_genes %>%
  mutate(
    Variants = if_else(str_starts(simplified_variants, "Copy"),
                       simplified_variants,
                       Original_Variants)
  ) %>%
  select(-simplified_variants, -Original_Variants)

top_20_full <- top_20_full %>%
  select(Freq,Variants,Genes,Patients) %>%
  rename(Counts = Freq)

filtered_df <- all_results[all_results$Gene == "KRAS",]
unique_rows <- filtered_df[!duplicated(filtered_df$Variant), ]

patients_with_variants <- all_results[all_results$Gene %in% c("KRAS", "ERBB2", "NRAS"), ]
unique_patients <- length(unique(patients_with_variants$Patient_ID))
