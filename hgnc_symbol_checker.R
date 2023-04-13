################################################################################
## Script: hgnc_symbol_checker.R
################################################################################

# required packages -------------------------------------------------------

if (!require("magrittr")) install.packages("magrittr")
library("magrittr")

if (!require("dplyr")) install.packages("dplyr")
library("dplyr")

if (!require("tidyr")) install.packages("tidyr")
library("tidyr")

if (!require("logr")) install.packages("logr")
library("logr")

# create log file ---------------------------------------------------------

## option to generate a log file to have a record of the 
## date the HGNC file was accessed through the FTP
## create a directory for log files

tmp <- file.path("./data/logs", "hgnc_checker.log")

lf <- log_open(tmp)

# hgnc gene file ----------------------------------------------------------

# If the file doesn't already exist, we will read it via FTP
filename <- "gene_with_protein_product.txt"
if (!file.exists(filename)) {
  filename <- paste("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types",
                    filename,
                    sep = "/")
}

protein_coding_genes <- readr::read_delim(filename,
                                          delim = "\t",
                                          col_names = TRUE)

## a given gene symbol can belong to one of different categories according to
## hgnc:
##   - symbol: The official gene symbol approved by the HGNC, which is
##             typically a short form of the gene name
##   - alias symbol: Alternative names for the gene. Aliases may be from
##                   literature, from other databases or may be added to
##                   represent membership of a gene group
##   - previous symbol: This field displays any symbols that were previously
##                      HGNC-approved nomenclature



# hgnc checker function ---------------------------------------------------


hgnc.checker <- function(gene_symbols, gene_file) {
  ## Purpose:
  ##   check gene symbols and retrieve hgnc id (stable identifiers)
  ##
  ## Description:
  ##   this function returns hgnc ids for protein coding genes symbols
  ##   however, the input file could be modified to include other genes
  ##
  ## Input:
  ##   gene_symbols: a vector of gene symbols to be checked
  ##   gene.file: a data frame 'gene_with_protein_product.txt' (see README)
  ##              file from the hgnc website
  ##
  ## Output:
  ##   A dataframe with 3 columns:
  ##     c1) 'hgnc_id': corresponding hgnc id ('-' , if no hgnc id was found)
  ##     c2) 'gene_symbol': gene symbol provided
  ##     c3) 'type': mapping type (approved_symbol,
  ##                               synonym_symbol,
  ##                               notfound_proteinpoding_symbol,...)
  
  
  check.approved <- function(input_genes, database) {
    ## A function to check if the input matches an approved gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, symbol) %>%
             mutate_if(is.factor, as.character) %>%
             filter(symbol != "") %>%
             filter(!is.na(symbol)) %>%
             filter(symbol %in%
                      input_genes) %>%
             dplyr::rename(gene_symbol = symbol) %>%
             mutate(type = "approved_symbol"))
  }
  
  
  check.synonyms <- function(input_genes, database) {
    ## A function to check if the input symbol corresponds to an
    ## alias / synonym gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, alias_symbol) %>%
             mutate_if(is.factor,
                       as.character) %>%
             filter(alias_symbol != "") %>%
             filter(!is.na(alias_symbol)) %>%
             separate_rows(alias_symbol, sep = "\\|") %>%
             mutate(gene_symbol = trimws(alias_symbol)) %>%
             dplyr::select(hgnc_id, gene_symbol) %>%
             filter(gene_symbol %in% input_genes) %>%
             mutate(type = "synonym_symbol"))
  }
  
  
  check.previous <- function(input_genes, database) {
    ## A function to check if the input symbol corresponds to a previous
    ## official gene symbol
    
    return(database %>%
             dplyr::select(hgnc_id, prev_symbol) %>%
             mutate_if(is.factor,
                       as.character) %>%
             filter(prev_symbol != "") %>%
             filter(!is.na(prev_symbol)) %>%
             separate_rows(prev_symbol, sep = "\\|") %>%
             mutate(gene_symbol = trimws(prev_symbol)) %>%
             dplyr::select(hgnc_id,
                           gene_symbol) %>%
             filter(gene_symbol %in% input_genes) %>%
             mutate(type = "previous_symbol"))
  }
  
  
  check.duplicates.symbol <- function(file_to_check_symbols,
                                      duplicates_symbol) {
    ## A function to check if we have any duplicated symbols
    
    if (!length(duplicates_symbol)) {
      return(file_to_check_symbols)
      
    } else {
      
      final_nodup_symbol <- file_to_check_symbols %>%
        filter(!gene_symbol %in%
                 duplicates_symbol)
      
      duplicate_symbols_df <- data.frame(hgnc_id = rep("-",
                                                       length(duplicates_symbol)),
                                         gene_symbol = duplicates_symbol,
                                         type = "ambiguous_symbol")
      
      final_dups_symbol <- rbind(final_nodup_symbol, duplicate_symbols_df)
      
      return(final_dups_symbol)
      
    }
  }
  
  
  check.duplicates.id <- function(file_to_check_ids, duplicates_id) {
    ## A function to check if we have any duplicated hgnc ids in the
    ## resulting file
    
    if (!length(duplicates_id)) {
      return(file_to_check_ids)
    } else {
      
      final_nodup_id <- file_to_check_ids %>%
        filter(!hgnc_id %in% duplicates_id)
      
      duplicate_ids_df <- file_to_check_ids %>%
        filter(hgnc_id %in% duplicates_id) %>%
        mutate(hgnc_id = "-", type = "ambiguous_symbol")
      
      final_dups_id <- rbind(final_nodup_id, duplicate_ids_df)
      
      return(final_dups_id)
      
    }
  }
  
  # we make sure we remove any leading or trailing white space in our vector
  # with gene names
  
  genes <- trimws(gene_symbols)
  
  hgnc <- gene_file
  
  # we first check with gene symbols match the official gene symbol
  
  hgnc_approved <- check.approved(genes, hgnc)
  
  
  # we next check with gene symbols not matching the official gene symbols are
  # synonyms or alias of an approved symbol
  
  hgnc_synonyms <- check.synonyms(genes[!genes %in% hgnc_approved$gene_symbol],
                                  hgnc)
  
  
  # we next check with gene symbols not matching the official gene symbols are
  # synonyms or alias of an approved symbol
  
  hgnc_previous <- check.previous(genes[!genes %in% c(hgnc_approved$gene_symbol,
                                                      hgnc_synonyms$gene_symbol)],
                                  hgnc)
  
  
  # we next identified with genes have not been identified in the hgnc file
  # and create a data frame with 3 columns
  
  genes_not_found <- genes[!genes %in% c(hgnc_approved$gene_symbol,
                                         hgnc_synonyms$gene_symbol,
                                         hgnc_previous$gene_symbol)]
  
  
  hgnc_notfound <- data.frame(hgnc_id = rep("-", length(genes_not_found)),
                              gene_symbol = genes_not_found) %>%
    mutate(type = "notfound_proteincoding_symbol") %>%
    mutate_if(is.factor, as.character)
  
  
  # we bind all the previous data frames
  
  hgnc_all <- hgnc_approved %>%
    bind_rows(hgnc_synonyms) %>%
    bind_rows(hgnc_previous) %>%
    bind_rows(hgnc_notfound)
  
  
  # we look for any potential duplicated gene symbols or hgnc ids in our
  # dataset
  
  duplicates_symbol <- hgnc_all %>%
    filter(duplicated(gene_symbol)) %>%
    pull(gene_symbol)
  
  
  results_noduplicated_symbol <- check.duplicates.symbol(hgnc_all,
                                                         duplicates_symbol)
  
  
  duplicates_id <- results_noduplicated_symbol %>%
    filter(duplicated(hgnc_id)) %>%
    filter(hgnc_id != "-") %>%
    pull(hgnc_id)
  
  results_noduplicated_id <- check.duplicates.id(results_noduplicated_symbol,
                                                 duplicates_id)
  
  
  # final dataframe
  
  results_final <- results_noduplicated_id
  
  return(results_final)
  
  # log file ----------------------------------------------------------------
  
  
  log_print("hgnc checker")
  
  log_close()
}
