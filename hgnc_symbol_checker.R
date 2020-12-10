###############################################################################################################################
###############################################################################################################################
###  Script: hgnc_symbol_checker.R                                                                                          ###
###  Purpose: check gene symbols and retrieve hgnc id (stable identifiers)                                                  ###
###  Description: this function returns hgnc ids for protein coding genes symbols                                           ###
###  however, the input file could be modified to include other genes                                                       ###
###  Input:                                                                                                                 ###
###  1) a data frame: "gene_with_protein_product.txt" (see README) file from the hgnc website and                           ###
###  2) a vector of gene symbols to be checked                                                                              ###
###  Output: dataframe with 3 columns:                                                                                      ###
###  c1) "HGNC.ID":corresponding hgnc id ("-" , if no hgnc id was found);                                                   ###
###  c2) "Gene.Symbol": gene symbol provided;                                                                               ###
###  c3) "Type": mapping type (Approved.Symbol,Synonym.Symbol,Notfound.ProteinCoding.Symbol,...)                            ###
###############################################################################################################################
###############################################################################################################################

### required packages


if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')


### import hgnc file from hgnc site: https://www.genenames.org/ (only protein coding genes)


protein.coding.genes <- read_delim("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt",
                                   delim = "\t", col_names = TRUE)


### vector of gene symbols to be checked:

mygenes <- c("PILAR","CPAMD9","GABAT","KIAA1511","ZNF323","RNF133","GNL3B","GNL3L","PHKG2","ADGRA1-AS1","ABHD17AP6","OR6C64P",
             "TMEM31","OR2T2","CPAMD6","PZP","BCAN","BCT1","BCAT1")

###############################################################################################################################
###############################################################################################################################


## a given gene symbol can belong to one of different categories according to hgnc:
## - symbol: The official gene symbol approved by the HGNC, which is typically a short form of the gene name
## - alias symbol: Alternative names for the gene. Aliases may be from literature, 
## from other databases or may be added to represent membership of a gene group
## - previous symbol: This field displays any symbols that were previously HGNC-approved nomenclature

##############################################################################################################################

### function

hgnc.checker <- function(gene.symbols,gene.file){
  
  ## we first check if the input matches an approved gene symbol
  
  check.approved <- function(input.genes,database) {
    
    return(database %>% 
             select(hgnc_id,symbol) %>% 
             mutate_if(is.factor,as.character) %>%
             filter(symbol!="") %>%
             filter(!is.na(symbol)) %>%
             filter(symbol %in% input.genes) %>% 
             rename(Gene.Symbol = symbol,HGNC.ID = hgnc_id) %>% 
             mutate(Type = "Approved.Symbol"))
    
  }
  
  ## next we check if the input symbol corresponds to an alias /synonym gene symbol
  
  check.synonyms <- function(input.genes,database){
    
    
    return(database %>% 
             select(hgnc_id,alias_symbol) %>% 
             mutate_if(is.factor,as.character) %>%
             filter(alias_symbol!="") %>%
             filter(!is.na(alias_symbol)) %>%
             separate_rows(alias_symbol,sep="\\|") %>% 
             rename(HGNC.ID = hgnc_id) %>%
             mutate(Gene.Symbol = trimws(alias_symbol)) %>%  
             select(HGNC.ID,Gene.Symbol)  %>% 
             filter(Gene.Symbol %in% input.genes) %>%
             mutate(Type = "Synonym.Symbol"))
    
  }
  
  ## next we check if the input symbol corresponds to a previous official gene symbol
  
  check.previous <- function(input.genes,database) {
    
    return(database %>% 
             select(hgnc_id,prev_symbol) %>% 
             mutate_if(is.factor,as.character) %>%
             filter(prev_symbol!="") %>%
             filter(!is.na(prev_symbol)) %>%
             separate_rows(prev_symbol,sep="\\|") %>% 
             rename(HGNC.ID = hgnc_id) %>%
             mutate(Gene.Symbol = trimws(prev_symbol)) %>%  
             dplyr::select(HGNC.ID,Gene.Symbol)  %>% 
             filter(Gene.Symbol %in% input.genes) %>%
             mutate(Type = "Previous.Symbol"))
    
  }
  
  ## next we check if we have any duplicated symbols
  
  check.duplicates.symbol <- function(file.to.check.symbols,duplicates.symbol){
    
    if(!length(duplicates.symbol)) { return(file.to.check.symbols)}
    
    else{ 
      
      final.nodup.symbol <- file.to.check.symbols %>% 
        filter(!Gene.Symbol %in% duplicates.symbol)
      
      duplicate.symbols.df <- data.frame(HGNC.ID = rep("-",length(duplicates.symbol)), 
                                         Gene.Symbol = duplicates.symbol,
                                         Type = "Ambiguous.Symbol")
      
      final.dups.symbol <- rbind(final.nodup.symbol,duplicate.symbols.df)
      
      return(final.dups.symbol)
      
    }
  }
  
  
  ## next we check if we have any duplicated hgnc ids in the resulting file
  
  check.duplicates.id <- function(file.to.check.ids,duplicates.id){
    
    if(!length(duplicates.id)) { return(file.to.check.ids)}
    
    else{ 
      
      final.nodup.id <- file.to.check.ids %>% 
        filter(!HGNC.ID %in% duplicates.id)
      
      duplicate.ids.df <- file.to.check.ids %>% 
        filter(HGNC.ID %in% duplicates.id) %>% 
        mutate(HGNC.ID = "-",Type="Ambiguous.Symbol")
      
      final.dups.id <- rbind(final.nodup.id,duplicate.ids.df)
      
      return(final.dups.id)
      
    }
  }
  
  # we make sure we remove any leading or trailing white space in our vector with gene names
  
  genes <- trimws(gene.symbols)
  
  hgnc <- gene.file
  
  #  genes <- trimws(mygenes )
  
  #hgnc = protein.coding.genes 
  
  # we first check with gene symbols match the official gene symbol
  
  hgnc.approved <- check.approved(genes,hgnc)
  
  # we next check with gene symbols not matching the official gene symbols are synonyms or alias of an approved symbol
  
  
  hgnc.synonyms <- check.synonyms(genes[!genes %in% hgnc.approved$Gene.Symbol],hgnc)
  
  # we next check with gene symbols not matching the official gene symbols are synonyms or alias of an approved symbol
  
  
  hgnc.previous <- check.previous(genes [!genes %in% c(hgnc.approved$Gene.Symbol,hgnc.synonyms$Gene.Symbol)],hgnc)
  
  
  # we next identified with genes have not been identified in the hgnc file and create a data frame with 3 columns
  
  genes.not.found <- genes[! genes %in% c(hgnc.approved$Gene.Symbol,hgnc.synonyms$Gene.Symbol,hgnc.previous$Gene.Symbol)]
  
  
  hgnc.notfound <- data.frame(HGNC.ID = rep("-",length(genes.not.found)),
                              Gene.Symbol = genes.not.found) %>%  
    mutate(Type = "Notfound.ProteinCoding.Symbol") %>%
    mutate_if(is.factor,as.character) 
  
  # we bind all the previous data frames
  
  hgnc.all <- hgnc.approved %>%
    bind_rows(hgnc.synonyms) %>%
    bind_rows(hgnc.previous) %>%
    bind_rows(hgnc.notfound)
  
  
  # we look for any potential duplicated gene symbols or  hgnc ids in our dataset
  
  duplicates.symbol <- hgnc.all %>%
    filter(duplicated(Gene.Symbol)) %>%
    pull(Gene.Symbol)
  
  
  results.noduplicated.symbol <- check.duplicates.symbol(hgnc.all,duplicates.symbol)
  
  
  duplicates.id <- results.noduplicated.symbol %>%
    filter(duplicated(HGNC.ID)) %>%
    filter(HGNC.ID !="-") %>%
    pull(HGNC.ID)
  
  results.noduplicated.id <- check.duplicates.id(results.noduplicated.symbol,duplicates.id)
  
  # final dataframe
  
  results.final <- results.noduplicated.id
  
  return(results.final)
  
} 

################################################################################################################################
################################################################################################################################

# run the function with vector of genes provided as en example

hgnc.checker(mygenes,protein.coding.genes)

# check in the orignal hgnc input file the problems with: GNL3L & GNL3B; PZP & CPAMD6; BCT1 & BCAT1

################################################################################################################################
################################################################################################################################


