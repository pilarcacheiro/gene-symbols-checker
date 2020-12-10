This repo contains the script "hgnc_symbol_checker.R"
###  Purpose
To check a set of gene symbols and retrieve the corresponding hgnc id (stable identifiers).  
###  Description
This function returns hgnc ids for protein coding genes symbols. However, the input file could be modified to include other genes                                                      
###  Input:                                                                                                               
1) a data frame: "gene_with_protein_product.txt" (see README) file from the hgnc website and                         
2) a vector of gene symbols to be checked                                                                          
###  Output: dataframe with 3 columns:                                                                             
c1) "HGNC.ID":corresponding hgnc id ("-" , if no hgnc id was found);                                             
c2) "Gene.Symbol": gene symbol provided;                                                                        
c3) "Type": mapping type (Approved.Symbol,Synonym.Symbol,Notfound.ProteinCoding.Symbol,...)   
