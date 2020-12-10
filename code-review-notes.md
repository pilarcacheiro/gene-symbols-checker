This repository contains the script "hgnc_symbol_checker.R"
####  Purpose
To check a set of gene symbols and retrieve the corresponding hgnc id (stable identifiers). 
Merging different files through gene symbols is a common and painful task for those of us working with genetic data. Gene symbols are not stable, they often get updated / removed, so it is advisable to use a set of stable identifier, in this case HGNC IDs, since they remain stable even if a name or symbol change.
Anyone wishing to retrieve the set of HGNC IDs from a list of gene symbols could use this script. This code is tailored to one specific input file.
####  Description
This function returns hgnc ids for protein coding genes symbols. However, the input file could be modified to include other genes. 
The user only needs to provide a vector with a set of genes that need to be checked. The function calls the require hgnc file from the ftp repository (This file could also be downloaded and be imported from a local directory. A vector of gene symbols is provided as an example.
####  Input                                                                                                             
1) a data frame: "gene_with_protein_product.txt" (see README) file from the hgnc website and                         
2) a vector of gene symbols to be checked                                                                          
####  Output
A dataframe with 3 columns:                                                                             
c1) "HGNC.ID":corresponding hgnc id ("-" , if no hgnc id was found);                                             
c2) "Gene.Symbol": gene symbol provided;                                                                        
c3) "Type": mapping type (Approved.Symbol,Synonym.Symbol,Notfound.ProteinCoding.Symbol,...)   
#### How to run the code
The code can be run as follows:
```
source("hgnc_symbol_checker.R")
hgnc.checker(gene.symbols,gene.file)
```
### Expected output from the example provided in the script
|HGNC.ID|Gene.Symbol|Type|
| ------------- |------------- |------------- |
|HGNC:6722|BCAM|Approved.Symbol|
|HGNC:14725|OR2T2|Approved.Symbol|
|HGNC:8931|PHKG2|Approved.Symbol|
|HGNC:21154| RNF133|Approved.Symbol|
|HGNC:28601|TMEM31|Approved.Symbol|
|HGNC:23|GABAT|Synonym.Symbol|
|HGNC:24191|PILAR|Synonym.Symbol|
|HGNC:29299|KIAA1511|Synonym.Symbol|
|HGNC:23336|CPAMD9|Previous.Symbol|
|HGNC:14097|ZNF323|Previous.Symbol|
|-|ADGRA1-AS1|Notfound.ProteinCoding.Symbol|
|-|ABHD17AP6|Notfound.ProteinCoding.Symbol|
|-|OR6C64P|Notfound.ProteinCoding.Symbol|
|-|BCAT1|Ambiguous.Symbol|
|-|GNL3L|Ambiguous.Symbol|
|-|PZP|Ambiguous.Symbol|
|-|GNL3B|Ambiguous.Symbol|
|-|CPAMD6|Ambiguous.Symbol|
|-|BCT1|Ambiguous.Symbol|

### Some notes

* This is work in progress. This function relies on a series on assumptions that may not necessarily be true (if the gene symbol provided matches the current approved gene symbol, the function stops there. This symbol could match an old gene symbol or an alias as well; it is also sequential: approved > alias > previous).
* The output does not provide solution for ambiguous calls.

### Desirable feedback

* Mainly coding style. I am aware the code could be simplified and optimised.
* Identify potential errors
* Suggestions for improving the output file
* Suggestions for improving potential ambiguous calls
* How useful do you find this function? Is it something that you think it could be converted into a package if improved/ considering some other scenarios (non protein coding genes)?
