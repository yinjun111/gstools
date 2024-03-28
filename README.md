# Gstools

Gstools suite for gene set analyses in AWS HPC environment.


gstools
version: $version
Usage: gstools [tool] [parameters]

Parameters:

    gs-report    Automatically generate Gene Set Analyses results for rnaseq-summary

    gs-fisher    Perform Fisher's Exact Test for gene lists
    gs-fisher-summary   Summarize gs-fisher results

    gsea-gen     Run GSEA analysis
    gsea-gen-summary   Summarize GSEA analysis results

    ipa-gen      Generate files for IPA analysis (not implemented)
    ipa-summary  Summary IPA analysis results

    gs-keggpathview/gs-kegg  Generate pathway figure for KEGG using Pathview
    gs-sbgnview  Generate pathway figure for Reactome/Wikipathways using Sbgnview

    #format conversion tools
    list2matrix  Convert gene lists by columns to matrix format
    matrix2ms    Convert matrix format to Metascape format

    list-dbs     List current databases for gs-fisher
