# Code to deconvolute CMC bulk tissue using psychENCODE derived signature matrix as reference

## File list

- `script_MAST.R`

    Script to run MAST and identify differentially expressed genes between inferred cell types using psychENCODE single cell data. Need to run on a computing cluster.

- `signature_genes.Rmd`

    Script to generate the signature matrix from psychENCODE single cell data.

    - `signature_genes.html`

    - `DE_gene_anno.rds`

        File summarizing differentially expressed genes identified by MAST.

- `deconvolution.Rmd`

    Script to deconvolute CMC data using ICeDT and CIBERSORT.

    - `deconvolution.html`

    - `deconvolution.pdf`

    - `CIBERSORT.Output_CMC_using_sig_genes_from_psychENCODE600_no_QN.csv`
    
        CIBERSORT estimates of cell type proportions of CMC data using the CIBERSORT online platform.

- `ExonicGeneLengths_GRCh37.RData`

    Exonic gene lengths in GRCh37. Used to calculate TPM of CMC data.

- `ExonicGeneLengths_v29.RData`

    Exonic gene lengths in GRCh38 and GENCODE v29. Used to calculate TPM of psychENCODE data.

- `tpm_signature_genes.rds`

    Signature matrix in TPM parsed from psychENCODE data.

- `CMC_ICeDT_fitw0_genes_from_psychENCODE.rds`

    ICeDT estimates of cell type proportions of CMC dataset.

- `signature_genes_brain_from_psychENCODE.txt`

    Signature matrix in TPM parsed from psychENCODE data. Input of CIBERSORT.

- `mixture_brain_from_CMC.txt`

    Bulk tissue expression in TPM parsed from CMC data. Input of CIBERSORT.
