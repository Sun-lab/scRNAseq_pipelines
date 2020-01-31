=================================
Human MTG
=================================

RNA sequencing data of single nuclei isolated from human middle temporal gyrus cortical area (MTG).
The data set includes 15928 single nuclei collected from six cortical layers of MTG.
The sequencing results were aligned to exons and introns in the GRCh38.p2 reference genome using the STAR algorithm, and aggregated counts at the gene level were calculated.
For more details, please see the Documentation tab in the Cell Types web application.

Gene expression data matrices
	human_MTG_2018-06-14_exon-matrix.csv
		Contains the (row, column) matrix of read counts obtained for each (gene, nucleus) based on alignment to exons.
		The first row contains the unique identifiers of the samples (nuclei)
		The first column contains the unique identifiers of the genes
	human_MTG_2018-06-14_intron-matrix.csv
		Contains the (row, column) matrix of read counts obtained for each (gene, nucleus) based on alignment to introns.
		The first row contains the unique identifiers of the samples (nuclei)
		The first column contains the unique identifiers of the genes
		

Sample information (human_MTG_2018-06-14_samples-columns.csv)
	sample_name
		Unique sample identifier
	sample_id
		Unique numeric sample identifier
	sample_type
		Indicates that the samples are nuclei
	organism
		Donor species
	donor_id
		Unique donor identifier
	sex
		Sex of donor
	age_days
		Age of donor in days (age in years multiplied by 365)
	brain_hemisphere
		Brain hemisphere of dissected cells
	brain_region
		Brain region targeted for sampling
	brain_subregion
		Brain subregion targeted for sampling	
	facs_date
		The date on which cells were collected by FACS. All dissections were performed on the same date as the FACS date for each sample
	facs_container
		FACS container unique identifier
	facs_sort_criteria
		Gating criteria used to select cells for sorting
	rna_amplification_set
		Amplification plate
	library_prep_set
		Library plate
	library_prep_avg_size_bp
		Average bp size of Library (Fragment Analyzer™ Automated CE)
	seq_name
		Unique identifier for sequencing
	seq_tube
		Sequencing lane
	seq_batch
		Sequencing batch
	total_reads
		Total number of sequencing reads
	percent_exon_reads
		% unique genomic reads aligned to exons (STAR)
	percent_intron_reads
		% unique genomic reads aligned to introns (STAR)
	percent_intergenic_reads
		% unique genomic reads aligned to intergenic sequence (STAR) [corrected 10/2018]
	percent_rrna_reads
		% total reads aligned to rRNA (STAR)
	percent_unique_reads
		% total reads aligned to genome and unique (STAR)
	percent_synth_reads
		% total reads aligned to ERCC synthetic mRNA (STAR)      
	percent_ecoli_reads
		% total reads aligned to E. coli (STAR)
	percent_aligned_reads_total
		% total reads aligned (STAR)
	complexity_cg
		Dinucleotide odds ratio (PRINSEQ)
	genes_detected_cpm_criterion
		# of genes with CPM values greater than 0, intron and exon counts included
	genes_detected_fpkm_criterion
		# of genes with FPKM values greater than 0, only exon counts included
	class
		Broad cell type class (e.g. GABAergic, Glutamatergic or Non-neuronal), or "No Class" if identified as low-quality or a doublet, or if originally assigned to an outlier or donor-specific cluster 
	cluster
		Cell type cluster name

	
Gene information (human_MTG_2018-06-14_genes-rows.csv)
	gene
		Gene symbol
	chromosome
		Chromosome location of gene
	entrez_id
		NCBI Entrez ID
	gene_name
		Gene name
	mouse_homologenes
		Mouse ortholog
