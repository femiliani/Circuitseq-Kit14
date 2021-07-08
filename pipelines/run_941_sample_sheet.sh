
./nextflow PlasmidSeq/pipelines/plasmid_assembly_sample_sheet.nf \
	   --samplesheet /analysis/2021_05_06_nanopore_pipeline/assess_assembly_101/fex_101_BC_to_ref_index.tsv \
	   --barcodes /analysis/2021_05_06_nanopore_pipeline/2021_05_25_nanopore_24_plasmids/guppy_plasmid_barcodes/ \
	   --fast5 /analysis/2021_05_06_nanopore_pipeline/fastq5_101/ \
	   --guppy_model dna_r9.4.1_450bps_sup.cfg \
	   --medaka_model r941_min_sup_g507 \
	   --gpu_slot cuda:0 \
	   --tn5proj /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa.prj \
	   --bcmat /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/bc.mat \
	   -with-report report.html \
	   --tn5ref /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa \
	   --barcode_min_score 40 \
	   -resume
