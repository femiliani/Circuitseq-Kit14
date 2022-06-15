#It is safest to use absolute paths  
NXF_VER=21.10.6 nextflow run <path to /pipelines/CircuitSeq.nf> \
           --GPU ON \
           -c <path to /pipelines/nextflow.config> \
           -with-singularity <path to .sif file> \
           --samplesheet <path to sample_sheet.tsv> \
           --use_existing_basecalls true \
           --fast5 "" \
           --basecalling_dir <path_to_fastq_dir> \
           --base_calling_summary_file <path_to_summary.txt> \
           --barcodes /plasmidseq/barcodes/v2/ \
           --guppy_model dna_r9.4.1_450bps_sup.cfg \
           --medaka_model r941_min_sup_g507 \
           --gpu_slot cuda:0 \
           --barcode_min_score 65 \
           --quality_control_processes true \
           -resume

