# rnaseqc_plus


TODO: remove the current fastq_polyA trimmer,
     replace with cutadapt like so:
     cutadapt -u 1 -a A{10} -m 50 -o SeracareFusion_rep1_S1_R1_001.cutadapt.fastq.gz SeracareFusion_rep1_S1_R1_001.fastq.gz | tee cutadapt.report
     according to the Sequoia documentation.
     and simply:   pip install cutadapt
     for installation.

     add fastqc for post-trimming


