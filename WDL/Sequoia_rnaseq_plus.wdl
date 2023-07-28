version 1.0

workflow Sequoia_rnaseq_plus_wf {

    input {
        String sample_name
        File fq_file
        File star_genome_idx_tar
        Boolean trim_polyA

        File collapsed_ref_annot_gtf

        Int min_quality = 20
        Int min_read_length = 50

        String docker="trinityctat/short_rnaseq_plus:latest"
        Int cpu = 10
        String memory="50G"
        Int preemptible = 0
        Int maxRetries = 0
        Float disk_space_multiplier = 3.0
    }



     call trim_polyA_from_fastqs_task {
           input:
           sample_name=sample_name,
           fq_file = fq_file,

           docker=docker,
           cpu=1,
           memory="4G",
           preemptible=preemptible,
           maxRetries=maxRetries,
           disk_space_multiplier=10

    }


    call quality_trim_task {
          input:
           sample_name=sample_name,
           fq_file = trim_polyA_from_fastqs_task.polyA_trimmed_fq_file,
           min_read_length = 50,
           min_quality = 20,
      
           docker=docker,
           cpu=4,
           memory="10G",
           preemptible=preemptible,
           maxRetries=maxRetries,
           disk_space_multiplier=10
    }


    call align_reads_task {
         input:
         sample_name=sample_name,
         fq_file = quality_trim_task.quality_trimmed_fq_file,
         star_genome_idx_tar = star_genome_idx_tar,

         docker=docker,
         cpu=cpu,
         memory=memory,
         preemptible=preemptible,
         maxRetries=maxRetries,
         disk_space_multiplier=disk_space_multiplier 
    }
    

    call rnaseqc_task {
        input:
           sample_name=sample_name,
           collapsed_ref_annot_gtf=collapsed_ref_annot_gtf,
           aligned_bam=align_reads_task.aligned_bam,
           aligned_bam_bai=align_reads_task.aligned_bam_bai,

           docker=docker,
           cpu=1,
           memory="10G",
           preemptible=preemptible,
           maxRetries=maxRetries,
           disk_space_multiplier=disk_space_multiplier 

     }
}


task trim_polyA_from_fastqs_task {

    input {
        String sample_name
        File fq_file

        String docker
        Int cpu
        String memory
        Int preemptible
        Int maxRetries
        Float disk_space_multiplier
    }

    Int disk_space = ceil(size(fq_file, "GB") * disk_space_multiplier)

    command <<<

      set -ex

      fastq_polyAT_trimmer.py --left_fq ~{fq_file} --out_prefix ~{sample_name} 

      gzip ~{sample_name}.polyA-trimmed.fastq

    >>>


    output {
      File polyA_trimmed_fq_file = "~{sample_name}.polyA-trimmed.fastq.gz"
    }


    runtime {
            docker: "~{docker}"
            disks: "local-disk " + disk_space + " HDD"
            memory: "~{memory}"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
    }
}



task quality_trim_task {
    input {
      String sample_name
      File fq_file
      Int min_read_length
      Int min_quality
      
      String docker
      Int cpu
      String memory
      Int preemptible
      Int maxRetries
      Float disk_space_multiplier
    }

    Int disk_space = ceil(size(fq_file, "GB") * disk_space_multiplier)

    command <<<
      set -ex

      java -jar /usr/local/bin/trimmomatic.jar SE -threads $cpu ~{fq_file} ~{sample_name}.qtrim_min~{min_quality}.fastq \
              SLIDINGWINDOW:4:~{min_quality} LEADING:~{min_quality} TRAILING:~{min_quality} MINLEN:~{min_read_length}


      gzip ~{sample_name}.qtrim_min~{min_quality}.fastq

    >>>

    output {
      File quality_trimmed_fq_file = "~{sample_name}.qtrim_min~{min_quality}.fastq.gz"
    }

    runtime {
            docker: "~{docker}"
            disks: "local-disk " + disk_space + " HDD"
            memory: "~{memory}"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
    }
}
 


task align_reads_task {
  input {
    String sample_name
    File fq_file
    File star_genome_idx_tar
    
    String docker
    Int cpu
    String memory
    Int preemptible
    Int maxRetries
    Float disk_space_multiplier
  }

  Int disk_space = ceil( (size(fq_file, "GB") + size(star_genome_idx_tar, "GB") ) * disk_space_multiplier)

  Int max_mate_dist = 100000
  
  command <<<
    set -ex

    # untar the genome lib
    tar xvf ~{star_genome_idx_tar}
    rm ~{star_genome_idx_tar}
    
    genomeDir="ctat_genome_lib_build_dir/ref_genome.fa.star.idx"


    
    readFilesCommand=""

    if [[ "~{fq_file}" == *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
    fi

    star_cmd="STAR \
            --runMode alignReads \
            --genomeDir ref_genome.fa.star.idx \
            --runThreadN $cpu \
            --readFilesIn ~{fq_file} \
            $readFilesCommand \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ~{sample_name}. \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --alignSJDBoverhangMin 10 \
            --limitBAMsortRAM 47271261705 \
            --alignInsertionFlush Right \
            --alignMatesGapMax ~{max_mate_dist} \
            --alignIntronMax  ~{max_mate_dist} "

    
    ${star_cmd}


    samtools index "~{sample_name}.Aligned.sortedByCoord.out.bam"
    
   >>>

   output {
     File aligned_bam = "~{sample_name}.Aligned.sortedByCoord.out.bam"
     File aligned_bam_bai = "~{sample_name}.Aligned.sortedByCoord.out.bam.bai"
   }

    runtime {
            docker: "~{docker}"
            disks: "local-disk " + disk_space + " HDD"
            memory: "~{memory}"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
    }
   
}



    

task rnaseqc_task {
  input {
    String sample_name
    File collapsed_ref_annot_gtf
    File aligned_bam
    File aligned_bam_bai
    
    String docker
    Int cpu
    String memory
    Int preemptible
    Int maxRetries
    Float disk_space_multiplier
  }


   Int disk_space = ceil( (size(aligned_bam, "GB") + size(collapsed_ref_annot_gtf, "GB") ) * disk_space_multiplier)

   command <<<
     set -ex

     rnaseqc ~{collapsed_ref_annot_gtf} ~{aligned_bam} . -u

     mv ~{aligned_bam}.gene_reads.gct ~{sample_name}.rnaseqc.gene_reads.gct
     mv ~{aligned_bam}.gene_fragments.gct ~{sample_name}.rnaseqc.gene_fragments.gct
     mv ~{aligned_bam}.gene_tpm.gct ~{sample_name}.rnaseqc.gene_tpm.gct
     mv ~{aligned_bam}.exon_reads.gct ~{sample_name}.rnaseqc.exon_reads.gct
     mv ~{aligned_bam}.exon_cv.tsv ~{sample_name}.rnaseqc.exon_cv.tsv
     mv ~{aligned_bam}.metrics.tsv ~{sample_name}.rnaseqc.metrics.tsv

     gzip ~{sample_name}.rnaseqc.*
     

   >>>

   output {
     File rnaseqc_gene_reads_gct = "~{sample_name}.rnaseqc.gene_reads.gct.gz"
     File rnaseqc_gene_fragments_gct = "~{sample_name}.rnaseqc.gene_fragments.gct.gz"
     File rnaseqc_gene_tpm_gct = "~{sample_name}.rnaseqc.gene_tpm.gct.gz"
     File rnaseqc_exon_reads_gct = "~{sample_name}.rnaseqc.exon_reads.gct.gz"
     File rnaseqc_exon_cv_tsv = "~{sample_name}.rnaseqc.exon_cv.tsv.gz"
     File rnaseqc_metrics_tsv = "~{sample_name}.rnaseqc.metrics.tsv.gz"
   }

    runtime {
            docker: "~{docker}"
            disks: "local-disk " + disk_space + " HDD"
            memory: "~{memory}"
            cpu: cpu
            preemptible: preemptible
            maxRetries: maxRetries
    }
}
   
