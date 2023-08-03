version 1.0

workflow Pbio_rnaseqc_plus_wf {

    input {
        String sample_name
        File fq_file
        File mm2_genome_idx
        File mm2_splice_bed
        File collapsed_ref_annot_gtf


        String docker = "trinityctat/rnaseqc_plus:latest"
        Int cpu = 10
        String memory="50G"
        Int preemptible = 0
        Int maxRetries = 0
        Float disk_space_multiplier = 3.0
    }


    call align_reads_mm2_task {
         input:
           sample_name=sample_name,
           fq_file = fq_file,

           mm2_genome_idx=mm2_genome_idx,
           mm2_splice_bed=mm2_splice_bed,

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
           aligned_bam=align_reads_mm2_task.aligned_bam,
           aligned_bam_bai=align_reads_mm2_task.aligned_bam_bai,

           docker=docker,
           cpu=1,
           memory="10G",
           preemptible=preemptible,
           maxRetries=maxRetries,
           disk_space_multiplier=disk_space_multiplier 

     }
}



task align_reads_mm2_task {
  input {
    String sample_name
    File fq_file

    File mm2_genome_idx
    File mm2_splice_bed
    
    String docker
    Int cpu
    String memory
    Int preemptible
    Int maxRetries
    Float disk_space_multiplier
  }

  Int disk_space = ceil( (size(fq_file, "GB") + size(mm2_genome_idx, "GB") ) * disk_space_multiplier)

  
  command <<<
    set -ex


    minimap2 --junc-bed ~{mm2_splice_bed} -ax splice -u b -t ~{cpu} {mm2_genome_idx} ~{fq_file} > mm2.sam

    samtools view -Sb -o mm2.unsorted.bam mm2.sam

    rm mm2.sam # free up disk

    samtools sort -@~{cpu} -o ~{sample_name}.mm2.bam mm2.unsorted.bam

    samtools index ~{sample_name}.mm2.bam
    
   >>>

   output {
     File aligned_bam = "~{sample_name}.mm2.bam"
     File aligned_bam_bai = "~{sample_name}.mm2.bam.bai"
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


     ln -s ~{aligned_bam} ~{sample_name}.bam
     ln -s ~{aligned_bam_bai} ~{sample_name}.bam.bai

     rnaseqc ~{collapsed_ref_annot_gtf} ~{sample_name}.bam . -u

     
     mv ~{sample_name}.bam.gene_reads.gct ~{sample_name}.rnaseqc.gene_reads.gct
     mv ~{sample_name}.bam.gene_fragments.gct ~{sample_name}.rnaseqc.gene_fragments.gct
     mv ~{sample_name}.bam.gene_tpm.gct ~{sample_name}.rnaseqc.gene_tpm.gct
     mv ~{sample_name}.bam.exon_reads.gct ~{sample_name}.rnaseqc.exon_reads.gct
     mv ~{sample_name}.bam.exon_cv.tsv ~{sample_name}.rnaseqc.exon_cv.tsv
     mv ~{sample_name}.bam.metrics.tsv ~{sample_name}.rnaseqc.metrics.tsv

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
   
