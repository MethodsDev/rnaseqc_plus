version 1.0

workflow Sequoia_rnaseq_plus_wf {

    input {
        String sample_name
	File fq_file
        File genome_lib_tar_with_STAR_idx
	Boolean trim_polyA

        File collapsed_ref_annot_gtf


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
	       docker=docker,
	       cpu=4,
               memory="10G",
               preemptible=preemptible,
               maxRetries=maxRetries,
               disk_space_multiplier=disk_space_multiplier
    }


    call align_reads_task {
         input:
	     sample_name=sample_name,
	     fq_file = quality_trim_task.quality_trimmed_fq_file,
	     genome_lib_tar = genome_lib_tar_with_STAR_idx,

             docker=docker,
             cpu=cpu,
             memory=memory,
             preemptible=preemptible,
             maxRetries=maxRetries,
             disk_space_multiplier=disk_space_multiplier 
    }
    

    call rnaseqc_task {
        input:
           collapsed_ref_annot_gtf=collapsed_ref_annot_gtf,
	   aligned_bam=align_reads_task.aligned_bam,
	   aligned_bam_bai=aligned_reads_task.aligned_bam_bai

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

      