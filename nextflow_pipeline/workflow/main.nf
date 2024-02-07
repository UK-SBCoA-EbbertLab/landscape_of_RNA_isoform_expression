// Make this pipeline a nextflow 2 implementation
nextflow.enable.dsl=2


log.info """
   OXFORD NANOPORE cDNA SEQUENCING PIPELINE - Bernardo Aguzzoli Heberle - EBBERT LAB - University of Kentucky
 ==============================================================================================================
 nanopore fastq files                                           : ${params.ont_reads_fq}
 nanopore sequencing summary files                              : ${params.ont_reads_txt}
 reference genome                                               : ${params.ref}
 reference annotation                                           : ${params.annotation}
 housekeeping genes 3' bias assessment                          : ${params.housekeeping}
 nanopore library prep kit                                      : ${params.cdna_kit}
 multiqc configuration file                                     : ${params.multiqc_config}

 reference genome is CHM13                                      : ${params.is_chm13}
 transcript discovery status                                    : ${params.is_discovery}

 nanopore fast5 files (basecall only)                           : ${params.fast5_dir}
 nanopore basecall config (basecall only)                       : ${params.basecall_config}
 nanopore basecall id (basecall only)                           : ${params.basecall_id}

 NDR Value for Bambu (Novel Discovery Rate)                     : ${params.NDR}
 Track read_ids with bambu?                                     : ${params.track_reads}

 MAPQ value for filtering bam file                              : ${params.mapq}

 Step: 1 = basecalling, 2 = mapping, 3 = quantification         : ${params.step}

 Path to pre-processed bambu RDS files                          : ${params.bambu_rds}
 Path to QC files that go into MultiQC report                   : ${params.multiqc_input}   

 Is this a direct RNAseq dataset?                               : ${params.is_dRNA}
 ==============================================================================================================
 """


// Import Workflows
include {NANOPORE_STEP_1} from '../sub_workflows/nanopore_workflow_STEP_1'
include {NANOPORE_cDNA_STEP_2} from '../sub_workflows/nanopore_cDNA_workflow_STEP_2'
include {NANOPORE_dRNA_STEP_2} from '../sub_workflows/nanopore_dRNA_workflow_STEP_2'
include {NANOPORE_STEP_2_BAM} from '../sub_workflows/nanopore_workflow_STEP_2_BAM'
include {NANOPORE_STEP_3} from '../sub_workflows/nanopore_workflow_STEP_3'


// Define initial files and channels
ont_reads_fq = Channel.fromPath(params.ont_reads_fq).map { file -> tuple(file.baseName, file) }
ont_reads_txt = Channel.fromPath(file(params.ont_reads_txt))
ref = file(params.ref)
housekeeping = file(params.housekeeping)
annotation = file(params.annotation)
fast5_dir = Channel.fromPath(params.fast5_dir)
basecall_config = Channel.from(params.basecall_config)
basecall_id = Channel.from(params.basecall_id)
cdna_kit = Channel.value(params.cdna_kit)
multiqc_config = Channel.fromPath(params.multiqc_config)
NDR = Channel.value(params.NDR)
track_reads = Channel.value(params.track_reads)
mapq = Channel.value(params.mapq)
bambu_rds = Channel.fromPath(params.bambu_rds)
multiqc_input = Channel.fromPath(params.multiqc_input, type: "file")
fai = file(params.fai)
bam = Channel.fromPath(params.bam).map { file -> tuple(file.baseName, file) }
bai = Channel.fromPath(params.bai)


if (params.ercc != "None") {
    ercc = Channel.fromPath(params.ercc)
    }
else {
    ercc = params.ercc
    }

if (params.ont_reads_txt == "None") {
    ont_reads_txt = Channel.value(params.ont_reads_txt)
    } else {
    // Make sure ONT sequencing summary and fastq files are in the same order
    ont_reads_txt = ont_reads_txt.toSortedList( { a, b -> a.baseName <=> b.baseName } ).flatten()
    }

if (params.ont_reads_fq != "None") {
    
    // Make sure files are in same order
    ont_reads_fq = ont_reads_fq.toSortedList( { a, b -> a[0] <=> b[0] } ).flatten().buffer(size:2)

    }

if ((params.bam != "None") && (params.bai != "None")) {

    // Make sure bam and bai files are in the correct order
    bam = bam.toSortedList( { a, b -> a[0] <=> b[0] } ).flatten().buffer(size:2)
    bai = bai.toSortedList( { a, b -> a.baseName <=> b.baseName } ).flatten()
    
    }

workflow {

    if (params.step == 1){
        NANOPORE_STEP_1(fast5_dir, basecall_config, basecall_id)
    }

    else if ((params.step == 2) && (params.bam == "None")){

        if (params.is_dRNA == "False") {
        
            NANOPORE_cDNA_STEP_2(ref, annotation, housekeeping, ont_reads_txt, ont_reads_fq, ercc, cdna_kit, track_reads, mapq)
        }

        else if (params.is_dRNA = "True") {

            NANOPORE_dRNA_STEP_2(ref, annotation, housekeeping, ont_reads_txt, ont_reads_fq, ercc, cdna_kit, track_reads, mapq)

        }
    }


    else if ((params.step == 2) && (params.bam != "None")) {
        
        NANOPORE_STEP_2_BAM(ref, annotation, bam, bai, ercc, track_reads, mapq)

    }

    else if(params.step == 3){
        
        NANOPORE_STEP_3(ref, fai, annotation, NDR, track_reads, bambu_rds, multiqc_input, multiqc_config)
    }

}
