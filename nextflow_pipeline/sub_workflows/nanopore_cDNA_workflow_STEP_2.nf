// Import Modules
include {MAKE_FAI} from '../modules/make_fai'
include {MAKE_INDEX_cDNA} from '../modules/make_index'
include {CHM13_GTF; CHM13_GTF_ERCC} from '../modules/chm13_gff3_to_gtf'
include {PYCHOPPER} from '../modules/pychopper'
include {PYCOQC} from '../modules/pycoqc'
include {MINIMAP2_cDNA; FILTER_BAM} from '../modules/minimap2'
include {RSEQC} from '../modules/rseqc'
include {BAMBU_PREP} from '../modules/bambu'


workflow NANOPORE_cDNA_STEP_2 {

    take:
        ref
        annotation
        housekeeping
        ont_reads_txt
        ont_reads_fq
        ercc
        cdna_kit
        track_reads
        mapq

    main:
        MAKE_FAI(ref)
        MAKE_INDEX_cDNA(ref)
        PYCHOPPER(ont_reads_fq, ont_reads_txt, cdna_kit)
        MINIMAP2_cDNA(PYCHOPPER.out.id, PYCHOPPER.out.fastq,  MAKE_INDEX_cDNA.out, PYCHOPPER.out.txt)
        FILTER_BAM(MINIMAP2_cDNA.out.id, mapq, MINIMAP2_cDNA.out.bam, MINIMAP2_cDNA.out.bai)
        
        if (params.ont_reads_txt != "None") {
            PYCOQC(MINIMAP2_cDNA.out.id, MINIMAP2_cDNA.out.fastq, MINIMAP2_cDNA.out.txt, MINIMAP2_cDNA.out.bam, MINIMAP2_cDNA.out.bai)
        }

        if (params.is_chm13 == true)
        {
            if (params.ercc == "None") 
            { 
                CHM13_GTF(annotation)
                annotation = CHM13_GTF.out.collect()
            }
            
            else 
            {
                CHM13_GTF_ERCC(annotation, ercc)
                annotation = CHM13_GTF_ERCC.out.collect()
            }
        }

        else
        {
            RSEQC(FILTER_BAM.out.id, FILTER_BAM.out.bam, FILTER_BAM.out.bai, housekeeping)
        }
        
        
        BAMBU_PREP(FILTER_BAM.out.id, mapq, FILTER_BAM.out.bam, FILTER_BAM.out.bai, ref, annotation, MAKE_FAI.out, track_reads)

}
