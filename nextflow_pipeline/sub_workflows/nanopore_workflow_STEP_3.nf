// Import Modules
include {BAMBU_DISCOVERY; BAMBU_QUANT} from '../modules/bambu'
include {GFFCOMPARE} from '../modules/gffcompare'
include {MAKE_TRANSCRIPTOME} from '../modules/make_transcriptome'
include {MULTIQC_GRCh38 ; MULTIQC_CHM13} from '../modules/multiqc'


workflow NANOPORE_STEP_3 {

    take:
        ref
        fai
        annotation
        NDR
        track_reads
        bambu_rds
        multiqc_input
        multiqc_config

    main:

        if (params.is_chm13 == true)
        {
            MULTIQC_CHM13(multiqc_input.collect(), multiqc_config)
        }

        else
        {
            MULTIQC_GRCh38(multiqc_input.collect(), multiqc_config)
        }
        
        if (params.is_discovery == true)
        {

            BAMBU_DISCOVERY(bambu_rds.collect(), ref, annotation, fai, NDR, track_reads)
            new_annotation = BAMBU_DISCOVERY.out.gtf
            GFFCOMPARE(new_annotation, annotation)

        }

        else
        {
            
            BAMBU_QUANT(bambu_rds.collect(), ref, annotation, fai)
            new_annotation = BAMBU_QUANT.out.gtf

        }

        MAKE_TRANSCRIPTOME(ref, fai, new_annotation)
        
}
