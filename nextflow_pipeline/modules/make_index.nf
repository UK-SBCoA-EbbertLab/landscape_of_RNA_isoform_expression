process MAKE_INDEX_cDNA {

    label 'small'

    input:
        path ref

    output:
        path 'ref.mmi'


    script:
        """
        minimap2 -t 8 -ax splice -uf -d ref.mmi $ref
        """
}

process MAKE_INDEX_dRNA {

    label 'small'

    input:
        path ref

    output:
        path 'ref.mmi'


    script:
        """
        minimap2 -t 8 -k14 -ax splice -uf -d ref.mmi $ref
        """
}
