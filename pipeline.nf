#!/usr/bin/env nextflow

/*
* Use echo para imprimir 'Hello World!' en un archivo
*/
process find_cases {

    output:
    path '.motu_gcc/'

    script:
    """
    Rscript workflow_interactions.r -m mOTUs_Species_Cluster -b gcf -s 5 -i /home/user/MAGs_BGCs_interactions/ -o /home/user/exp/2025-interacions/
    """
}

process graphs2 {
    output: 
    path '' 

}

workflow {

    main:
    // emite un saludo
    sayHello()
}
