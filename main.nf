#!/usr/bin/env nextflow

/* 

QC
 - Fastp
 - FastQC
 - MultiQC
 - Quast

Mapping and Variant calling
 - Snippy
 
Assembly
 - Shovill (spades) 
*/

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
Modules
#==============================================
*/

include { ASSEMBLY } from './modules/processes.nf'
include { QUAST_FROM_READS; QUAST_FROM_CONTIGS } from './modules/processes.nf'
include { SNIPPYFASTQ } from './modules/processes.nf'
include { SNIPPYFASTA } from './modules/processes.nf' 
include { SNIPPYCORE } from './modules/processes.nf'
include { AMR_PLM_FROM_READS; AMR_PLM_FROM_CONTIGS; PLATON_READS; PLATON_CONTIGS; MOBTYPER; SUMMARY_BLASTN; AMR_ABRFORMAT} from './modules/processes.nf'
include { MLST_FROM_READS; MLST_FROM_CONTIGS } from './modules/processes.nf'
include { GUBBINS; SNP_SITES; SNP_DISTS; PHYLOTREE } from './modules/processes.nf'

/*
#==============================================
Parameters
#==============================================
*/

params.ref = "test_data/EC958.fasta"
params.reads = ""
params.outdir = " "
params.contigs = " "
params.mlstdb = "ecoli"
params.prefix = "core"
params.blastn = " "
params.segment = ""

/*
#==============================================
Channels
#==============================================
*/

//Channel.fromFilePairs(params.reads, checkIfExists: true)
       //.map{it}
       //.view()
       //.set{reads}

//Channel.fromPath(params.ref, checkIfExists:true)
        //.view()       
       //.set{refFasta}

//Channel.fromPath(params.contigs, checkIfExists:true)
       //.set{assembly}
       //.view()

//return

/*
#==============================================
Workflows
#==============================================
*/


workflow qc_contigs {
    Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}
    main:       
    QUAST(assembly)
    MULTIQC(QUAST.out.collect())
}

workflow mlst_amr_plm_reads {
    Channel.fromPath(params.reads, checkIfExists:true)
           //.view()
           .set{reads}
    main:
    FASTP(reads)
    ASSEMBLY(FASTP.out.reads)
    MLST_FROM_READS(ASSEMBLY.out.assembly)
    AMR_PLM_FROM_READS(ASSEMBLY.assembly)
}

workflow mlst_amr_plm_contigs {
    Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}
    main:
    MLST_FROM_CONTIGS(assembly)
    AMR_PLM_FROM_CONTIGS(assembly)
    
}
      

workflow snippy_fasta {
    Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}

    Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}

    main:
    SNIPPYFASTA(assembly.combine(refFasta))
    
    emit:
    SNIPPYFASTA.out // results
}  

//return

workflow  snippy_core_tree {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       Channel.fromPath(params.ref, checkIfExists:true)
           //.view()       
           .set{refFasta}
       FASTP(reads)
       SNIPPYFASTQ(FASTP.out.reads.combine(refFasta))
       SNIPPYCORE(SNIPPYFASTQ.out, refFasta)
       GUBBINS(SNIPPYCORE.out.for_gubbins)
       SNP_SITES(GUBBINS.out.polymorphsites) 
       SNP_DISTS(SNP_SITES.out.nonrec)
       PHYLOTREE(SNP_SITES.out)      
}


workflow assembly_plmcharac_amr_ns {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       
       main:
       RAWFASTQC(reads)
       FASTP(reads)
       FASTP_SINGLE(reads)
       //CLEANFASTQC(FASTP.out.reads)
       CLEANFASTQC_SINGLE(FASTP_SINGLE.out.cat_fastq)
       //MULTIQC_READS(RAWFASTQC.out.mix(CLEANFASTQC.out).collect())
       //MULTIQC_READS(CLEANFASTQC.out.collect())
       MULTIQC_READS(CLEANFASTQC_SINGLE.out.collect())
       ASSEMBLY(FASTP.out.reads)
       QUAST_FROM_READS(ASSEMBLY.out.assembly)
       MULTIQC_CONTIGS(QUAST_FROM_READS.out.collect())
       MLST_FROM_READS(ASSEMBLY.out.assembly)
       AMR_PLM_FROM_READS(ASSEMBLY.out.assembly)
       PLATON_READS(ASSEMBLY.out.assembly)
       //MOBTYPER(ASSEMBLY.out.assembly)
       
}

workflow rawqc {
       Channel.fromFilePairs(params.reads, checkIfExists: true)
           .map{it}
           //.view()
           .set{reads}
       
       main:
       RAWFASTQC_SINGLE(reads)
}

workflow hpru_annot {
       Channel.fromPath(params.contigs, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{assembly}

       main:
       MLST_FROM_CONTIGS(assembly)
       AMR_PLM_FROM_CONTIGS(assembly)
       PLATON_CONTIGS(assembly)
    
}

workflow flanker {
       Channel.fromPath(params.segment, checkIfExists:true)
           .map({it -> tuple(it.baseName, it)})
           //.view()
           .set{segment}
       main:
       FLANKER(segment)
}

workflow extractgenes {
       Channel.fromPath(params.genbank_dir)
              .map({it -> tuple(it.baseName, it)})
              .set {genbank}

       main:
       EXTRACTGENE(genbank)
}