#!/usr/bin/env nextflow

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
QC of raw reads
#==============================================
*/

process RAWFASTQC {
	cpus 4

    conda './conda/fastqc.yaml'
	
    tag {"FastQC raw ${id} reads"}
	
	publishDir "$params.outdir/raw_fastqc", mode: 'copy'

	input:
    tuple val(id), path(reads) 
	
	output:
	path("*")
	
	script:
	
	"""
	fastqc --threads ${task.cpus} ${reads[0]}
	fastqc --threads ${task.cpus} ${reads[1]}
	"""

}

process RAWFASTQC_SINGLE {
	cpus 4

    conda './conda/fastqc.yaml'
	
    tag {"FastQC raw ${id} reads"}
	
	publishDir "$params.outdir/raw_fastqc_single", mode: 'copy'

	input:
    tuple val(id), path(reads) 
    
	
	output:
    path("${id}_raw_reads.fastq.gz"), emit: cat_raw
	path("*")
	
	script:
	
	"""
    cat ${reads[0]} ${reads[1]} > ${id}_raw_reads.fastq.gz
	fastqc --threads ${task.cpus} ${id}_raw_reads.fastq.gz
    unzip ${id}_raw_reads_fastqc.zip
    rm ${id}_raw_reads_fastqc.zip
    mv ${id}_raw_reads_fastqc/fastqc_data.txt ${id}_raw_reads_fastqc/${id}.txt
	
	"""
}

/*
#==============================================
Read cleaning with Fastp
#==============================================
*/

process FASTP {
	cpus 8

    conda './conda/fastp.yaml'
    

    tag {"filter $id reads"}
    publishDir "$params.outdir/clean_fastqs/", mode: "copy"

    input:
    tuple val(id), path(reads) 
    
    output:
    tuple val(id), path("${id}_clean_R*.fastq.gz"), emit: reads
    path("${id}.fastp.json"), emit: json
    

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${id}_clean_R1.fastq.gz -O ${id}_clean_R2.fastq.gz -w ${task.cpus} -j ${id}.fastp.json
    
    """
}

process FASTP_SINGLE {
    cpus 8

    conda './conda/fastp.yaml'
    

    tag {"filter $id reads"}
    publishDir "$params.outdir/clean_fastqs/", mode: "copy"

    input:
    tuple val(id), path(reads) 
    
    output:
    tuple val(id), path("${id}_clean_R*.fastq.gz"), emit: reads
    path("${id}.fastp.json"), emit: json
    tuple val(id), path("${id}.fastq.gz"), emit: cat_fastq

    script:
    """
    
    fastp -i ${reads[0]} -I ${reads[1]} -o ${id}_clean_R1.fastq.gz -O ${id}_clean_R2.fastq.gz -w ${task.cpus} -j ${id}.fastp.json
    cat ${id}_clean_R*.fastq.gz > ${id}.fastq.gz

    """

}

//return

/*
#==============================================
QC clean reads
#==============================================
*/

process CLEANFASTQC {
	cpus 4
	
    conda './conda/fastqc.yaml'
	
	tag {"FastQC clean ${id} reads"}

    publishDir "$params.outdir/clean_fastqc", mode: 'copy', pattern: "${id}*"

	input:
    tuple val(id), path(reads)
	
	output:
	path ("*")
    
	
	script:
	"""
    fastqc --threads ${task.cpus} ${reads}
    
	"""

}

process CLEANFASTQC_SINGLE {
	cpus 4
	
    conda './conda/fastqc.yaml'
	
	tag {"FastQC clean ${id} reads"}

    publishDir "$params.outdir/clean_fastqc_single", mode: 'copy', pattern: "${id}*"

	input:
    tuple val(id), path(cat_fastq)
    
	output:
	path ("*")
    
	script:
	"""
    fastqc --threads ${task.cpus} ${cat_fastq}
    unzip ${id}_fastqc.zip
    cp ${id}_fastqc/fastqc_data.txt ${id}_fastqc/${id}.txt
	"""

}

//return

/*
#==============================================
Collate and summarize all read QC files
#==============================================
*/

process MULTIQC_READS {
	
	conda './conda/multiqc.yaml'

	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc_reads/", mode: "copy"

	input:
    path ('*')
    
    output:
    path ("*multiqc_report.html")
    path ("*_data")

    script:

    """
    multiqc . 
    """
}

/*
#==============================================
De novo assembly
#==============================================
*/

process ASSEMBLY {
  	cpus 8

	tag { "assemble ${id}" }

	conda './conda/shovill.yaml'
  
  	publishDir "$params.outdir/assemblies/", mode: "copy"

  	input:
  	tuple val(id), path(reads)

  	output:
  	tuple val(id), path("${id}_contigs.fa"), emit: assembly

  	script:
 	"""
  	shovill --R1 ${reads[0]} --R2 ${reads[1]} --outdir shovill_${id} --cpus ${task.cpus}
	mv shovill_${id}/contigs.fa ${id}_contigs.fa
  	"""
}

/*
#==============================================
QC assembled genomes
#==============================================
*/

process QUAST_FROM_READS  {
    
    tag { "QC assembly using Quast" }

    conda './conda/quast.yaml'
    
    publishDir "$params.outdir/quast", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)  
    
    output:
    path("quast_${id}")

    script:
    """
    quast.py  ${assembly} -o quast_${id}
    cp quast_${id}/report.tsv quast_${id}/${id}_Quastreport.tsv
    """
}

process QUAST_FROM_CONTIGS  {
    
    tag { "QC assembly using Quast" }

    conda './conda/quast.yaml'
    
    publishDir "$params.outdir/quast", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)  
    
    output:
    path("quast_${assembly}")

    script:
    """
    quast.py --threads ${task.cpus} ${assembly} -o quast_${assembly}
    """
}

process MULTIQC_CONTIGS {
	
	conda './conda/multiqc.yaml'


	tag {"Collate and summarize QC files"}

	publishDir "$params.outdir/multiqc_contigs/", mode: "copy"

	input:
    path ('*')
    
    output:
    file "*multiqc_report.html"
    file "*_data"

    script:

    """
    multiqc . 
    """
}


/*
#===============================================
AMRG and Plasmid Type Profiling
#===============================================
*/

process AMR_PLM_FROM_READS {
    
    tag { "${id}: AMR finding with Abricate" }
    
    conda './conda/abricate.yaml'

    publishDir "$params.outdir/amr_plasmid/", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)  
    
    
    output:
    path("${id}_card.tab")
    path("${id}_resfinder.tab")
    path("${id}_plasmidfinder.tab")
    path("${id}_card_summary.tsv")
    path("${id}_resfinder_summary.tsv")
    path("${id}_plasmidfinder_summary.tsv")
    path("*")

    script:
    """
    abricate  --db card ${assembly} > ${id}_card.tab
    abricate --summary ${id}_card.tab > ${id}_card_summary.tsv
    abricate  --db resfinder ${assembly} > ${id}_resfinder.tab
    abricate --summary ${id}_resfinder.tab > ${id}_resfinder_summary.tsv
    abricate  --db plasmidfinder ${assembly} > ${id}_plasmidfinder.tab
    abricate --summary ${id}_plasmidfinder.tab > ${id}_plasmidfinder_summary.tsv
    """
}

process AMR_PLM_FROM_CONTIGS {
    
    tag { "AMR finding with Abricate" }
    
    conda './conda/abricate.yaml'

    publishDir "$params.outdir/amr_plasmid/", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)  
    
    
    output:
    path("*")

    script:
    """
    abricate  --db card ${assembly} > ${id}_card.tab
    abricate --summary ${id}_card.tab > ${id}_card_amr_summary.tsv
    abricate  --db resfinder ${assembly} > ${id}_resfinder_amr.tab
    abricate --summary ${id}_resfinder_amr.tab > ${id}_resfinder_summary.tsv
    abricate  --db plasmidfinder ${assembly} > ${id}_plasmidfinder.tab
    abricate --summary ${id}_plasmidfinder.tab > ${id}_plasmidfinder_summary.tsv
    """
}

/*
#===============================================
Determine MLST
#===============================================
*/

process MLST_FROM_READS {
    cpus 4

    tag {"MLST: ${id}"}

    conda 'conda/mlst.yaml'

    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)  
    
    
    output:
    path("${id}_ST.tsv")

    """
    mlst --scheme $params.mlstdb ${assembly} > ${id}_ST.tsv
    """

}

process MLST_FROM_CONTIGS {
    cpus 4

    tag {"MLST: ${assembly}"}

    conda 'conda/mlst.yaml'

    publishDir "$params.outdir/mlst/", mode: 'copy'
    
    input:
    tuple val(id), path(assembly)
    
    
    output:
    path("${id}_ST.tsv")

    """
    mlst --scheme $params.mlstdb ${assembly} > ${id}_ST.tsv
    """

}

/*
#===============================================
Read mapping and SNP calling from Illumina reads
#===============================================
*/

process SNIPPYFASTQ {
	cpus 4

	tag { "call snps from FQs: ${id}" }

	conda './conda/snippy.yaml'
	    
	publishDir "$params.outdir/snps/", mode: "copy"

    input:
    tuple val(id), path(reads), path(refFasta)
    

    output:
	path("${id}_snippy/*"), emit: snippy // whole output folder

    """
    snippy --cpus $task.cpus --outdir ${id}_snippy --prefix ${id} --reference ${refFasta} --R1 ${reads[0]} --R2 ${reads[1]} 
    """

}


/*
#==================================================
Read mapping and SNP calling from assembled contigs
#==================================================
*/

process SNIPPYFASTA {
	cpus 4

	tag { "call snps from contigs: ${id}" }

	conda './conda/snippy.yaml'
        
	publishDir "$params.outdir/snps/", mode: "copy"

    input:
    tuple val(id), path(fasta)
    path(refFasta)

    output:
    file("${id}_snippy/*"), emit: snippy // whole output folder

    """
    snippy --cpus $task.cpus --report --outdir $id --prefix $id --reference ${refFasta} --ctgs $fasta 
    """
}

/*
#==================================================
SNP alignment
#==================================================
*/

process SNIPPYCORE {
    
    tag { "create snp alignments using snippy-core" }

    conda './conda/snippy.yaml'
    

    publishDir "${params.outdir}/snippy_core", mode: "copy"
    

    input:
    path(snippy)
    path(refFasta)  // collected list of snippy output directories + ref genome
    

    output:
    //path("*")
    path("snp.core.fasta")
    path("snp.core.vcf")
    path("wgs.core.fasta"), emit: for_gubbins

    """
    snippy-core --ref ${refFasta} --prefix core ${params.outdir}/snps/*_snippy/
    mv core.aln snp.core.fasta
    mv core.vcf snp.core.vcf
    snippy-clean_full_aln core.full.aln > wgs.core.fasta
    
    """

}

process GUBBINS {
    cpus 16

    tag { "remove recombinant segments with Gubbins" }

    
    publishDir "${params.outdir}/gubbins", mode: "copy"
    

    input:
    path(for_gubbins)

    output:
    path("wgs.core.filtered_polymorphic_sites.fasta"), emit: polymorphsites
    path("*")

    script:
    """
    run_gubbins.py ${for_gubbins} --threads  $task.cpus --remove-identical-sequences --tree-builder fasttree 

    """
}

/*
#==================================================
Create a SNP-only (non-recombinant) alignment
#==================================================
*/

process SNP_SITES {

    tag { "create a SNP-only, non-rec alignment" }

    
    publishDir "${params.outdir}/snp-sites", mode: "copy"
    
    input:
    path(polymorphsites)

    output:
    path("core.full.nonrec.aln"), emit: nonrec

    script:
    """
    snp-sites -c ${polymorphsites} > core.full.nonrec.aln 

    """
    
}

/*
#==================================================
Create a pairwise SNP matrix
#==================================================
*/

process SNP_DISTS {
    tag { "create a pairwise SNP matrix" }

    publishDir "${params.outdir}/snp-dists", mode: "copy"

    input:
    path(nonrec)

    output:
    path("*")
    //path("${nonrec}.snp-dists.tsv"), snp_matrix_tsv
    //path("${nonrec}.snp-dists.csv"), snp_matrix_csv
    //path("${nonrec}.snp-dists_molten.csv"), snp_matrix_molten

    script:
    """
    snp-dists ${nonrec} > ${nonrec}.snp-dists.tsv
    snp-dists -c ${nonrec} > ${nonrec}.snp-dists.csv
    snp-dists -c ${nonrec} -m > ${nonrec}.snp-dists_molten.csv

    """
}

/*
#===================================================================
Generate a non-recombinat SNP-based phylogenomic tree using FastTree
#===================================================================
*/

process PHYLOTREE {

    tag { "generate a non-rec tree using FastTree" }

    
    publishDir "${params.outdir}/tree", mode: "copy"
    
    input:
    path(nonrec)

    output:
    path("*")

    script:
    """
    FastTree -nt ${nonrec} > core.full.nonrec.fasttree.tre
    """

}

/*
#=================================================================
Create a multi-fasta alignment using MAFFT
#=================================================================
*/

process MSA {
    conda './conda/mafft.yaml'

    publishDir "${params.outdir}/msa", mode: 'copy'

    input:
    path("*")
    
    output:
    path("mafft_consensus_msa.aln")

    script:
    """
    cat *.consensus.fa > consensus_msa.aln
    mafft consensus_msa.aln > mafft_consensus_msa.aln
    """

}


process PLATON_READS {

    cpus 8

    tag { "ID plasmids in short-read draft assemblies" }

    conda './conda/platon.yaml'

    publishDir "${params.outdir}/platon", mode: "copy"

    input:
    tuple val(id), path(assembly)

    output:
    path("*")

    script:
    """
    platon  ${assembly} --prefix ${id} --threads ${task.cpus}
    """ 

}

process PLATON_CONTIGS {
    cpus 8

    tag { "${id}: ID plasmids in short-read draft assemblies" }

    conda './conda/platon.yaml'

    publishDir "${params.outdir}/platon", mode: "copy"

    input:
    tuple val(id), path(assembly)

    output:
    path("*")

    script:
    """
    platon  ${assembly} --prefix ${id} --threads ${task.cpus} --db /mnt/microbio/HOMES/arund/HPRU/Bugflow_DSL2/bin/platon/db/
    """ 

}

process MOBTYPER {
    tag { "MOB-typer and MOB-recon: ${id}" }

    conda './conda/mobsuite.yaml'

    publishDir "${params.outdir}/mobsuite", mode: "copy"

    input:
    tuple val(id), path(assembly)

    output:
    path("*")

    script:
    """
    mob_typer -x -i ${assembly} -o ${id}_mobtyper_results.txt
    #mob_recon --infile ${assembly} --outdir ${id}_mobrecon_outdir
    """ 

}

process FLANKER {
    conda './conda/flanker.yaml'

    publishDir "${params.outdir}/flanker", mode: "copy"

    input:
    tuple val(id), path(segment)

    output:
    path("*"), emit: flank

    script:
    """
    flanker --flank both --window 0 -wstop 5000 -wstep 100 --gene blaCTX-M-15 --fasta_file ${id}.fasta --include_gene
    """

}

process EXTRACTGENE {

    publishDir "${params.outdir}/extracted_genes", mode: "copy"

    input:
    tuple val(id), path(genbank)

    output:
    path("*"), emit: extracted_gene

    script:
    """
    python3 bin/extract_startgene_endgene.py ${id}.gbk -s "IS1380 family transposase ISEcp1" -e "Beta-lactamase CTX-M-1"
    """
}
