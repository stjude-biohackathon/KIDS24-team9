#!/usr/bin/env nextflow
/***********************************************************************************
 This is a Nextflow pipeline-BALL Classification Container  
 Input: Fastq filelist   
 Output: 
 Written by: Wenchao Zhang
 The Center for Applied Bioinformatics,St. Jude Children Research Hosptial
 Date: 09/20/2022
***********************************************************************************/


RANK_Training_RDS  =""
RF_Training_RDS =""
def parse_config_parameters() {
    if( params.Lineage_Type == 'BALL' )
    {
  
      RANK_Training_RDS = params.ML_Classifier.BALL.RANK_Training_RDS
      RF_Training_RDS =  params.ML_Classifier.BALL.RF_Training_RDS
  
     }
     else if( params.Lineage_Type  == 'TALL' )
     {

       RANK_Training_RDS = params.ML_Classifier.TALL.RANK_Training_RDS
       RF_Training_RDS =  params.ML_Classifier.TALL.RF_Training_RDS
     }
     else if( params.Lineage_Type  == 'AML' )
     {
        RANK_Training_RDS= params.ML_Classifier.AML.RANK_Training_RDS
        RF_Training_RDS= params.ML_Classifier.AML.RF_Training_RDS
     }
     else
     {
        error "Invalid Lineage_Type: ${params.Lineage_Type}"
        exit 1
     }
}

def DispConfig() {
 log.info """
Welocme to run Nextflow Pipeline DIVIA.nf  
Your configuration are the following:
  project               : ${params.project}
  fastq_filelist        : ${params.fastq_filelist}
  outdir                : ${params.outdir}
  
  Lineage_Type          : ${params.Lineage_Type}
  ChrName_withChr       : ${params.ChrName_withChr}

  Select_Trim_Galore    : ${params.Select_Trim_Galore}
  Select_Fusion_Catcher : ${params.Select_Fusion_Catcher}
  Select_STAR_Fusion    : ${params.Select_STAR_Fusion}
  Select_Arriba_Fusion  : ${params.Select_Arriba_Fusion}

  Select_GATK           : ${params.Select_GATK}
  Select_Haplo_PerChrom : ${params.Select_Haplo_PerChrom}
  Select_VarDict        : ${params.Select_VarDict}

  Select_Variant_ANNOVAR: ${params.Select_Variant_ANNOVAR}

  Select_RSEM           : ${params.Select_RSEM} 
  Select_ML_RANK        : ${params.Select_ML_Classifier_RANK}
  Select_ML_tSNE        : ${params.Select_ML_Classifier_tSNE}
  Select_ML_RF          : ${params.Select_ML_Classifier_RF}

  Select_RNASeqCNV      : ${params.Select_RNASeqCNV}  
  Select_RNAIndel       : ${params.Select_RNAIndel}  
  
  Select_Pindel         : ${params.Select_Pindel}
  Select_iAdmix         : ${params.Select_iAdmix}

  The following detail information is based your configuration: 
 
 """  
//  exit 0    
}

def helpMessage() {
  log.info """
        Welocme to run Nextflow Pipeline DIVIA.nf 
        Usage:
        A typical command for running the pipeline is as follows:
        nextflow run DIVIA.nf -profile cluster_singularity 
        A command with more configurable arguments can be           
        nextflow run DIVIA.nf -profile cluster -w ~/Nextflow_work -- ./

        Configurable arguments:
        --project                     Project name/atlas, a folder that used to distinguish other analysis 
        --Fastq_filelist              Fastq file list that you want to do BALL classification  
        --Lineage_Type                ALL Lineage TYpe. BALL |TALL |AML   
        --ChrName_withChr             Whether ChrName has prefix chr. Y|N

        --outdir                      The directory for the pipeline output           
        --Select_Trim_Galore          Whether select Trim_Galore for trimming reads and FastQC. Y(Default) | N 
        --Select_Fusion_Catcher       Whether select Fusion_Catcher for Fusion Detection. Y(Default) | N 
        --Select_STAR_Fusion          Whether select STAR_Fusion for Fusion Detection. Y(Default) | N
        --Select_Arriba_Fusion        Whether select Arriba for Fusion Detection. Y(Default) | N 
             
        --Select_GATK                 Whether select GATK for rna variant calling. Y(Default) | N
        --Select_Haplo_PerChrom       Whether select Chrom Per chrom to run HaplotypeCaller. Y (Default) |N
        --Select_VarDict              Whether select to VarDict for Variant calling.  Y (Default) |N.
        --Select_Variant_ANNOVAR      Whether select ANNOVAR for the called RNA variant annotation. Y(Default) |N
       
        --Select_RSEM                 Whether select RSEM for rnaseq quantification. Y(Default) | N
        --Select_ML_Classifier_RANK   Whether select RANK for GE-ML Classifier. Y(Default) | N
        --Select_ML_Classifier_tSNE   Whether select tSNE for GE-ML Classifier. Y(Default) | N
        --Select_ML_Classifier_RF     Whether select RandomForest for GE-ML Classifier. Y(Default) | N
        --Select_RNASeqCNV            Whether select RNASeqCNV for CNV calling. Y(Default) | N 
        --Select_RNAIndel             Whether select RNAIndel for calling coding indels from tumor RNA-Seq data. Y(Default) | N
        --Select_Pindel               Whether select Pindel for detecting the breakpoints of large deletions, insertions, inversions, and tandem duplications Y(Default) | N
        --Select_iAdmix               Whether select iAdmix to estimate admixture coefficients from aligned BAM. Y(Default) | N. 
       --help | h                    This usage statement.
       Note:
         All the above arguments can be configured by the command line interface or in the nextflow.config (default)  
        """
}


// Show help message
if ((params.help) || (params.h)) {
    helpMessage()
    exit 0
}

//Parse the input Parameters and configs. 
parse_config_parameters() 
// Display the configuration
DispConfig() 

// Set up the two external input channels for samples (fastq files)  
Read_Fastq_Ch  = Channel
            .fromPath(params.fastq_filelist)
            .splitText()
            .splitCsv(sep: '\t')
            //.splitCsv(sep: ',')

// HaplotypeCalling can be implemented chromosome by chromosome. Need to consider two scenario with/without prefixing "Chr"     
chromosomes_ch = (params.ChrName_withChr=="Y")? Channel
    .from( "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM" ): Channel
    .from( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M" )


/*This process is responsible for zcat the multiple lanes' fastq.gz into one mergred fastq.gz  */
process Zcat_MergeFastq {
   publishDir "${params.outdir}/${params.project}/Zcat_MergeFastq/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, Read_Fastq, Strandness from Read_Fastq_Ch
   
   output:
   set SampleName, file("${SampleName}.Merge_Read_R1.fastq.gz"), file("${SampleName}.Merge_Read_R2.fastq.gz"), Strandness into MergeFastq_Ch
   
   // scratch ${params.TMP_DIR}

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   zcat `echo $Read_Fastq | cut -d' ' -f1 |tr ',' ' '` >${SampleName}.Merge_Read_R1.fastq && gzip ${SampleName}.Merge_Read_R1.fastq
   zcat `echo $Read_Fastq | cut -d' ' -f2 |tr ',' ' '` >${SampleName}.Merge_Read_R2.fastq && gzip ${SampleName}.Merge_Read_R2.fastq
   """  
}

MergeFastq_Ch.into {MergeFastq_Ch1; MergeFastq_Ch2 } 
/**********************************************************************************************
* This Process is responsible for trimming reads by Trim_Galore
***********************************************************************************************/
Trim_Galore_In_Ch = (params.Select_Trim_Galore == "N" ) ? Channel.empty(): MergeFastq_Ch1
process Trim_Galore {
   publishDir "${params.outdir}/${params.project}/Trim_Galore/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(Merge_Read_Fastq_R1), file(Merge_Read_Fastq_R2), Strandness from Trim_Galore_In_Ch
   
   output:
   set SampleName, Strandness, file("${SampleName}_val_1.fq.gz"), file("${SampleName}_val_2.fq.gz"), file("${SampleName}_R1_unpaired_1.fq.gz"), file("${SampleName}_R2_unpaired_2.fq.gz"), file("${SampleName}_val_1_fastqc.html"), file("${SampleName}_val_2_fastqc.html"),  file("${Merge_Read_Fastq_R1}_trimming_report.txt"), file("${Merge_Read_Fastq_R2}_trimming_report.txt"), file("${SampleName}_val_1_fastqc.zip"), file("${SampleName}_val_2_fastqc.zip") into Trim_Galore_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   trim_galore \
     --paired  \
     --retain_unpaired \
     --cores 4 \
     --output_dir . \
     --basename ${SampleName} \
     --fastqc_args "--outdir . " \
     ${Merge_Read_Fastq_R1} \
     ${Merge_Read_Fastq_R2}
   """
}

/***********************************************************************************************
*  This process is responsible for adapting the trim_galore output to STAR mapping
************************************************************************************************/
process Trim_Galore_Adaptor_STARMapping {
   input:
   set SampleName, Strandness, file(val_1_fq_gz), file(val_2_fq_gz), file(unpaired_1_fq_gz), file(unpaired_2_fq_gz), file(va1_1_fastqc_html), file(va1_2_fastqc_html),  file(R1_fastq_gz_trimming_report), file(R2_fastq_gz_trimming_report), file(val_1_fastqc_zip), file(val_2_fastqc_zip) from Trim_Galore_Ch

   output:
   set SampleName, file("${SampleName}.Trim_Read_R1.fastq.gz"), file("${SampleName}.Trim_Read_R2.fastq.gz"), Strandness into Trim_Galore_Adaptor_STARMapping_Ch
   
   script:
   """
   ln -s ${val_1_fq_gz} ${SampleName}.Trim_Read_R1.fastq.gz
   ln -s ${val_2_fq_gz} ${SampleName}.Trim_Read_R2.fastq.gz
   """
}


STAR_Mapping_Upstream_Ch = (params.Select_Trim_Galore == "N" ) ? MergeFastq_Ch2 : Trim_Galore_Adaptor_STARMapping_Ch
STAR_Mapping_Upstream_Ch.into { STAR_Mapping_In_Ch; Fusion_Catcher_Adaptor_Ch}
/*************************************************************************************************************
* This Process is be responsible for Fusion Catcher, which start from the fastq file instead of BAM as input
*************************************************************************************************************/
Fusion_Catcher_In_Ch= (params.Select_Fusion_Catcher == "N" ) ? Channel.empty() :Fusion_Catcher_Adaptor_Ch
process Fusion_Catcher {
   publishDir "${params.outdir}/${params.project}/Fusion_Catcher/${SampleName}", mode: 'copy', overwrite: true
   memory { (64.GB + (32.GB * task.attempt)) } // First attempt 64GB, second 96GB, etc
   errorStrategy 'retry'
   maxRetries 3
   
   input:
   set SampleName, file(R1_Read_Fastq), file(R2_Read_Fastq), Strandness from Fusion_Catcher_In_Ch   
   
   output:
   set SampleName, file("info.txt"), file("junk-chimeras.txt"), file("summary_candidate_fusions.txt"), file("viruses_bacteria_phages.txt"), file("final-list_candidate-fusion*.*"), file("supporting-reads_gene-fusions*.*") into Fusion_Catcher_Ch 
   
   script:
   """ 
   export TMPDIR=${params.TMP_DIR}
   fusioncatcher \
   -d ${params.Fusion_Catcher.data_dir} \
   -i ${R1_Read_Fastq},${R2_Read_Fastq} \
   --threads ${params.Fusion_Catcher.ThreadN} \
   --o . \
   --aligners blat,star,bowtie2 \
   --skip-star \
   --skip-blat
   """ 
}

/***********************************************************************************************
* This Process is be responsible for STAR Mapping that will used in the downstream STAR_Fusion
************************************************************************************************/
process STAR_Mapping {
   
   publishDir "${params.outdir}/${params.project}/STAR_Mapping/${SampleName}", mode: 'copy', overwrite: true

    memory { (64.GB + (32.GB * task.attempt)) } // First attempt 64GB, second 96GB, etc
    errorStrategy 'retry'
    maxRetries 3

   input:
// set SampleName, Read_Fastq from Read_Fastq_Ch
   set SampleName, file(R1_Read_Fastq), file(R2_Read_Fastq), Strandness from STAR_Mapping_In_Ch   

   output:
   set SampleName, Strandness, file("${SampleName}.STAR.Aligned.sortedByCoord.out.bam"), file("${SampleName}.STAR.Aligned.toTranscriptome.out.bam"), file("${SampleName}.STAR.ReadsPerGene.out.tab"),  file("${SampleName}.STAR.SJ.out.tab"), file("${SampleName}.STAR.Chimeric.out.junction"), file("${SampleName}.STAR.Log.final.out") into STAR_Mapping_Ch  
  
   // --sjdbOverhang ${params.read_length - 1}
   // --outStd BAM_Unsorted \
   // --readFilesIn $Read_Fastq 
   // def STRAND_TYPE=(Strandness=="Unstranded")? "Unstranded" : "Stranded"  

   script:
   def STRAND_TYPE=(Strandness =~ /Stranded_/)? "Stranded" : "Unstranded"
   """ 
   export TMPDIR=${params.TMP_DIR}
   STAR \
 --limitOutSJcollapsed 4000000 \
 --limitSjdbInsertNsj 3000000  \
 --limitBAMsortRAM 67108864000 \
 --genomeDir ${params.STAR_Mapping.STAR_Index} \
 --runThreadN ${params.STAR_Mapping.ThreadN} \
 --readFilesIn ${R1_Read_Fastq} ${R2_Read_Fastq} \
 --readFilesCommand zcat \
 --outFilterType BySJout \
 --outFilterMultimapNmax 20 \
 --alignSJoverhangMin 8 \
 --alignSJstitchMismatchNmax 5 -1 5 5 \
 --alignSJDBoverhangMin 10 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 100000 \
 --alignMatesGapMax 100000 \
 --genomeLoad NoSharedMemory \
 --outFileNamePrefix ${SampleName}.STAR. \
 --outSAMmapqUnique 60 \
 --outSAMmultNmax 1 \
 --outSAMstrandField intronMotif \
 --outSAMattributes NH HI AS nM NM MD \
 --outSAMunmapped Within \
 --outSAMtype BAM SortedByCoordinate \
 --outReadsUnmapped None \
 --outSAMattrRGline ID:${SampleName} LB:LIB01 PL:ILLUMINA SM:${SampleName} PU:H3MHFDMXX \
 --chimSegmentMin 12 \
 --chimJunctionOverhangMin 12 \
 --chimSegmentReadGapMax 3 \
 --chimMultimapNmax 10 \
 --chimMultimapScoreRange 10 \
 --chimNonchimScoreDropMin 10 \
 --chimOutJunctionFormat 1 \
 --chimOutType Junctions WithinBAM SoftClip \
 --quantMode TranscriptomeSAM GeneCounts \
 --twopassMode Basic \
 --peOverlapNbasesMin 12 \
 --peOverlapMMp 0.1 \
 --outWigType wiggle \
 --outWigStrand ${STRAND_TYPE} \
 --outWigNorm RPM
 """  
}

/*Fork the STAR Mapping stream into several downstream channaels*/
STAR_Mapping_Ch.into{ STAR_Mapping_Ch1; STAR_Mapping_Ch2; STAR_Mapping_Ch3; STAR_Mapping_Ch4; STAR_Mapping_Ch5; STAR_Mapping_Ch6 }

/***********************************************************************************************
* This Process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream STAR Fusion
************************************************************************************************/
Adaptor_STAR_Fusion_In_Ch = (params.Select_STAR_Fusion == "N" ) ? Channel.empty() : STAR_Mapping_Ch1
process STARMapping_Adaptor_STAR_Fusion {
// publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_STAR_Fusion/$SampleName", mode: 'copy', overwrite: true

   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(Chimeric_out_junction), file(Log_final_out) from Adaptor_STAR_Fusion_In_Ch
   
   output:
   set SampleName, file("${SampleName}.Chimeric.out.junction") into Adaptor_STAR_Fusion_Ch
   //ln -s ${Chimeric_out_junction}  ${SampleName}.Chimeric.out.junction

   script:
   """
   echo "need to convert 21-column Chimeric junction file to 16-column Chimeric junction file" 
   cat ${Chimeric_out_junction}| grep "^#" > Last_twolines_comment.txt
   cat ${Chimeric_out_junction}| grep -v "^#" | awk 'BEGIN {OFS="\t"}; NR>1 { print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14,\$15,\$21}' > ${SampleName}.Chimeric.out.junction
cat Last_twolines_comment.txt >> ${SampleName}.Chimeric.out.junction
   """
}


/***********************************************************************************************
* This Process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream Arriba_Fusion
************************************************************************************************/
Adaptor_Arriba_Fusion_In_Ch = (params.Select_Arriba_Fusion == "N" ) ? Channel.empty() : STAR_Mapping_Ch2
process STARMapping_Adaptor_Arriba_Fusion {
 //publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_Arriba_Fusion/$SampleName", mode: 'copy', overwrite: true

   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(Chimeric_out_junction), file(Log_final_out) from Adaptor_Arriba_Fusion_In_Ch
   
   output:
   set SampleName, file("${SampleName}.aligned.sorted.bam") into Adaptor_Arriba_Fusion_Ch

   script:
   """
   ln -s ${Align_sortedByCoord_Bam}  ${SampleName}.aligned.sorted.bam
   """
}

/***********************************************************************************************
* This Process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream RNASeq_RSEM
**********************************************************************************************/
Adaptor_RNAseq_RSEM_In_Ch = (params.Select_RSEM == "N") ? Channel.empty() : STAR_Mapping_Ch3
process STARMapping_Adaptor_RNAseq_RSEM {
 //publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_RNAseq_RSEM/$SampleName", mode: 'copy', overwrite: true

   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from Adaptor_RNAseq_RSEM_In_Ch
   
   output:
   set SampleName, Strandness, file("${SampleName}.aligned.toTranscriptome.bam") into Adaptor_RNAseq_RSEM_Ch

   script:
   """
   ln -s $Align_toTranscriptome_Bam ${SampleName}.aligned.toTranscriptome.bam
   """
}

/*****************************************************************************************************************************************
* This Process is be responsible for generating the mark duplicated bam, which can be used in the downstream GATK based Variant Calling and RNAIndel calling
*******************************************************************************************************************************************/
process STARMapping_MarkDuplicate {
    publishDir "${params.outdir}/${params.project}/STARMapping_MarkDuplicate/$SampleName", mode: 'copy', overwrite: true
   
   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from STAR_Mapping_Ch4
   
   output:
   set SampleName, file("${SampleName}.marked_dup.bam"), file("${SampleName}.marked_dup.bam.bai") into STARMapping_MarkDuplicate_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}

   samtools view -f 0x2 \
   -b ${Align_sortedByCoord_Bam} \
   -o ${SampleName}_remove_sameread.bam
 
   gatk MarkDuplicates \
   -I ${SampleName}_remove_sameread.bam \
   -O ${SampleName}.marked_dup.bam \
   -M ${SampleName}.marked_dup.metrics.txt \
   --TMP_DIR ${params.TMP_DIR}
 
   samtools index ${SampleName}.marked_dup.bam  
   """    
}

STARMapping_MarkDuplicate_Ch.into { STARMapping_MarkDuplicate_Ch1; STARMapping_MarkDuplicate_Ch2 }
/***********************************************************************************************
* This Process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream GATK based Variant Calling
************************************************************************************************/
// Adaptor_RNAvar_GATK_In_Ch = (params.Select_GATK == "N") ? Channel.empty() : STAR_Mapping_Ch4
Adaptor_RNAvar_GATK_In_Ch = (params.Select_GATK == "N") ? Channel.empty() : STARMapping_MarkDuplicate_Ch1
process STARMapping_Adaptor_RNAvar_GATK {
// publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_RNAvar_GATK/$SampleName", mode: 'copy', overwrite: true

   input: 
   set SampleName, file(marked_dup_bam), file(marked_dup_bam_bai) from Adaptor_RNAvar_GATK_In_Ch
   // set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from Adaptor_RNAvar_GATK_In_Ch
   
   output:
   set SampleName, file("${SampleName}.GATK.marked_dup.bam"), file("${SampleName}.GATK.marked_dup.bam.bai") into Adaptor_RNAvar_GATK_Ch
   
   // export TMPDIR=${params.TMP_DIR} 
   // ln -s ${Align_sortedByCoord_Bam}  ${SampleName}.aligned.sorted.bam
   // samtools index ${SampleName}.aligned.sorted.bam
   script:
   """ 
   ln -s ${marked_dup_bam}  ${SampleName}.GATK.marked_dup.bam
   ln -s ${marked_dup_bam_bai}  ${SampleName}.GATK.marked_dup.bam.bai
   """
   
}

/****************************************************************************************************************************
* This process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream moudle- Pindel
*****************************************************************************************************************************/
Adaptor_Pindel_In_Ch = (params.Select_Pindel == "N") ? Channel.empty() : STAR_Mapping_Ch5
process STARMapping_Adaptor_Pindel {
 //publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_Pindel/$SampleName", mode: 'copy', overwrite: true
   
   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from Adaptor_Pindel_In_Ch
   
   output:
   set SampleName, file("${SampleName}.aligned.sorted.bam"), file("${SampleName}.aligned.sorted.bam.bai") into Adaptor_Pindel_Ch
   
   script:
   """
   export TMPDIR=${params.TMP_DIR}
   ln -s ${Align_sortedByCoord_Bam}  ${SampleName}.aligned.sorted.bam  
   samtools index ${SampleName}.aligned.sorted.bam
   """
     
}


/****************************************************************************************************************************
* This process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream moudle- iAdmix
*****************************************************************************************************************************/
Adaptor_iAdmix_In_Ch = (params.Select_iAdmix == "N") ? Channel.empty() : STAR_Mapping_Ch6
process STARMapping_Adaptor_iAdmix {
 // publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_iAdmix/$SampleName", mode: 'copy', overwrite: true
   
   input:
   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from Adaptor_iAdmix_In_Ch
   
   output:
   set SampleName, file("${SampleName}.aligned.sorted.bam") into Adaptor_iAdmix_Ch
   
   script:
   """
   ln -s ${Align_sortedByCoord_Bam}  ${SampleName}.aligned.sorted.bam  
   """ 
     
}

/****************************************************************************************************************************
* This process is be responsible for adjusting/editting STAR Mapping's outputs to adapt to the downstream moudle- RNAIndel
*****************************************************************************************************************************/
Adaptor_RNAIndel_In_Ch = (params.Select_RNAIndel == "N") ? Channel.empty() : STARMapping_MarkDuplicate_Ch2
process STARMapping_Adaptor_RNAIndel {
   publishDir "${params.outdir}/${params.project}/STARMapping_Adaptor_RNAIndel/$SampleName", mode: 'copy', overwrite: true
   
   input:
   set SampleName, file(marked_dup_bam), file(marked_dup_bam_bai) from Adaptor_RNAIndel_In_Ch
//   set SampleName, Strandness, file(Align_sortedByCoord_Bam), file(Align_toTranscriptome_Bam), file(ReadsPerGene_out_tab), file(SJ_Out_Tab),  file(STAR_Chimeric_out_junction), file(Log_final_out) from Adaptor_RNAIndel_In_Ch
   
   output:
   set SampleName, file("${SampleName}.RNAIndel.marked_dup.bam"), file("${SampleName}.RNAIndel.marked_dup.bam.bai") into Adaptor_RNAIndel_Ch
   
   script:
   """
   ln -s ${marked_dup_bam} ${SampleName}.RNAIndel.marked_dup.bam
   ln -s ${marked_dup_bam_bai} ${SampleName}.RNAIndel.marked_dup.bam.bai
   """ 
}


/********************************************************************************************************
* This Process is be responsible for STAR_Fusion
**********************************************************************************************************/
process STAR_Fusion {
   publishDir "${params.outdir}/${params.project}/STAR_Fusion/$SampleName", mode: 'copy', overwrite: true
  // export SINGULARITYENV_APPEND_PATH="/usr/local/arriba_v2.3.0:/usr/local/STAR-Fusion/ctat-genome-lib-builder:/usr/local/STAR-Fusion:/usr/local/RSEM" 

   memory { (64.GB + (32.GB * task.attempt)) } // First attempt 64GB, second 96GB, etc
   errorStrategy 'retry'
   maxRetries 3

   input:
   set SampleName, file(Chimeric_out_junction) from Adaptor_STAR_Fusion_Ch
   
   output:
   set SampleName, file("${SampleName}_STAR_Fusion_Result.tar.gz") into STAR_Fusion_Out_Ch
  // set file("star-fusion.fusion_predictions.abridged.coding_effect.tsv"), file("star-fusion.fusion_predictions.tsv"), file("star-fusion.fusion_predictions.abridged.coding_effect.tsv") into STAR_Fusion_Out_Ch
  
 // /usr/local/src/STAR-Fusion/STAR-Fusion   ${SampleName}_STAR_Fusion_Result
   script:
   """       
    export TMPDIR=${params.TMP_DIR}
    STAR-Fusion \
    --genome_lib_dir ${params.STAR_Fusion.Ref} \
    -J ${Chimeric_out_junction} \
    --CPU ${params.STAR_Fusion.ThreadN} \
    --examine_coding_effect \
    --output_dir ./${SampleName}_STAR_Fusion_Result
       
    tar -zcvf ${SampleName}_STAR_Fusion_Result.tar.gz ${SampleName}_STAR_Fusion_Result 
   """  
}

/********************************************************************************
* This Process is be responsible for Arriba_Fusion 
*********************************************************************************/ 
process Arriba_Fusion {
   publishDir "${params.outdir}/${params.project}/Arriba_Fusion/$SampleName", mode: 'copy', overwrite: true
 //export SINGULARITYENV_APPEND_PATH="/usr/local/arriba_v2.3.0:/usr/local/STAR-Fusion/ctat-genome-lib-builder:/usr/local/STAR-Fusion:/usr/local/RSEM" 
   
   memory { (100.GB + (16.GB * task.attempt)) } // First attempt 100GB, second 116GB, etc
   errorStrategy 'retry'
   maxRetries 3

   input:
   set SampleName, file(STAR_aligned_sorted_bam) from Adaptor_Arriba_Fusion_Ch
    
   output:
   set file("${SampleName}.arriba.fusions.tsv"), file("${SampleName}.arriba.fusions.discarded.tsv") into Arriba_Fusion_Out_Ch
   
   script:
   """
   export TMPDIR=${params.TMP_DIR}
   arriba \
    -a ${params.Arriba_Fusion.fa} \
    -g ${params.Arriba_Fusion.gtf} \
    -b ${params.Arriba_Fusion.blacklist} \
    -k ${params.Arriba_Fusion.known_fusion_k} \
    -t ${params.Arriba_Fusion.known_fusion_t} \
    -p ${params.Arriba_Fusion.protein_domains} \
    -x ${STAR_aligned_sorted_bam} \
    -o ${SampleName}.arriba.fusions.tsv \
    -O ${SampleName}.arriba.fusions.discarded.tsv
   """  
}


/***********************************************************************************************
* This Process is be responsible for RNAseq Quantification RSEM 
************************************************************************************************/
process RNAseq_RSEM {
   publishDir "${params.outdir}/${params.project}/RNAseq_RSEM/${SampleName}", mode: 'copy', overwrite: true
// export SINGULARITYENV_APPEND_PATH="/usr/local/arriba_v2.3.0:/usr/local/STAR-Fusion/ctat-genome-lib-builder:/usr/local/STAR-Fusion:/usr/local/RSEM" 
   
   input:
   set SampleName, Strandness, file(STAR_Aligned_toTranscriptome_Bam) from Adaptor_RNAseq_RSEM_Ch
   output:
   set SampleName, file("${SampleName}.RSEM.genes.results"), file("${SampleName}.RSEM.isoforms.results"), file("${SampleName}.RSEM.stat.tar.gz") into RNAseq_RSEM_Ch

   // --paired-end  \
   // def STRAND_TYPE=(Strandness=="Unstranded")? "none" : "reverse"

   script:
   def STRAND_TYPE=(Strandness =~ /Stranded_/)? "reverse" : "none"  
   """
   export TMPDIR=${params.TMP_DIR}
   rsem-calculate-expression \
     --num-threads ${params.RSEM.ThreadN} \
     --no-bam-output  \
     --alignments  \
     --paired-end  \
     --strandedness ${STRAND_TYPE} \
     ${STAR_Aligned_toTranscriptome_Bam} \
     ${params.RSEM.Ref} \
     ${SampleName}.RSEM

   tar zcvf ${SampleName}.RSEM.stat.tar.gz ${SampleName}.RSEM.stat
   """
}

/***********************************************************************************************
* The following Process modules take the RSEM's output as input and are responsible for Gene 
expression base ML-Classifier (RANK, tSNE, and RF|RBF|Polynomial ) 
************************************************************************************************/
RNAseq_RSEM_Ch.into{ RNAseq_RSEM_Ch1; RNAseq_RSEM_Ch2; RNAseq_RSEM_Ch3; RNAseq_RSEM_Ch4; RNAseq_RSEM_Ch5 }
process ExpectedCount_GenesResults_Aadptor {
  //  publishDir "${params.outdir}/${params.project}/ExpectedCount_GenesResults_Aadptor/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(RSEM_genes_results), file(RSEM_isoforms_results), file(RSEM_stat_tar_gz) from RNAseq_RSEM_Ch1
   output:
   file("${SampleName}") into ExpectedCount_GenesResult_Ch
 //file("${SampleName}.ExpectedCount") into ExpectedCount_GenesResult_Ch    
 
   //ln -s $RSEM_genes_results ${SampleName}.ExpectedCount
   script:
   """
   ln -s $RSEM_genes_results ${SampleName}
   """   
}

/****************************************************************************************************
* This process use grep to filter some unnecessary rows in the RSEM/genes.results, which can be used used by rsem-generate-data-matrix to generate a modified Count Matrix.   
****************************************************************************************************/
process ModifiedCount_GenesResults_Aadptor {
   input:
   set SampleName, file(RSEM_genes_results), file(RSEM_isoforms_results), file(RSEM_stat_tar_gz) from RNAseq_RSEM_Ch2
   output:
   file("${SampleName}.ModifiedCount") into ModifiedCount_GenesResult_Ch    

   script:
   """
   cat $RSEM_genes_results|grep -v '^ERCC'| grep -v 'gene_id' > ${SampleName}.ModifiedCount
   """  
 /*
   echo -e "__no_feature\t0\t0\t0\t16681862\t0\t0" >> ${SampleName}.ModifiedCount
   echo -e "__ambiguous\t0\t0\t0\t1639711\t0\t0" >> ${SampleName}.ModifiedCount
   echo -e "__too_low_aQual\t0\t0\t0\t0\t0\t0" >> ${SampleName}.ModifiedCount
   echo -e "__not_aligned\t0\t0\t0\t0\t0\t0" >> ${SampleName}.ModifiedCount
   echo -e "__alignment_not_unique\t0\t0\t0\t7068675\t0\t0" >> ${SampleName}.ModifiedCount
   echo -e "__DUX4\t0\t0\t0\t19\t0\t0" >> ${SampleName}.ModifiedCount
 */ 

}

/****************************************************************************************************
* This process adopt a trick to switch the column 5 and 6 of the RSEM/genes.results, which can be used used by rsem-generate-data-matrix to generate the corresponding TPM Matrix   
****************************************************************************************************/
process TPM_GenesResults_Aadptor {
  //  publishDir "${params.outdir}/${params.project}/TPM_GenesResults_Aadptor/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(RSEM_genes_results), file(RSEM_isoforms_results), file(RSEM_stat_tar_gz) from RNAseq_RSEM_Ch3
   output:
   file("${SampleName}") into TPM_GenesResult_Ch
// file("${SampleName}.TPM") into TPM_GenesResult_Ch    

// cat $RSEM_genes_results| awk 'BEGIN {OFS="\t"}; {print \$1, \$2, \$3, \$4, \$6, \$5, \$7}'> ${SampleName}.TPM
   script:
   """
   cat $RSEM_genes_results| awk 'BEGIN {OFS="\t"}; {print \$1, \$2, \$3, \$4, \$6, \$5, \$7}'> ${SampleName}
   """   
}

/****************************************************************************************************
* This process adopt a trick to switch the column 5 and 7 of the RSEM/genes.results, which can be used used by rsem-generate-data-matrix to generate the corresponding FPKM Matrix   
****************************************************************************************************/
process FPKM_GenesResults_Aadptor {
  //  publishDir "${params.outdir}/${params.project}/FPKM_GenesResults_Aadptor/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(RSEM_genes_results), file(RSEM_isoforms_results), file(RSEM_stat_tar_gz) from RNAseq_RSEM_Ch4
   output:
   file("${SampleName}") into FPKM_GenesResult_Ch
// file("${SampleName}.FPKM") into FPKM_GenesResult_Ch    

// cat $RSEM_genes_results| awk 'BEGIN {OFS="\t"}; {print \$1, \$2, \$3, \$4, \$7, \$6, \$5}'> ${SampleName}.FPKM  
   script:
   """
   cat $RSEM_genes_results| awk 'BEGIN {OFS="\t"}; {print \$1, \$2, \$3, \$4, \$7, \$6, \$5}'> ${SampleName}
   """   
}

/*****************************************************************************************************
* This process is responsible for generateing a two-colum reformatted count, which can be used in the following RNAseqCNV analysis 
*****************************************************************************************************/
process RSEM_Reformatted_Count {
   publishDir "${params.outdir}/${params.project}/RSEM_Reformatted_Count/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(RSEM_genes_results), file(RSEM_isoforms_results), file(RSEM_stat_tar_gz) from RNAseq_RSEM_Ch5
   
   output:
   set SampleName, file("${SampleName}.reformatted.counts") into RESM_Reformatted_Count_Ch

//cat $RSEM_genes_results| grep -v 'ERCC' |grep -v 'gene_id' |awk 'BEGIN {OFS="\t"}; {print \$1, \$5}'> ${SampleName}.reformatted.counts 
// cat $RSEM_genes_results| grep -v 'ERCC' |grep -v 'gene_id'  |awk '{printf "%s\\t%.0f\\n",\$1,\$5}'> ${SampleName}.reformatted.counts
// cat $RSEM_genes_results| grep -v 'ERCC' |grep -v 'gene_id' |awk 'BEGIN {OFS="\t"}; {print \$1, \$5}'| tr '.' '\t' |cut -f 1,3 > ${SampleName}.reformatted.counts
  
   script:
   """
   cat $RSEM_genes_results| grep -v 'ERCC' |grep -v 'gene_id'  |awk '{printf "%s\\t%.0f\\n",\$1,\$5}'|tr '.' '\\t' |cut -f 1,3> ${SampleName}.reformatted.counts
   """   
}

/*****************************************************************************************************
* This process is responsible for generating the Expected CountMatrix across samples, which can be used for GE-ML (Gene Expression based Machine Learning Classifier) 
******************************************************************************************************/
process Get_ExpectedCount_Matrix {
  publishDir "${params.outdir}/${params.project}/ExpectedCount_Matrix/", mode: 'copy', overwrite: true
  input:
  file(ExpectedCount_genes_results_list) from ExpectedCount_GenesResult_Ch.toSortedList() //.collect()  
  
  output:
  file("${params.project}.ExpectedCount.matrix") into ExpectedCount_Matrix_Ch
  
  script:
  """
  rsem-generate-data-matrix \
  ${ ExpectedCount_genes_results_list.collect { " $it" }.join()} >${params.project}.ExpectedCount.matrix
  """
}

/*****************************************************************************************************
* This process is responsible for generating the Modified CountMatrix across samples, which can be used for GE-ML (Gene Expression based Machine Learning Classifier) 
******************************************************************************************************/
process Get_ModifiedCount_Matrix {
  publishDir "${params.outdir}/${params.project}/ModifiedCount_Matrix/", mode: 'copy', overwrite: true
  input:
  file(ModifiedCount_genes_results_list) from ModifiedCount_GenesResult_Ch.toSortedList() //.collect()  
  
  output:
  file("${params.project}.ModifiedCount.matrix") into ModifiedCount_Matrix_Ch
  
  script:
  """
  rsem-generate-data-matrix \
  ${ ModifiedCount_genes_results_list.collect { " $it" }.join()} >${params.project}.ModifiedCount.matrix
  """
}

/*****************************************************************************************************
* This process is responsible for generating the TPM Matrix across samples, which can be used for GE-ML (Gene Expression based Machine Learning Classifier) 
******************************************************************************************************/
process Get_TPM_Matrix {
  publishDir "${params.outdir}/${params.project}/TPM_Matrix/", mode: 'copy', overwrite: true
  input:
  file(TPM_genes_results_list) from TPM_GenesResult_Ch.toSortedList()     //collect()  
  
  output:
  file("${params.project}.TPM.matrix") into TPM_Matrix_Ch
  
  script:
  """
  rsem-generate-data-matrix \
  ${ TPM_genes_results_list.collect { " $it" }.join()} > ${params.project}.TPM.matrix
  """
}

/*****************************************************************************************************
* This process is responsible for generating the FPKM Matrix across samples, which can be used for GE-ML (Gene Expression based Machine Learning Classifier) 
******************************************************************************************************/
process Get_FPKM_Matrix {
  publishDir "${params.outdir}/${params.project}/FPKM_Matrix/", mode: 'copy', overwrite: true
  input:
  file(FPKM_genes_results_list) from FPKM_GenesResult_Ch.toSortedList()    //collect()  
  
  output:
  file("${params.project}.FPKM.matrix") into FPKM_Matrix
  
  script:
  """
  rsem-generate-data-matrix \
  ${ FPKM_genes_results_list.collect { " $it" }.join()} > ${params.project}.FPKM.matrix
  """
}

ExpectedCount_Matrix_Ch.into { ExpectedCount_Matrix_Ch1; ExpectedCount_Matrix_Ch2 }
tSNE_Count_Matrix_Ch = (params.Select_ML_Classifier_tSNE == "N") ? Channel.empty() : ModifiedCount_Matrix_Ch // ExpectedCount_Matrix_Ch1
RF_Count_Matrix_Ch = (params.Select_ML_Classifier_RF == "N") ? Channel.empty() : ExpectedCount_Matrix_Ch2
RANK_TPM_Matrix_Ch = (params.Select_ML_Classifier_RANK == "N") ? Channel.empty() : TPM_Matrix_Ch 

/* This Process module is responsible for Gene-Expression Machine Learning Classifier -RANK*/
process ML_Classifier_RANK {
  publishDir "${params.outdir}/${params.project}/ML_Classifier_RANK/", mode: 'copy', overwrite: true
  input:
  file(TPM_Matrix) from RANK_TPM_Matrix_Ch
  val RANK_RDS from RANK_Training_RDS

  output:file("${params.project}*.*")
  // file("${params.project}_RANK_codingGene_classification.txt")
  // ${params.ML_Classifier.RANK_Training_RDS}
  
  script: 
  """
  export R_LIBS=${params.R_4_2_LIBS}
  Rscript ${params.ML_Classifier.RScript_Path}/ML_Classifier_Rank.R \
  ${params.project} \
  ${RANK_RDS} \
  $TPM_Matrix 
  """ 
}

/* This Process module is responsible for Gene-Expression Machine Learning Classifier -tSNE */
process ML_Classifier_tSNE {
  publishDir "${params.outdir}/${params.project}/ML_Classifier_tSNE/", mode: 'copy', overwrite: true
  input:
  file(Count_Matrix) from tSNE_Count_Matrix_Ch
  
  output:
  set file("${params.project}_top4_PCs*"), file("${params.project}_tSNE_2D.*"), file("${params.project}_tSNE_3D.*"), file("tsne_color_map.csv"), file("images.tar.gz") into ML_Classifier_tSNE_Ch 
  
  //template 'ML_Classifier_tSNE.sh' ML_Classifier_tSNE.R generate_tSNE.R
  //def cmd = (workflow.profile != "singularity") ? "Rscript" : ""
  // ${params.ML_Classifier.tSNE.LibraryType} \
  script: 
  """
  export R_LIBS=${params.R_4_2_LIBS}
  Rscript ${params.ML_Classifier.RScript_Path}/ML_Classifier_tSNE.R \
  ${params.project} \
  ${Count_Matrix} \
  ${params.ML_Classifier.tSNE.Metadata} \
  ${params.ML_Classifier.tSNE.ReferenceCountsFile} \
  ${params.ML_Classifier.tSNE.ReferenceMetaData} \
  ${params.ML_Classifier.tSNE.Top1kGenesFile}  

  tar zcvf images.tar.gz images
  """
}

/* This Process module is responsible for Gene-Expression Machine Learning Classifier -RF(Random Forest) */
process ML_Classifier_RF {
  publishDir "${params.outdir}/${params.project}/ML_Classifier_RF/", mode: 'copy', overwrite: true
  input:
  file(Count_Matrix) from RF_Count_Matrix_Ch
  val RF_RDS from RF_Training_RDS  

  output:
  file("${params.project}_RF_predictions.tsv")
  
  //def cmd = (workflow.profile == "singularity") ? "Rscript" : ""
  // echo "ML_Classifier_RF"  ${params.ML_Classifier.RF_Training_RDS} \

  script:
  """
  export R_LIBS=${params.R_4_2_LIBS}
  Rscript ${params.ML_Classifier.RScript_Path}/ML_Classifier_RF.R \
  ${params.project} \
  ${RF_RDS} \
  $Count_Matrix
  """ 
}


/***********************************************************************************************
* The following Process modules are responsible for GATK based RNA Variant Calling
************************************************************************************************/
/***********************************************************************************************
* This Process, as the 1st process of GATK Variant Calling, is responsible for SplitNCigarReads 
************************************************************************************************/
process GATK_SplitNCigarReads {
   publishDir "${params.outdir}/${params.project}/GATK_SplitNCigarReads", mode: 'copy', overwrite: true

   input:
   set SampleName, file(STAR_Aligned_out_bam), file(STAR_Aligned_out_bam_bai) from Adaptor_RNAvar_GATK_Ch

   output:
   set SampleName, file("${SampleName}.splitN.bam"), file("${SampleName}.splitN.bai") into SplitNCigarRead_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk SplitNCigarReads \
    --java-options "-Xms8g" \
    --reference ${params.GATK.Ref} \
    --input ${STAR_Aligned_out_bam} \
    --output ${SampleName}.splitN.bam
   """
}

/***********************************************************************************************
* This Process, as the 2nd process of GATK Variant Calling, is responsible for BaseRecalibration 
************************************************************************************************/
process GATK_BaseRecalibrator {
   publishDir "${params.outdir}/${params.project}/GATK_BaseRecalibrator", mode: 'copy', overwrite: true

   input:
   set SampleName, file(splitN_bam), file(splitN_bam_bai) from SplitNCigarRead_Ch

   output:
   set SampleName, file("${SampleName}.recal_data.table"), file("${SampleName}_cached.splitN.bam"), file("${SampleName}_cached.splitN.bai") into BaseRecalibrator_Ch   

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk BaseRecalibrator \
    --java-options "-Xms32g" \
    --use-original-qualities \
    --reference ${params.GATK.Ref} \
    --known-sites ${params.GATK.KNOWN_SITE1} \
    --known-sites ${params.GATK.KNOWN_SITE2} \
    --known-sites ${params.GATK.KNOWN_SITE3}  \
    --input $splitN_bam \
    --output ${SampleName}.recal_data.table

    ln -s $splitN_bam ${SampleName}_cached.splitN.bam
    ln -s ${splitN_bam_bai} ${SampleName}_cached.splitN.bai
   """
}

/***********************************************************************************************
* This Process, as the 3nd process of GATK Variant Calling, is responsible for Applying the recalibrated data table to BQSR. 
************************************************************************************************/
process GATK_ApplyBQSR {
   publishDir "${params.outdir}/${params.project}/GATK_ApplyBQSR", mode: 'copy', overwrite: true

   input:
   set SampleName, file(recal_data_table), file(cached_splitN_bam), file(cached_splitN_bam_bai) from BaseRecalibrator_Ch
   output:
   set SampleName, file("${SampleName}.recal.bam"), file("${SampleName}.recal.bai") into ApplyBQSR_Ch   

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk ApplyBQSR \
    --java-options "-Xms32g" \
    --use-original-qualities \
    --add-output-sam-program-record \
    --reference ${params.GATK.Ref} \
    --bqsr-recal-file ${recal_data_table} \
    --input ${cached_splitN_bam} \
    --output ${SampleName}.recal.bam
   """
}

ApplyBQSR_Ch.into {ApplyBQSR_Ch1; ApplyBQSR_Ch2; ApplyBQSR_Ch3}

/***********************************************************************************************
* This Process is responsible for using Vardict  for Variant Calling.
************************************************************************************************/
VarDict_InCh= (params.Select_VarDict=="N")? Channel.empty() : ApplyBQSR_Ch3
process VarDict {
   publishDir "${params.outdir}/${params.project}/VarDict", mode: 'copy', overwrite: true
   input:
   set SampleName, file(BQSR_recal_bam), file(BQSR_recal_bai) from VarDict_InCh
   output:
   set SampleName, file("${SampleName}.VarDict.vcf.gz"), file("${SampleName}.VarDict.vcf.gz.tbi") into VarDict_Ch
   
   // ${VARDICT_DIR}/vardict-java -f 0.05 -c 1 -S 2 -E 3 -g 4 -r 2 -t -th $THREADS -v -G ${REFERENCE} -b $BAM $CALLBED | ${VARDICT_DIR}/teststrandbias.R | ${VARDICT_DIR}/var2vcf_valid.pl -N ${SAMPLE} -E -f 0.02 > ${SAMPLE}.${CALLER}.vcf

   script:
   """
    export R_LIBS=${params.R_4_1_LIBS}

    vardict \
    -f 0.05 \
    -c 1 \
    -S 2 \
    -E 3 \
    -g 4 \
    -r 2 \
    -t \
    -th ${params.VarDict.ThreadN} \
    -v \
    -G ${params.VarDict.REFERENCE} \
    -b ${BQSR_recal_bam} \
    ${params.VarDict.CALLBED} | \
    teststrandbias.R | \
    var2vcf_valid.pl -N ${SampleName} \
    -E \
    -f 0.02 > ${SampleName}.VarDict.vcf
   
    bgzip \
    ${SampleName}.VarDict.vcf
 
    tabix \
    ${SampleName}.VarDict.vcf.gz
   """    
}

/***********************************************************************************************
* This Process, as the 4th process of GATK Variant Calling and also the key module for Haplotype Variant Calling.
************************************************************************************************/
GATK_HaplotypeCaller_InCh = (params.Select_Haplo_PerChrom=="Y")? Channel.empty() : ApplyBQSR_Ch1
process GATK_HaplotypeCaller {
   publishDir "${params.outdir}/${params.project}/GATK_HaplotypeCaller", mode: 'copy', overwrite: true

   memory { (64.GB + (32.GB * task.attempt)) } // First attempt 64GB, second 96GB, etc
   errorStrategy 'retry'
   maxRetries 3

   input:
   set SampleName, file(BQSR_recal_bam), file(BQSR_recal_bai) from GATK_HaplotypeCaller_InCh
   output:
   set SampleName, file("${SampleName}.GATK.vcf.gz"), file("${SampleName}.GATK.vcf.gz.tbi")  into HaplotypeCaller_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk HaplotypeCaller \
    --java-options "-Xms32g" \
    --annotation-group StandardAnnotation \
    --annotation-group StandardHCAnnotation \
    --standard-min-confidence-threshold-for-calling 20 \
    --dont-use-soft-clipped-bases \
    --reference ${params.GATK.Ref} \
    --input ${BQSR_recal_bam} \
    --output ${SampleName}.GATK.vcf.gz
   """
}

/************************************************************************************************************
*  This process is 4th process of GATK Variant Calling, but it is a speed up module, by which the haplotype calling can be executed chrom by chrom.
*  The chrom level parallel structure can greatly speed up the time-consuming module-HaplotypeCaller
************************************************************************************************************/

HaplotypeCaller_perChrom_inch = (params.Select_Haplo_PerChrom=="Y")? ApplyBQSR_Ch2.spread(chromosomes_ch) : Channel.empty()
process GATK_HaplotypeCaller_perChrom {
     
    input:   
    set SampleName, file(BQSR_recal_bam), file(BQSR_recal_bai), Chrosome_Interval from HaplotypeCaller_perChrom_inch
  
	output:
    set SampleName, file("${SampleName}.GATK.${Chrosome_Interval}.vcf.gz") , file("${SampleName}.GATK.${Chrosome_Interval}.vcf.gz.tbi") into Chrosome_Interval_vcf_Ch

    script:
    """  
    export TMPDIR=${params.TMP_DIR}
    gatk HaplotypeCaller \
    --java-options "-Xms32g" \
    --annotation-group StandardAnnotation \
    --annotation-group StandardHCAnnotation \
    --standard-min-confidence-threshold-for-calling 20 \
    --dont-use-soft-clipped-bases \
    --reference ${params.GATK.Ref} \
    --input ${BQSR_recal_bam} \
    --intervals ${Chrosome_Interval} \
    --output ${SampleName}.GATK.${Chrosome_Interval}.vcf.gz 
    
    """
} 

/************************************************************************************************************
/* This Process is be responsible for Merging the interval (eg. chrosome level) VCFs into sample level VCF 
*************************************************************************************************************/
process MergeVCF {
    publishDir "${params.outdir}/${params.project}/MergeVCF", mode: 'copy', overwrite: true
	
	input:
	set SampleName, file (vcfs), file (vcftbis) from Chrosome_Interval_vcf_Ch.groupTuple()

	output:
    set SampleName, file("${SampleName}.GATK.vcf.gz") , file("${SampleName}.GATK.vcf.gz.tbi") into HaplotypeCaller_merge_vcf_Ch

    script:    
    """
    export TMPDIR=${params.TMP_DIR}
    gatk --java-options '-Xmx64g' \
    MergeVcfs \
    ${vcfs.collect{"--INPUT $it " }.join()} \
    -O ${SampleName}.GATK.vcf.gz
    """
} 

/***********************************************************************************************
* This Process, as the last process module of GATK Variant Calling , is responsible for filtering the  called Variant from Haplotype_Caller.
************************************************************************************************/
GATK_VariantFiltration_InCh=(params.Select_Haplo_PerChrom=="Y")? HaplotypeCaller_merge_vcf_Ch : HaplotypeCaller_Ch
process GATK_VariantFiltration {
   publishDir "${params.outdir}/${params.project}/GATK_VariantFiltration", mode: 'copy', overwrite: true

   input:
   set SampleName, file(GATK_vcf), file(GATK_vcf_tbi) from GATK_VariantFiltration_InCh
   output:
   set SampleName, file("${SampleName}.GATK.hardfiltered.vcf"), file("${SampleName}.GATK.hardfiltered.vcf.idx") into RNAvar_GATK_Ch
 
   script:
   """
   export TMPDIR=${params.TMP_DIR}
   gatk VariantFiltration \
    --java-options "-Xms32g" \
    --cluster-window-size 35 \
    --cluster-size 3 \
    --filter-name "FS" \
    --filter-expression "FS > 30.0" \
    --filter-name "QD" \
    --filter-expression "QD < 2.0" \
    --reference ${params.GATK.Ref} \
    --variant ${GATK_vcf} \
    --output ${SampleName}.GATK.hardfiltered.vcf
   """
}


RNAvar_GATK_Ch.into {RNAvar_GATK_Ch1; RNAvar_GATK_Ch2; RNAvar_GATK_Ch3}
/********************************************************************************************
*  This process is responsible for annotating the RNA-seq SNV varaints called by GATK
*********************************************************************************************/
Variant_ANNOVAR_In_Ch= (params.Select_Variant_ANNOVAR == "N") ? Channel.empty() :RNAvar_GATK_Ch1
process Variant_ANNOVAR {
   publishDir "${params.outdir}/${params.project}/Variant_ANNOVAR/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(GATK_vcf), file(GATK_vcf_idx) from Variant_ANNOVAR_In_Ch 
   
   output:
   set file("${SampleName}.vcf.ann.avinput"), file("${SampleName}.*.annovar.*.tab"), file("${SampleName}.vcf.ann.*.txt"), file("${SampleName}.vcf.ann.*.vcf.gz"), file("${SampleName}.vcf.ann.*.vcf.gz.tbi") into Variant_ANNOVAR_Ch

   script:
   """
   export TMPDIR=${params.TMP_DIR}
   ${params.ANNOVAR.ANNOVAR_PATH}/table_annovar.pl \
   ${GATK_vcf} \
   ${params.ANNOVAR.ANNOVAR_DB} \
   -buildver ${params.ANNOVAR.REF_BUILD} \
   -out ${SampleName}.vcf.ann -remove -protocol refGene,avsnp150,1000g2015aug_all,exac03,exac03nontcga,esp6500siv2_all,gnomad211_exome,gnomad30_genome,dbnsfp35a,revel,intervar_20180118,dbscsnv11,cosmic70,nci60,clinvar_20180603 -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
   
   bgzip ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf
   tabix ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf.gz
   
   gatk VariantsToTable \
   -V ${SampleName}.vcf.ann.${params.ANNOVAR.REF_BUILD}_multianno.vcf.gz \
   -O ${SampleName}.GATK.annovar.${params.ANNOVAR.REF_BUILD}.tab \
   ${params.ANNOVAR.CONFIG_STANDARD} ${params.ANNOVAR.CONFIG_INFO} ${params.ANNOVAR.CONFIG_FORMAT} ${params.ANNOVAR.CONFIG_ANNOVAR} 
   """
}


/**********************************************************************************************
* This process is responsible for RNAseqCNV calling, which is based on the two-column reformatted count file from RESM_Reformatted_Count_Ch and hard-filttered vcf from RNAvar_GATK_Ch2
**********************************************************************************************/
RNASeqCNV_Ch_In= (params.Select_RNASeqCNV == "N") ? Channel.empty() : RESM_Reformatted_Count_Ch.join(RNAvar_GATK_Ch2)
process RNASeqCNV {
   publishDir "${params.outdir}/${params.project}/RNASeqCNV/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(reformatted_count), file(GATK_hardfiltered_vcf), file(GATK_hardfiltered_vcf_idx) from RNASeqCNV_Ch_In
  
   output:
   set SampleName, file("alteration_matrix.tsv"), file("estimation_table.tsv"), file("manual_an_table.tsv"), file("${SampleName}.tar.gz") into RNASeqCNV_Ch  
   
   script:
   """
   echo -e out_dir=\\"./\\"> config
   echo -e count_dir=\\"./\\" >>config
   echo -e snv_dir=\\"./\\" >>config
   
   echo -e "${SampleName}\\t${reformatted_count}\\t${GATK_hardfiltered_vcf}"> metadata 
   
   export R_LIBS=${params.R_4_1_LIBS}
   Rscript ${params.RNASeq_CNV.RScript_Path}/RNASeqCNV.R >RNASeq_CNV.log 2>RNASeq_CNV.err
   tar zcvf ${SampleName}.tar.gz ${SampleName}

   """
}


/**********************************************************************************************
* This process is responsible for RNAIndel calling, which is based on the STAR mapping bam and hard-filttered vcf from RNAvar_GATK_Ch3
**********************************************************************************************/
RNAIndel_Ch_In= (params.Select_RNAIndel == "N") ? Channel.empty() : Adaptor_RNAIndel_Ch.join(RNAvar_GATK_Ch3)
process RNAIndel {
  publishDir "${params.outdir}/${params.project}/RNAIndel/${SampleName}", mode: 'copy', overwrite: true

  input:
  set SampleName, file(STAR_marked_dup_bam), file(STAR_marked_dup_bam_bai), file(GATK_hardfiltered_vcf), file(GATK_hardfiltered_vcf_idx) from RNAIndel_Ch_In
  
  output:
  set SampleName, file("${SampleName}_RNAIndel_gatk.vcf.gz"), file("${SampleName}_RNAIndel_gatk.vcf.gz.tbi"), file("${SampleName}_RNAIndel_gatk.tab") into RNAIndel_Ch

  // -p ${params.RNAIndel.ThreadN} \  
  script:
  """
  export TMPDIR=${params.TMP_DIR}

  bgzip ${GATK_hardfiltered_vcf}
  tabix ${GATK_hardfiltered_vcf}.gz
  
  rnaindel PredictIndels \
  -i ${STAR_marked_dup_bam} \
  -o ${SampleName}_RNAIndel_gatk.vcf \
  -r ${params.RNAIndel.Ref_Fa} \
  -m 16000m \
  -d ${params.RNAIndel.Database_model} \
  -v ${GATK_hardfiltered_vcf}.gz \
  -p ${params.RNAIndel.ThreadN}

  gatk VariantsToTable \
   -V ${SampleName}_RNAIndel_gatk.vcf.gz \
   -O ${SampleName}_RNAIndel_gatk.tab \
   ${params.RNAIndel.CONFIG_STANDARD} ${params.RNAIndel.CONFIG_INFO} ${params.RNAIndel.CONFIG_FORMAT} 
  """
}


/**********************************************************************************************
* This process is responsible for Pindel, which is based on the aligned bam for detecting breakpoints of large deletions, mdeium size insertions,
inversions, tandem duplications and other SVs at signle-based resolution. 
**********************************************************************************************/
process Pindel {
   publishDir "${params.outdir}/${params.project}/Pindel/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(STAR_Aligned_out_bam), file(STAR_Aligned_out_bam_bai) from Adaptor_Pindel_Ch

   output:
   set SampleName, file("${SampleName}_Pindel.out_BP"), file("${SampleName}_Pindel.out_D"), file("${SampleName}_Pindel.out_INV"), file("${SampleName}_Pindel.out_LI"), file("${SampleName}_Pindel.out_SI"), file("${SampleName}_Pindel.out_TD"), file("${SampleName}_Pindel.out_RP") into Pindel_Ch
   
   script:
   """
   export TMPDIR=${params.TMP_DIR}
   echo -e "${STAR_Aligned_out_bam}\\t${params.Pindel.InsertSize}\\t${SampleName} "> ${SampleName}.PindelConfig
   pindel \
   -f ${params.Pindel.Ref} \
   -i ${SampleName}.PindelConfig \
   --number_of_threads ${params.Pindel.ThreadN} \
   --min_num_matched_bases 3 \
   --report_inversions TRUE \
   --min_inversion_size 20 \
   --report_duplications TRUE \
   --report_long_insertions TRUE \
   --report_breakpoints TRUE \
   --report_interchromosomal_events TRUE \
   --report_breakpoints \
   --minimum_support_for_event 2 \
   --maximum_allowed_mismatch_rate 0.3 \
   --sensitivity 0.99 \
   --min_distance_to_the_end 5 \
   --include ${params.Pindel.ITD_BED} \
   -o ${SampleName}_Pindel.out  
   """
}

/**********************************************************************************************
* This process is responsible for converting the Pindel output into VCF format
**********************************************************************************************/
process Pindel2VCF {
    
   publishDir "${params.outdir}/${params.project}/Pindel2VCF/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(Pindel_out_BP), file(Pindel_out_D), file(Pindel_out_INV), file(Pindel_out_LI), file(Pindel_out_SI), file(Pindel_out_TD), file(Pindel_out_RP) from Pindel_Ch

   output:
   set SampleName, file("${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf.gz"), file("${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf.gz.tbi")  into Pindel2VCF_Ch
   
   script:
   """
   export TMPDIR=${params.TMP_DIR}
   pindel2vcf \
   -P ${SampleName}_Pindel.out \
   -r ${params.Pindel.Ref} \
   -R ${params.Pindel.REF_BUILD} \
   -d ${params.Pindel.REF_Date} \
   --min_coverage ${params.Pindel.Min_Coverage}  \
   --het_cutoff ${params.Pindel.Het_Cutoff} \
   --hom_cutoff ${params.Pindel.Hom_Cutoff} \
   --vcf ${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf

   bgzip -f ${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf
   tabix -f -p vcf ${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf.gz
   bcftools view -Ov --exclude-uncalled --min-ac=1 ${SampleName}.Pindel.${params.Pindel.REF_BUILD}.raw.vcf.gz 
   """
    
}

/*************************************************************************************************************************
* This process iAdmix is responsible for estimating admixture coefficients, which is based on the aligned bam and population 
* allele frequencies for common SNPs. The output is admixture coefficients for each reference population
*
* // --path ${which ANCESTRY|sed -e 's/\/ANCESTRY//g'}, `which ANCESTRY|sed -e 's/\/ANCESTRY//g'`
* --path "/hpcf/authorized_apps/rhel7_apps/iadmix/vendor/iAdmix"
***************************************************************************************************************************/
process iAdmix {
   publishDir "${params.outdir}/${params.project}/iAdmix/${SampleName}", mode: 'copy', overwrite: true
   input:
   set SampleName, file(STAR_Aligned_out_bam) from Adaptor_iAdmix_Ch

   output:
   set SampleName, file("${SampleName}.ancestry.input"), file("${SampleName}.ancestry.out"), file("${SampleName}.forGLL"), file("${SampleName}.GLL") into iAdmix_Ch

   //ANCESTRY_PATH="\${which ANCESTRY|sed -e 's/\/ANCESTRY//g'}"
      
   script:  
   """   
   export TMPDIR=${params.TMP_DIR}
   runancestry.py \
   -f ${params.iAdmix.AFFILE} \
   --bam ${STAR_Aligned_out_bam} \
   --addchr TRUE -o ${SampleName} \
   --path ${params.iAdmix.iAdmix_PATH}
   """
}

/* process iAdmix {
   publishDir "${params.outdir}/${params.project}/iAdmix/${SampleName}", mode: 'copy', overwrite: true

   input:
   set SampleName, file(STAR_Aligned_out_bam) from Adaptor_iAdmix_Ch

   output:
   set SampleName, file("${SampleName}.ancestry.input"), file("${SampleName}.ancestry.out"), file("${SampleName}.forGLL"), (file("${SampleName}.GLL") into iAdmix_Ch
   
   script:
   """   
   runancestry.py \
   -f ${params.iAdmix.AFFILE} \
   --bam ${STAR_Aligned_out_bam} \
   --addchr TRUE \
   -o ${SampleName} \
   --path ${params.iAdmix.iAdmix_PATH} 
   """
}
*/

// Local Variables:
// mode: groovy
// End:
