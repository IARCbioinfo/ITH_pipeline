#!/usr/bin/env nextflow

// Copyright (C) 2018 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.help = null
params.config	= null
params.cpu  = 2
params.cpu2 = 2
params.mem  = 8
params.samtools_folder = "/usr/bin/"
params.bcftools_folder = "/usr/bin/"
params.bnpy_folder = "/usr/bin/"
params.hatchet_folder = "/usr/bin/"
params.bam_folder = null
params.ref = null
params.output_folder = "results"
params.bed = "NO_BED"
params.known_snps = "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz"


params.cluBB_d  = 0.08
params.cluBB_tR_min = 0.1
params.cluBB_tR_max = 1
params.cluBB_tB_min = 0.02
params.cluBB_tB_max = 0.2

log.info ""
log.info "----------------------------------------------------------------"
log.info "            Intratumoral heterogeneity with HATCHet             "
log.info "----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "------------------- ITH_pipeline----------------------------"
    log.info ""
    log.info "nextflow run script.nf --input_folder path/to/input/ --ref path/to/ref/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--bam_folder         PATH        Folder containing bam files"
    log.info "--ref                  PATH        WHOLE Path to reference fasta file (should be indexed)"
    log.info "--correspondance       FILE        File containing correspondence between path to normal and path to tumor bams for each patient  "
    log.info ""
    log.info "Optional arguments:"
    log.info "--cpu                  INTEGER     Number of cpu to use for pileup, RDR and BAF computation (default=2)"
    log.info "--cpu                  INTEGER     Number of cpu to use for downstream computations (default=2)"
    //log.info "--config               FILE        Use custom configuration file"
    log.info "--known_snps           FILE        VCF with known SNPs for calling germline heterozygous (default=https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-All.vcf.gz)"
    log.info "--bed                  FILE        Restrict to genomic regions in bed file (default=null)"
    log.info "--mem                  INTEGER     Size of memory used in GB (default=8)"
    log.info "--output_folder				 PATH				 Path to output folder (default=.)"
    log.info "--samtools_folder      PATH        samtools installation dir (default=/usr/bin/)"
    log.info "--bcftools_folder      PATH        bcftools installation dir (default=/usr/bin/)"
    log.info "--hatchet_folder       PATH        hatchet installation dir (default=/usr/bin/)"
    log.info "--bnpy_folder          PATH        bnpy-dev installation dir (default=/usr/bin/)"
    log.info "--cluBB_d              NUMERIC     d parameter for cluBB (default=0.08)"
    log.info "--cluBB_tB             NUMERIC     tB parameter for cluBB (default=0.15)"
    log.info "--cluBB_tR             NUMERIC     tR parameter for cluBB (default=0.03)"
    log.info ""
    log.info "Flags:"
    log.info "--help                             Display this message"
    log.info ""
    exit 0
}

assert (params.bam_folder != null) : "please provide the --bam_folder option"
assert (params.ref != null) : "please provide the --ref option"
assert (params.correspondance != null) : "please provide the --correspondance option"

correspondance = file(params.correspondance)
bams = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
                  .map{ row -> [ row.sample , file(params.bam_folder + "/" + row.tumor), file(params.bam_folder + "/" + row.tumor+'.bai'),
                                 file(params.bam_folder + "/" + row.normal), file(params.bam_folder + "/" + row.normal+'.bai') ] }

tn_bambai = bams.groupTuple(by: 0)
                  .map { row -> tuple(row[0] , row[1], row[2] , row[3][0] , row[4][0]  ) }


tn_bambai.into{ tn_bambai_4binBAM ; tn_bambai_4bafBAM ; tn_bambai_4SNPcallingBAM}

bed  = file(params.bed)
ref  = file(params.ref)
ref_fai  = file(params.ref+".fai")
dict = file(params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict')

process binBAM {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bin/", mode: 'copy'

     input:
     set val(sampleID), file(bamT), file(baiT), file(bamN), file(baiN) from tn_bambai_4binBAM
     file bed
     file ref
     file ref_fai
     file dict

     output:
     set val(sampleID), file('*_RDRnormal.1bed'), file('*_RDRtumor.1bed'), file('*total.tsv') into bins
     file('*.log') into bin_logs

     shell :
      if( bed.name=='NO_BED' ){
        opt = ""
     }else{
        opt = "-r ${bed}"
     }
     '''
     ALLNAMES="!{bamN} !{bamT}"
     ALLNAMES="${ALLNAMES//.bam/}"
     python3 -m hatchet binBAM -N !{bamN} -T !{bamT} -S ${ALLNAMES} -b 50kb -g !{ref} -j !{params.cpu} -q 20 !{opt}\
                               -O !{sampleID}_RDRnormal.1bed -o !{sampleID}_RDRtumor.1bed -t !{sampleID}_total.tsv -v &> !{sampleID}_bins.log
     '''
}

process SNPcallingBAM {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bin/", mode: 'copy'

     input:
     set val(sampleID), file(bamT), file(baiT), file(bamN), file(baiN) from tn_bambai_4SNPcallingBAM
     file ref
     file ref_fai
     file dict

     output:
     set val(sampleID), file('*.vcf.gz') into SNPcalls
     file('*.log') into SNPcalling_logs

     shell :
      if( bed.name=='NO_BED' ){
        opt = ""
     }else{
        opt = "-r ${bed}"
     }
     '''
     python3 -m hatchet SNPCaller -N !{bamN} -r !{ref} -j !{params.cpu} -c 8 -C 300 -q 20 -Q 20 \
                                  -R !{params.known_snps} -o . |& tee !{sampleID}_SNPcalling.log
     '''
}

bamSNPs = tn_bambai_4bafBAM.join(SNPcalls)

process deBAF {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/baf/", mode: 'copy'

     input:
     set val(sampleID), file(bamT), file(baiT), file(bamN), file(baiN), file(SNPs) from bamSNPs
     file ref
     file ref_fai
     file dict

     output:
     set val(sampleID), file('*_BAFnormal.1bed'), file('*_BAFtumor.1bed') into bafs
     file('*.log') into baf_logs

     shell :
     '''
     ALLNAMES="!{bamN} !{bamT}"
     ALLNAMES="${ALLNAMES//.bam/}"
     python3 -m hatchet deBAF -N !{bamN} -T !{bamT} -S ${ALLNAMES} -r !{ref} -j !{params.cpu} -q 20 -Q 20 -U 20 -c 8 -C 300 \
                              -L *.vcf.gz -O !{sampleID}_BAFnormal.1bed -o !{sampleID}_BAFtumor.1bed -v &> !{sampleID}_bafs.log
     '''
}
binbafs = bins.join(bafs)

process comBBo {
     cpus params.cpu2
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bb/", mode: 'copy'

     input:
     set val(sampleID), file(bin_normal), file(bin_tumor), file(RDRtotal), file(baf_normal), file(baf_tumor) from binbafs

     output:
     set val(sampleID), file('*.bb') into bbs

     shell :
     '''
     python3 -m hatchet comBBo -c !{bin_normal} -C !{bin_tumor} -B !{baf_tumor} -t !{RDRtotal} > !{sampleID}_bulk.bb
     '''
}

tR_list = [0.1,0.15,0.2,0.3,0.4,0.5]//(1..10)*0.1//params.cluBB_tR_min.step(params.cluBB_tR_max, 0.1 )
tB_list = [0.02,0.04,0.06,0.1,0.15]//(1..10)*0.02//params.cluBB_tB_min.step(params.cluBB_tB_max, 0.05)

process cluBB {
     cpus params.cpu2
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bbc/", mode: 'copy'

     input:
     set val(sampleID), file(bb) from bbs
     each tB from tB_list
     each tR from tR_list

     output:
     set val(sampleID), file('*.seg'), file('*.bbc'), val(tB) , val(tR) into bbcs

     shell :
     '''
     python3 -m hatchet cluBB !{bb} -o !{sampleID}_tR!{tR}_tB!{tB}_tumor.seg -O !{sampleID}_tR!{tR}_tB!{tB}_tumor.bbc -tB !{tB} -tR !{tR} -d !{params.cluBB_d}
     '''
}
bbcs.into{bbcs4plots; bbcs4hatchet }

process plots {
     //cpus params.cpu2
     //memory params.mem+'G'
     //tag { sampleID }
     publishDir "${params.output_folder}/plots/${sampleID}", mode: 'copy'

     input:
     set sampleID, file(seg), file(bbc), val(tB) , val(tR) from bbcs4plots

     output:
     set sampleID, file("*.pdf"), file("*.png") into plots

     shell :
     '''
     python3 -m hatchet BBot -c RD --figsize 6,3 !{bbc} &
     python3 -m hatchet BBot -c CRD --figsize 6,3 !{bbc} &
     python3 -m hatchet BBot -c BAF --figsize 6,3 !{bbc} &
     python3 -m hatchet BBot -c BB !{bbc} &
     python3 -m hatchet BBot -c CBB !{bbc} &

     wait
     for f in *.pdf ; do mv -- $f "!{sampleID}_tB!{tB}_tR!{tR}_$f" ; done
     for f in *.png ; do mv -- $f "!{sampleID}_tB!{tB}_tR!{tR}_$f" ; done
     '''
}

process hatchet {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/results/", mode: 'copy'

     input:
     set val(sampleID), file(seg), file(bbc), val(tB) , val(tR) from bbcs4hatchet

     output:
     set val(sampleID), file('*ucn*') into ucns
     file '*.log' into logs_hatchet

     shell :
     prefix = bbc.baseName
     '''
     python3 -m hatchet solve -i !{prefix} -n2,10 -p 1000 -u 0.02 -j !{params.cpu} -eD 6 -eT 12 -g 0.35 \
                          -l 0.6 &> >(tee >(grep -v Progress > !{sampleID}_hatchet.log))
     mkdir !{sampleID}
     for f in *ucn* ; do mv -- $f "!{sampleID}/!{sampleID}_tB!{tB}_tR!{tR}_$f" ; done
     '''
}

process plot_final {
     cpus params.cpu2
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/plots/", mode: 'copy'

     input:
     set val(sampleID), file(ucn) from ucns

     output:
     file('*.pdf')

     shell :
     '''
     python3 -m hatchet BBeval !{sampleID}_best.bbc.ucn -u 0.02
     mkdir !{sampleID}
     for f in *.pdf ; do mv -- $f "!{sampleID}/!{sampleID}_$f" ; done
     '''
}
