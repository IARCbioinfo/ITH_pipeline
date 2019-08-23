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
params.cpu = 1
params.mem = 4
params.samtools_folder = "/usr/bin/"
params.bcftools_folder = "/usr/bin/"
params.bnpy_folder = "/usr/bin/"
params.hatchet_folder = "/usr/bin/"
params.bam_folder = null
params.ref = null
params.output_folder = "results"

params.cluBB_d = 0.08

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
    log.info "--cpu                  INTEGER     Number of cpu to use (default=1)"
    log.info "--config               FILE        Use custom configuration file"
    log.info "--mem                  INTEGER     Size of memory used in GB (default=4)"
    log.info "--output_folder				 PATH				 Path to output folder (default=.)"
    log.info "--samtools_folder      PATH        samtools installation dir (default=/usr/bin/)"
    log.info "--bcftools_folder      PATH        bcftools installation dir (default=/usr/bin/)"
    log.info "--hatchet_folder       PATH        hatchet installation dir (default=/usr/bin/)"
    log.info "--bnpy_folder          PATH        bnpy-dev installation dir (default=/usr/bin/)"
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


tn_bambai.into{ tn_bambai_4binBAM ; tn_bambai_4bafBAM }

process binBAM {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bin/", mode: 'copy'

     input:
     set val(sampleID), file(bamT), file(baiT), file(bamN), file(baiN) from tn_bambai_4binBAM

     output:
     set val(sampleID), file('*_normal.bin'), file('*_bulk.bin') into bins
     file('*.log') into bin_logs
     file('*total.bin')

     shell :
     '''
     ALLNAMES="!{bamN} !{bamT}"
     ALLNAMES="${ALLNAMES//.bam/}"
     python2 !{params.hatchet_folder}/utils/binBAM.py -N !{bamN} -T !{bamT} -S ${ALLNAMES} -b 50kb -g hg38 -j !{params.cpu} -q 20 -O !{sampleID}_normal.bin -o !{sampleID}_bulk.bin -t !{sampleID}_total.bin -v &> !{sampleID}_bins.log
     '''
}

process deBAF {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/baf/", mode: 'copy'

     input:
     set val(sampleID), file(bamT), file(baiT), file(bamN), file(baiN) from tn_bambai_4bafBAM

     output:
     set val(sampleID), file('*_normal.baf'), file('*_bulk.baf') into bafs
     file('*.log') into baf_logs

     shell :
     '''
     ALLNAMES="!{bamN} !{bamT}"
     ALLNAMES="${ALLNAMES//.bam/}"
     python2 !{params.hatchet_folder}/utils/deBAF.py -N !{bamN} -T !{bamT} -S ${ALLNAMES} -r !{params.ref} -j !{params.cpu} -q 20 -Q 20 -U 20 -c 4 -C 300 -O !{sampleID}_normal.baf -o !{sampleID}_bulk.baf -v &> !{sampleID}_bafs.log
     '''
}
binbafs = bins.join(bafs)

process comBBo {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bb/", mode: 'copy'

     input:
     set val(sampleID), file(bin_normal), file(bin_bulk), file(baf_normal), file(baf_bulk) from binbafs

     output:
     set val(sampleID), file('*.bb') into bbs

     shell :
     '''
     python2 !{params.hatchet_folder}/utils/comBBo.py -c !{bin_normal} -C !{bin_bulk} -B !{baf_bulk} -m MIRROR -e 12 > !{sampleID}_bulk.bb
     '''
}

process cluBB {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/bbc/", mode: 'copy'

     input:
     set val(sampleID), file(bb) from bbs

     output:
     set val(sampleID), file('*.seg'), file('*.bbc') into bbcs

     shell :
     '''
     python2 !{params.hatchet_folder}/utils/cluBB.py !{bb} -by !{params.bnpy_folder} -o !{sampleID}_bulk.seg -O !{sampleID}_bulk.bbc -e 12 -tB 0.03 -tR 0.15 -d !{params.cluBB_d}
     '''
}
bbcs.into{bbcs4plots; bbcs4hatchet }

process plots {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/plots/", mode: 'copy'

     input:
     set val(sampleID), file(seg), file(bbc) from bbcs4plots

     output:
     file('*.pdf') into plots

     shell :
     '''
     python2 !{params.hatchet_folder}/utils/BBot.py -c RD  --figsize 6,3 !{bbc} &
     python2 !{params.hatchet_folder}/utils/BBot.py -c CRD --figsize 6,3 !{bbc} &
     python2 !{params.hatchet_folder}/utils/BBot.py -c BAF --figsize 6,3 !{bbc} &
     python2 !{params.hatchet_folder}/utils/BBot.py -c BB  !{bbc} &
     python2 !{params.hatchet_folder}/utils/BBot.py -c CBB !{bbc} &
     wait
     for f in *.pdf ; do mv -- "$f" "!{sampleID}_$f" ; done
     '''
}

process hatchet {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/results/", mode: 'copy'

     input:
     set val(sampleID), file(seg), file(bbc) from bbcs4hatchet

     output:
     set val(sampleID), file('*ucn*') into ucns
     file '*.log' into logs_hatchet

     shell :
     prefix = bbc.baseName
     '''
     python2 !{params.hatchet_folder}/bin/HATCHet.py !{params.hatchet_folder}/build/solve -i !{prefix} -n2,6 -p 100 -v 2 -u 0.1 -r 12 -j !{params.cpu} -eD 6 -eT 12 -l 0.5 &> >(tee >(grep -v Progress > !{sampleID}_hatchet.log))
     for f in *ucn* ; do mv -- "$f" "!{sampleID}_$f" ; done
     '''
}

process plot_final {
     cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder+"/plots/", mode: 'copy'

     input:
     set val(sampleID), file(ucn) from ucns

     output:
     file('*.pdf')

     shell :
     '''
     python2 !{params.hatchet_folder}/utils/BBeval.py !{sampleID}_best.bbc.ucn
     for f in *.pdf ; do mv -- "$f" "!{sampleID}_$f" ; done
     '''
}
