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
params.input_folder = null
params.ref = null
params.output_folder = "."

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
    log.info "--input_folder         PATH        Folder containing bam files"
    log.info "--ref                  PATH        WHOLE Path to reference fasta file (should be indexed)"
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

assert (params.input_folder != null) : "please provide the --input_folder option"
assert (params.ref != null) : "please provide the --ref option"

bams = Channel.fromPath( params.input_folder+'/*.bam')
              .ifEmpty { error "Cannot find any BAM file in: ${params.input_folder}" }

process hatchet {
		 cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder, mode: 'copy'

     input :
     file bam from bams

     output:
     file "${sampleID}*.vcf" into vcf
     file "${sampleID}.alignments.txt.gz" into alignments

     shell :
     sampleID=bam.baseName.replace("bam","")
     '''
     !{baseDir}/bin/run_HATCHet.sh ${params.cpu} !{params.ref} !{params.samtools_folder} !{params.bcftools_folder} !{params.bnpy_folder} !{params.hatchet_folder} !{params.input_folder}
     '''
}
