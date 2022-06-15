#!/usr/bin/env nextflow

def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --bams sample.bam [Options]
    
    Inputs Options:
    --input         Input file

    Resource Options:
    --max_cpus      Maximum number of CPUs (int)
                    (default: $params.max_cpus)  
    --max_memory    Maximum memory (memory unit)
                    (default: $params.max_memory)
    --max_time      Maximum time (time unit)
                    (default: $params.max_time)
    See here for more info: https://github.com/lifebit-ai/hla/blob/master/docs/usage.md
    """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

// Define channels from repository files
projectDir = workflow.projectDir
ch_run_sh_script = Channel.fromPath("${projectDir}/bin/run.sh")
ch_template_config_yaml = Channel.fromPath("${projectDir}/bin/template_config.yaml")
ch_template_application_properties = Channel.fromPath("${projectDir}/bin/application.properties")

// Define Channels from input
// runID, Proband ID, VCF_path, VCF_index_path,Proband sex, Mother ID, Father ID
Channel
    .fromPath(params.families_file)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1, sep:'\t')
    .map { run_id, proband_id, vcf_path, vcf_index_path, proband_sex, mother_id, father_id -> [ run_id, proband_id, file(vcf_path), file(vcf_index_path), proband_sex, mother_id, father_id ] }
    .set { ch_input }

	ch_input.view()
//
// // Define Process
// process run_exomiser {
//     tag "$sample_name"
//     label 'low_memory'
//     publishDir "${params.outdir}", mode: 'copy'
//     container "exomiser/exomiser-cli:sha@2f0d869de8b0"
//
//     input:
//     set val(sample_name), file(input_file) from ch_input
//     file(template_config_yaml) from ch_template_config_yaml
//
//     output:
//     file "${output_file}*" into ch_out
//
//     script:
//     """
//     #-v "/data/exomiser-data:/exomiser-data" \
//     # -v "/opt/exomiser/exomiser-config/:/exomiser"  \
//     # -v "/opt/exomiser/exomiser-cli-${project.version}/results:/results"  \
//     cp ${template_config_yaml} exomiser_analysis.yml
//     sed -i  "s/assembly_placeholder/${assembly}/" exomiser_analysis.yml
//     sed -i  "s/vcf_placeholder/${vcf}/" exomiser_analysis.yml
//     sed -i  "s/ped_placeholder/${ped}/" exomiser_analysis.yml
//     sed -i  "s/proband_placeholder/${proband}/" exomiser_analysis.yml
//     sed -i  "s/hpo_placeholder/${hpo}/" exomiser_analysis.yml
//     sed -i  "s/output_file_placeholder/${output_file}/" exomiser_analysis.yml
//
//      exomiser-cli  \
//      --analysis exomiser_analysis.yml  \
//      --spring.config.location=/exomiser/application.properties
//     """
//   }
