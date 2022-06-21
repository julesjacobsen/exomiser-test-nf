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
ch_genome_data = Channel.fromPath("${params.genome_data_path}")
ch_phenotype_data = Channel.fromPath("${params.phenotype_data_path}")

// Define Channels from input
// runID, Proband ID, VCF_path, VCF_index_path,Proband sex, Mother ID, Father ID
Channel
    .fromPath(params.families_file)
    .ifEmpty { exit 1, "Cannot find input file : ${params.input}" }
    .splitCsv(skip:1, sep:'\t')
    .map { run_id, proband_id, hpo, vcf_path, vcf_index_path, proband_sex, mother_id, father_id -> [ run_id, proband_id, hpo, file(vcf_path), file(vcf_index_path), proband_sex, mother_id, father_id ] }
    .into { ch_input; ch_input_view }

    ch_input_view.view()

// Define Process
process run_exomiser {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'
    container 'docker.io/exomiser/exomiser-cli@sha256:2f0d869de8b06feb0abf8ac913f52937771ec947f8bdf956167925ad78b273e2'

    input:
    set val(run_id), val(proband_id), val(hpo), file(vcf_path), file(vcf_index_path), val(proband_sex), val(mother_id), val(father_id) from ch_input
    file(template_config_yaml) from ch_template_config_yaml
    file(template_application_properties) from ch_template_application_properties

    output:
    file "${output_file}*" into ch_out

    script:
    // TODO this needs adding to the input TSV file
    def hpo_ids = formatHpoIds("${hpo}")
    // Limitation - it will only cope with a trio max. Singletons will need 0 in place of the missing mother/father ids
    // Family_ID	Individual_ID	Paternal_ID	Maternal_ID	sex	Phenotype (1=unaffected, 2=affected)
    // need mother and father rows!
    def ped = createPed("${run_id}", "${proband_id}", "${father_id}", "${mother_id}", "${proband_sex}")
    def ped_path = new File("${proband_id}.ped")
    ped_path.write(ped)

    """
    cp ${template_config_yaml} exomiser_analysis.yml
    sed -i  "s/assembly_placeholder/${params.assembly}/" exomiser_analysis.yml
    sed -i  "s/vcf_placeholder/${vcf_path}/" exomiser_analysis.yml
    sed -i  "s/ped_placeholder/${ped_path}/" exomiser_analysis.yml
    sed -i  "s/proband_placeholder/${proband_id}/" exomiser_analysis.yml
    sed -i  "s/hpo_placeholder/${hpo_ids}/" exomiser_analysis.yml
    sed -i  "s/output_file_placeholder/${proband_id}/" exomiser_analysis.yml

    # volume information - should have been handled by nextflow?
    # -v "/data/exomiser-data:/exomiser-data" \
    # -v "/opt/exomiser/exomiser-config/:/exomiser"  \
    # -v "/results:/results"  \

    exomiser-cli  \
     --analysis exomiser_analysis.yml  \
     --spring.config.location=${template_application_properties}
     --exomiser.data-directory='.'
     --exomiser.${params.assembly}.data-version=${params.assembly_data_version}
     --exomiser.phenotype.data-version=${params.phenotype_data_version}
    """
  }

def createPed(runId, probandId, fatherId, motherId, probandSex) {
    def motherLine = personLine(runId, motherId, 0, 0, toPedSex('F'), 1)
    def fatherLine = personLine(runId, fatherId, 0, 0, toPedSex('M'), 1)
    def probandLine = personLine(runId, probandId, fatherId, motherId, toPedSex(probandSex), 2)
    "${motherLine}${fatherLine}${probandLine}"
}

def personLine(runId, probandId, fatherId, motherId, sex, affected) {
    "${runId}\t${probandId}\t${fatherId}\t${motherId}\t${sex}\t${affected}\n"
}

def toPedSex(sex) {
    switch (sex) {
        case 'M': return 1
        case 'F': return 2
        default: return 0
    }
}

def formatHpoIds(hpo) {
    if (!hpo || hpo == '.') {
        return ''
    }
    hpo.split(',').each(hp -> "'${hp}'").join(', ')
}