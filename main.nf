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

Channel
    .fromPath("${params.exomiser_data_path}", type: 'dir')
    .map { exomiser_data_path -> [file(exomiser_data_path)] }
    .set { ch_exomiser_data_path }


// Define Process
process run_exomiser {
    tag "$sample_name"
    label 'low_memory'
    publishDir "${params.outdir}", mode: 'copy'
    // container 'docker.io/exomiser/exomiser-cli@sha256:2f0d869de8b06feb0abf8ac913f52937771ec947f8bdf956167925ad78b273e2'
    container 'quay.io/lifebitai/exomiser:12.1.0'
//     containerOptions "-v ${params.exomiser_data_path}:/exomiser-data"
    containerOptions "-v ${params.exomiser_data_path}:/exomiser-data"
    // TODO: this bit is broken - the ${params.exomiser_data_path} works locally using an absolute path
    // https://github.com/julesjacobsen/exomiser-test-nf/issues/3
//     containerOptions "-v exomiser-data:/exomiser-data"

    input:
    set val(run_id), val(proband_id), val(hpo), file(vcf_path), file(vcf_index_path), val(proband_sex), val(mother_id), val(father_id) from ch_input
    file(template_config_yaml) from ch_template_config_yaml
    file(template_application_properties) from ch_template_application_properties
    file(exomiser_data_path) from ch_exomiser_data_path

    output:
    file "${output_file}*" into ch_out

    script:
    // TODO this needs adding to the input TSV file
    def hpo_ids = formatHpoIds("${hpo}")
    // Limitation - it will only cope with a trio max. Singletons will need 0 in place of the missing mother/father ids
    def ped = createPed("${run_id}", "${proband_id}", "${father_id}", "${mother_id}", "${proband_sex}")
    def ped_path = "${proband_id}.ped"
    """
    # hacky stuff to get pedigree into working directory instead of root
    printf '${ped}' >> ${ped_path}
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

    java -jar /exomiser/exomiser-cli-12.1.0.jar  \
     --analysis exomiser_analysis.yml  \
     --exomiser.data-directory=/exomiser-data \
     --exomiser.hg19.data-version=2202 \
     --exomiser.phenotype.data-version=2202
    """
  }

def createPed(familyId, probandId, fatherId, motherId, probandSex) {
    // Family_ID	Individual_ID	Paternal_ID	Maternal_ID	sex	Phenotype (1=unaffected, 2=affected)
    def motherLine = personLine(familyId, motherId, 0, 0, 'F', 1)
    def fatherLine = personLine(familyId, fatherId, 0, 0, 'M', 1)
    def probandLine = personLine(familyId, probandId, fatherId, motherId, probandSex, 2)
    if (!motherLine && !fatherLine) {
        return personLine(familyId, probandId, 0, 0, probandSex, 2)
    }
    "${motherLine}${fatherLine}${probandLine}"
}

def personLine(familyId, individualId, fatherId, motherId, individualSex, affected) {
    if (!individualId || individualId == '.') {
        return ''
    }
    "${familyId}\t${individualId}\t${fatherId}\t${motherId}\t${-> toPedSex(individualSex)}\t${affected}\n"
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