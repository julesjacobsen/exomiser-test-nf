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
process prepare_exomiser_input_files {
    tag "$sample_name"
    label 'low_memory'

    input:
    set val(run_id), val(proband_id), val(hpo), path(vcf_path), path(vcf_index_path), val(proband_sex), val(mother_id), val(father_id) from ch_input
    path(template_config_yaml) from ch_template_config_yaml

    output:
    file "${proband_id}*.{ped,yml,gz}" into ch_input_files
    val(proband_id) into ch_proband_id

    script:
    // TODO this needs adding to the input TSV file
//     def hpo_ids = formatHpoIds("${hpo}")
    // Limitation - it will only cope with a trio max. Singletons will need 0 in place of the missing mother/father ids
    def ped = createPed("${run_id}", "${proband_id}", "${father_id}", "${mother_id}", "${proband_sex}")
    def ped_path = "${proband_id}.ped"
    """
    # hacky stuff to get pedigree into working directory instead of root
    printf '${ped}' >> ${ped_path}
    cp ${template_config_yaml} ${proband_id}-analysis.yml
    mv ${vcf_path} ${proband_id}.vcf.gz
    sed -i  "s/assembly_placeholder/${params.assembly}/" ${proband_id}-analysis.yml
    sed -i  "s/vcf_placeholder/${proband_id}.vcf.gz/" ${proband_id}-analysis.yml
    sed -i  "s/ped_placeholder/${ped_path}/" ${proband_id}-analysis.yml
    sed -i  "s/proband_placeholder/${proband_id}/" ${proband_id}-analysis.yml
    sed -i  "s/hpo_placeholder/${hpo}/" ${proband_id}-analysis.yml
    sed -i  "s/output_file_placeholder/${proband_id}/" ${proband_id}-analysis.yml
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

// def formatHpoIds(hpo) {
//     if (!hpo || hpo == '.') {
//         return ''
//     }
//     hpo.split(',').each(hp -> '${hp}').join(', ')
// }

process run_exomiser {
    disk = '2.GB'
    cpus = 4
    memory = '12.GB'
    // TODO: how to enable NextFlow to run this using the container specified entry-point rather than /bin/bash?
//     container = 'docker.io/exomiser/exomiser-cli@sha256:2f0d869de8b06feb0abf8ac913f52937771ec947f8bdf956167925ad78b273e2'
    container 'quay.io/lifebitai/exomiser:12.1.0'
    containerOptions = "-v ${params.exomiser_data_path}:/data/exomiser-data"
//     publishDir "${params.report_dir}", mode: 'copy'

    input:
        set val(proband_id) from ch_proband_id
        path(x) from ch_input_files
        path(exomiser_data_path) from ch_exomiser_data_path
        path(application_properties) from ch_template_application_properties

//     output:
//         file "${proband_id}*.{html,json,tsv}" into ch_out

    script:
    """
    # Create symlink between absolut path of the staged folder and the folder defined in application.properties "/data/exomiser-data"
    # see here: https://github.com/julesjacobsen/exomiser-test-nf/blob/0df9324df65a358e226e3898f92d1783e9702990/bin/application.properties#L26
    ln -s "\$PWD/$exomiser_data_path/" /data/exomiser-data
    java -jar /exomiser/exomiser-cli-12.1.0.jar  \
     --analysis "${proband_id}"-analysis.yml  \
     --spring.config.location=$application_properties \
     --exomiser.data-directory='.' \
     --exomiser.${params.assembly}.data-version=${params.assembly_data_version} \
     --exomiser.phenotype.data-version=${params.phenotype_data_version}
    """
}