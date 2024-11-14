import java.nio.file.*


params.db = "${baseDir}/db"
params.output = "${launchDir}"


def pathDb = Paths.get(params.db).toAbsolutePath().normalize()
def pathOutput = Paths.get(params.output).toAbsolutePath().normalize()


chSamples = Channel.fromPath( "${baseDir}/data/test-genomes.tsv" )
    .splitCsv( sep: '\t', skip: 1  )
    .view()
    .map( {
        def accession = it[0]
        def genus = it[1]
        def species = it[2]
        def strain = it[3]
        def complete = it[4] == 'yes' ? true : false
        def url = it[5]
        return [accession, genus, species, strain, complete, url]
    } )


process download {

    tag "${accession}"

    executor 'local'
    maxForks 4
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(accession), val(genus), val(species), val(strain), val(complete), val(url) from chSamples

    output:
    tuple path('assembly.fna'), val(accession), val(genus), val(species), val(strain), val(complete) into chBakta

    script:
    """
    wget -O assembly.fna.gz ${url}/${url.split('/').last()}_genomic.fna.gz
    gunzip assembly.fna.gz
    """
}

process bakta {

    tag "${accession}"

    errorStrategy 'retry'
    maxRetries 3
    cpus 8
    memory { 16.GB * task.attempt }
    conda "${baseDir}/../environment.yml"

    input:
    tuple path('assembly.fna'), val(accession), val(genus), val(species), val(strain), val(complete) from chBakta

    output:
    tuple val(accession), path("${accession}.json.gz") into chBaktaIO, chBaktaPlot
    path("${accession}.*") into chFinalBakta
    publishDir pattern: "${accession}.*", path: "${pathOutput}/${accession}/annotation/", mode: 'copy'

    script:
    completeOption = complete ? '--complete' : ''
    """
    ${baseDir}/../bin/bakta --db ${pathDb} --verbose --prefix ${accession} --genus ${genus} --species "${species}" --strain "${strain}" --keep-contig-headers --threads ${task.cpus} --force ${completeOption} assembly.fna
    gzip -k ${accession}.json
    """
}

process bakta_io {

    tag "${accession}"

    errorStrategy 'retry'
    maxRetries 3
    cpus 1
    memory 1.GB 
    conda "${baseDir}/../environment.yml"

    input:
    tuple val(accession), path('genome.json.gz') from chBaktaIO

    output:
    path("${accession}.*") into chFinalBaktaIO
    publishDir pattern: "${accession}.*", path: "${pathOutput}/${accession}/io", mode: 'copy'

    script:
    """
    ${baseDir}/../bin/bakta_io --verbose --prefix ${accession} --force genome.json.gz
    """
}

process bakta_plot {

    tag "${accession}"

    errorStrategy 'retry'
    maxRetries 3
    cpus 1
    memory 1.GB 
    conda "${baseDir}/../environment.yml"

    input:
    tuple val(accession), path('genome.json.gz') from chBaktaPlot

    output:
    path("${accession}.*") into chFinalBaktaPlot
    publishDir pattern: "${accession}.*", path: "${pathOutput}/${accession}/plot", mode: 'copy'

    script:
    """
    ${baseDir}/../bin/bakta_plot --verbose --sequences all --type cog --prefix ${accession} --force genome.json.gz
    """
}
