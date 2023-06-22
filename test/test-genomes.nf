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
    memory { 8.GB * task.attempt }
    conda "${baseDir}/../environment.yml"

    input:
    tuple path('assembly.fna'), val(accession), val(genus), val(species), val(strain), val(complete) from chBakta

    output:
    path("${accession}.*") into chFinalBakta
    publishDir pattern: "${accession}.*", path: "${pathOutput}/${accession}/", mode: 'copy'

    script:
    completeOption = complete ? '--complete' : ''
    """
    ${baseDir}/../bin/bakta --db ${pathDb} --verbose --prefix ${accession} --genus ${genus} --species "${species}" --strain "${strain}" --keep-contig-headers --threads ${task.cpus} --force ${completeOption} assembly.fna
    """
}
