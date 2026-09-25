
nextflow.enable.dsl=2

params.in = 'psc.faa'
params.out = 'diamond.tsv'

params.id = 90
params.qcov = 80
params.scov = 80
params.block = 1000

process diamond {
    errorStrategy 'retry'
    maxRetries 5
    cpus 8
    memory '32 GB'
    time { 12.h * task.attempt }

    input:
    path('input.faa')

    output:
    path('diamond.tsv')

    script:
    def pathDb = java.nio.file.Paths.get(params.db).toAbsolutePath().normalize()
    """
    diamond blastp \
        --query input.faa \
        --db "${pathDb}" \
        --id ${params.id} \
        --query-cover ${params.qcov} \
        --subject-cover ${params.scov} \
        --max-target-seqs 1 \
        -b4 \
        --threads ${task.cpus} \
        --load-threads ${task.cpus} \
        --out diamond.tsv \
        --outfmt 6 qseqid sseqid stitle length pident qlen slen evalue \
        --fast
    """
}

workflow {
    def pathInput = java.nio.file.Paths.get(params.in).toAbsolutePath().normalize()
    def pathDb = java.nio.file.Paths.get(params.db).toAbsolutePath().normalize()
    def pathOutput = java.nio.file.Paths.get(params.out).toAbsolutePath().normalize()

    println("run diamond")
    println("query: ${pathInput}")
    println("DB: ${pathDb}")
    println("Output: ${pathOutput}")

    chAAs = Channel.fromPath( pathInput )
        .splitFasta( by: params.block, file: true )

    diamond(chAAs)

    diamond.out.collectFile( sort: false, name: pathOutput, storeDir: '.')
}
