
params.in = 'psc.faa'
params.out = 'diamond.tsv'

params.id = 90
params.qcov = 80
params.scov = 80
params.block = 1000

Channel.fromPath( params.in )
    .splitFasta( by: params.block, file: true )
    .set( { chAAs } )

process diamond {
    errorStrategy 'finish'
    maxRetries 3
    cpus 8
    memory '32 GB'
    conda 'diamond=2.0.11'

    input:
    file('input.faa') from chAAs

    output:
    file('diamond.tsv') into chDiamondResults

    script:
    """
    diamond blastp \
        --query input.faa \
        --db ${params.db} \
        --id ${params.id} \
        --query-cover ${params.qcov} \
        --subject-cover ${params.scov} \
        --max-target-seqs 1 \
        -b4 \
        --threads ${task.cpus} \
        --out diamond.tsv \
        --outfmt 6 qseqid sseqid stitle length pident qlen slen evalue \
        --fast
    """
}

chDiamondResults.collectFile( sort: false, name: params.out, storeDir: '.')