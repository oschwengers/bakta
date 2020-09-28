
params.psc = 'psc.faa'

Channel.fromPath( params.psc )
    .splitFasta( by: 100000, file: true )
    .set( { chAAs } )

process hmmsearch {
    errorStrategy 'finish'
    maxRetries 3
    cpus 4
    memory '1 GB'
    conda 'hmmer=3.3.1'

    input:
    file('input.faa') from chAAs

    output:
    file('hmm.out') into chHmmResults

    script:
    """
    hmmsearch --cut_tc --noali --tblout hmm.out --cpu ${task.cpus} ${params.db} input.faa
    """
}

chHmmResults.collectFile( sort: false, name: 'hmm-amr.out', storeDir: '.')