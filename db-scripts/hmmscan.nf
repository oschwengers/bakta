
params.fasta = 'psc.faa'

Channel.fromPath( params.fasta )
    .splitFasta( by: 40000, file: true )
    .set( { chAAs } )

process hmmscan {
    errorStrategy 'finish'
    maxRetries 3
    cpus 1
    memory '1 GB'
    conda 'hmmer=3.3.1'

    input:
    file('input.faa') from chAAs

    output:
    file('hmm.out') into chHmmResults

    script:
    """
    hmmscan --cut_tc --noali --tblout hmm.out --cpu ${task.cpus} ${params.db} input.faa
    """
}

chHmmResults.collectFile( sort: false, name: 'hmmscan.tblout', storeDir: '.')