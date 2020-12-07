
import java.nio.file.Paths

params.in = 'psc.faa'
params.out = 'hmmsearch'

params.block = 10000

chAAs = Channel.fromPath( params.in )
    .splitFasta( by: params.block, file: true )

process hmmsearch {
    errorStrategy 'ignore'
    maxRetries 3
    cpus 2
    memory '1 GB'
    clusterOptions '-l virtual_free=1G'
    conda 'hmmer=3.3.1'

    input:
    path('input.faa') from chAAs

    output:
    path('hmm.tblout') into chHmmResults

    script:
    """
    hmmsearch -E 1E-10 --noali --tblout hmm.tblout --cpu ${task.cpus} ${params.db} input.faa
    """
}

chHmmResults.collectFile( sort: false, name: "${params.out}.tblout", storeDir: '.', skip: 3, keepHeader: true)