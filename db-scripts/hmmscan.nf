
import java.nio.file.Paths

params.in = 'psc.faa'
params.out = 'hmmscan.tblout'

params.block = 1000

Channel.fromPath( params.in )
    .splitFasta( by: params.block, file: true )
    .set( { chAAs } )

process hmmscan {
    errorStrategy 'finish'
    maxRetries 3
    cpus 2
    memory '1 GB'
    clusterOptions '-l virtual_free=1G'
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

chHmmResults.collectFile( sort: false, name: params.out, storeDir: '.')