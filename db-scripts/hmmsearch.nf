
import java.nio.file.Paths

params.in = 'psc.faa'
params.out = 'hmmsearch.tblout'

params.block = 100000

def pathInput = Paths.get(params.in).toAbsolutePath().normalize()
def pathDb = Paths.get(params.db).toAbsolutePath().normalize()
def pathOutput = Paths.get(params.out).toAbsolutePath().normalize()

print("run Diamond")
print("query: ${pathInput}")
print("DB: ${pathDb}")
print("Output: ${pathOutput}")

chAAs = Channel.fromPath( pathInput )
    .splitFasta( by: params.block, file: true )

process hmmsearch {
    errorStrategy 'ignore'
    maxRetries 3
    cpus 1
    memory '1 GB'
    conda 'hmmer=3.3.2'

    input:
    path('input.faa') from chAAs

    output:
    path('hmm.tblout') into chHmmResults

    script:
    """
    hmmsearch --cut_tc --noali --tblout hmm.tblout --cpu ${task.cpus} ${pathDb} input.faa
    """
}

chHmmResults.collectFile( sort: false, name: pathOutput, storeDir: '.', skip: 3, keepHeader: true)