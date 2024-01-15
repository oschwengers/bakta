
import java.nio.file.*

params.in = 'psc.faa'
params.out = 'diamond.tsv'

params.id = 90
params.qcov = 80
params.scov = 80
params.block = 1000

def pathInput = Paths.get(params.in).toAbsolutePath().normalize()
def pathDb = Paths.get(params.db).toAbsolutePath().normalize()
def pathOutput = Paths.get(params.out).toAbsolutePath().normalize()

print("run diamond")
print("query: ${pathInput}")
print("DB: ${pathDb}")
print("Output: ${pathOutput}")

Channel.fromPath( pathInput )
    .splitFasta( by: params.block, file: true )
    .set( { chAAs } )

process diamond {
    errorStrategy 'finish'
    maxRetries 3
    cpus 8
    memory '32 GB'
    conda 'diamond=2.1.8'

    input:
    file('input.faa') from chAAs

    output:
    file('diamond.tsv') into chDiamondResults

    script:
    """
    diamond blastp \
        --query input.faa \
        --db ${pathDb} \
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

chDiamondResults.collectFile( sort: false, name: pathOutput, storeDir: '.')