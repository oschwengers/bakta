
import java.nio.file.Paths

params.in = 'psc.faa'
params.out = 'hmmsearch.tblout'
params.dom_out = 'hmmsearch.domtblout'
params.no_tc = false
params.dom = false
params.block = 100000

def pathInput = Paths.get(params.in).toAbsolutePath().normalize()
def pathDb = Paths.get(params.db).toAbsolutePath().normalize()
def pathOutput = Paths.get(params.out).toAbsolutePath().normalize()
def pathDomOutput = Paths.get(params.dom_out).toAbsolutePath().normalize()
def useTC = params.no_tc ? false : true
def useDom = params.dom ? true : false

print("run hmmsearch")
print("query: ${pathInput}")
print("DB: ${pathDb}")
print("Output: ${pathOutput}")
print("Output Domain: ${pathDomOutput}")
print("TC: ${useTC}")

chAAs = Channel.fromPath( pathInput )
    .splitFasta( by: params.block, file: true )

process hmmsearch {
    errorStrategy 'ignore'
    maxRetries 3
    cpus 1
    memory { 1.GB * task.attempt }
    conda 'hmmer=3.3.2'

    input:
    path('input.faa') from chAAs

    output:
    path('hmm.tblout') into chHmmResults
    path('hmm.dom.tblout') optional true  into chHmmDomResults

    String paramTC = useTC ? "--cut_tc" : "-E 1E-10"
    String paramDom = useDom ? "--domtblout hmm.dom.tblout" : ""
    script:
    """
    hmmsearch ${paramTC} -o /dev/null --noali --tblout hmm.tblout ${paramDom} --cpu ${task.cpus} ${pathDb} input.faa
    """
}

chHmmResults.collectFile( sort: false, name: pathOutput, skip: 3, keepHeader: true, storeDir: '.', tempDir: "${workDir}/tblout")
chHmmDomResults.collectFile( sort: false, name: pathDomOutput, skip: 3, keepHeader: true, storeDir: '.', tempDir: "${workDir}/domtblout")