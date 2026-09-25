
nextflow.enable.dsl=2

params.in = 'psc.faa'
params.out = 'hmmsearch.tblout'
params.dom_out = 'hmmsearch.domtblout'
params.no_tc = false
params.dom = false
params.block = 100000

process hmmsearch {
    errorStrategy 'retry'
    maxRetries 5
    cpus 1
    memory { 1.GB * task.attempt }
    time { 12.h * task.attempt }

    input:
    path('input.faa')

    output:
    path('hmm.tblout'), emit: tblout
    path('hmm.dom.tblout'), optional: true, emit: domtblout

    script:
    def pathDb = java.nio.file.Paths.get(params.db).toAbsolutePath().normalize()
    def useTC = params.no_tc ? false : true
    def useDom = params.dom ? true : false
    def paramTC = useTC ? "--cut_tc" : "-E 1E-10"
    def paramDom = useDom ? "--domtblout hmm.dom.tblout" : ""
    """
    hmmsearch ${paramTC} -o /dev/null --noali --tblout hmm.tblout ${paramDom} --cpu ${task.cpus} "${pathDb}" input.faa
    """
}

workflow {
    def pathInput = java.nio.file.Paths.get(params.in).toAbsolutePath().normalize()
    def pathDb = java.nio.file.Paths.get(params.db).toAbsolutePath().normalize()
    def pathOutput = java.nio.file.Paths.get(params.out).toAbsolutePath().normalize()
    def pathDomOutput = java.nio.file.Paths.get(params.dom_out).toAbsolutePath().normalize()
    def useTC = params.no_tc ? false : true

    println("run hmmsearch")
    println("query: ${pathInput}")
    println("DB: ${pathDb}")
    println("Output: ${pathOutput}")
    println("Output Domain: ${pathDomOutput}")
    println("TC: ${useTC}")

    chAAs = Channel.fromPath( pathInput )
        .splitFasta( by: params.block, file: true )

    hmmsearch(chAAs)

    hmmsearch.out.tblout.collectFile( sort: false, name: pathOutput, skip: 3, keepHeader: true, storeDir: '.')
    hmmsearch.out.domtblout.collectFile( sort: false, name: pathDomOutput, skip: 3, keepHeader: true, storeDir: '.')
}
