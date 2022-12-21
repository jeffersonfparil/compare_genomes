///////////////////////////////////
// SETUP THE WORKING ENVIRONMENT //
///////////////////////////////////
// Assumed extension names:
//  - genome: *.fna
//  - annotation: *.gff
//  - coding DNA: *.cds
//  - proteome: *.faa

process TEST {
    label "LOW_MEM_LOW_CPU"
    input:
        val dir
    output:
        val dir
    shell:
    '''
    cd !{dir}
    echo !{params.urls} > test-1.tmp
    bash !{projectDir}/../scripts/test.sh !{dir}
    '''
}


workflow {
    TEST(params.dir)
}
