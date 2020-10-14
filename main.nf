$HOSTNAME = ""
params.outdir = 'results'  


if (!params.bam){params.bam = ""} 
if (!params.tab){params.tab = ""} 
if (!params.bai){params.bai = ""} 

Channel.fromPath(params.bam, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_bam_file_g_0}
Channel.fromPath(params.tab, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_2_outputFileTab_g_0}
Channel.fromPath(params.bai, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_3_bam_index_g_0}


process next {

input:
 set val(name), file(bam) from g_1_bam_file_g_0
 set val(name), file(tab) from g_2_outputFileTab_g_0
 set val(name), file(bai) from g_3_bam_index_g_0


container "onuryukselen/psi_sigma_pipeline"

script:
"""
mkdir afolder
for file in *.bam ; do mv \$file \${file//_sorted\\.bam/\\.Aligned\\.sortedByCoord\\.out\\.bam} ; done
for file in *.bai ; do mv \$file \${file//_sorted\\.bam\\.bai/\\.Aligned\\.sortedByCoord\\.out\\.bam\\.bai} ; done
for file in *SJ.out.tab ; do mv \$file \${file//_Merged_SJ\\.out\\.tab/\\.SJ\\.out\\.tab} ; done
mv *.bam *.bai *SJ.out.tab afolder
cd afolder
ln -s ${params.gtf} genes.gtf

# 0. Create groupa.txt and groupb.txt files
total=\$(ls *.SJ.out.tab|wc -l)
counta=\$(printf "%.0f" \$((\$total/2)))
countb=\$((total-counta))
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n \$counta > groupa.txt
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|tail -n \$countb > groupb.txt

"""

}

params.PSIsigma_db_path =  ""  //* @input

process create_db {



"""
perl ~/tools/PSI-Sigma-1.9j/PSIsigma-db-v.1.0.pl Homo_sapiens.GRCh38.99.sorted.gtf 1 5 1 1  &
perl ${params.PSIsigma_db_path} ${params.gtf} 1 5 1 1  && echo 0 &
echo 1
perl ${params.PSIsigma_db_path} ${params.gtf} 2 5 1 1  && echo 2 &
echo 3
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
