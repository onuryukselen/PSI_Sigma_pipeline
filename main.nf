$HOSTNAME = ""
params.outdir = 'results'  


if (!params.bam){params.bam = ""} 
if (!params.tab){params.tab = ""} 
if (!params.bai){params.bai = ""} 
if (!params.all_groups){params.all_groups = ""} 

Channel.fromPath(params.bam, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_bam_file_g_0}
Channel.fromPath(params.tab, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_2_outputFileTab_g_0}
Channel.fromPath(params.bai, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_3_bam_index_g_0}
g_11_all_groups_g_10 = file(params.all_groups, type: 'any') 


process rename {

input:
 set val(name), file(bam) from g_1_bam_file_g_0
 set val(name), file(tab) from g_2_outputFileTab_g_0
 set val(name), file(bai) from g_3_bam_index_g_0

output:
 file "*.bam"  into g_0_bam_file_g_7, g_0_bam_file_g_9, g_0_bam_file_g_12
 file "*.bai"  into g_0_bam_index_g_7, g_0_bam_index_g_9, g_0_bam_index_g_12
 file "*.tab"  into g_0_tab_file_g_7, g_0_tab_file_g_12

script:
"""
for file in *.bam ; do mv \$file \${file//_sorted\\.bam/\\.Aligned\\.sortedByCoord\\.out\\.bam} ; done
for file in *.bai ; do mv \$file \${file//_sorted\\.bam\\.bai/\\.Aligned\\.sortedByCoord\\.out\\.bam\\.bai} ; done
for file in *SJ.out.tab ; do mv \$file \${file//_Merged_SJ\\.out\\.tab/\\.SJ\\.out\\.tab} ; done

"""

}

chromosome_list = params.PSI_Sigma_prep.chromosome_list
chromosome_list = chromosome_list.split(',')

process PSI_Sigma_prep {


output:
 val chromosome_list  into g_5_name_g_7

script:
"""
echo ${chromosome_list}
"""
}

params.PSIsigma_db_path =  ""  //* @input

process create_db {

input:
 file bam from g_0_bam_file_g_7.collect()
 file bai from g_0_bam_index_g_7.collect()
 file tab from g_0_tab_file_g_7.collect()
 val chr_name from g_5_name_g_7.unique().flatten()

output:
 file "*.txt"  into g_7_groups
 file "*.tmp"  into g_7_tmp_file_g_8

container "onuryukselen/psi_sigma_pipeline"

script:
"""
# 0. Create groupa.txt and groupb.txt files
total=\$(ls *.SJ.out.tab|wc -l)
counta=\$(printf "%.0f" \$((\$total/2)))
countb=\$((total-counta))
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n \$counta > groupa.txt
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|tail -n \$countb > groupb.txt
perl ${params.PSIsigma_db_path} ${params.gtf} ${chr_name} 5 1 1
"""
}


process merge_db {

input:
 file tmp from g_7_tmp_file_g_8.collect()

output:
 file "PSIsigma1d9j.bed"  into g_8_bed_file
 file "PSIsigma1d9j.db"  into g_8_db_file_g_9

"""
cat *.db.tmp > PSIsigma1d9j.db
cat *.bed.tmp > PSIsigma1d9j.bed

"""
}

params.PSIsigma_ir_path =  ""  //* @input

process create_intronic_read {

input:
 file bam from g_0_bam_file_g_9
 file bai from g_0_bam_index_g_9
 file db from g_8_db_file_g_9

output:
 file "*.IR.out.tab"  into g_9_tab_file_g_12

container "onuryukselen/psi_sigma_pipeline"

script:
"""
perl ${params.PSIsigma_ir_path} $db $bam 1 
"""
}


process prepare_groups {

input:
 file groups from g_11_all_groups_g_10

output:
 file "*"  into g_10_all_groups_g_12

// # DMSO_1	A0	DMSO
// # DMSO_2	A0	DMSO
// # EXP1_10nM_1	A1	EXP1_10nM
// # EXP1_10nM_2	A1	EXP1_10nM
// # EXP2_10nM_1	A2	EXP2_10nM
// # EXP2_10nM_2	A2	EXP2_10nM
// # DMSO_CHX_1	B0	DMSO_CHX
// # DMSO_CHX_2	B0	DMSO_CHX
// # EXP3_10nM_1	B1	EXP3_10nM
// # EXP3_10nM_2	B1	EXP3_10nM
// # EXP4_10nM_1	B2	EXP4_10nM
// # EXP4_10nM_2	B2	EXP4_10nM

// # create directory of third column
// # insert groupa and groupb into directory.
shell:
'''
#!/usr/bin/env perl
use Data::Dumper;
use strict;

### Parse group file
my $group_file = "!{groups}";
open(IN, "<$group_file") or die "Can't open $group_file.\\n";
#all_groups: 'A' => { '1' => {'cond' => 'EXP1_10nM' file' => ['EXP1_10nM_1','EXP1_10nM_2']}, '0' => {'cond' => 'DMSO','file' => ['DMSO_1','DMSO_2']}
my %all_groups;
while(<IN>){
  chomp;
  my ($file, $group, $cond) = split("\\t");
  my ($groupName, $groupId) = split /(\\d+)/, $group;

  if(!exists($all_groups{$groupName}{$groupId})){
    $all_groups{$groupName}{$groupId} = {"cond" => $cond, "file" => $file};
  } else{
    my $files =  $all_groups{$groupName}{$groupId}{"file"};
    $all_groups{$groupName}{$groupId}{"file"} = $files."\\n".$file;
  }
}


foreach my $groupName (keys %all_groups){
  my $contFiles = $all_groups{$groupName}{"0"}{"file"};
  foreach my $groupId (keys %{ $all_groups{$groupName}}){
	if ($groupId ne "0"){
		my $dir = $all_groups{$groupName}{$groupId}{"cond"};
		my $expfiles = $all_groups{$groupName}{$groupId}{"file"};

		mkdir $dir unless -d $dir; # Check if dir exists. If not create it.
		open my $out, ">", "$dir/groupa.txt" or die "Can't open '$dir/groupa.txt'\n";
		print $out $contFiles;
		close $out;
		
		open my $out, ">", "$dir/groupb.txt" or die "Can't open '$dir/groupb.txt'\n";
		print $out $expfiles;
		close $out;
    	
	} 
  }
}


print Dumper \\%all_groups;
'''
}

params.dummyai_path =  ""  //* @input

process run_PSI_Sigma_ {

input:
 file cond from g_10_all_groups_g_12.flatten()
 file ir_out_tab from g_9_tab_file_g_12.collect()
 file out_tab from g_0_tab_file_g_12.collect()
 file bam from g_0_bam_file_g_12.collect()
 file bai from g_0_bam_index_g_12.collect()


container "onuryukselen/psi_sigma_pipeline"

script:
"""
mv $ir_out_tab $out_tab $bam $bai $cond/.
cd $cond
perl ${params.dummyai_path} --gtf ${params.gtf} --nread 10 --type 1 --name ${cond}_PSIsigma1d9j --irmode 1 --skipratio 0.05 --denominator 0 --irrange 10 --adjp 0
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
