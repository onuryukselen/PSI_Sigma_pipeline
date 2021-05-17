$HOSTNAME = ""
params.outdir = 'results'  

$HOSTNAME = "default"
params.DOWNDIR = (params.DOWNDIR) ? params.DOWNDIR : ""
//* params.genome_build =  ""  //* @dropdown @options:"human_hg19, human_hg38, macaque_macFas5, rat_rn6, rat_rn6ens, mousetest_mm10, custom"
//* params.run_PSI_Sigma =  "yes"  //* @dropdown @options:"yes","no" @show_settings:"PSI_Sigma_prep", "prepare_groups"


def _species;
def _build;
def _share;
//* autofill
if (params.genome_build == "mousetest_mm10"){
    _species = "mousetest"
    _build = "mm10"
} else if (params.genome_build == "human_hg19"){
    _species = "human"
    _build = "hg19"
} else if (params.genome_build == "human_hg38"){
    _species = "human"
    _build = "hg38"
} else if (params.genome_build == "mouse_mm10"){
    _species = "mouse"
    _build = "mm10"
} else if (params.genome_build == "macaque_macFas5"){
    _species = "macaque"
    _build = "macFas5"
} else if (params.genome_build == "rat_rn6"){
    _species = "rat"
    _build = "rn6"
} else if (params.genome_build == "rat_rn6ens"){
    _species = "rat"
    _build = "rn6ens"
}
if ($HOSTNAME == "default"){
    _share = "${params.DOWNDIR}/genome_data"
    $SINGULARITY_IMAGE = "https://galaxyweb.umassmed.edu/pub/dnext_data/singularity/UMMS-Biocore-rna-seq-2.0.simg"
}

if ($HOSTNAME == "fs-bb7510f0"){
    _share = "/mnt/efs/share/genome_data"
    $SINGULARITY_IMAGE = "https://galaxyweb.umassmed.edu/pub/dnext_data/singularity/UMMS-Biocore-rna-seq-2.0.simg"
	$SINGULARITY_OPTIONS = "--bind /mnt"
} else if ($HOSTNAME == "192.168.20.150"){
    _share = "/home/botaoliu/share/genome_data"
    $SINGULARITY_IMAGE = "https://galaxyweb.umassmed.edu/pub/dnext_data/singularity/UMMS-Biocore-rna-seq-2.0.simg"
} else if ($HOSTNAME == "50.228.141.2"){
    _share = "/share/genome_data"
    $SINGULARITY_IMAGE = "https://galaxyweb.umassmed.edu/pub/dnext_data/singularity/UMMS-Biocore-rna-seq-2.0.simg"
    $CPU  = 1
    $MEMORY = 10
} else if ($HOSTNAME == "ghpcc06.umassrc.org"){
    _share = "/share/data/umw_biocore/genome_data"
    $SINGULARITY_IMAGE = "/project/umw_biocore/singularity/UMMS-Biocore-rna-seq-2.0.simg"
    $TIME = 500
    $CPU  = 1
    $MEMORY = 32 
    $QUEUE = "long"
}
if (params.genome_build && $HOSTNAME){
    params.genomeDir ="${_share}/${_species}/${_build}/"
    params.genome ="${_share}/${_species}/${_build}/${_build}.fa"
    params.gtf ="${_share}/${_species}/${_build}/ucsc.gtf"
    params.gff ="${_share}/${_species}/${_build}/genes.gff3"
    params.exon_file  = "${_share}/${_species}/${_build}/exons.txt.gz"
    params.bed_file_genome ="${_share}/${_species}/${_build}/${_build}.bed"
    params.ref_flat ="${_share}/${_species}/${_build}/ref_flat"
    params.genomeSizePath ="${_share}/${_species}/${_build}/${_build}.chrom.sizes"

}
if ($HOSTNAME){
    params.genomeCoverageBed_path = "genomeCoverageBed"
    params.wigToBigWig_path = "wigToBigWig"
    params.samtools_path = "samtools"
    params.pdfbox_path = "/usr/local/bin/dolphin-tools/pdfbox-app-2.0.0-RC2.jar"
    params.gtf2bed_path = "/usr/local/bin/dolphin-tools/gtf2bed"
    params.PSIsigma_db_path = "/usr/local/bin/PSI-Sigma-1.9l/PSIsigma-db-v.1.0.pl"
    params.PSIsigma_ir_path = "/usr/local/bin/PSI-Sigma-1.9l/PSIsigma-ir-v.1.2.pl"
    params.dummyai_path = "/usr/local/bin/PSI-Sigma-1.9l/dummyai.pl"
    $CPU  = 1
    $MEMORY = 10
}
//*


if (!params.bam){params.bam = ""} 
if (!params.tab){params.tab = ""} 
if (!params.custom_gtf){params.custom_gtf = ""} 
if (!params.gtf){params.gtf = ""} 

Channel.fromPath(params.bam, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_bam_file_g_18}
Channel.fromPath(params.tab, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_2_outputFileTab_g_18}
Channel.value(params.custom_gtf).into{g_20_gtfFilePath_g_26;g_20_gtfFilePath_g_46}
Channel.value(params.gtf).into{g_28_gtfFilePath_g_27;g_28_gtfFilePath_g_47}


process rename {

input:
 set val(name_bam), file(bam) from g_1_bam_file_g_18
 set val(name_tab), file(tab) from g_2_outputFileTab_g_18

output:
 file "*.bam"  into g_18_bam_file_g_26, g_18_bam_file_g_27, g_18_bam_file_g_44, g_18_bam_file_g_45, g_18_bam_file_g_46, g_18_bam_file_g_47
 file "*.bai"  into g_18_bam_index_g_26, g_18_bam_index_g_27, g_18_bam_index_g_44, g_18_bam_index_g_45, g_18_bam_index_g_46, g_18_bam_index_g_47
 file "*.tab"  into g_18_tab_file_g_26, g_18_tab_file_g_27, g_18_tab_file_g_46, g_18_tab_file_g_47

when:
(params.run_PSI_Sigma && (params.run_PSI_Sigma == "yes")) || !params.run_PSI_Sigma

script:
"""
if [ "$bam" != "${name_bam}.bam" ]; then
    mv $bam ${name_bam}.bam
fi
if [ "$tab" != "${name_tab}.tab" ]; then
    mv $tab ${name_tab}.tab
fi
samtools index ${name_bam}.bam
mv ${name_bam}.bam ${name_bam}.Aligned.sortedByCoord.out.bam
mv ${name_bam}.bam.bai ${name_bam}.Aligned.sortedByCoord.out.bam.bai
mv ${name_tab}.tab ${name_tab}.SJ.out.tab

"""

}


process PSI_Sigma_prep {


output:
 val chromosome_list  into g_40_name_g_26, g_40_name_g_27
 val PSI_sigma_parameters  into g_40_parameters_g_46, g_40_parameters_g_47

when:
(params.run_PSI_Sigma && (params.run_PSI_Sigma == "yes")) || !params.run_PSI_Sigma

script:
PSI_sigma_parameters = params.PSI_Sigma_prep.PSI_sigma_parameters
chromosome_list = params.PSI_Sigma_prep.chromosome_list
chromosome_list = chromosome_list.split(',')
chromosome_list*.trim()
"""
echo ${chromosome_list}
"""

}

//* params.PSIsigma_db_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process create_db {

input:
 file bam from g_18_bam_file_g_27.collect()
 file bai from g_18_bam_index_g_27.collect()
 file tab from g_18_tab_file_g_27.collect()
 val chr_name from g_40_name_g_27.flatten()
 val custom_gtf from g_28_gtfFilePath_g_27

output:
 file "*.txt"  into g_27_groups
 file "*.tmp"  into g_27_tmp_file_g_29
 val db_name  into g_27_db_name_g_29

container "onuryukselen/psi_sigma_pipeline:3.0"

when:
custom_gtf.toString().indexOf("/") > -1

script:
println custom_gtf.toString().indexOf("/") > -1
println custom_gtf
custom_gtf = custom_gtf.toString()
db_name = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "PSIsigma1d9l_StringTie": "PSIsigma1d9l" 
println db_name
"""
if [ -e "${custom_gtf}" ]; then
    gtfPath="${custom_gtf} "
elif [ -e "${params.gtf}" ]; then
    gtfPath="${params.gtf} "
fi
echo "gtf path: \${gtfPath}"
# 0. Create groupa.txt and groupb.txt files
total=\$(ls *.SJ.out.tab|wc -l)
counta=\$(printf "%.0f" \$((\$total/2)))
countb=\$((total-counta))
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n \$counta > groupa.txt
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|tail -n \$countb > groupb.txt
perl ${params.PSIsigma_db_path} \$gtfPath $chr_name 5 1 1
"""
}


process merge_db {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${db_name}.bed$/) "PSIsigma/$filename"
	else if (filename =~ /${db_name}.db$/) "PSIsigma/$filename"
}

input:
 file tmp from g_27_tmp_file_g_29.collect()
 val db_name from g_27_db_name_g_29.collect()

output:
 file "${db_name}.bed"  into g_29_bed_file
 file "${db_name}.db"  into g_29_db_file_g_45, g_29_db_file_g_47

script:
db_name = db_name[0]
"""	
cat *.db.tmp > ${db_name}.db
cat *.bed.tmp > ${db_name}.bed
"""

}

//* params.PSIsigma_ir_path =  ""  //* @input

process create_intronic_read {

input:
 file bam from g_18_bam_file_g_45
 file bai from g_18_bam_index_g_45
 file db from g_29_db_file_g_45

output:
 file "*.IR.out.tab"  into g_45_tab_file_g_47

container "onuryukselen/psi_sigma_pipeline:3.0"

script:
"""
perl ${params.PSIsigma_ir_path} $db $bam 1 
"""
}

//* params.PSIsigma_db_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process create_db_stringtie {

input:
 file bam from g_18_bam_file_g_26.collect()
 file bai from g_18_bam_index_g_26.collect()
 file tab from g_18_tab_file_g_26.collect()
 val chr_name from g_40_name_g_26.flatten()
 val custom_gtf from g_20_gtfFilePath_g_26

output:
 file "*.txt"  into g_26_groups
 file "*.tmp"  into g_26_tmp_file_g_8
 val db_name  into g_26_db_name_g_8

container "onuryukselen/psi_sigma_pipeline:3.0"

when:
custom_gtf.toString().indexOf("/") > -1

script:
println custom_gtf.toString().indexOf("/") > -1
println custom_gtf
custom_gtf = custom_gtf.toString()
db_name = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "PSIsigma1d9l_StringTie": "PSIsigma1d9l" 
println db_name
"""
if [ -e "${custom_gtf}" ]; then
    gtfPath="${custom_gtf} "
elif [ -e "${params.gtf}" ]; then
    gtfPath="${params.gtf} "
fi
echo "gtf path: \${gtfPath}"
# 0. Create groupa.txt and groupb.txt files
total=\$(ls *.SJ.out.tab|wc -l)
counta=\$(printf "%.0f" \$((\$total/2)))
countb=\$((total-counta))
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|head -n \$counta > groupa.txt
ls *.SJ.out.tab|sed 's/.SJ.out.tab//'|tail -n \$countb > groupb.txt
perl ${params.PSIsigma_db_path} \$gtfPath $chr_name 5 1 1
"""
}


process merge_db_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${db_name}.bed$/) "PSIsigma_StringTie/$filename"
	else if (filename =~ /${db_name}.db$/) "PSIsigma_StringTie/$filename"
}

input:
 file tmp from g_26_tmp_file_g_8.collect()
 val db_name from g_26_db_name_g_8.collect()

output:
 file "${db_name}.bed"  into g_8_bed_file
 file "${db_name}.db"  into g_8_db_file_g_44, g_8_db_file_g_46

script:
db_name = db_name[0]
"""	
cat *.db.tmp > ${db_name}.db
cat *.bed.tmp > ${db_name}.bed
"""

}

def downFile(path){
    if (path.take(1).indexOf("/") == 0){
      target=path
    } else {
      a=file(path)
      fname = a.getName().toString()
      target = "${baseDir}/work/${fname}"
      a.copyTo(target) 
    }
    return target
}
process prepare_groups {


output:
 file "*"  into g_10_all_groups_g_46, g_10_all_groups_g_47

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
when:
(params.run_PSI_Sigma && (params.run_PSI_Sigma == "yes")) || !params.run_PSI_Sigma

shell:
all_groups_path = params.prepare_groups.all_groups_path

'''
#!/usr/bin/env perl
use Data::Dumper;
use strict;

### Parse group file
my $group_file = "!{all_groups_path}";
open(IN, "<$group_file") or die "Can't open $group_file.\\n";
#all_groups: 'A' => { '1' => {'cond' => 'EXP1_10nM' file' => ['EXP1_10nM_1','EXP1_10nM_2']}, '0' => {'cond' => 'DMSO','file' => ['DMSO_1','DMSO_2']}
my %all_groups;
while(<IN>){
  chomp;
  my ($file, $group, $cond) = split("\\t");
  $cond =~ s/[\\r\\n]+$//;
  my ($groupName, $groupId) = split /(\\d+)/, $group;
  if ($groupName ne "" && $groupId ne ""){
	if(!exists($all_groups{$groupName}{$groupId})){
    	$all_groups{$groupName}{$groupId} = {"cond" => $cond, "file" => $file};
	} else{
    	my $files =  $all_groups{$groupName}{$groupId}{"file"};
    	$all_groups{$groupName}{$groupId}{"file"} = $files."\\n".$file;
	}
  }
}


foreach my $groupName (keys %all_groups){
  my $contFiles = $all_groups{$groupName}{"0"}{"file"};
  foreach my $groupId (keys %{ $all_groups{$groupName}}){
	if ($groupId ne "0"){
		my $dir = $all_groups{$groupName}{$groupId}{"cond"};
		$dir = "_conds_".$dir;
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

g_28_gtfFilePath_g_47= g_28_gtfFilePath_g_47.ifEmpty([""]) 

//* params.dummyai_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process PSI_Sigma_run {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "PSIsigma/$filename"
}

input:
 file cond from g_10_all_groups_g_47.flatten()
 file ir_out_tab from g_45_tab_file_g_47.collect()
 file out_tab from g_18_tab_file_g_47.collect()
 file bam from g_18_bam_file_g_47.collect()
 file bai from g_18_bam_index_g_47.collect()
 val custom_gtf from g_28_gtfFilePath_g_47
 file db from g_29_db_file_g_47
 val PSI_sigma_parameters from g_40_parameters_g_47

output:
 file "*"  into g_47_outputDir

container "onuryukselen/psi_sigma_pipeline:3.0"

script:
mainGtf = params.gtf.toString()
custom_gtf = custom_gtf.toString()
nameCustomGtf = custom_gtf.substring(custom_gtf.lastIndexOf('/')+1,custom_gtf.length())
nameMainGtf = mainGtf.substring(mainGtf.lastIndexOf('/')+1,mainGtf.length())
dbName = db.baseName
"""
if [ -e "${custom_gtf}" ]; then
	cp ${custom_gtf} ${nameCustomGtf}
    gtfFile="${nameCustomGtf} "
elif [ -e "${params.gtf}" ]; then
	cp ${params.gtf} ${nameMainGtf}
    gtfFile="${nameMainGtf} "
fi
echo "gtf name: \${gtfFile}"

## merge groupa.txt and groupb.txt and mv related files into cond folder
# find ${cond}/group*.txt | xargs -I{} sh -c "cat {}; echo ''" > allgroups.txt
# while IFS="" read -r line || [ -n "\$line" ]
# do
#   mv \${line}* ${cond}/.
# done < allgroups.txt
# mv  $db ${cond}/.

mv $ir_out_tab $out_tab $bam $bai $db ${cond}/.
cp \$gtfFile ${cond}/.

## clean "_conds_" prefix from condition name
cond="${cond}"
realcond=\${cond#"_conds_"}
mv ${cond} \$realcond

cd \$realcond
tmpdir="${baseDir}/work"
export TMPDIR=\$tmpdir
perl ${params.dummyai_path} --gtf \$gtfFile  --name ${dbName} ${PSI_sigma_parameters} > log.txt 2>&1
rm $bam $bai \$gtfFile $db
"""
}

//* params.PSIsigma_ir_path =  ""  //* @input

process create_intronic_read {

input:
 file bam from g_18_bam_file_g_44
 file bai from g_18_bam_index_g_44
 file db from g_8_db_file_g_44

output:
 file "*.IR.out.tab"  into g_44_tab_file_g_46

container "onuryukselen/psi_sigma_pipeline:3.0"

script:
"""
perl ${params.PSIsigma_ir_path} $db $bam 1 
"""
}

g_20_gtfFilePath_g_46= g_20_gtfFilePath_g_46.ifEmpty([""]) 

//* params.dummyai_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process PSI_Sigma_run_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "PSIsigma_StringTie/$filename"
}

input:
 file cond from g_10_all_groups_g_46.flatten()
 file ir_out_tab from g_44_tab_file_g_46.collect()
 file out_tab from g_18_tab_file_g_46.collect()
 file bam from g_18_bam_file_g_46.collect()
 file bai from g_18_bam_index_g_46.collect()
 val custom_gtf from g_20_gtfFilePath_g_46
 file db from g_8_db_file_g_46
 val PSI_sigma_parameters from g_40_parameters_g_46

output:
 file "*"  into g_46_outputDir

container "onuryukselen/psi_sigma_pipeline:3.0"

script:
mainGtf = params.gtf.toString()
custom_gtf = custom_gtf.toString()
nameCustomGtf = custom_gtf.substring(custom_gtf.lastIndexOf('/')+1,custom_gtf.length())
nameMainGtf = mainGtf.substring(mainGtf.lastIndexOf('/')+1,mainGtf.length())
dbName = db.baseName
"""
if [ -e "${custom_gtf}" ]; then
	cp ${custom_gtf} ${nameCustomGtf}
    gtfFile="${nameCustomGtf} "
elif [ -e "${params.gtf}" ]; then
	cp ${params.gtf} ${nameMainGtf}
    gtfFile="${nameMainGtf} "
fi
echo "gtf name: \${gtfFile}"

## merge groupa.txt and groupb.txt and mv related files into cond folder
# find ${cond}/group*.txt | xargs -I{} sh -c "cat {}; echo ''" > allgroups.txt
# while IFS="" read -r line || [ -n "\$line" ]
# do
#   mv \${line}* ${cond}/.
# done < allgroups.txt
# mv  $db ${cond}/.

mv $ir_out_tab $out_tab $bam $bai $db ${cond}/.
cp \$gtfFile ${cond}/.

## clean "_conds_" prefix from condition name
cond="${cond}"
realcond=\${cond#"_conds_"}
mv ${cond} \$realcond

cd \$realcond
tmpdir="${baseDir}/work"
export TMPDIR=\$tmpdir
perl ${params.dummyai_path} --gtf \$gtfFile  --name ${dbName} ${PSI_sigma_parameters} > log.txt 2>&1
rm $bam $bai \$gtfFile $db
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
