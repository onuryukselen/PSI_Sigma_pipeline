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
Channel.value(params.custom_gtf).into{g_20_gtfFilePath_g_26;g_20_gtfFilePath_g_59;g_20_gtfFilePath_g_109;g_20_gtfFilePath_g_113}
Channel.value(params.gtf).into{g_28_gtfFilePath_g_27;g_28_gtfFilePath_g_60;g_28_gtfFilePath_g_110;g_28_gtfFilePath_g_114}


process rename {

input:
 set val(name_bam), file(bam) from g_1_bam_file_g_18
 set val(name_tab), file(tab) from g_2_outputFileTab_g_18

output:
 file "*.bam"  into g_18_bam_file_g_26, g_18_bam_file_g_27, g_18_bam_file_g_44, g_18_bam_file_g_45, g_18_bam_file_g_59, g_18_bam_file_g_60
 file "*.bai"  into g_18_bam_index_g_26, g_18_bam_index_g_27, g_18_bam_index_g_44, g_18_bam_index_g_45, g_18_bam_index_g_59, g_18_bam_index_g_60
 file "*.tab"  into g_18_tab_file_g_26, g_18_tab_file_g_27, g_18_tab_file_g_59, g_18_tab_file_g_60

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
 file "*"  into g_54_all_groups_g_60
 val all_groups_path  into g_54_all_groups_path_g_108, g_54_all_groups_path_g_113, g_54_all_groups_path_g_114

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


process prepare_groups_stringtie {

input:
 val all_groups_path from g_54_all_groups_path_g_108

output:
 file "*"  into g_108_all_groups_g_59

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


process PSI_Sigma_prep {


output:
 val chromosome_list  into g_77_name_g_26, g_77_name_g_27
 val PSI_sigma_parameters  into g_77_parameters_g_59, g_77_parameters_g_60
 val gct_parameters  into g_77_gct_parameters_g_111, g_77_gct_parameters_g_112, g_77_gct_parameters_g_113, g_77_gct_parameters_g_114

when:
(params.run_PSI_Sigma && (params.run_PSI_Sigma == "yes")) || !params.run_PSI_Sigma

script:
PSI_sigma_parameters = params.PSI_Sigma_prep.PSI_sigma_parameters
gct_parameters = params.PSI_Sigma_prep.gct_parameters

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
 val chr_name from g_77_name_g_27.flatten()
 val custom_gtf from g_28_gtfFilePath_g_27

output:
 file "*.txt"  into g_27_groups
 file "*.tmp"  into g_27_tmp_file_g_29
 val db_name  into g_27_db_name_g_29

container "dolphinnext/psi_sigma_pipeline:3.0"

when:
custom_gtf.toString().indexOf("/") > -1

script:
custom_gtf = custom_gtf.toString()
db_name = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "PSIsigma1d9l_StringTie": "PSIsigma1d9l" 
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
	if (filename =~ /${db_name}.bed$/) "PSI_sigma_alone/$filename"
	else if (filename =~ /${db_name}.db$/) "PSI_sigma_alone/$filename"
}

input:
 file tmp from g_27_tmp_file_g_29.collect()
 val db_name from g_27_db_name_g_29.collect()

output:
 file "${db_name}.bed"  into g_29_bed_file
 file "${db_name}.db"  into g_29_db_file_g_45, g_29_db_file_g_60, g_29_db_file_g_114

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
 file "*.IR.out.tab"  into g_45_tab_file_g_60

container "dolphinnext/psi_sigma_pipeline:3.0"

script:
"""
perl ${params.PSIsigma_ir_path} $db $bam 1 
"""
}

g_28_gtfFilePath_g_60= g_28_gtfFilePath_g_60.ifEmpty([""]) 

//* params.dummyai_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process PSI_Sigma_run {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${realcond}$/) "PSI_sigma_alone/$filename"
}

input:
 file cond from g_54_all_groups_g_60.flatten()
 file ir_out_tab from g_45_tab_file_g_60.collect()
 file out_tab from g_18_tab_file_g_60.collect()
 file bam from g_18_bam_file_g_60.collect()
 file bai from g_18_bam_index_g_60.collect()
 val custom_gtf from g_28_gtfFilePath_g_60
 file db from g_29_db_file_g_60
 val PSI_sigma_parameters from g_77_parameters_g_60

output:
 file "${realcond}"  into g_60_outputDir_g_110, g_60_outputDir_g_114

container "dolphinnext/psi_sigma_pipeline:3.0"

script:
mainGtf = params.gtf.toString()
custom_gtf = custom_gtf.toString()
nameCustomGtf = custom_gtf.substring(custom_gtf.lastIndexOf('/')+1,custom_gtf.length())
nameMainGtf = mainGtf.substring(mainGtf.lastIndexOf('/')+1,mainGtf.length())
realcond = cond.toString() - '_conds_' 
dbName = db.baseName 
realDBname = realcond + "_" + dbName
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
mv ${dbName}.db ${realDBname}.db
tmpdir="${baseDir}/work"
export TMPDIR=\$tmpdir
perl ${params.dummyai_path} --gtf \$gtfFile  --name ${realDBname} ${PSI_sigma_parameters} > log.txt 2>&1
rm $bam $bai \$gtfFile ${realDBname}.db
"""
}


process generate_gct {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_summary.*$/) "PSI_sigma_alone/$filename"
}

input:
 file psi_sigma from g_60_outputDir_g_114.collect()
 val all_groups from g_54_all_groups_path_g_114
 val custom_gtf from g_28_gtfFilePath_g_114
 file db_file from g_29_db_file_g_114
 val gct_parameters from g_77_gct_parameters_g_114

output:
 val "yes"  into g_114_run_process_g_110
 file "*_summary*"  into g_114_outputDir_g_110, g_114_outputDir_g_112
 val suffix  into g_114_suffix_g_112

script:
nread=gct_parameters.split()[0]
dPSI=gct_parameters.split()[1]
adjp=gct_parameters.split()[2]
direction=gct_parameters.split()[3]
suffix=db_file.getBaseName()
dbgtf = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? custom_gtf : params.gtf
gtfName = file(params.gtf).getName().toString()
dbgtfName = file(custom_gtf).getName().toString()
all_groupsName = file(all_groups).getName().toString()

"""
ln -s ${params.gtf} $gtfName
if [ "${gtfName}" != "${dbgtfName}" ]; then
	ln -s ${custom_gtf} $dbgtfName
fi
ln -s $all_groups $all_groupsName
# mv $db_file ${psi_sigma}/.
perl /home/share/tools/DolphinNext/rnaseq/src/gct_v5.1.pl $all_groupsName ${gtfName} $dbgtfName ${suffix}.db $suffix $nread $dPSI $adjp $direction

"""
}


process generate_sequence_logos {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "PSI_sigma_alone/$filename"
}

input:
 file gct from g_114_outputDir_g_112
 val gct_parameters from g_77_gct_parameters_g_112
 val suffix from g_114_suffix_g_112

output:
 file "*"  into g_112_outputDir_g_110

script:
nread=gct_parameters.split()[0]
dPSI=gct_parameters.split()[1]
adjp=gct_parameters.split()[2]
direction=gct_parameters.split()[3]

"""
perl /home/share/tools/DolphinNext/rnaseq/src/gct2fasta_v2.pl ${params.genome} ${params.gtf} ${suffix}_r${nread}_${direction}_dPSI${dPSI}_adjp${adjp}_summary 6 18 inclusion ${dPSI} ${adjp}
perl /home/share/tools/DolphinNext/rnaseq/src/gct2fasta_v2.pl ${params.genome} ${params.gtf} ${suffix}_r${nread}_${direction}_dPSI${dPSI}_adjp${adjp}_summary 6 18 exclusion ${dPSI} ${adjp}

"""
}


process generate_gene_and_cluster_tables {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /output\/.*$/) "PSI_sigma_alone/$filename"
}

input:
 file psi_sigma from g_60_outputDir_g_110.collect()
 val run_process from g_114_run_process_g_110
 file gctFiles from g_114_outputDir_g_110
 file out from g_112_outputDir_g_110
 val custom_gtf from g_28_gtfFilePath_g_110

output:
 file "output/*"  into g_110_outputDir

script:
suffix = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "sorted.annotated.txt" : "sorted.txt"
suffix2 = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? ".annotated" : ""
"""
mkdir output
cd output
find ../*/* -name '*ir3.${suffix}' -exec ln -s {} . \\;
dPSI=0
cp="0.05"
for fn in `ls *ir3.${suffix}`; do
	perl /home/share/tools/DolphinNext/rnaseq/src/gene_v1.pl \$fn \$dPSI \$cp 3
	perl /home/share/tools/DolphinNext/rnaseq/src/cluster_v1.pl \$fn \$dPSI \$cp 3
done
tar zcvf gene_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz gene.table.dPSI\$dPSI.adjp\$cp.*
tar zcvf cluster_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz cluster.table.dPSI\$dPSI.adjp\$cp.*

dPSI=0
cp="1"
for fn in `ls *ir3.${suffix}`; do
	perl /home/share/tools/DolphinNext/rnaseq/src/gene_v1.pl \$fn \$dPSI \$cp 3
	perl /home/share/tools/DolphinNext/rnaseq/src/cluster_v1.pl \$fn \$dPSI \$cp 3
done
tar zcvf gene_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz gene.table.dPSI\$dPSI.adjp\$cp.*
tar zcvf cluster_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz cluster.table.dPSI\$dPSI.adjp\$cp.*

find -type l -exec rm {} \\;
"""
}

//* params.PSIsigma_db_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process create_db_stringtie {

input:
 file bam from g_18_bam_file_g_26.collect()
 file bai from g_18_bam_index_g_26.collect()
 file tab from g_18_tab_file_g_26.collect()
 val chr_name from g_77_name_g_26.flatten()
 val custom_gtf from g_20_gtfFilePath_g_26

output:
 file "*.txt"  into g_26_groups
 file "*.tmp"  into g_26_tmp_file_g_8
 val db_name  into g_26_db_name_g_8

container "dolphinnext/psi_sigma_pipeline:3.0"

when:
custom_gtf.toString().indexOf("/") > -1

script:
custom_gtf = custom_gtf.toString()
db_name = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "PSIsigma1d9l_StringTie": "PSIsigma1d9l" 
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
	if (filename =~ /${db_name}.bed$/) "PSI_Sigma/$filename"
	else if (filename =~ /${db_name}.db$/) "PSI_Sigma/$filename"
}

input:
 file tmp from g_26_tmp_file_g_8.collect()
 val db_name from g_26_db_name_g_8.collect()

output:
 file "${db_name}.bed"  into g_8_bed_file
 file "${db_name}.db"  into g_8_db_file_g_44, g_8_db_file_g_59, g_8_db_file_g_113

script:
db_name = db_name[0]
"""	
cat *.db.tmp > ${db_name}.db
cat *.bed.tmp > ${db_name}.bed
"""

}

//* params.PSIsigma_ir_path =  ""  //* @input

process create_intronic_read_stringtie {

input:
 file bam from g_18_bam_file_g_44
 file bai from g_18_bam_index_g_44
 file db from g_8_db_file_g_44

output:
 file "*.IR.out.tab"  into g_44_tab_file_g_59

container "dolphinnext/psi_sigma_pipeline:3.0"

script:
"""
perl ${params.PSIsigma_ir_path} $db $bam 1 
"""
}

g_20_gtfFilePath_g_59= g_20_gtfFilePath_g_59.ifEmpty([""]) 

//* params.dummyai_path =  ""  //* @input
//* params.gtf =  ""  //* @input

process PSI_Sigma_run_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${realcond}$/) "PSI_Sigma/$filename"
}

input:
 file cond from g_108_all_groups_g_59.flatten()
 file ir_out_tab from g_44_tab_file_g_59.collect()
 file out_tab from g_18_tab_file_g_59.collect()
 file bam from g_18_bam_file_g_59.collect()
 file bai from g_18_bam_index_g_59.collect()
 val custom_gtf from g_20_gtfFilePath_g_59
 file db from g_8_db_file_g_59
 val PSI_sigma_parameters from g_77_parameters_g_59

output:
 file "${realcond}"  into g_59_outputDir_g_109, g_59_outputDir_g_113

container "dolphinnext/psi_sigma_pipeline:3.0"

script:
mainGtf = params.gtf.toString()
custom_gtf = custom_gtf.toString()
nameCustomGtf = custom_gtf.substring(custom_gtf.lastIndexOf('/')+1,custom_gtf.length())
nameMainGtf = mainGtf.substring(mainGtf.lastIndexOf('/')+1,mainGtf.length())
realcond = cond.toString() - '_conds_' 
dbName = db.baseName 
realDBname = realcond + "_" + dbName
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
mv ${dbName}.db ${realDBname}.db
tmpdir="${baseDir}/work"
export TMPDIR=\$tmpdir
perl ${params.dummyai_path} --gtf \$gtfFile  --name ${realDBname} ${PSI_sigma_parameters} > log.txt 2>&1
rm $bam $bai \$gtfFile ${realDBname}.db
"""
}


process generate_gct_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*_summary.*$/) "PSI_Sigma/$filename"
}

input:
 file psi_sigma from g_59_outputDir_g_113.collect()
 val all_groups from g_54_all_groups_path_g_113
 val custom_gtf from g_20_gtfFilePath_g_113
 file db_file from g_8_db_file_g_113
 val gct_parameters from g_77_gct_parameters_g_113

output:
 val "yes"  into g_113_run_process_g_109
 file "*_summary*"  into g_113_outputDir_g_109, g_113_outputDir_g_111
 val suffix  into g_113_suffix_g_111

script:
nread=gct_parameters.split()[0]
dPSI=gct_parameters.split()[1]
adjp=gct_parameters.split()[2]
direction=gct_parameters.split()[3]
suffix=db_file.getBaseName()
dbgtf = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? custom_gtf : params.gtf
gtfName = file(params.gtf).getName().toString()
dbgtfName = file(custom_gtf).getName().toString()
all_groupsName = file(all_groups).getName().toString()

"""
ln -s ${params.gtf} $gtfName
if [ "${gtfName}" != "${dbgtfName}" ]; then
	ln -s ${custom_gtf} $dbgtfName
fi
ln -s $all_groups $all_groupsName
# mv $db_file ${psi_sigma}/.
perl /home/share/tools/DolphinNext/rnaseq/src/gct_v5.1.pl $all_groupsName ${gtfName} $dbgtfName ${suffix}.db $suffix $nread $dPSI $adjp $direction

"""
}


process generate_sequence_logos_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*$/) "PSI_Sigma/$filename"
}

input:
 file gct from g_113_outputDir_g_111
 val gct_parameters from g_77_gct_parameters_g_111
 val suffix from g_113_suffix_g_111

output:
 file "*"  into g_111_outputDir_g_109

script:
nread=gct_parameters.split()[0]
dPSI=gct_parameters.split()[1]
adjp=gct_parameters.split()[2]
direction=gct_parameters.split()[3]

"""
perl /home/share/tools/DolphinNext/rnaseq/src/gct2fasta_v2.pl ${params.genome} ${params.gtf} ${suffix}_r${nread}_${direction}_dPSI${dPSI}_adjp${adjp}_summary 6 18 inclusion ${dPSI} ${adjp}
perl /home/share/tools/DolphinNext/rnaseq/src/gct2fasta_v2.pl ${params.genome} ${params.gtf} ${suffix}_r${nread}_${direction}_dPSI${dPSI}_adjp${adjp}_summary 6 18 exclusion ${dPSI} ${adjp}

"""
}


process generate_gene_and_cluster_tables_stringtie {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /output\/.*$/) "PSI_Sigma/$filename"
}

input:
 file psi_sigma from g_59_outputDir_g_109.collect()
 val run_process from g_113_run_process_g_109
 file gctFiles from g_113_outputDir_g_109
 file out from g_111_outputDir_g_109
 val custom_gtf from g_20_gtfFilePath_g_109

output:
 file "output/*"  into g_109_outputDir

script:
suffix = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? "sorted.annotated.txt" : "sorted.txt"
suffix2 = (custom_gtf.indexOf("StringTie.sorted.gtf") > -1) ? ".annotated" : ""
"""
mkdir output
cd output
find ../*/* -name '*ir3.${suffix}' -exec ln -s {} . \\;
dPSI=0
cp="0.05"
for fn in `ls *ir3.${suffix}`; do
	perl /home/share/tools/DolphinNext/rnaseq/src/gene_v1.pl \$fn \$dPSI \$cp 3
	perl /home/share/tools/DolphinNext/rnaseq/src/cluster_v1.pl \$fn \$dPSI \$cp 3
done
tar zcvf gene_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz gene.table.dPSI\$dPSI.adjp\$cp.*
tar zcvf cluster_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz cluster.table.dPSI\$dPSI.adjp\$cp.*

dPSI=0
cp="1"
for fn in `ls *ir3.${suffix}`; do
	perl /home/share/tools/DolphinNext/rnaseq/src/gene_v1.pl \$fn \$dPSI \$cp 3
	perl /home/share/tools/DolphinNext/rnaseq/src/cluster_v1.pl \$fn \$dPSI \$cp 3
done
tar zcvf gene_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz gene.table.dPSI\$dPSI.adjp\$cp.*
tar zcvf cluster_level${suffix2}.dPSI\$dPSI.adjp\$cp.tar.gz cluster.table.dPSI\$dPSI.adjp\$cp.*

find -type l -exec rm {} \\;
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
