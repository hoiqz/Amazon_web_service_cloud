<?php
####################################################################
# Author: Hoi Qiangze
# This script is used for jobscheduler
# Contains all the functions to create the queue message for 
# individual applications
# Editied:12/2/2011
###################################################################3
function exome_nssnvs($s3location)
{
	//Take in the post variables
	$ORGANISM=$_POST['organism'];
	$EMAIL=$_POST['address'];
	$MULTI_TRANSCRIPTS=$_POST['multi_transcripts'];
	$EBL_GENE_ID=$_POST['oo1'];
	$GENE_NAME=$_POST['oo2'];
	$GENE_DESC=$_POST['oo3'];
	$EBL_TRANS_STATUS=$_POST['oo4'];
	$EBL_PROTEIN_FAMILY_ID=$_POST['oo5'];
	$EBL_PROTEIN_FAMILY_ID_DESC=$_POST['oo6'];
	$PROTEIN_FAMILY_SIZE=$_POST['oo7'];
	$KA_KS_HUMAN_MOUSE=$_POST['oo8'];
	$KA_KS_HUMAN_MACAQUE=$_POST['oo9'];
	$OMIM_DISEASE=$_POST['oo10'];
	$ALL_WEIGHTED_AVERAGE=$_POST['oo11'];
	$CEU_POPULATION=$_POST['oo12'];
	$HCB_POPULATION=$_POST['oo13'];
	$JPT_POPULATION=$_POST['oo14'];
	$YRI_POPULATION=$_POST['oo15'];
	$ALL_AVAIL_POPULATIONS=$_POST['oo16'];
	$EUROPEAN_POPULATION=$_POST['oo17'];
	$EAST_ASIAN_POPULATION=$_POST['oo18'];
	$AFRICAN_POPULATION=$_POST['oo19'];

	//create the message now + location of s3
	return $message="command=curl --max-time 600 --connect-timeout 600 -F \"isCLI=true\" -F \"address=$EMAIL\" -F \"organism=$ORGANISM\" -F \"multi_transcripts=$MULTI_TRANSCRIPTS\" -F \"oo1=$EBL_GENE_ID\" -F \"oo2=$GENE_NAME\" -F \"oo3=$GENE_DESC\" -F \"oo4=$EBL_TRANS_STATUS\" -F \"oo5=$EBL_PROTEIN_FAMILY_ID\" -F \"oo6=$EBL_PROTEIN_FAMILY_ID_DESC\" -F \"oo7=$PROTEIN_FAMILY_SIZE\" -F \"oo8=$KA_KS_HUMAN_MOUSE\" -F \"oo9=$KA_KS_HUMAN_MACAQUE\" -F \"oo10=$OMIM_DISEASE\" -F \"oo11=$ALL_WEIGHTED_AVERAGE\" -F \"oo12=$CEU_POPULATION\" -F \"oo13=$HCB_POPULATION\" -F \"oo14=$JPT_POPULATION\" -F \"oo15=$YRI_POPULATION\" -F \"oo16=$ALL_AVAIL_POPULATIONS\" -F \"oo17=$EUROPEAN_POPULATION\" -F \"oo18=$EAST_ASIAN_POPULATION\" -F \"oo19=$AFRICAN_POPULATION\" -F \"CHR_file=@\nfile=$s3location\napp=/sift-bin/Extended_SIFT_feed_to_chr_coords.pl";
	#return $message="curl --max-time 600 --connect-timeout 600 -d \"address=$EMAIL&CHR_file=$s3location&organism=$ORGANISM&multi_transcripts=$MULTI_TRANSCRIPTS&oo1=$EBL_GENE_ID&oo2=$GENE_NAME&oo3=$GENE_DESC&oo4=$EBL_TRANS_STATUS&oo5=$EBL_PROTEIN_FAMILY_ID&oo6=$EBL_PROTEIN_FAMILY_ID_DESC&oo7=$PROTEIN_FAMILY_SIZE&oo8=$KA_KS_HUMAN_MOUSE&oo9=$KA_KS_HUMAN_MACAQUE&oo10=$OMIM_DISEASE&oo11=$ALL_WEIGHTED_AVERAGE&oo12=$CEU_POPULATION&oo13=$HCB_POPULATION&oo14=$JPT_POPULATION&oo15=$YRI_POPULATION&oo16=$ALL_AVAIL_POPULATIONS&oo17=$EUROPEAN_POPULATION&oo18=$EAST_ASIAN_POPULATION&oo19=$AFRICAN_POPULATION\" ";
}

function blink($s3location){}
function exome_indels($s3location){
	$ORGANISM=$_POST['organism'];
	$EMAIL=$_POST['address'];
	$ALL_TRANSCRIPTS=$_POST['ALL_TRANSCRIPTS'];
	$GENE_NAME=$_POST['oo1'];
	$GENE_DESC=$_POST['oo2'];
	$EBL_PROTEIN_FAMILY_ID=$_POST['oo3'];
	$EBL_PROTEIN_FAMILY_ID_DESC=$_POST['oo4'];
	$EBL_TRANS_STATUS=$_POST['oo5'];
	$PROTEIN_FAMILY_SIZE=$_POST['oo6'];
	$KA_KS_HUMAN_MOUSE=$_POST['oo7'];
	$KA_KS_HUMAN_MACAQUE=$_POST['oo8'];
	$OMIM_DISEASE=$_POST['oo9'];
	$OVERLAP_W_1000_GENOME=$_POST['oo10'];
	return $message="command=curl --max-time 600 --connect-timeout 600 -F \"isCLI=true\" -F \"address=$EMAIL\" -F \"ALL_TRANSCRIPTS=$ALL_TRANSCRIPTS\" -F \"organism=$ORGANISM\" -F \"oo1=$GENE_NAME\" -F \"oo2=$GENE_DESC\" -F \"oo3=$EBL_PROTEIN_FAMILY_ID\" -F \"oo4=$EBL_PROTEIN_FAMILY_ID_DESC\" -F \"oo5=$EBL_TRANS_STATUS\" -F \"oo6=$PROTEIN_FAMILY_SIZE\" -F \"oo7=$KA_KS_HUMAN_MOUSE\" -F \"oo8=$KA_KS_HUMAN_MACAQUE\" -F \"oo9=$OMIM_DISEASE\" -F \"oo10=$OVERLAP_W_1000_GENOME\" -F \"CHR_file=@\nfile=$s3location\napp=/sift-bin/SIFT_chr_coords_indels_submit.pl";
	#return $message="curl --max-time 600 --connect-timeout 600 -d \"address=$EMAIL&CHR_file=$s3location&ALL_TRANSCRIPTS=$ALL_TRANSCRIPTS&organism=$ORGANISM&oo1=$GENE_NAME&oo2=$GENE_DESC&oo3=$EBL_PROTEIN_FAMILY_ID&oo4=$EBL_PROTEIN_FAMILY_ID_DESC&oo5=$EBL_TRANS_STATUS&oo6=$PROTEIN_FAMILY_SIZE&oo7=$KA_KS_HUMAN_MOUSE&oo8=$KA_KS_HUMAN_MACAQUE&oo9=$OMIM_DISEASE&oo10=$OVERLAP_W_1000_GENOME\" ";
}
function snp_classifier($s3location){
$BUILD=$_POST['build'];
$DATALINE=$_POST['data'];
$FILEFORMAT=$_POST['fileformat'];
return $message="command=curl --max-time 600 --connect-timeout 600 -F \"isCLI=true\" -F \"build=$BUILD\" -F \"fileformat=$FILEFORMAT\" -F \"data=$DATALINE\" -F \"file=@\nfile=$s3location\napp=/sift-bin/SIFT_intersect_cds_submit.pl";
#return $message="curl --max-time 600 --connect-timeout 600 -d \"build=$BUILD&file=${INPUT_DATA1}&fileformat=$FILEFORMAT&data=$DATALINE\" ";
}
function sequence($s3location){}
function batch($s3location){}
function dbsnp($s3location){}
function snpedia_pileup($s3location){}

?>
