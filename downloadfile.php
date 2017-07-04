<?php
require_once '/var/www/gis/config_jobsch.php';
require_once '/var/www/gis/sdk.class.php';
$bucket='user_upload_bucket';
$db_bucket='pauline.c.ng.lab';

$s3 = new AmazonS3(awsAccessKey, awsSecretKey);
$my_downloadlist=$s3->get_object_list($db_bucket);
print_r($my_downloadlist);
$thousand_genome_path="/mnt1/db/hg19_2011_10_27/1000genomes";
foreach ($my_downloadlist as $file){
	print "file: $file\n";
	$filepath="mnt1/db/hg19_2011_10_27/$file";
	$handle=fopen($filepath,'w+');
	$s3 = new AmazonS3(awsAccessKey, awsSecretKey);
	$dlfile=$s3->get_object(
			$db_bucket,
			$file,
			array('fileDownload'=> $handle));
	echo " downloaded file name is $file at $filepath\n";
	close($handle);
}
?>
