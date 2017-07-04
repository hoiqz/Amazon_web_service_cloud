<?php
if (!class_exists('S3')) require_once 'S3.php';
require_once '/var/www/gis/config_jobsch.php';
require_once '/var/www/gis/sdk.class.php';

#################################
# Setup some que and nucket variables
#################################
//$fh=fopen('/mnt3/tmp/log.txt','w');
$msg_handle=$argv[1];
#$user_email=$argv[2]; // for sending email notification
$uploadfile=$argv[2]; // function takes in the name of the file
$filepid=$argv[3];
/*
   fwrite($fh, "msg handle is this $msg_handle\n");
   fwrite($fh, "upload file is this $uploadfile\n");
   fwrite($fh, "file id is $filepid\n");

   if(array_count_values($argv) >3){
   $printedvalues=print_r($argv);
//	fwrite($fh,$printedvalues);
$email=$argv[4];
}

foreach ($argv as $key => $value) {
$toFile = "Key: $key; Value: $value \n";
// write to file 
fwrite($fh, "$toFile") or die('Could not write to file'); 
// close file 
}
fclose($fh);
 */
#$bucket='pauline.c.ng.lab';
$bucketName='user_upload_bucket';
$queue_url='https://queue.amazonaws.com/665323014563/sift_queue';

// Delete the completed job
$sqs = new AmazonSQS(awsAccessKey, awsSecretKey);
$delete_msg=$sqs->delete_message($queue_url,$msg_handle);


// Delete the file in S3
// Instantiate the class
$s3 = new S3(awsAccessKey, awsSecretKey);
if ($s3->deleteObject($bucketName, $uploadFile)) {
	echo "S3::deleteObject(): Deleted file\n";
}
else{
	echo "S3::deleteObject(): Failed to delete file\n";
}
#######################################
# Upload the fileto S3 bucket for user
# to download 
######################################
$uploadfile2="/mnt3/tmp/$filepid\_nssnv_$uploadfile\_predictions.tsv";
// chedck if there are messges in the que
if ($s3->putObjectFile($uploadFile2, $bucketName, baseName($uploadFile2), S3::ACL_PUBLIC_READ_WRITE)) {
	echo "S3::putObjectFile(): File copied to {$bucketName}/".baseName($uploadFile2).PHP_EOL;
} else {
	echo "S3::putObjectFile(): Failed to copy file\n";
}

#####################################
# Poll for Message
#####################################
while(1){
	$message_in_queue=`php /var/www/gis/receive_job.php`;
	if(preg_match("/There is nothing in the queue/",$message_in_queue)){
		sleep(5); #try again
	}else{
		break;
	}
}

?>
