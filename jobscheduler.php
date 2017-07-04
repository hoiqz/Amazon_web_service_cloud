<html>
<body>
<?php
if (!class_exists('S3')) require_once 'S3.php';
require_once '/var/www/gis/config_jobsch.php';
require_once '/var/www/gis/sdk.class.php';
include 'siftapp.php';
#######################
# Capturing the upload
######################

$bucket = 'user_upload_bucket';
$queue_URL='https://queue.amazonaws.com/665323014563/sift_queue';
$original_name_no_uniq=$_FILES['CHR_file']['name'];
$original_name=uniqid($original_name_no_uniq);
$fileholding_area="/var/www/gis/temp/$original_name";
#print_r(array_count_values($_POST);
#$fileholding_area="/mnt3/tmp/upload/$original_name";
if(isset($_FILES['CHR_file'])){
	echo "file upload detected\n";
	$filename=$_FILES['CHR_file'];
	$movefile=move_uploaded_file($filename['tmp_name'],$fileholding_area);
}
if($movefile){
	echo " hurray successfully moved file \n";
}
else{
	echo "problem moving file ".$_FILES['CHR_file']['error'];
}

#######################
# upload to S3 and return the URL
########################

$uploadFile = $fileholding_area; // must be a full path
// Check if our upload file exists
if (!file_exists($uploadFile) || !is_file($uploadFile))
	exit("\nERROR: No such file: $uploadFile\n\n");

	// Check for CURL
	if (!extension_loaded('curl') && !@dl(PHP_SHLIB_SUFFIX == 'so' ? 'curl.so' : 'php_curl.dll'))
	exit("\nERROR: CURL extension not loaded\n\n");

	$s3 = new S3(awsAccessKey, awsSecretKey);

	// List your buckets:
	echo "S3::listBuckets(): ".print_r($s3->listBuckets(), 1)."\n";

	//upload to our bucket
	if ($s3->putObjectFile($uploadFile, $bucket, baseName($uploadFile), S3::ACL_PUBLIC_READ_WRITE)) {
		echo "S3::putObjectFile(): File copied to {$bucket}/".baseName($uploadFile).PHP_EOL;
	shell_exec("rm /var/www/gis/temp/$original_name");
	} else {
		echo "S3::putObjectFile(): Failed to copy file\n";
	}
// Get object info
$info = $s3->getObjectInfo($bucket, baseName($uploadFile));
echo "S3::getObjecInfo(): Info for {$bucket}/".baseName($uploadFile).': '.print_r($info, 1);
$fileURL="https://$bucket.s3.amazomaws.com/".baseName($uploadFile);
unlink($fileholding_area);
########################
# Determine the application
#######################
if(isset($_POST['app']) && !empty($_POST['app'])){


	// read determine what app to run
	$APP=$_POST['app'];
	$Q_msg=generate_message($APP,$fileURL);
	echo "the message is $Q_msg \n";
}
else{
	die(" please choose the app you want to run");
}

###########################
# SEND A MESSAGE TO THE QUERE
###########################

$sqs = new AmazonSQS(awsAccessKey, awsSecretKey);
$sendmsg=$sqs->send_message($queue_URL,$Q_msg);
if($sendmsg->isOK()){
	echo "successfully queued job!";
}
else{
	echo "queue failed!";
}



function generate_message($app,$s3location){
	switch(strtoupper($app)){
		case 'EXOME_NSSNVS':
			return exome_nssnvs($s3location);
			break;
		case 'BLINK':
			return blink($s3location);
			break;
		case 'EXOME_INDELS':
			return exome_indels($s3location);
			break;
		case 'SNP_CLASSIFIER':
			return snp_classifier($s3location);
			break;
		case 'SEQUENCE':
			return sequence($s3location);
			break;
		case 'BATCH':
			return batch($s3location);
			break;
		case 'DBSNP':
			return dbsnp($s3location);
			break;
		case 'SNPEDIA_PILEUP':
			return snpedia_pileup($s3location);
			break;
	}
}
/*
   echo " processing the POST request< /br>";
//print_r($_POST);
//print_r($_FILES);

echo "using array_key() to store variables";

foreach(array_keys($_POST) as $key){
echo "key name is $key";
echo "value is $_POST[$key]";
}
 */
?>
</body>
</html>
