<?php
if (!class_exists('S3')) require_once 'S3.php';
require_once '/var/www/gis/config_jobsch.php';
require_once '/var/www/gis/sdk.class.php';
$bucket = 'pauline.c.ng.lab';
#$subfolder= '/mnt1/db/hg19_2011_10_27/';
//$queue_URL='https://queue.amazonaws.com/665323014563/sift_queue';
$uploadFile = $argv[1]; // must be a full path
#$fullpath=$subfolder.baseName($uploadFile);
echo "full path file name is $fullpath\n";
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
        #if ($s3->putObjectFile($uploadFile, $bucket, $fullpath, S3::ACL_PUBLIC_READ_WRITE)) {
        if ($s3->putObjectFile($uploadFile, $bucket, baseName($uploadFile), S3::ACL_PUBLIC_READ_WRITE)) {
                echo "S3::putObjectFile(): File copied to {$bucket}/".baseName($uploadFile).PHP_EOL;
        //shell_exec("rm /var/www/gis/temp/$original_name");
        } else {
                echo "S3::putObjectFile(): Failed to copy file\n";
        }
// Get object info
$info = $s3->getObjectInfo($bucket, baseName($uploadFile));
echo "S3::getObjecInfo(): Info for {$bucket}/".baseName($uploadFile).': '.print_r($info, 1);
$fileURL="https://$bucket.s3.amazomaws.com/".baseName($uploadFile);
?>
