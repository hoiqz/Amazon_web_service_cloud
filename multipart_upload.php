<?php
require_once '/var/www/gis/sdk.class_for_loadbalancer.php';
$bucket='pauline.c.ng.lab';
$filename='sift_hg19_db.tar.gz';
/*
Initiate a new multipart upload using initiate_multipart_upload().
Upload the parts using upload_part().
Complete the upload using complete_multipart_upload().
*/

// Define a megabyte
define('MB', 1024 * 1024);

// Instantiate the class
$s3 = new AmazonS3();


 // Initiate a new multipart upload
$response = $s3->initiate_multipart_upload($bucket, $filename, array(
        'contentType' => application/octet-stream,
    'acl' => AmazonS3::ACL_PUBLIC,
    'storage' => AmazonS3::STORAGE_STANDARD,
));

// Get the Upload ID
$upload_id = (string) $response->body->UploadId;
echo $upload_id;

// Get the list of pieces
$parts = $s3->get_multipart_counts(filesize($filename), 1000*MB);
 print_r($parts);

// Queue batch requests
foreach ($parts as $i => $part)
{
    $s3->batch()->upload_part($bucket, $filename, $upload_id, array(
        'expect' => '100-continue',
        'fileUpload' => "$filename",
        'partNumber' => ($i + 1),
        'seekTo' => (integer) $part['seekTo'],
        'length' => (integer) $part['length'],
    ));
}
$batch_responses = $s3->batch()->send();
print_r($batch_responses);
echo " now we list parts";
$parts2 = $s3->list_parts($bucket, $filename, $upload_id);
 print_r($parts2);
// Initiate a new multipart upload
$response = $s3->complete_multipart_upload($bucket, $filename, $upload_id, $parts2);

#print_r($response);
// Send batch requests

?>
