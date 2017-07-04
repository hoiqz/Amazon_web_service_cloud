<?php
require_once '/var/www/gis/config_jobsch.php';
require_once '/var/www/gis/sdk.class.php';
$bucket='user_upload_bucket';
$db_bucket='pauline.c.ng.lab';
$queue_url='https://queue.amazonaws.com/665323014563/sift_queue';
$handle=fopen("/var/www/gis/temp/debug_log",'w+');


##################################
# download the databse sqlite
# file from S3 bucket
##################################
`sh /home/hoi/mountdrive.sh`;

#################################
# Receive MEssage from the queue
#################################

$sqs = new AmazonSQS(awsAccessKey, awsSecretKey);
$response= $sqs->receive_message($queue_url,array('VisibilityTimeout'=>10));
$job_command=$response->body->ReceiveMessageResult->Message->Body;
//print_r($response);
print_r($job_command);

if(isset($job_command)&& !empty($job_command)){

	$receipt_handle=$response->body->ReceiptHandle(0);
	echo " i have grabbed a message".$job_command." with this receipt handle: $receipt_handle";
	$prtstr= " i have grabbed a message".$job_command." with this receipt handle: $receipt_handle";
	fwrite($handle,$prtstr);
}
else{
	fwrite($handle,"nothing in queue\n");
	die("There is nothing in the queue. Quitting program. Cya!");
}
##################################
# download the user uploaded 
# file from S3 bucket
##################################

// extract the name of the file
$file_url=extracturl($job_command);
$filename=basename($file_url);
// $wgetcommand="wget --no-check-certificate $file_url";
#$wgetcommand="wget --no-check-certificate -O /var/www/gis/temp/$fileuniqid $file_url";
$filepath="/mnt3/tmp/$filename";
$handle=fopen($filepath,'w+');
$s3 = new AmazonS3(awsAccessKey, awsSecretKey);
//$response = $s3->get_object_acl($bucket, 'S25.nssnvs.test');
// $dlfile=$s3->get_object('pauline.c.ng.lab',$filename);
$dlfile=$s3->get_object(
		$bucket,
		$filename,
		array('fileDownload'=> $handle));
echo " downloaded file name is $filename at $filepath\n";
$prtstr=" downloaded file name is $filename at $filepath\n";
fwrite=($handle,$prtstr);

###########################################
# RUN THE JOB WITH THE MESSAGE AND THE FILE
###########################################

// formate the URl address to curl to
echo "Executing job command now!\n\n";
$execute_cmd=format_msg($job_command,$filepath,$filename,$receipt_handle);
$output=`$execute_cmd`;
echo "Job Output : ";
echo $output;
$prtstr="has job output:$output";
fwrite($handle,$prtstr);







########################################
# FUNCTIONS
########################################
function format_msg($msg_str,$file,$fn,$rh){
	$cmd=extractcmd($msg_str);
	preg_match('/(.*9=0\" )(.*)/',$cmd,$matches);
	$cmd1=$matches[1];
	$cmd2=$matches[2];

	$app=extractapp($msg_str);
	$instanceIP=`curl http://169.254.169.254/latest/meta-data/public-ipv4`; // to get the public ip
		//echo " MY INSTANCE HAS IP: $instanceIP\n";
		$instanceIP=preg_replace("/\./","-",$instanceIP);
	$instance_add="ec2-$instanceIP.compute-1.amazonaws.com";

	// msg_str contains message of the format:
	// command=xxxxx
	// file= url to file
	// app= sift-bin/program.pl

	$formatted_job_cmd=$cmd1."-F \"msg_handle=$rh\" -F \"objectname=$fn\" ".$cmd2.$file."\" ".$instance_add.$app;
	echo "$formatted_job_cmd\n";
	return $formatted_job_cmd;
}

function extracturl($msg_str){
	$msg_str_array=explode("\n",$msg_str);
	$msg_str_new=str_replace("file=","",$msg_str_array[1]);
	echo "url of file is : $msg_str_new\n";
	return $msg_str_new;
}
function extractcmd($msg_str){
	$msg_str_array=explode("\n",$msg_str);
	$msg_str_new=str_replace("command=","",$msg_str_array[0]);
	//        echo "url of file is : $msg_str_new\n";
	return $msg_str_new;
}
function extractapp($msg_str){
	$msg_str_array=explode("\n",$msg_str);
	$msg_str_new=str_replace("app=","",$msg_str_array[2]);
	//        echo "url of file is : $msg_str_new\n";
	return $msg_str_new;
}

?>
