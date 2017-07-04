#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

require 'config.pl';

#my $bin             = "/usr/local/web/sift-bin";  #handled by config.pl
#my $tmp             = "/opt/www/sift/tmp";

$batch_file = @ARGV[0];

## Changes due to md5
my $list_of_pids = @ARGV[1];
my $pid_key = $ARGV[2];
my $multi_transcripts = $ARGV[3];
if (!defined $multi_transcripts) {
    $multi_transcripts = 0;
}

my $message_handle = "no";
my $s3filename = "none";
my $isCLI = $ARGV[4];
if ($isCLI && $isCLI eq "true") {
    $message_handle = $ARGV[5];
    $s3filename = $ARGV[6];
}

=pod
my $debug = "$tmp/batch.debug";
open(DEBUG, ">$debug");
print DEBUG "Just before Extended_SIFT_chr_coords_submit.pl\n";
=cut

open (BATCH_FILE ,"$batch_file");
while (<BATCH_FILE>){
	chomp;
	@elts = split /\t/, $_;
	$master_pid = $elts[0];
	$pid = $elts[1];
	$chrfile = $elts[2];
	$org = $elts[3];
	$seq_identity_filter = $elts[4];
	$COORD_SYSTEM = $elts[5];
	$last_partition = $elts[6];
	$output_options = $elts[7];
	$address = $elts[8];	
#	    system ("perl $bin/SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $list_of_pids $pid_key $address"); 

#	    system ("perl $bin/Extended_SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $list_of_pids $pid_key $multi_transcripts $address"); # multi_transcripts must come before address because users may not provide email address.
	    system ("perl $bin/Extended_SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $list_of_pids $pid_key $multi_transcripts $isCLI $message_handle $s3filename $address"); # multi_transcripts must come before address because users may not provide email address.


#	print DEBUG "perl $bin/Extended_SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $list_of_pids $pid_key $multi_transcripts $isCLI $message_handle $s3filename $address\n";


}
close(BATCH_FILE);

#close(DEBUG);
