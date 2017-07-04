#!usr/bin/perl
#########################################################
# subroutine for users using command line interface (CLI)
# Used to handle the message handle and uniq file id to 
# do with AWS que msg deleteion + S3 file deletion
# author: Hoi
#########################################################

# usage: my($msg_handle,$objectname)=&webservice(%names);

use strict;
sub webservice {
    chomp @_;
    my ($form_field_href) = @_;
    my %form_field = %{$form_field_href};
    my $message_handle='';
    my $s3filename='';
    
    $message_handle=$form_field{msg_handle};
    $s3filename=$form_field{objectname};
chomp($message_handle);
chomp($s3filename);
$message_handle=~s/^\n//g;
$s3filename=~s/^\n//g;
    return ($message_handle,$s3filename);
    
}

1;
