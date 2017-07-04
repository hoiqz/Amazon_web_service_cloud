#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
#use lib '/usr/local/packages/sift/site_perl/5.8.8'; # for Template 
#use Template;
use List::Util qw[min max];
use Digest::MD5 qw (md5 md5_hex md5_base64);


$| = 1;
require 'SIFT_subroutines.pm';
require 'config.pl';
require 'utr_filtering.pl';
require 'CLI.pl';

system("umask 006");
#my $bin             = "/usr/local/web/sift-bin"; #handled by config.pl
#my $tmp             = "/opt/www/sift/tmp";

#my $pid             = $$;

# Change to MD5
my $date = `date`; chomp ($date);
my $toolname = "nssnv";
my $key = substr (md5_hex($$ . $date), 0,10);
$key =~ s/\///g;
my $pid = $key . "_" . $toolname;
my $FILENAME = "";


my $num_coords_per_split = 1000;
my $content = "";

print "Content-type: text/html\n\n";

# SIMNL DEBUG
my $debugfile = "$tmp/$pid\.feed.debug";
open(DEBUG, ">$debugfile");


if($ENV{"REQUEST_METHOD"} ne "POST") {
  $content = "This script should be referenced with a METHOD of POST\n";
  &finish_script($content, -1);
} # end if($ENV{"REQUEST_METHOD"} ne "POST")

read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
%names = general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );

print DEBUG "qstring: $QUERY_STRING\n";


### Modify pid if there is a filepath
if ($FILENAME ne "") {
    $pid .= "_" . $FILENAME;
}

my $outpage_url = "/sift-bin/catfile.csh?$tmp/$pid.outpage.html";


my $address;

chomp($names{address});

$names{address} =~ s/\s+//;

if($names{address} ne "") {
  $address = $names{address};
} # end if($names{address} ne "")

## Check for validity of user inputs
my $all_chr_file = $tmp . "/$pid.allchrfile";
#print "iiiii$organism<br>";

## $names{CHR_file} provides the contents of the file passed in by user
if($names{CHR} eq "" && $names{CHR_file} eq "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please enter some chromosome coordinates with substitutions.</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR} eq "" && $names{CHR_file} eq "")

my $organism = $names{organism};
$organism =~ s/^\s+//;
$organism =~ s/\s+$//;


# $organism is database directory (eg. Human_db_36)
if($organism =~ /Select Organism/i) {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please select organism after pressing back button on your browser.</p>\n";
  &finish_script($content, -1);
} # end if($organism =~ /Select Organism/i)


# SIMNL<2011-12-07>: Need to read in organism_map.txt to get location of directory
my @organism_dirs = `cat $bin/organism_map.txt`;
my %organism_dir_map = ();
foreach my $org (@organism_dirs) {
    chomp $org;
    my ($key, $dir) = split(/\t/, $org);
    $organism_dir_map{$key} = $dir;
}
my $database_loc = $organism_dir_map{$organism};
if (!defined($database_loc) || $database_loc eq "") {
    $database_loc = "hg19_2011_10_27";
}
# SIMNL<2011-12-07>: Cater to Hoi's CLI
my $isCLI = "false";
if (defined($names{"isCLI"})) {
    my $cli = $names{"isCLI"};
    $cli =~ s/^\n//;
    chomp $cli;
    $isCLI = $cli;
}
my $message_handle = "";
my $s3filename = "";
if ($isCLI ne "false") {
    ($message_handle, $s3filename) = webservice(\%names);
    my @keys = keys(%names);
    foreach my $k (@keys) {
	my $val = $names{$k};
	print DEBUG "$key -> $val\n";
    }
}


$seq_identity_filter = "90";
#Read input list of chromosome coordinates and add to all_chr_string
my $all_chr_string;

# This appears to be an extra check that there are at least coordinates provided?
if($names{CHR_file} !~ /\d/) {
  $names{CHR_file} = "";
} # end if($names{CHR_file} !~ /\d/)

if($names{CHR} !~ /\d/) {
  $names{CHR} = "";
} # end if($names{CHR} !~ /\d/)


my $oo1 = $names{oo1}==1 ? 1 : 0;
my $oo2 = $names{oo2}==1 ? 1 : 0;
my $oo3 = $names{oo3}==1 ? 1 : 0;
my $oo4 = $names{oo4}==1 ? 1 : 0;
my $oo5 = $names{oo5}==1 ? 1 : 0;
my $oo6 = $names{oo6}==1 ? 1 : 0;
my $oo7 = $names{oo7}==1 ? 1 : 0;
my $oo8 = $names{oo8}==1 ? 1 : 0;
my $oo9 = $names{oo9}==1 ? 1 : 0;
my $oo10 = $names{oo10}==1 ? 1 : 0;
my $oo11 = $names{oo11}==1 ? 1 : 0;
my $oo12 = $names{oo12}==1 ? 1 : 0;

# SIM changed
my $oo13 = $names{oo13}==1 ? 1 : 0; # HCB
my $oo14 = $names{oo14}==1 ? 1 : 0; # JPT
my $oo15 = $names{oo15}==1 ? 1 : 0; # YRI
my $oo16 = $names{oo16}==1 ? 1 : 0; # 1000 Genome All
my $oo17 = $names{oo17}==1 ? 1 : 0; # 1000 Genome European
my $oo18 = $names{oo18}==1 ? 1 : 0; # 1000 Genome E. Asian
my $oo19 = $names{oo19}==1 ? 1 : 0; # 1000 Genome W. African
my $oo20 = $names{oo20}==1 ? 1 : 0; # 1000 Genome S. Asian
my $oo21 = $names{oo21}==1 ? 1 : 0; # 1000 Genome American
# end SIM changed


my $show_multi_transcripts = $names{multi_transcripts} ? 1 : 0; # Display multiple transcripts if exists

my $output_options = "$oo1,$oo2,$oo3,$oo4,$oo5,$oo6,$oo7,$oo8,$oo9,$oo10,$oo11,$oo12";

# SIM changed
$output_options = $output_options . ",$oo13,$oo14,$oo15,$oo16,$oo17,$oo18,$oo19,$oo20,$oo21";
#$output_options = $output_options . ",$oo20,$oo21"; # For future when 1000 genome has S. Asian and American data.
# end SIM changed

# Allow only direct inputs or by file, but not both
if($names{CHR_file} ne "" && $names{CHR} ne "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please choose only one of the following input methods after clicking the back button on your browser</p>\n";
  $content .= "<p>1. Paste the input in the relevant textbox</p>\n";
  $content .= "<p>2. Upload the text file containing input data</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR_file} ne "" && $names{CHR} ne "")

if($names{CHR_file} eq "" && $names{CHR} eq "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please choose one of the following input methods after clicking the back button on your browser</p>\n";
  $content .= "<p>1. Paste the input in the relevant textbox</p>\n";
  $content .= "<p>2. Upload the text file containing input data</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR_file} eq "" && $names{CHR} eq "")

# Write all data into all chr file, changing all to uppercase.
open( CHR_FILE, ">$all_chr_file" );

if($names{CHR_file} ne "") {
  $names{CHR_file} =~ s/\r/\n/g;
  $names{CHR_file} =~ tr/A-Z/a-z/;

  if($names{CHR_file} =~ /\d/) {
    print CHR_FILE uc( $names{CHR_file} );
  } # end if($names{CHR_file} =~ /\d/)

  $all_chr_string = uc($names{CHR_file}), "\t";
  $input_method = "FILE";
} # end if($names{CHR_file} ne "")
else {
  $names{CHR} =~ tr/A-Z/a-z/;
  $names{CHR} =~ s/^\n//;
  $names{CHR} =~ s/\r/\n/g;
  print CHR_FILE uc($names{CHR});
  $all_chr_string .= uc($names{CHR}), "\t";
  $input_method = "TEXT";
} # end else

close(CHR_FILE);


### 
### If user input is > 30K, the website throws out a message that the file is too large and to do 
### (1) Intersect With Coding webpage before submitting to the SIFT_chr_coords_submit. 
###
my $number_of_inputs = 0;
open(CHR_FILE, "<$all_chr_file");
while(my $l = <CHR_FILE>) {
    $l =~ s/^\s+//;
    $l =~ s/\s+$//;
    if ($l ne "") {
	$number_of_inputs++;
    }
}
close(CHR_FILE);

# SIMNL<2011-12-07>: Changes to cater to Hoi's CLI
if ($isCLI eq "false" && $number_of_inputs > 100010) {
#if ($number_of_inputs > 100010) {
    print "The file has over 100000 inputs.<BR>";
    print "Please partition into separate files of not more than 10000 each and submit separately.<BR>";
    print "You might also wish to run your inputs here to restrict to coding regions: <a href=\"/www/SIFT_intersect_coding_submit.html\">Restrict to Coding Regions</a>.<BR>";
    exit(0);
}

# Take first line and determine if it is RESIDUE or SPACE COORDS
open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file for validation");
while(<CHR_FILE>) {
  if($_ =~ /\d/ && $_ =~ /\,/) {
    $first_line = $_;
    last;
  } # end if($_ =~ /\d/ && $_ =~ /\,/)
} # end while(<CHR_FILE>)
close(CHR_FILE);

if($first_line =~ /[\d+,X,x,Y,y],\d+,\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i) {
  $COORD_SYSTEM = "SPACE";
} # end if($first_line =~ /[\d+,X,x,Y,y],\d+,\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i)
elsif($first_line =~ /[\d+,x,X,y,Y],\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i) {
  $COORD_SYSTEM = "RESIDUE";
} # end elsif($first_line =~ /[\d+,x,X,y,Y],\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i)
else {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Incorrect input format. Please see <a href=\"\/chr_coords_example.html\">Sample format</a></p>";
  last;	
}

$content .= "<p>Your input data has been recognized to use <a href=\"\/chr_coords_example.php\" rel=\"external\">$COORD_SYSTEM based coordinate system</a>. Your job id is $key and is currently running.  Your job has been partitioned into datasets of $num_coords_per_split positions and the status of each job can be viewed in the <a href=\"$outpage_url\" rel=\"external\">SIFT results status page</a>. Once the status of a job is <span class=\"green\">'Complete'</span>, you may view or download the results. A partitioned job typically takes 6-7 min to complete.</p><p>Proceed to <a href=\"$outpage_url\" rel=\"external\">SIFT results status page.</a><p>Problems? Contact <a href=\"contact.pl\">us</a> with your job id.</p>\n";

=pod
&finish_script($content, 0); 

# Close the I/O handles
#close(STDin);
close(STDOUT);
#close(STDERR); 
=cut

#prepare output status page
my $outpage_content = '';
my $outpage_header = 'SIFT Results Status';
my $outpage_counter = 0;

open (OUTPAGE,">>$tmp/$pid.outpage.html") || die ("cannot open outpage");
$outpage_content .= "<p><a href=\"$outpage_url\">Refresh page</a></p>";
$outpage_content .= '<table>';
$heading = "<tr class=\"tableHeader\"><td><p>Job</p></td><td><p>Job size</p></td><td><p>Job ID</p></td><td><p>Job status</p></td><td><p>View results</p></td><td><p>Download results</p></td></tr>";
$outpage_content .= "$heading\n";


#if COORD SYSTEM is ABSOLUTE then convert chr file to space based.
if ($COORD_SYSTEM eq "RESIDUE"){
	open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
	open (CHR_FILE_NEW,">$all_chr_file.new") || die ("Cannot open new all chr file");
	while (<CHR_FILE>){
		chomp;
		if ($_ !~ /\d/){next}
		@elts = split /\,/, $_;
		$chr = @elts[0];
		$coord2 = @elts[1];
		$coord1 = $coord2-1;
		$orn = @elts[2];
		$alleles = @elts[3];
		$comment = @elts[4];

		# SIM Bug fix (2011-03-08): Comments may contain commas 
		if (scalar @elts >= 5) {
		    my @comments = splice(@elts, 4);
		    $comment = join(',', @comments);
		}
		# SIM end bug fix

		print CHR_FILE_NEW "$chr,$coord1,$coord2,$orn,$alleles,$comment\n";
	}
	close(CHR_FILE);
	close (CHR_FILE_NEW);
	system("mv $all_chr_file.new $all_chr_file");
}



##### Pauline add UTR and filter coding
$snpfile_to_gff = "perl $bin/snpfile_to_gff.pl"; # takes SIFT format file and puts it in gff format
$intersect_locations = "java -Xmx1000m -jar $bin/IntersectFeatures.jar "; # SIM, please modify;
#my $Variation_db_dir = "/usr/local/web/packages/db/$organism";

# SIMNL<2011-12-07>: Cater to Hoi's CLI
#my $Variation_db_dir = "/mnt1/db/$organism";
my $Variation_db_dir = "/mnt1/db/$database_loc";

print "DEBUG: $Variation_db_dir\n";

$utr_gff = $Variation_db_dir . "/non_coding/" . "GENE_INFO_v1.0_UTR_build37.gff";
$coding_gff = $Variation_db_dir . "/non_coding/" . "GENE_INFO_v1.0_CODING_build37.gff";
$flank_utr_gff = $Variation_db_dir . "/non_coding/" . "GENE_INFO_v1.0_flankUTR_build37.gff";
# create gff file
system("$snpfile_to_gff $all_chr_file");
$all_chr_gff_file = "$all_chr_file.gff";
# Intersect flanking UTR. order matters so that utr-gene_index will
# overwrite (UTR annotation priority over flanking UTR)
my $noncoding_annotation_file = $all_chr_file . ".noncoding_annot";
my $utr_gene_index_href = intersectFlankingUTR($all_chr_file, $intersect_locations, $all_chr_gff_file,
					       $flank_utr_gff, $flankutr_intersect_file,
					       $noncoding_annotation_file);
my %utr_gene_index = %{$utr_gene_index_href};

# Intersect UTR
my $utr_intersect_file= $all_chr_file . ".intersect_utr_variants";

intersectUTR($intersect_locations, $all_chr_gff_file, $utr_gff, $utr_interset_file,
             $utr_gene_index_href, $utr_intersect_file, $noncoding_annotation_file);

# Intersect Coding
$coding_intersect_file = $all_chr_file . ".intersect_cds_variants";
$all_chr_file_cds_only = $all_chr_file . ".cds_only";
intersectCoding($intersect_locations, $all_chr_gff_file,
                $coding_gff, $coding_intersect_file,
                $all_chr_file_cds_only, $all_chr_file);

##### End Pauline's additions on Sept 29, 2011

# SIMNL<2011-11-11>: Here we create a file containing user queries that are NOT coding
#grep -vhFxf b690bce271_nssnv_segmentaa.allchrfile.cds_only b690bce271_nssnv_segmentaa.allchrfile > b690bce271_nssnv_segmentaa.non_coding_list
my $non_coding_list = "$all_chr_file.noncoding.list";
`grep -vhFxf $all_chr_file_cds_only $all_chr_file > $non_coding_list`;



&finish_script($content, 0); 

# Close the I/O handles
#close(STDin);
close(STDOUT);
#close(STDERR); 



open (CHR_FILE,"$all_chr_file_cds_only") || die ("Cannot open all chr file cds only\n");
#open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
$count = 0;

# CHANGE TO MD5
#$new_pid = $pid+1;
my $tmp_pid = $$ + 1;
my $new_date = `date`; chomp ($new_date);
my $new_key = substr (md5_hex($tmp_pid . $new_date), 0,10);
$new_pid = $new_key . "_" . $toolname;
if ($FILENAME ne "") {
    $new_pid .= "_" . $FILENAME;
}
$tmp_pid++;
$input_set = 1;
open (OUTFILE,">$tmp/$new_pid.chrfile");

while (<CHR_FILE>){
	chomp;
	if ($_ !~ /\d/){next}
	$count++;
	if ($count % $num_coords_per_split== 0){
		push @pid_list,$new_pid;
		print OUTFILE "$_\n";
                $job_size = $num_coords_per_split;
		$start = $count-$job_size+1;
		$end = $count ;
                $outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Partitioned set $input_set</p></td>";
                $outpage_content .= "<td><p>Input rows $start to $end</p></td>";
                $outpage_content .= "<td><p>$new_pid</p></td>";
                $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
                $outpage_content .= "</tr>\n";
		close(OUTFILE);
		$input_set++;

		$tmp_pid++;
		my $new_date = `date`; chomp ($new_date);
		my $new_key = substr (md5_hex($tmp_pid . $new_date), 0,10);
		$new_pid = $new_key . "_" . $toolname;
		if ($FILENAME ne "") {
		    $new_pid .= "_" . $FILENAME
		}
		open (OUTFILE,">$tmp/$new_pid.chrfile");
	}
	else{
		print OUTFILE "$_\n";
		
	}
}
close(CHR_FILE);
close(OUTFILE);
$job_size = $count % $num_coords_per_split;
$start = $end+1;
$end = $start + $job_size-1;
if ($job_size == 0){
	system("rm -f $tmp/$new_pid.chrfile");
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td>Complete set</td>";
        $outpage_content .= "<td><p>Input rows 1 to $end</p></td>";
        $outpage_content .= "<td><p>$pid</p></td>";
        $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
        $outpage_content .= "</tr></table>\n";

}
else{
	push @pid_list, $new_pid;
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Partitioned set $input_set</p></td>";
	$outpage_content .= "<td><p>Input rows $start to $end</p></td>";
	$outpage_content .= "<td><p>$new_pid</p></td>";
	$outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
	$outpage_content .= "</tr>\n";
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td>Complete set</td>";
        $outpage_content .= "<td><p>Input rows 1 to $end</p></td>";
        $outpage_content .= "<td><p>$pid</p></td>";
        $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
        $outpage_content .= "</tr></table>\n";

}

$outpage_content .= "<p class=\"header1\">Batch Report</p>\n";
$outpage_content .= "<p>Number of input (non-intronic) variants: <br />\n";
$outpage_content .= "Coding variants: <br />\n";
$outpage_content .= "Coding variants predicted: <br />\n";
$outpage_content .= "Tolerated: <br />\n";
$outpage_content .= "Damaging: <br />\n";
$outpage_content .= "Nonsynonymous: <br />\n";
$outpage_content .= "Synonymous: <br />\n";
$outpage_content .= "Novel: <br />\n";



### Changes due to change to md5
my $list = join(":", @pid_list);


#Pauline Sept 7 2010 commented out beacuse didn't have template file
#  my $outpage_title = "SIFT: Feed to CHR Coords";
#  my @stylesheets = qw(/stylesheets/main.css);
 
#  my $outpage_tt = Template->new({
#    ABSOLUTE => 1,
#  });

#  my $template_file = ''; #'/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/templates/3_column_fixed_width.tpl';
#  my $vars = {
#    main_content => $outpage_content,
#    page_header => $outpage_header,
#    title => $outpage_title,
#    stylesheets => \@stylesheets,
#  };

#  my $template_output = '';
  
#  $outpage_tt->process($template_file, $vars, \$template_output) || die $oupage_tt->error();
#  print OUTPAGE $template_output;
  print OUTPAGE $outpage_content; # Pauline Sept 7 2010 added in ;

close (OUTPAGE);
#system("cp $tmp/$pid.outpage.html $tmp/$pid.outpage.swap.html");
#direct to the outpage url
#print "Location: $outpage_url\n\n";
#print_outpage();






open (BATCH_FILE,">$tmp/$pid.batchfile");
for ($i = 0; $i < scalar @pid_list; $i++){
	$new_pid = @pid_list[$i];
	chomp $new_pid;
	if ($i == scalar @pid_list -1) {
	    # SIMNL<2011-12-07>: Cater to other organisms
	    print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$database_loc\t$seq_identity_filter\t$COORD_SYSTEM\tLAST\t$output_options\t$address\n";
#	    print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tLAST\t$output_options\t$address\n";
	}
	else{
	    # SIMNL<2011-12-07>: Cater to other organisms
	    print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$database_loc\t$seq_identity_filter\t$COORD_SYSTEM\tNOT_LAST\t$output_options\t$address\n";
#	    print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tNOT_LAST\t$output_options\t$address\n";
	}
}

# Changes due to md5 support
# SIMNL<2011-12-07>: Cater to Hoi's CLI
system("perl $bin/Extended_SIFT_feed_to_chr_coords_batch.pl $tmp/$pid.batchfile $list $key $show_multi_transcripts $isCLI $message_handle $s3filename");
#system("perl $bin/Extended_SIFT_feed_to_chr_coords_batch.pl $tmp/$pid.batchfile $list $key $show_multi_transcripts");

print DEBUG "perl $bin/Extended_SIFT_feed_to_chr_coords_batch.pl $tmp/$pid.batchfile $list $key $show_multi_transcripts $isCLI $message_handle $s3filename\n";

close (BATCH_FILE);

sub print_outpage{
	open(OUTPAGE, "$tmp/$pid.outpage.html") || die("cannot open outpage");
	while (<OUTPAGE>){
        	print;
	}
	close (OUTPAGE);
}








sub general_parse {
        local ( $content_type, $query_string ) = @_;
        local ( %ans, @q, $pair, $loc, $boundary, $temp, $temp1 );

	# SIM TEST
	#print "$query_string<BR>";
	# END SIM TEST

        if ( $content_type eq "application/x-www-form-urlencoded" ) {

                # break up into individual name/value lines
                @q = split( /&/, $query_string );

		### NOTE: To add support for adding filename to pid here

                foreach $pair (@q) {

                        # break the name/value pairs up
                        # use split rather than regular expressions because the value
                        # may have
                        #  newlines in it
                        split( /=/, $pair, 2 );

                        # change '+' to ' '
                        $_[1] =~ s/\+/ /g;

                        # change the escaped characters (must be after the split on '&'
                        # and '=')
                        $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

                        $ans{ $_[0] } = $_[1];
                }    #end of foreach $pair

        }    #end of if ($content_type)
        else {
                $loc = index( $content_type, "boundary=" );
                if ( $loc > 0 ) {
                        $temp = substr( $content_type, $loc + 9 );

                 #               Why is this necessary? (boundary= doesn't match actual)
                        $boundary = "--" . $temp;

                        # break up into individual name/value lines
                        @q = split( /$boundary/, $query_string );

                        foreach $pair (@q) {

                                # break the name/value pairs up
                                $loc = index( $pair, "name=" );
                                $temp = substr( $pair, $loc + 5 );

                                #         $loc = index($temp, "\n\n");
                                $loc = index( $temp, "\n" );
                                $temp1 = substr( $temp, $loc + 2 );


				# Store file names
				if ($temp =~ /CHR_file.*filename=\"(.*)\"/) {
				    my $filepath = $1;
				    my $separator = '\/';
				    if ($filepath =~ /\\/) {
					$separator = '\\';
				    }
				    my $last_index = rindex($filepath, $separator);
				    $FILENAME = substr($filepath, $last_index + 1);
				}


                                #   Get rid of stuff after the name; including semicolon if any
                                $loc_semi = index( $temp, ";" );
                                $loc_eol  = index( $temp, "\n" );
                                $loc      = $loc_eol;
                                if ( $loc_semi > 0 && $loc_semi < $loc ) {
                                        $loc = $loc_semi;
                                }
                                if ( $loc > 0 ) { $temp = substr( $temp, 0, $loc ); }

                                #               Get rid of quotes around the name
                                $temp =~ s/\"//g;

                                #               Still has a trailing whitespace character ...
                                $temp =~ s/\s//g;

                  #               Substitute spaces with nothing
                  #               Need to strip leading/ending whitespace off of $temp1,
                  #               but be careful not to strip off internal CRs
                  #               MAC file lines end in just \r, no \n, so makelis won't
                  # find all
                  #               of the sequences; DOS file lines end in \r\n, UNIX in
                  #\n.
                  #               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
                 #In my case, I want to keep the newlines in fields which have "file" or
                 # 'seq"
                 # and remove newlines everywhere else.
                 #if ( $temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/ ) {
                                $temp1 =~ s/\r\n/\n/g;
                                $temp1 =~ s/\r/\n/g;

                                #}

                         # for other variables that are NOT files or file-like, remove extra
                         #whitespace
                         #else { $temp1 =~ s/\s//g; }
                                if ( $temp ne "" ) { $ans{$temp} = $temp1; }
                        }    # end of foreach
                }    #end of if loc > 0
                else {
                        print "Cannot parse\n";
                        print "content_type=$content_type\n";
                        print "query_string=$query_string\n";
                }
        }
        return %ans;

        #print "</PRE>";
}    # end of general_parse
sub hextodec {
        unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

  ###########################################################################
  # This function calls the template object and injects the HTML variables
  # into the template. 
  #
  # @param : the HTML content
  # @param : the exit code
  # @return void
  ###########################################################################
  sub finish_script {
    my($content, $code) = @_;

    print $content;

    my $title = "SIFT: Genome Center";
    my $header = "SIFT Genome Center";
    my @stylesheets = qw(/stylesheets/main.css);

#    my $tt = Template->new({
#      ABSOLUTE => 1,
#    });
# Pauline Sept 7 2010 commented out 
  #  my $template_file = '/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/templates/3_column_fixed_width.tpl';
#    my $vars = {
#      main_content => $content,
#      title => $title,
#      page_header => $header,
#      stylesheets => \@stylesheets,
#    };

  #  $tt->process($template_file, $vars) || die $tt->error();

    if($code < 0) {  
      exit($code);
    } # end if($code < 0)
  } # end finish_script
  ###########################################################################

