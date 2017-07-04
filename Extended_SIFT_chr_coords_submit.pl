#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use Class::Struct;
#
# Nov 28, 2000 added option to email results
# 7/4/01 gap option is turned off permanently
#11/11/03 Added SIFT_queue stuff  JGH
# 3/8/04 Don't wait for search to finish, send them format.pl URL  JGH
# 9/3/10 Commented out sequence extractor and datadir
#-------------------------------------------------------------------------
#use lib '/usr/local/packages/sift/site_perl/5.8.8';
use DBI;
#use Tie::IxHash;
use Template;

#########################################################################
#
# SNL <2011-08-01>:
# From SIFT_feed_to_chr_coords_batch.pl
# If users checked multiple-rows per transcripts checkbox
# This script will be invoked, otherwise, the original SIFT_chr_coords_submit.pl
# will be called.
#
#########################################################################

$| = 1;
require 'SIFT_subroutines.pm';
require 'config.pl';

#       Set file permissions to rw-rw----
system("umask 006");
#my $tmp             = "/opt/www/sift/tmp";
my $url_address = $base;
my $master_pid = $ARGV[0];
my $pid = $ARGV[1];
my $all_chr_file = $ARGV[2];
my $organism = $ARGV[3];

my $last_partition = 0;
if ($ARGV[5] !~ /NOT_LAST/i){
	$last_partition = 1;
}

my $debugfile = "$tmp/Extended_SIFT_chr_coords_submit.debug";
open(DEBUG, ">$debugfile");

# SIMNL<2011-10-27>: Need to start working on extending this to other than Human
#my $Human_db_dir    = "/usr/local/web/packages/db/$organism"; #"/opt/www/sift/db/$organism";
#my $database_directory = "/usr/local/web/packages/db/$organism"; #"/opt/www/sift/db/$organism";
my $database_directory = "/mnt1/db/$organism"; 


print DEBUG "database_directory: $database_directory\n";


my $output_options = $ARGV[7];


## Changes due to 
my $pids_list_as_str = $ARGV[8];
my @list_of_pids = split(":", $pids_list_as_str);
my $pid_key = $ARGV[9];


# SNL <2011-08-03> Add support for multi-row per transcripts
my $multi_transcripts = $ARGV[10];

#my $address = $ARGV[11];

# SIMNL<2011-12-06>: To cater to Hoi's CLI application
my $isCLI = $ARGV[11];
my $message_handler = $ARGV[12];
my $s3filename = $ARGV[13];
my $address = $ARGV[14];


my $seq_identity_filter = $ARGV[4];
my @stylesheets = qw(/stylesheets/main.css);

check_ip_counts();

my $COORD_SYSTEM = $ARGV[6];

struct VariantInfo => {

        rsid => '$',
        orn_db => '$',
        enst => '$',
        ensp => '$',
        region => '$',
        snp_type => '$',
        codon1 => '$',
        codon2 => '$',
        nt1 => '$',
        nt2 => '$',
        alleles_db => '$',
        AA1 => '$',
        AA2 => '$',
        AAPOS => '$',
        score => '$',
        median => '$',
        seqs_rep => '$',
        freq_av => '$',
        freq_ceu => '$',

	# Sim
	freq_hcb => '$',
	freq_jpt => '$',
	freq_yri => '$',
	# Sim end

	# SIMNL<2011-10-31> (empty if no allele frequencies)
	avg_1k_af => '$',
	european_1k_af => '$',
	east_asian_1k_af => '$',
	west_african_1k_af => '$',
	south_asian_1k_af => '$',
	american_1k_af => '$',


	# SIM 1000 genome
	all_1K_genome => '$',
	european_1K_genome => '$',
	east_asian_1K_genome => '$',
	west_african_1K_genome => '$',
	south_asian_1K_genome => '$',
	american_1K_genome => '$',
	# SIM end


# gene info
        ensg => '$',
        enst => '$',
        ensp => '$',
        gene_name => '$',
        gene_desc => '$',
        ensfm => '$',
        fam_desc => '$',
        gene_status => '$',
        fam_size => '$',
        kaks_mouse => '$',
        kaks_macaque => '$',
        mim_status => '$',


};

update_status($master_pid,$pid,"Running");

# SIMNL<2011-10-27>: bins.list is a soft link to bins_xxx.list
#$bin_file = "$Human_db_dir/bins.list";
$bin_file = "$database_directory/bins.list";


#Create binned SNP files for querying variation chr database- the script also creates snp_chr_map_file to map bin file to chr database and table
system("perl $bin/Extended_map_choords_to_bin.pl $bin_file $all_chr_file $pid");


print DEBUG "perl $bin/Extended_map_choords_to_bin.pl $bin_file $all_chr_file $pid\n";

close(DEBUG);

my @query_result = ();

# SIMNL <2011-10-31>: Now we need to find a different way to get 1000 genomes data
# it is now in a different database, with distinct key:
# CHR, COORD1, COORD2, UNIQ_ID, NT1, NT2
my %one_thousand_genome_hash = ();

# SIMNL <2011-10-27>: This array stores the uniq key id for each query
my @uniq_key_id_array = ();


# SIMNL <2011-11-09>: Optimize - query 1000 genomes only when user wants that data
my $needs_1000_genome = 0; # if > 0, then need to run query
my @allOptions = split(/,/, $output_options);
for(my $idx = 15; $idx < scalar(@allOptions); $idx++) {
    $needs_1000_genome += $allOptions[$idx];
}

open (SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt");
while (<SNP_CHR_MAP_FILE>){
        chomp;

        @elts = split /\t/, $_;
        $chr = @elts [0];
        $bin = @elts[1];
        $snp_chr_bin_file = @elts[2];
        $db_chr_str = @elts[3];
        $table_chr = @elts[4];
	
	# SIMNL<2011-10-27>: Change the way we view databases, no longer just Humans
	$db_chr =DBI->connect( "dbi:SQLite:dbname=$database_directory/Variation_CHR$chr.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );

#        $db_chr =DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/Human_CHR$chr.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );

	$db_chr->do('PRAGMA synchronous=1');
	$db_chr->do('PRAGMA cache_size=4000');
	#$sth_db_chr = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD1 = ? AND COORD2 = ? ");


        open (BIN_FILE, "$tmp/$snp_chr_bin_file") || die ("Cannot open snp_chr_bin_file");
        $coord2_str = "";
	$coord1_str = "";
	while (<BIN_FILE>){
		$count ++;
                chomp;
                @elts2 = split '\t',$_;
                $coord1 = @elts2[1];
                $coord2 = @elts2[2];
                $orn_usr = @elts2[3];
                $alleles = @elts2[4];
		$usr_comment = @elts2[5];
                $alleles =~ /(.+?)\/(.+?)/;
                $nt2a = $1;
                $nt2b = $2;

	# PCN if given in - strand, then reverse complement
                if ($orn_usr == -1) {
                        $nt2a =~ tr/ACGTacgt/TGCAtgca/;
                        $nt2b =~ tr/ACGTacgt/TGCAtgca/;
                        $orn_usr = 1;
                }

		$coord2_str.="$coord2,";
		$coord1_str.="$coord1,";
		### comment this block if using in operator see below.
		#$sth_db_chr->execute($coord1,$coord2);
		#while (@rows = $sth_db_chr->fetchrow_array()){
        	#	push @query_result, join("\t",@rows);
		
	        #}
		#system("echo $snp_chr_bin_file $count >> $tmp/$pid.query_status.txt");		
		###
		
	} #end while (<BIN_FILE>)
	close (BIN_FILE);
	chop $coord1_str;
	chop $coord2_str;

	### uncomment following block if using in operator - see above
#	$sth_db_chr = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD2 in ($coord2_str) AND COORD1 in ($coord1_str)");
	$sth_db_chr = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD2 in ($coord2_str)");
	
	#system ("echo \"select * from $table_chr where chr = \'chr$chr\' AND COORD2 in \($coord2_str\) AND COORD1 in \($coord1_str\)\" > $tmp/remove");
	$sth_db_chr->execute();

        while (my @rows = $sth_db_chr->fetchrow_array()){
	    push @query_result, join("\t",@rows);	    

	    # SIMNL <2011-10-24> The previous "ENST" has now been replaced with the uniq key id.  
	    my $uniq_id = $rows[6];
	    push @uniq_key_id_array, $uniq_id;

	} #end while

	system("echo $snp_chr_bin_file $count $time>> $tmp/$pid.query_status.txt"); 
	###

	############################################################################
	# SIMNL<2011-10-31>: Get 1000 genome data
	# At this point, let's find the corresponding 1000 Genomes

	if ($needs_1000_genome > 0) {

	    my $genome_1K_db = "$database_directory/1000Genomes/chr$chr" . "_1kg.sqlite";
	    
	    my $db_1000 =DBI->connect( "dbi:SQLite:dbname=$genome_1K_db","", "", { RaiseError => 1, AutoCommit => 1 } );
	    my $g1k_query = "select * from $table_chr where chr = \'chr$chr\' AND COORD2 in ($coord2_str) AND COORD1 in ($coord1_str)";
	    my $sth_db_1000 = $db_1000->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD2 in ($coord2_str)");
	    $sth_db_1000->execute();
	    
	    while (my @g1K_rows = $sth_db_1000->fetchrow_array()) {
		my ($g1k_chr, $g1k_coord1, $g1k_coord2, $g1k_rsid, $g1k_uniq_id, $g1k_nt1, $g1k_nt2,
		    $avg, $european, $east_asian, $west_african, $south_asian, $american,
		    $avg_ref, $avg_alt, $european_ref, $european_alt, $east_asian_ref, $east_asian_alt,
		    $west_african_ref, $west_african_alt, $south_asian_ref, $south_asian_alt, $american_ref, $american_alt) = @g1K_rows;	    
		
		my $g1k_result = join("\t", @g1K_rows);
		
		$avg_ref = &roundDecimals($avg_ref);
		$european_ref = &roundDecimals($european_ref);
		$east_asian_ref = &roundDecimals($east_asian_ref);
		$west_african_ref = &roundDecimals($west_african_ref);
		$south_asian_ref = &roundDecimals($south_asian_ref);
		$american_ref = &roundDecimals($american_ref);
		
		$avg_alt = &roundDecimals($avg_alt);
		$european_alt = &roundDecimals($european_alt);
		$east_asian_alt = &roundDecimals($east_asian_alt);
		$west_african_alt = &roundDecimals($west_african_alt);
		$south_asian_alt = &roundDecimals($south_asian_alt);
		$american_alt = &roundDecimals($american_alt);
		
		my $avg_af = "";
		if ($avg_ref ne "EMPTY" && $avg_alt ne "EMPTY") { $avg_af = "$avg_ref,$avg_alt"; }
		
		my $european_af = "";
		if ($european_ref ne "EMPTY" && $european_alt ne "EMPTY") { $european_af = "$european_ref,$european_alt"; }
		
		my $east_asian_af = "";
		if ($east_asian_ref ne "EMPTY" && $east_asian_alt ne "EMPTY") { $east_asian_af = "$east_asian_ref,$east_asian_alt"; }
		
		my $west_african_af = "";
		if ($west_african_ref ne "EMPTY" && $west_african_alt ne "EMPTY") { $west_african_af = "$west_african_ref,$west_african_alt"; }
		
		my $south_asian_af = "";
		if ($south_asian_ref ne "EMPTY" && $south_asian_alt ne "EMPTY") { $south_asian_af = "$south_asian_ref,$south_asian_alt"; }
		
		my $american_af = "";
		if ($american_ref ne "EMPTY" && $american_alt ne "EMPTY") { $american_af = "$american_ref,$american_alt"; }
		
		my $genome_1000 = "$avg_af\t$european_af\t$east_asian_af\t$west_african_af\t$south_asian_af\t$american_af";
		
		my $key = "$g1k_chr:$g1k_coord1:$g1k_coord2:$g1k_nt1:$g1k_nt2";
		$one_thousand_genome_hash{$key} = $genome_1000;

	    } #end while
	    $db_1000->disconnect;
	} #end if ($needs_1000_genome > 0)
	############################################################################

} #end while (<SNP_CHR_MAP_FILE>)
close (SNP_CHR_MAP_FILE);

#index all results obtained into index_query_result. hash of arrays.



# SIMNL<2011-10-27>: Human_Supp.sqlite to Variation_Supp.sqlite
#$db_supp = DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp = DBI->connect( "dbi:SQLite:dbname=$database_directory/Supplement_DB.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp->do('PRAGMA synchronous=1');

# SIMNL <2011-10-24> New database that maps the unique key id to the corresponding transcript id
#$db_map = DBI->connect( "dbi:SQLite:dbname=$Variation_db_dir/GeneAnnotationMap.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_map = DBI->connect( "dbi:SQLite:dbname=$database_directory/GeneAnnoDB.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_map->do('PRAGMA synchronous=1');

my $db_map_stm = $db_map->prepare("SELECT * FROM MAPPING_TABLE WHERE UNIQ_ID = ?");
# SIMNL <2011-10-24> The fields ENST/ENSP/ENSG in GENE_INFO table
# should have been renamed, but decided to keep it to 
# The ENST field can return NM_XXX as well
# The prepare part for gene_info is moved down 
#$sth_db_supp_geneinfo = $db_supp->prepare("select * from GENE_INFO where ENST = ?"); replaced by the following 2 queries
my $query_using_enst = $db_supp->prepare("select * from GENE_INFO WHERE ENST = ?");
my $query_using_nm = $db_supp->prepare("select * from GENE_INFO WHERE ENST = ?");


$sth_db_supp_allelefreq = $db_supp->prepare("select * from ALLELE_FREQ where RSID = ?"); 



### NGAKLENG <2011-10-24>: Orn now comes from Variation_CHRZ.sqlite
### I checked and verified that my SIFT_db.sqlite database's orn is the same as the 
### transcript orn found in ensGene.txt.gz, refGene.txt.gz, etc. 
### That is, I didn't do any flipping of the strand from xxxGene.txt.gz when 
### transferring to my database.
### The original TRANSCRIPT_DB.sqlite for the most part matches, but there were
### some that were different. To look into this later.
=pod
## PCN - Jan 25, 2010 add transcript info
#$db_tx = DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/TRANSCRIPT_DB.sqlite",
         "",  "", { RaiseError => 1, AutoCommit => 1} );
$db_tx->do('PRAGMA synchronous=1');
$sth_db_strand_tx = $db_tx->prepare ("select orn from ENST_REGION where ENST = ?");
=cut

%tx_enst_hash;


##### SNL <2011-08-01> to support multiple rows
# This stores all results returned by the database using chr:coord1:coord2 as key
my $multi_href = &prepareMultiTranscriptHash(\@query_result);
my %multi_hash = %{$multi_href}; # <k,v> = <chr:coord1:coord2, @datum>

# SIMNL <2011-10-27>: Now we need to try to get one row that has Ensembl or RefSeq Transcript IDs to the forefront of the array
# This is because we only have Ensembl/RefSeq Gene Information, 
my %non_ensembl_refseq_rows = (); # This stores <$key, @rows_without_enst_refseq>, we push this into     
                                  # $index_query_result{$key} pushing those WITH Ensembl/RefSeq first

 
### SNL <2011-08-01> At this point in time, the %multi_hash contains all returned values from the db
### with chr:coord1:coord2 as key
# SIMNL <2011-10-27>: Now need to get the uniq_key_id
my $num_results = scalar(@query_result);

foreach (my $index = 0; $index < $num_results; $index++) {  
#foreach $row (@query_result){

    my $row = $query_result[$index];

	chomp $row;

    my ($key, $row2, $transcript_identities) = &getRowData($row, $index, "YES");

    # SIMNL<2011-10-27>: We want those with Ensembl and/or RefSeq transcripts
    # to be in front, this is to cater to the default case where the user
    # does not want multiple rows for transcripts cutting across the same
    # chr,coord1,coord2,nt1,nt2.
    # We need ENS and/or NM to get Gene Info if user selects Gene Info options
    if ($transcript_identities =~ /ENS/) {
	push @ {$index_query_result{$key}}, $row2;
    } elsif ($transcript_identities =~ /NM/) {
	push @ {$index_query_result{$key}}, $row2; # Priority given to Ensembl gene information since we had those in the old database.
    } else {	
	push @ {$non_ensembl_refseq_rows{$key}}, $row2;
    }

# Original
#    push @ {$index_query_result{$key}}, $row2;	
# End Original

} #end large for loop



# SIMNL<2011-10-27>: Now, put those without Ensemble/RefSeq back into the data
# This is to push rows with ens transcripts or refseq transcripts to the forefront
my @no_ensembl_refseq_keys = keys(%non_ensembl_refseq_rows);
foreach my $key (@no_ensembl_refseq_keys) {
    my @data = @ {$non_ensembl_refseq_rows{$key}};
    foreach my $datum (@data) {
	push @ {$index_query_result{$key}}, $datum;
    }
}


$rsid = "";$enst = "";$ensp="";$coord1="";$coord2 = "";$chr = "";$nt1 = "";$nt2="";


#parse snp_chr_map_file to fetch data from previously created $index_query_result.
open (SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt");
open (BAD_ENST_FILE, ">$tmp/$pid.bad_enst_file.txt");


# SIM
#my $thousand_genome_hashref = &create_thousand_genome_hash();
#my %thousand_genome_hash = %{$thousand_genome_hashref};
# SIM end


while (<SNP_CHR_MAP_FILE>){
	chomp;
	@elts = split /\t/, $_;
	$chr = @elts [0];	
	$bin = @elts[1];
	$snp_chr_bin_file = @elts[2];
	$db_chr_str = @elts[3];
	$table_chr = @elts[4];

	my $aa_invalid_nt1 = "";
	my $aa_invalid_nt2 = "";

	# SIMNL <2011-09-09> Dr. Azim bug: If after iterating, and not found, need to use the alleles matching the original user's nt
	my $orig_user_nt1 = ""; # to handle the block with 'not a reference allele, don't know how to process'
	my $orig_user_nt2 = ""; # because the nt1,nt2 has already been iterated over.
	my $orig_user_orn = 1;
	open (BIN_FILE, "$tmp/$snp_chr_bin_file") || die ("Cannot open snp_chr_bin_file");
	while (<BIN_FILE>){
		chomp;

		@elts2 = split '\t',$_;
		$coord1 = @elts2[1];
		$coord2 = @elts2[2];
		$orn_usr = @elts2[3];
## PCN if orn_usr is -1, reverse complement alleles, because in database
# it's referenced to +1 orientation
		$alleles = @elts2[4];
		$alleles =~ /(.+?)\/(.+?)/;
		$usr_comment = @elts2[5];

		$nt2a = $1; 
		$nt2b = $2; 

		if ($orn_usr == -1) {
                        $nt2a =~ tr/ACGTacgt/TGCAtgca/;
                        $nt2b =~ tr/ACGTacgt/TGCAtgca/;
                        $orn_usr = 1;
                }


		$orig_user_nt1 = $nt2a; # SIMNL <2011-09-09> Dr. Azim case, store user alleles oriented to positive orientation
		$orig_user_nt2 = $nt2b;

		$key_ensp_CDS = "chr$chr\t$coord1\t$coord2\t1\t1"; 			#CDS is true, ensp is true,
 
		@rows_ensp_CDS = @{ $index_query_result{$key_ensp_CDS} };
		$key_no_ensp_CDS = "chr$chr\t$coord1\t$coord2\t1\t0";                   #CDS is true, ensp is false,
                @rows_no_ensp_CDS = @{ $index_query_result{$key_no_ensp_CDS} };
#Prateek had below commented out, only looking at rows_ensp_CDS
# Pauline allowed rows_no_ensp_CDS to be seen
		@rows_CDS = (@rows_ensp_CDS,@rows_no_ensp_CDS);
		#replaced above with below - not using no-ensp rows any more
		#will also not output anything for not found in db
		#@rows_CDS = @rows_ensp_CDS;
# Pauline commented out above, print out rows with non-ensp id's
		$key_ensp_no_CDS = "chr$chr\t$coord1\t$coord2\t0\t1";			#CDS is false , ensp is true
                @rows_ensp_no_CDS = @{ $index_query_result{$key_ensp_no_CDS} };
		$key_no_ensp_no_CDS = "chr$chr\t$coord1\t$coord2\t0\t0";    		#CDS is false , ensp is false
                @rows_no_ensp_no_CDS = @{ $index_query_result{$key_no_ensp_no_CDS} };
		@rows_no_CDS = (@rows_ensp_no_CDS,@rows_no_ensp_no_CDS); 		# final combined - this and $rows_CDS will now be used.

			
		if (scalar @rows_CDS == 0){						# no CDS found for this coords


			if (@rows_no_CDS == 0){ 					#coord not in db
		#		next;
				#added above for not outputting anything not in db
				my $m = make_new_variant ();
                                $m->region ("NON-GENIC");
				# SIMNL <2011-09-09> Found this bug while fixing Dr. Azim bug 
				# by this time, $orn_user is set to 1, but $alleles is still user's allele
				# In the results page, we show the users query with alleles oriented to positive wrt ref genome.
				my $reoriented_alleles = $nt2a . "\/" . $nt2b; 
                                push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$reoriented_alleles\t$usr_comment";
                                $coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$reoriented_alleles\t$usr_comment"} = $m;
                                #original: push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment";
                                #original: $coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment"} = $m; # ORIGINAL 

			}
			else{			# rows outside of CDS found (promotor, Downstream or UTR)
				my $v = make_variant_with_elts (@rows_no_CDS[0]);

# PCN 012610 changed to print orientation of user, which should be +1 strand
# because alleles will be printed out with respect to +1 strand, so it makes more sense
                                push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment";
				$coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment"} = $v; # original

			} #end else from if (@rows_no_CDS == 0)

		} # end if (scalar @rows_CDS == 0) 
		else {			#CDS found - iterate through rows- choose one with right nt2 - if nt2 is ref then modify subst

			$selected_row = -1;
			$ref_allele_row = -1;
			$AA1_VALID_found = 0;
			for ($i = 0; $i < scalar @rows_CDS; $i++){

				$row = @rows_CDS[$i];
				chomp $row;

				#print "$row<BR>";
				#system("echo $row >> $tmp/remove");
				@elts_CDS = split /\t/, $row;
				$enst = @elts_CDS[6];
                                $ensg = @elts_CDS[26];
				$orn_db = @elts_CDS[3];
				$nt1 = @elts_CDS[10];
				$nt2 = @elts_CDS[11];
				$AA1_VALID = @elts_CDS[21];
				if ($AA1_VALID == 1){
					$AA1_VALID_found = 1;
					# SIMNL <2011-09-09>: Dr Azim bug: if the if block below is never accessed, then
					# ref_allele_row ends up being the last row of the results from database, in most cases
					# the alleles from that last row will not match the user's query alleles.
					# original: $ref_allele_row = $i;  #save last row so if seclected row==-1, this row is used assuming  ref allele analysis reqd.
					my $tmp_nt1 = $nt1; my $tmp_nt2 = $nt2; my $tmp_orn = $orn_db; # do not want to mess up original variables
					if ($tmp_orn == -1) {
					    $tmp_nt1 =~ tr/ACGTacgt/TGCAtgca/;
					    $tmp_nt2 =~ tr/ACGTacgt/TGCAtgca/;
					}	   			
					if ( (($tmp_nt2 eq $nt2b && $tmp_nt1 eq $nt2a)
							 ||
					     ($tmp_nt2 eq $nt2a && $tmp_nt1 eq $nt2b)) 
						&& $tmp_orn eq $orn_db) {
#						&& $tmp_orn eq 
#						 $tx_enst_hash{$enst}) {
					   $ref_allele_row = $i; 
                                	}





# Pauline commented out Prateek's code 01-22-10
# homozygous variants not processed, heterozygous only
# don't care which orientation, reference doesn't have to be first
# add check for orientation (orn)
#print "nucleotides $nt2 $nt1  and $nt2a $nt2b\n";
#print "transcript orientation is $tx_enst_hash{$enst} $enst $ensg orn db $orn_db\n";
					if ( (($nt2 eq $nt2b && $nt1 eq $nt2a)
							 ||
					     ($nt2 eq $nt2a && $nt1 eq $nt2b)) 
					     && $tmp_orn eq $orn_db ) { # Temporarily set this as such, to fix
#						&& $orn_db eq  
#						 $tx_enst_hash{$enst}) { 
#print "orn_tx $orn_tx orn_db $orn_db\n";
                                               	$selected_row = $i;
                                               	last;
                                	}

				}
				else{				#this CDS has bad enst. go to next 				    
				    # SIMNL <2011-08-16> fixing inconsistent allele1/allele2 issue (Yunlei)
				    # Check nt1, nt2 against database's nt1, nt
				    if ( (($nt2 eq $nt2b && $nt1 eq $nt2a) ||
					  ($nt2 eq $nt2a && $nt1 eq $nt2b)) && 
					     $tmp_orn eq $orn_db ) { # Temporarily set this as such, to fix 
#					 $orn_db eq $tx_enst_hash{$enst}) { 
					$ref_allele_row = $i;
				    }
				    $aa_invalid_nt1 = $nt2a;
				    $aa_invalid_nt2 = $nt2b;
				    next;
				    # end SIMNL <2011-08-16> fixing inconsistent allele1/allele2 issue (Yunlei)

=pod original
					 $ref_allele_row = $i;
					#save this coord in  BAD_ENST_FILE for later
					next;
=cut

				}
			} #end for ($i = 0; $i < scalar @rows_CDS; $i++)


			if($selected_row != -1){			#valid CDS
                                $row = @rows_CDS[$selected_row];        
                                chomp $row;
				my $v = make_variant_with_elts ($row);
				my $alleles_db = $v->alleles_db;

                                push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
                                $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;

                        } #end if($selected_row != -1)
	
			elsif($selected_row == -1 && $AA1_VALID_found ==1){		#this means user input requires pred on ref allele

				if ($nt2a eq $nt2b) {
#this means user input requires pred on ref allele
#print  "chr $chr $coord1 $coord2 in ab98ua\n";
				    $row = @rows_CDS[$ref_allele_row];
				    chomp $row;
				    my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
				    $v->alleles_db ($v->nt1 . "\/" . $v->nt1);
				    $v->snp_type ("Synonymous");
				    $v->codon2 ($v->codon1);
				    $v->AA2 ($v->AA1);
				    $v->score ("N/A");
				    $v->median ("N/A");
				    $v->seqs_rep ("N/A");
				    my $s = $v->nt1 . "\/" . $v->nt1;
				    push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$s\t$usr_comment";
				    $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$s\t$usr_comment"} = $v;
                                } else {
# not a reference allele, don't know how to process
# Pauline added error message on Sept 15 , 2011
				    $row = @rows_CDS[$ref_allele_row];
				    chomp $row;
                                 my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
		my $error_message = "<a href=\"/www/SIFT_help.html#SIFT_INPUT_GENOME_Error\" rel=\"external\">Input Error -Possible Strand Issue</a>";

				    $v->snp_type ($error_message);
				    $v->codon1 ("");
				    $v->codon2 ("");
				    $v->AA1 ("NA");
				    $v->AA2 ("");
				    $v->score ("N/A");
				    $v->median ("N/A");
				    $v->seqs_rep ("N/A");

				    # SIMNL <2011-09-09>: Fix Dr. Azim issue, need to display the original user's query
				    # original: $v->alleles_db ($nt1 . "\/" . $nt2);
				    $v->alleles_db ($orig_user_nt1 . "\/" . $orig_user_nt2);
				    $v->orn_db($orig_user_orn);
				    my $alleles_db = $v->alleles_db;
				    push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
				    $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;
				    
                                        print BAD_ENST_FILE "$chr\t$coord1\t$coord2\t$orn_usr\t$nt1\/$nt2\t$usr_comment\n";
                                next;

                                }

			} #end elsif($selected_row == -1 && $AA1_VALID_found ==1)

			elsif($AA1_VALID_found == 0){		#All ensts found bad. save this coord in a different file.

				 $row = @rows_CDS[$ref_allele_row];
                                        chomp $row;
                                 my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
				 my $error_message = "<a href=\"/www/SIFT_help.html#SIFT_OUTPUT_GENOME_Error\" rel=\"external\">Gene Annotation Error</a>";
				 $v->snp_type ("Unknown");
				 $v->codon1 ("");
				 $v->codon2 ("");
				 $v->AA1 ("NA");
				 $v->AA2 ("");
				 $v->ensp ($error_message);
				 $v->score ("N/A");
				 $v->median ("N/A");
				 $v->seqs_rep ("N/A");

# SIMNL <2011-08-16> fix to inconsistent allele1/allele2 (Yunlei)
#                                 $v->alleles_db ($nt1 . "\/" . $nt2);
				 $v->alleles_db ($aa_invalid_nt1 . "\/" .$aa_invalid_nt2);

				 my $alleles_db = $v->alleles_db;
				 push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
				 $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;

				print BAD_ENST_FILE "$chr\t$coord1\t$coord2\t$orn_usr\t$nt1\/$nt2\t$usr_comment\n";
				next;
			} #end elsif($AA1_VALID_found == 0)
		}
	}
	close (BIN_FILE)
	} #end while <SNP_CHR_MAP_FILE>
close (SNP_CHR_MAP_FILE);
close (BAD_ENST_FILE);


# SIMNL<2011-11-01>
# We want to ensure that the original row (the one that would have been printed out if multi-transcript not selected)
# is not re-printed out
my %duplicate_mgt_set = ();


#create enst_pred_list for coords in CDS .
open( ENSTFILE, ">>$tmp/$pid.enstfile" );
foreach $coord (@coord_arr_CDS) {

	my $variant_result    = $coord_result_CDS_hash{$coord};
        my $line = return_variant_line ($variant_result);
#print "line in here ENST is $line for coord $coord\n";
        my $alleles = (split /\t/, $coord)[4];
#print "coord after is $coord alleles is $alleles\n";
        my $key = $variant_result->enst . "\t" . $variant_result->AAPOS;

# Pauline commented out, so another prediction printed out even ifvariant
# is at same location

#	if (exists ($hash_enst_aapos_seen{$key})){
#		next;	
#	}	
#	else{
		$hash_enst_aapos_seen{$key} = "$alleles";
		$enst_line          = "$coord\t$line";
		push( @enst_pred_list, $enst_line );
		print ENSTFILE "$enst_line\n";	
#	}

	$duplicate_mgt_set{$enst_line} = 0;


	# SNL <2011-08-03>: Support for multi-row per transcripts
	if ($multi_transcripts != 0) {
	    my ($chr, $coord1, $coord2, $orn, $alleles, @others) = split(/\t/, $coord); # 1	39352270	39352271	1	G/T	
	    my $key = "chr$chr:$coord1:$coord2";
	    my $aref = $multi_hash{$key}; # If users did not select multi row per transcript, the multi_hash would have been empty anyways.
	    if (defined($aref)) {
		my @array = @{$aref};
		foreach my $data (@array) { # $data is simply the output from the database query (every field)
		    # check if nt1,nt2 oriented correctly is equal to $coord's variant
		    # if so, add to output
		    my ($chromo, $coordinate1, $coordinate2, $orientation, $rsid, $ensg, $enst, $ensp, $region,
			$snp, $nucleo1, $nucleo2, @others) = split(/\t/, $data);
#		# Note: nucleo1 and nucleo2 are from database, which is always + wrt to human genome
#		# $orientation is transcript orientation (orn from database), so we shouldn't do this
		    
		    if ($orientation == -1) {
			# translate to positive wrt reference genome
			$nucleo1 =~ tr/ACGTacgt/TGCAtgca/;
			$nucleo2 =~ tr/ACGTacgt/TGCAtgca/;
			$orientation = 1;
		    } #end if ($orientation == -1)
		    
		    my $ref_nt1 = $variant_result->nt1;
		    my $ref_nt2 = $variant_result->nt2;
		    if ($orn == -1) {
			$nt1 =~ tr/ACGTacgt/TGCAtgca/;
			$nt2 =~ tr/ACGTacgt/TGCAtgca/;
			$orn = 1;		    
		    } # if ($orn == -1)
		    
		    # Now compare
		    if (((uc($nucleo1) eq uc($nt1)) && (uc($nucleo2) eq uc($nt2))) || 
			((uc($nucleo1) eq uc($nt2)) && (uc($nucleo2) eq uc($nt1)))) {
			# Ensure that the original ENST is not re-added in
			if ($variant_result->enst ne $enst) {
			    # Create variant_result from $data
			    # Create enst_line using $coord, $variant_line return
			    my $v = make_variant_with_elts ($data);
			    my $alleles_db = $v->alleles_db;			
			    # Get 1000 genomes presence data
			    my ($chrom, $rs) = &get1000KeyElements($data);

			    my $new_line = return_variant_line($v);
			    my $new_enst_line = "$coord\t$new_line";
			    my $orig_codon1 = $variant_result->codon1;
			    my $orig_codon2 = $variant_result->codon2;
			    my $multi_codon1 = $v->codon1;
			    my $multi_codon2 = $v->codon2;
			    if (($orig_codon1 eq $multi_codon1) && ($orig_codon2 eq $multi_codon2)) { # codons have to be the same.
				push( @enst_pred_list, $new_enst_line ); # This is what gives us multiple rows per transcript



				# SIMNL<2011-11-01> To prevent duplication of rows.
				if (!defined($duplicate_mgt_set{$new_enst_line})) {
				    print ENSTFILE "$new_enst_line\n";

				}
# original				print ENSTFILE "$new_enst_line\n";
			    }
			}
		    }		    
		} #end foreach
	    } #end if (defined $aref)
	} #end if ($multi_transcripts != 0)
	# SNL <2011-08-03>: end change to support multi-row per line

} #end foreach $coord (@coord_arr_CDS)
close(ENSTFILE);


=pod ORIGINAL 2011-08-01
#create enst_pred_list for coords in CDS .
open( ENSTFILE, ">>$tmp/$pid.enstfile" );
foreach $coord (@coord_arr_CDS) {
	my $variant_result    = $coord_result_CDS_hash{$coord};
        my $line = return_variant_line ($variant_result);
#print "line in here ENST is $line for coord $coord\n";
        my $alleles = (split /\t/, $coord)[4];
#print "coord after is $coord alleles is $alleles\n";
        my $key = $variant_result->enst . "\t" . $variant_result->AAPOS;

# Pauline commented out, so another prediction printed out even ifvariant
# is at same location

#	if (exists ($hash_enst_aapos_seen{$key})){
#		next;	
#	}	
#	else{
		$hash_enst_aapos_seen{$key} = "$alleles";
		$enst_line          = "$coord\t$line";
		push( @enst_pred_list, $enst_line );
		print ENSTFILE "$enst_line\n";	
#	}
}
close(ENSTFILE);
=cut

#create enst_subst_list for coords not in CDS - these will not be predicted by SIFT.
open( ENSTFILE, ">>$tmp/$pid.enstfile" );
foreach $coord (@coord_arr_no_CDS) {
        $variant_result    = $coord_result_no_CDS_hash{$coord};
	my $line = return_variant_line ($variant_result);


	my ($none1, $na1, $non2, $na2, $non_genic, $dash, $na3) = split("\t", $line);
	$line = join("\t","", "", "", "", "", "", "",$na2,$non_genic, $dash, $na3);
#	$line = join("\t", $chr,$coord1,$coord2,$orn,$allele,$codons,$na," "," "," ",$non_genic,@others);
#	$line = join("\t", $chr,$coord1,$coord2,$orn,$allele," ", " ", @others);

        $enst_line          = "$coord\t$line";
        push( @enst_pred_list, $enst_line );

        print ENSTFILE "$enst_line\n";

}
close(ENSTFILE);






###################################################################


#read from combined file and print table
open( ENST_PRED_FILE, "$tmp/$pid.enstfile" )
  || die("Cannot open predictions file");

$heading_tsv = get_output_heading($output_options);
$heading_html = $heading_tsv;

my $outfile_table = '';
my $outfile_table_counter = 0;

open( OUTFILETABLE, ">$tmp/$pid\_predictions.html" );
open( OUTFILETSV,   ">$tmp/$pid\_predictions.tsv" );
print OUTFILETSV "$heading_tsv\n";

$outfile_table .= '<HTML><HEAD><link rel="stylesheet" type="text/css" href="/www/css/tables.css" /></HEAD>';


$outfile_table .= '<table border="1">';
$outfile_table .='<tr class="tableHeader"><td><p>';
$heading_html =~ s?\t?</p></td><td><p>?g;
$outfile_table .="$heading_html</p></td></tr>\n";
my $warningflag = 0;

while (<ENST_PRED_FILE>) {
	chomp;

	# SIM
	my @test = split("\t", $_);
	my $chr_1K = "chr" . @test[0];
	my $coord1_1K = @test[1];
	my $coord2_1K = @test[2];	
	my $nt1nt2 = @test[4];
	my ($nt1_1K, $nt2_1K) = split("\/", $nt1nt2);
	my $rsid_1K = @test[9];

#	$table_row = get_output_row($output_options, $_);

    # SIMNL<2011-10-31>: Now that 1K Genome is in $_, we get it from there
#	    $table_row = get_output_row($output_options, $_, $genome1KData);
	    $table_row = get_output_row($output_options, $_);

	$outfile_table .= "<tr class=\"tableRow".($outfile_table_counter++ % 2 == 0 ? "Even" : "Odd")."\">\n";
	if ($table_row =~ /Synonymous/ || $table_row =~ /\w\d+\*/){
		$is_synonymous = 1;
	}
	else{
		$is_synonymous = 0;
	}
	my @fields = split( '\t', $table_row );
	$num_cells = scalar @fields;
	$count = 0;
	for $cell (@fields) {
		chomp $cell;
		$count ++;
		$cell_tsv = $cell;
		$cell_html = $cell;

		if (defined $cell) {
		    $cell =~ s/^\s+//;
		    $cell =~ s/\s+$//;
		    if ($cell eq "") {
			$cell_html = "&nbsp;";
		    }
		}

		if (!defined $cell_html) {
		    $cell_html = "&nbsp;";
		}

		if ( $cell =~ /DAMAGING/ || $cell =~ /Warning/ ) {
			$cell_html =  "<span class=\"red\">$cell</span>";
			$cell_tsv = $cell;
		}
		elsif ( $cell =~ /NOT PREDICTED/i || $cell =~ /not found/i ) {
			$cell_html = "<span class=\"blue\">$cell</span>";
			$cell_tsv = $cell;
		}
		elsif ($cell =~ /(rs\d+):(.+?)/i){
			$dbsnp_id = $1;$mutation = $2;
			$cell_html =  "<a href=\"http:\/\/www.ncbi.nlm.nih.gov\/sites\/entrez?db=snp&cmd=search&term=\+$dbsnp_id\" rel=\"external\">$dbsnp_id:$mutation</a>";
			$cell_tsv = $cell;
		}


		# SIMNL<2011-10-31>: We added more columns
		if($count >= 13 && $count <= 15 && $is_synonymous == 1){
			$cell_tsv = "N\/A";			
		}
=pod ORIGINAL
		if($count >= 9 && $count <=12 && $is_synonymous == 1){
			$cell_html =  "N\/A";
			$cell_tsv = "N\/A";
		}
=cut

		$outfile_table .= "<td><p>$cell_html</p></td>";
		if ($count == $num_cells){
			print OUTFILETSV "$cell_tsv";
		}
		else{
			print OUTFILETSV "$cell_tsv\t";
		}

	}
	$outfile_table .= "</tr>\n";
	print OUTFILETSV "\n";

}

# print errors so people know what to look up
#open (BAD_ENST_FILE, "$tmp/$pid.bad_enst_file.txt");
#my $line;
my $error_message_invalid = "Error! Ensembl gene annotation did not match genome sequence. Please try ncbi37, or annotate by hand.";
#while ($line = <BAD_ENST_FILE>) {
#        chomp ($line);
#        my @fields = split (/\t/, $line);
#        if ($COORD_SYSTEM eq "SPACE"){
#                # join all coordinates, do nothing
#        } else {
                # splice out beg coordinate
#                splice (@fields,1,1);
#        }
#        my $variant = join (",",@fields);
#        $variant =~ s/\,$//;
#        print OUTFILETSV $variant . "\t" . $error_message_invalid . "\n";
#        print OUTFILETABLE "<tr><td>$variant</td><td colspan=\"11\">$error_message_invalid</td></tr>\n";
#}
#close (BAD_ENST_FILE);

$outfile_table .= "</table>\n<BR>";


$outfile_table .= "</HTML>\n";

if ( $warningflag == 1 ) {

	$outfile_table .=
"<p class=\"red\">* Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.</p>";
}

 # my $outpage_table_tt = Template->new({
 #   ABSOLUTE => 1,
 # });


 # my $template_file = '/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/perl_templates/1_column_fixed_width.tpl';
 # my $vars = {
 #   main_content => $outfile_table,
 #   page_header => 'SIFT: Predictions',
 #   title => 'SIFT: Predictions',
 #   stylesheets => \@stylesheets,
 # };

 # my $outfile_table_template = '';
  
 # $outpage_table_tt->process($template_file, $vars, \$outfile_table_template) || die $outpage_table_tt->error();
  
#  print OUTFILETABLE $outfile_table_template;
print OUTFILETABLE $outfile_table; # Pauline Sept 2010 because I don't know the template;
  
#close(FILE);
close(OUTFILETABLE);
close(OUTFILETSV);
#print "before update status\n";



update_status($master_pid,$pid,"Complete");

#print "after update status\n";

system("rm -f $tmp/$pid.*.siftresults.predictions");

#my $html_combined = '';
my $html_combined = '<HTML><HEAD><link rel="stylesheet" type="text/css" href="css/tables.css" /></HEAD>';
my $html_combined_counter = 0;

#if final batch - combine all prediction files and update status and stats
open(TSV_COMBINED,">$tmp/$master_pid\_predictions.tsv")|| die ("Cannot open prediction tsv files for combining: $tmp/$master_pid\_predictions.tsv\n");
open(HTML_COMBINED,">$tmp/$master_pid\_predictions.html")||die("Cannot open prediction tsv files for combining");
print TSV_COMBINED "$heading_tsv\n";
$html_combined .= '<table border="1">';
$html_combined .= '<tr class="tableHeader"><td>';

$html_combined .= "$heading_html</td></tr>\n";


if ($last_partition == 1){
#	for ($i = $master_pid + 1; $i <= $pid; $i ++){

    foreach my $i (@list_of_pids) {

		open(TSV,"$tmp/$i\_predictions.tsv")|| die ("Cannot open prediction tsv files for combining");
		open(HTML,"$tmp/$i\_predictions.html")|| die ("Cannot open prediction tsv files for combining");
		while (<TSV>){
			#combined TSV file
			chomp;
			if ($_ =~ /Coordinates/i){
				next;
			}
			else{
				print TSV_COMBINED "$_\n";
			}
			$total_num++;		
			if ($_ =~ /CDS/){
				$CDS_num ++;
			}
			if ($_ =~ /TOLERATED/){
                                $tolerated_num ++;
                        }
			if ($_ =~ /DAMAGING/){
                                $damaging_num ++;
                        }
			if ($_ =~ /Nonsynonymous/){
                                $nonsynonymous_num ++;
                        }
			if ($_ =~ /Synonymous/){
                                $synonymous_num ++;
                        }

			if ($_ =~ /novel/){
                                $novel_num ++;
                        }

			#combined html file
			@fields = split( '\t', $_ );
			if ($_ =~ /Synonymous/ || $_ =~ /\w\d+\*/){
				$is_synonymous = 1;
			}
			else{
				$is_synonymous = 0;
			}
			$html_combined .= "<tr class=\"tableRow".($html_combined_counter++ % 2 == 0 ? "Even" : "Odd")."\">\n";
			$count = 0;
                	for $cell (@fields) {
				$count ++;
                        	if ( $cell =~ /DAMAGING/ || $cell =~ /Warning/ ) {
                                	$cell =  "<span class=\"red\">$cell</span>";
                        	}
                        	elsif ( $cell =~ /NOT PREDICTED/i || $cell =~ /not found/i ) {
                                	$cell  = "<span class=\"blue\">$cell</span>";
                        	}
				elsif ($cell =~ /(rs\d+):(.+?)/i){
                                	$dbsnp_id = $1;$mutation = $2;
                                	$cell = "<a href=\"http:\/\/www.ncbi.nlm.nih.gov\/sites\/entrez?db=snp&cmd=search&term=\+$dbsnp_id\" rel=\"external\">$dbsnp_id:$mutation</a>";
                        	}
				if($count >= 9 && $count <=12 && $is_synonymous == 1){
		                        $cell =  "N\/A";
                		}
#				if ($cell =~ /Error/) {
#					print HTML_COMBINED "<td colspan=\"11\">$cell</td>";
#				} else {
				$cell =~ s/^\s+//;
				$cell =~ s/\s+$//;
				if ($cell eq "") {
				    $cell = "&nbsp;";
				}
				$html_combined .= "<td><p>$cell</p></td>";
#				}			
                	}
		}
                $html_combined .= "</tr>\n";
	    }


    ### SIMNL<2011-11-11>: Now add non-coding annotations if exists

    my $noncoding_hash = ();
    my $noncoding_annot_file = "$tmp/$master_pid.allchrfile.noncoding_annot";
    if (-e $noncoding_annot_file && -s $noncoding_annot_file > 0) {


	# The noncoding_annot_file has many duplicate coords, we just take 1
	my %uniq_non_coding_hash = ();
	open(NONCODING_ANNOT, "<$noncoding_annot_file");
	while(my $line = <NONCODING_ANNOT>) {
	    chomp $line;
	    my ($chr, $coord, $noncoding_region, $transcript_ids) = split(/\t/, $line);
	    my $key = "$chr:$coord";  
	    $uniq_non_coding_hash{$key} = $line;
	}
	close(NONCODING_ANNOT);

	# Now get user query from the non coding list, get annotation from uniq_non_coding_hash based on user query
	my $noncoding_list = "$tmp/$master_pid.allchrfile.noncoding.list";
	open(NONCODE, "<$noncoding_list");
	while(my $line = <NONCODE>) {
	    chomp $line;
	    if ($line =~ /ONTENT-TYPE/ || $line =~ /^$/) { next; }
	    my ($user_query, @others) = split(/\t/, $line);
	    my $chr = "", my $coord1 = ""; my $coord2 = "";
	    if ($COORD_SYSTEM eq "RESIDUE") {
		($chr, $coord2, $strand, $alleles) = split(/,/, $user_query);
	    } else {
		($chr, $coord1, $coord2, $strand, $alleles) = split(/,/, $user_query);
	    }
	    my $coord = $coord1;

	    if ($COORD_SYSTEM eq "RESIDUE") {
		$coord = $coord2;
	    }

	    my $key = "$chr:$coord";

	    my $non_coding_info = $uniq_non_coding_hash{$key};
	    if (defined($non_coding_info)) {
		my ($chr, $coord, $noncoding_region, $transcript_ids) = split(/\t/, $non_coding_info);

		my @transcripts = split(/;/, $transcript_ids);
		my $first_transcript = "";
		if (scalar(@transcripts) > 0) {
		    $first_transcript = $transcripts[0];
		}

#		my $non_coding_line = "$user_query\t-\t$first_transcript\t\t\t\t\t\t\t\t$noncoding_region\t";
		my $non_coding_line = "$user_query\t-\t";
		if ($first_transcript =~ /ENS/) {
		    $non_coding_line .= "$first_transcript\t\t\t\t\t\t\t\t$noncoding_region\t";
		} elsif ($first_transcript =~ /NM/) {
		    $non_coding_line .= "\t$first_transcript\t\t\t\t\t\t\t$noncoding_region\t";
		} elsif ($first_transcript =~ /CCDS/) {
		    $non_coding_line .= "\t\t\t$first_transcript\t\t\t\t\t$noncoding_region\t";
		} else {
		    $non_coding_line .= "\t\t$first_transcript\t\t\t\t\t\t$noncoding_region\t";
		}

		print TSV_COMBINED "$non_coding_line\n";
		$html_combined .= "<TR>";		
		my $non_coding_html_line = "<TD>$user_query</TD><TD>-</TD>";
		if ($first_transcript =~ /ENS/) {
		    $non_coding_html_line .= "<TD>$first_transcript</TD><TD></TD><TD></TD><TD></TD>";
		} elsif ($first_transcript =~ /NM/) {
		    $non_coding_html_line .= "<TD></TD><TD>$first_transcript</TD><TD></TD><TD></TD>";
		} elsif ($first_transcript =~ /CCDS/) {
		    $non_coding_html_line .= "<TD></TD><TD></TD><TD></TD><TD>$first_transcript</TD>";
		} else { # UCSC 
		    $non_coding_html_line .= "<TD></TD><TD></TD><TD>$first_transcript</TD><TD></TD>";
		}

		for (my $i = 0; $i < 4; $i++) {
		    $non_coding_html_line .= "<TD></TD>";
		}
		$non_coding_html_line .= "<TD>$noncoding_region</TD>";
		for (my $i = 0; $i < 8; $i++) {
		    $non_coding_html_line .= "<TD></TD>";
		}
		$html_combined .= $non_coding_html_line;
		$html_combined .= "</TR>";
	    } #end if (defined($user_query))


	    
	} #end while(my $line = <NONCODE>)
	close(NONCODE);
    } #end if (-e $noncoding_annot_file && -s $noncoding_annot_file > 0)

    ### END SIMNL<2011-11-11>: Now add non-coding annotations if exists


	$html_combined .= "</table>\n";

    $html_combined .= "</html>\n";

	update_status($master_pid,$master_pid,"Complete");
	update_stats($master_pid,$total_num,$CDS_num,$tolerated_num,$damaging_num,$nonsynonymous_num,$synonymous_num,$novel_num);
}


 # Pauline commentd out because Template is unknown as of Sept 16,2010.

#  my $html_combined_tt = Template->new({
#    ABSOLUTE => 1,
#  });

#  $vars = {
#    main_content => $html_combined,
#    page_header => 'SIFT: Predictions',
#    title => 'SIFT: Predictions',
#    stylesheets => \@stylesheets,
#  };

#  my $html_combined_template = '';
  
#  $html_combined_tt->process($template_file, $vars, \$html_combined_template) || die $html_combined_tt->error();

# Pauline commentd out because Template is unknown as of Sept 16,2010.
# when known, can put it back in  
#  print HTML_COMBINED $html_combined_template;

  print HTML_COMBINED $html_combined;
 
close (TSV_COMBINED);
close (HTML_COMBINED);


#email the results
if ( $address ne "" && $last_partition == 1) {
	open( MESSAGE, ">$tmp/$master_pid.email_message.txt" );

        print MESSAGE "Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent submission at\n$url_address/www/sift/tmp/$master_pid\_predictions.html\nOR your can download the tab-delimited file at\n$url_address/www/sift/tmp/$master_pid\_predictions.tsv.\n" . 
	    "(Right click on the link above to save the file)\n\n" .
	    "Respond with this job id \"$pid_key\" for any future correspondance.\nThank you for using SIFT.\n\nSincerely\nSIFT Team\nGenome Institute of Singapore\n60 Biopolis Street\nSingapore 138672\n(Do not reply to this email)\n";

	close(MESSAGE);

	my $email_results = system("mutt -F .muttrc -s \"SIFT Results Results for Job ID $pid_key\" $address <$tmp/$master_pid.email_message.txt");

}



sub update_stats{
#	$tmp = "/opt/www/sift/tmp";
	$master_pid = @_[0];
	$total_num = @_[1];
	$CDS_num = @_[2];
	$tolerated_num = @_[3];
	$damaging_num = @_[4];
	$nonsynonymous_num = @_[5];
	$synonymous_num = @_[6];
	$novel_num = @_[7];
	$predicted_num = $tolerated_num + $damaging_num;
	
	$CDS_pct = 0;
	if ($total_num != 0) {
	    $CDS_pct = int($CDS_num * 100 / $total_num);
	}

	$predicted_pct = 0;
	if ($CDS_num != 0) {
	    $predicted_pct = int($predicted_num * 100 / $CDS_num) ;
	}

	$tolerated_pct = 0;
	if ($predicted_num != 0) {
	    $tolerated_pct = int($tolerated_num *100 / $predicted_num);
	}

	$damaging_pct = 100 - $tolerated_pct;

	$nonsynonymous_pct = 0;
	if ($CDS_num != 0) {
	    $nonsynonymous_pct = int($nonsynonymous_num*100/$CDS_num);
	}	
	$synonymous_pct = 100 - $nonsynonymous_pct;

	$novel_pct = 0;
	if ($total_num != 0) {
	    $novel_pct = int($novel_num * 100 / $total_num);
	}

	open (OUTPAGE,"$tmp/$master_pid.outpage.html") || die ("Cannot open outpage for updating status");
	open (OUTPAGE_MODIFIED,">$tmp/$master_pid.outpage.modified.html") || die ("Cannot open swap outpage file");
	while (<OUTPAGE>){
		chomp;
		if ($_ =~ /Number of input \(non-intronic\) variants/){
                        print OUTPAGE_MODIFIED "<p>Number of input (non-intronic) variants: $total_num <br />\n";
                }
		elsif ($_ =~ /Coding variants/ && $_ !~ /predicted/){
			print OUTPAGE_MODIFIED "Coding variants: $CDS_pct% ($CDS_num out of $total_num) <br />\n";
		}
		elsif($_ =~ /Coding variants predicted/){
			print OUTPAGE_MODIFIED "Coding variants predicted: $predicted_pct% ($predicted_num out of $CDS_num) <br />\n";
		}
		elsif($_ =~ /Tolerated/){
                        print OUTPAGE_MODIFIED "Tolerated: $tolerated_pct% ($tolerated_num out of $predicted_num) <br />\n";
                }
		elsif($_ =~ /Damaging/){
                        print OUTPAGE_MODIFIED "Damaging: $damaging_pct% ($damaging_num out of $predicted_num) <br />\n";
                }
		elsif($_ =~ /Nonsynonymous/){
                        print OUTPAGE_MODIFIED "Nonsynonymous: $nonsynonymous_pct% ($nonsynonymous_num out of $CDS_num) <br />\n";
                }
		elsif($_ =~ /Synonymous/){
                        print OUTPAGE_MODIFIED "Synonymous: $synonymous_pct% ($synonymous_num out of $CDS_num) <br />\n";
                }
		elsif($_ =~ /Novel/){
                        print OUTPAGE_MODIFIED "Novel: $novel_pct% ($novel_num out of $total_num) <br />\n";
                }
		else{
			print OUTPAGE_MODIFIED "$_\n";
		}

	}
	close (OUTPAGE_MODIFIED);
        close (OUTPAGE);
        system ("chmod 755 $tmp/$master_pid.outpage.modified.html");
        rename "$tmp/$master_pid.outpage.modified.html", "$tmp/$master_pid.outpage.html"

}

sub update_status{
#	$tmp = "/opt/www/sift/tmp";
	$master_pid = @_[0];
	$slave_pid = @_[1];

	my $slave_md5 = $slave_pid;
	if ($slave_pid =~ m/^(.*)\_nssnv/) {
	    $slave_md5 = $1;
	} else {
	    $slave_md5 = $slave_pid;
	}

	$status = @_[2];
	my $outpage_counter = 0;
	#print "$master_pid $slave_pid $status<BR>";
	open (OUTPAGE,"$tmp/$master_pid.outpage.html") || die ("Cannot open outpage for updating status");
	open (OUTPAGE_MODIFIED,">$tmp/$master_pid.outpage.modified.html") || die ("Cannot open swap outpage file");
	while (<OUTPAGE>){
        	chomp;
        	if ($_ =~ /.+?(Partitioned set .+?)\<\/p>\<\/td\>\<td\>\<p\>(.+?)\<.+?\>\<.+?\>\<.+?\>\<.+?\>$slave_pid\<.+?/){
                	$partition = $1;
                	$job_size = $2;

	                if ($status =~ /running/i){
	                        print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td><p>$job_size</p></td><td><p>$slave_pid</p></td><td><p class=\"blue\">Running..</p></td><td><p>Not available</p></td><td><p>Not available</p></td></tr>\n";
	                }
	                elsif($status =~ /complete/i){
#	                        print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td>$job_size</p></td><td><p>$slave_pid</p></td><td><p class=\"green\">Complete</p></td><td><p><a href=\"/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html\" rel=\"external\">$slave_pid table</a></p></td><td><p><A HREF=\"$url_address\/www\/sift\/tmp\/$slave_pid\_predictions.tsv\">$slave_pid results</A></p></td></tr>\n";
	                        print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td>$job_size</p></td><td><p>$slave_md5</p></td><td><p class=\"green\">Complete</p></td><td><p><a href=\"/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html\" rel=\"external\">$slave_md5 table</a></p></td><td><p><A HREF=\"$url_address\/www\/sift\/tmp\/$slave_pid\_predictions.tsv\">$slave_md5 results</A></p></td></tr>\n";

	                }
	        }
		elsif($_ =~ /.+?Complete set\<\/td\>\<td\>(.+?)\<.+?\>$slave_pid\<.+?/){

			$partition = "Complete set";
                        $job_size = $1;
                        if ($status =~ /running/i){
                                print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td><p>$job_size</p></td><td><p>$slave_pid</p></td><td><p class=\"blue\">Running..</p></td><td><p>Not available</p></td><td><p>Not available</p></td></tr></table>\n";
                        }
                        elsif($status =~ /complete/i){
#                                print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td>$job_size</p></td><td><p>$slave_pid</p></td><td><p class=\"green\">Complete</p></td><td><p><a href=\"/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html\" rel=\"external\">$slave_pid table</a></p></td><td><p><a href=\"$url_address\/www\/sift\/tmp\/$slave_pid\_predictions.tsv\">$slave_pid results</a></p></td></tr></table>\n";
                                print OUTPAGE_MODIFIED "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? "Even" : "Odd")."\"><td><p>$partition</p></td><td>$job_size</p></td><td><p>$slave_md5</p></td><td><p class=\"green\">Complete</p></td><td><p><a href=\"/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html\" rel=\"external\">$slave_md5 table</a></p></td><td><p><a href=\"$url_address\/www\/sift\/tmp\/$slave_pid\_predictions.tsv\">$slave_md5 results</a></p></td></tr></table>\n";

                        }
                }

	        else{
	                print OUTPAGE_MODIFIED "$_\n";
	        }
	}
	close (OUTPAGE_MODIFIED);
	close (OUTPAGE);
	system ("chmod 755 $tmp/$master_pid.outpage.modified.html");
	rename "$tmp/$master_pid.outpage.modified.html", "$tmp/$master_pid.outpage.html"
}
#-------------------------------------------------------------------------


# SIMNL<2011-12-06>: To cater to Hoi's CLI application
# The hardcoded path needs to be fixed later
if ($last_partition == 1 && $isCLI eq "true") { 
    if ($address) {
	system("php /var/www/gis/complete_job.php $message_handler $s3filename $master_pid\_predictions.tsv $address");
    } else {
	system("php /var/www/gis/complete_job.php $message_handler $s3filename $master_pid\_predictions.tsv");
    }
}



exit(0);


sub get_output_heading{
        $oo = @_[0] ;
        @elts = split /,/, $oo;


	# SIMNL <2011-10-27>: To cater to user's request to separate transcript and protein ids (enst from nm, ensp from np, etc.)
        $heading = "Coordinates\tCodons\t" . 
	    "Ensembl Transcript ID\tRefSeq Transcript ID\tKnown Transcript ID\tCCDS Transcript ID\t" . 
	    "Ensembl Protein ID\tRefSeq Protein ID\tKnown Protein ID\t" .
	    "Substitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t\# Seqs at position";

=pod original
        $heading =
"Coordinates\tCodons\tTranscript ID\tProtein ID\tSubstitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t\# Seqs at position";
=cut

        @options = ("Gene ID","Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease","Average Allele Freqs","CEU Allele Freqs");


	# Sim
	push (@options, "HCB Allele Freqs", "JPT Allele Freqs", "YRI Allele Freqs", "Avg. Allele Freq. (1000 Genomes)", "European Allele Freq. (1000 Genome)", "East Asian Allele Freq. (1000 Genome)", "West African Allele Freq. (1000 Genome)");
	push (@options, "South Asian Allele Freq. (1000 Genome)", "American Allele Freq. (1000 Genome)");
	# Sim end

        for ($i = 0 ; $i < scalar @elts; $i++){
                if (@elts[$i] eq "1"){
                        $heading.= "\t@options[$i]";
                }
        }
        $heading.="\tUser Comment";
        return $heading;
}


sub get_output_row {
        $oo = @_[0];
        $pred_line = @_[1];
        @elts = split /\t/, $pred_line;
        $chrom= @elts [0];
        $coord_begin = @elts[1];
        $coord_end = @elts[2];
        $orn = @elts[3];
        $alleles = @elts[4];
        $usr_comment = @elts[5];
        $proteinid = @elts[6];
        $subst = @elts[7];
        $ensp = @elts[8];
        $rsid = @elts[9];
        $region = @elts[10];
        if ($region =~ /CDS/){
                $region = "EXON CDS";
        }
        elsif($region =~ /3\'UTR/){
                $region = "3\' UTR";
        }
        elsif($region =~ /5\'UTR/){
                $region = "5\' UTR";
        }
        elsif($region =~ /3\'UTR/){
                $region = "3\' UTR";
        }

        $codons = @elts[11];
        $snp_type = @elts[12];
        $score = @elts[13];
        $median = @elts[14];
        $seqs_rep = @elts[15];

        # SIMNL <2011-10-27>: Because of the way we store rows with no sift scores
        # we need to add this check
        if ($score == -1 || $median == -1 || $seqs_rep == -1) {
            $score = "";
            $median = "";
            $seqs_rep = "";
        }

	if ($score eq "" || $score == -1) {
	    $pred = "NA";
	}
        elsif ($score <= 0.05){
                if ($median > 3.25){
                        $pred = "DAMAGING *Warning! Low confidence.";
                }
                else{
                        $pred = "DAMAGING";
                }
        }
        else{
                $pred = "TOLERATED";
        }
        if ($score eq ""){
                $pred = "Not scored";
                $score = "NA";
                $median = "NA";
                $seqs_rep = "NA";
        }
        if ($COORD_SYSTEM eq "SPACE"){
                $coords = "$coord_begin-$coord_end";
        }
        else{
                $coords = "$coord_end";
        }

	# SIMNL <2011-10-27>: To cater to user's request to separate transcript ids and protein ids into separate columns.
	# For some reason, transcript ids was set as $proteinid. Perhaps an artifact of older code
	$proteinid =~ s/,/\t/g;
	$ensp =~ s/,/\t/g;
	$proteinid =~ s/:/,/g;
	$ensp =~ s/:/,/g;
	# END SIM <2011-10-27>

        $table_row =
"$chrom,$coords,$orn,$alleles\t$codons\t$proteinid\t$ensp\t$subst\t$region\t$rsid\t$snp_type\t$pred\t$score\t$median\t$seqs_rep";

	
	if ($table_row =~ /NON-GENIC/) { $table_row =~ s/DAMAGING//g; }

        @elts2 = split /,/, $oo;


	# Now we add Gene Info elements if user wants them
	for(my $i = 0; $i < 10; $i++) {
	    if ($elts2[$i] == 1) {
		$table_row .= "\t" . $elts[16+$i];
	    }
	}

	# Next columns are HapMap allele frequencies
	for (my $i = 0; $i < 5; $i++) {
	    if ($elts2[$i + 10] == 1) {
		my $index = $i + 10;
		$table_row .= "\t" . $elts[26 + $i];
	    }
	}


	# And finally, the remaining columns are 1000 Genome allele frequencies
	for (my $i = 0; $i < 6; $i++) {
	    if ($elts2[$i + 15] == 1) {
		# We add 1000 Genomes allelel frequencies here
		$table_row .= "\t" . $elts[31+$i];
	    }
	}

	# Take note of $user_comment position.
        $table_row.="\t$usr_comment\n";

        return $table_row;
}



#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment
#variable
# returns: an associative array of name/value pairs.  The name is the
#key.
sub parse_query {
	local ($query_string) = @_;
	local ( %ans, @q, $pair );

	#print $query_string;
	# break up into individual name/value lines
	@q = split( /&/, $query_string );

	foreach $pair (@q) {

		# break the name/value pairs up
		# use split rather than regular expressions because the value may
		# have
		#  newlines in it
		split( /=/, $pair, 2 );

		# change '+' to ' '
		$_[1] =~ s/\+/ /g;

		# change the escaped characters (has to be after the split on '&'
		# and '=')
		$_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

		$ans{ $_[0] } = $_[1];
	}

	return %ans;
}

#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a
# string)
# returns: the decimal representation of the number
sub hextodec {
	unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:   CONTENT_TYPE
#               QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the
# key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#               Content-Disposition: form-data; name="key1"
#               <blank line>
#               val1
#               <boundary>
#               Content-Disposition: form-data; name="key2"
#               <blank line>
#               val2
#               <boundary>

sub general_parse {
	local ( $content_type, $query_string ) = @_;
	local ( %ans, @q, $pair, $loc, $boundary, $temp, $temp1 );

	if ( $content_type eq "application/x-www-form-urlencoded" ) {

		# break up into individual name/value lines
		@q = split( /&/, $query_string );

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
			my $content = "<p>Cannot parse</p>\n";
			$content .= "<p>content_type=$content_type</p>\n";
			$content .= "<p>query_string=$query_string</p>\n";
			&finish_script($content, -1);
		}
	}
	return %ans;

	#print "</PRE>";
}    # end of general_parse

# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub make_hash {
	my ($file) = @_;
	my %hash;
	open( HASH, $file ) || die "can't open $file";
	my $line;
	while ( $line = <HASH> ) {
		chomp($line);
		if ( exists( $hash{$line} ) ) {
			$hash{$line}++;
		}
		else {
			$hash{$line} = 1;
		}
	}
	close(HASH);
	return (%hash);
}

sub update_IP_logfile {
	my ( $queuefile, $IP_address ) = @_;

	$lockqueuefile = "$queuefile.lock";

	# lockfile will wait until it can lock the file
	`./lockfile $lockqueuefile`;

	# append the address and command to the queue file
	open( FILE, ">>$queuefile" );
	print FILE "$IP_address\n";
	close(FILE);

	chmod( 0664, $queuefile );

	# remove the lock file
	unlink($lockqueuefile);

}

sub round {
    my($number) = shift;
    return int($number + .5);
}

sub
make_new_variant
{
        my $v = VariantInfo->new();

        $v->rsid ("NA");
        $v->orn_db ("");
        $v->enst ("");
        $v->ensp ("");
        $v->region ("");
        $v->snp_type ("NA");
        $v->codon1 ("");
        $v->codon2 ("");
        $v->nt1 ("");
        $v->nt2 ("");
        $v->alleles_db => '$';
        $v->AA1 ("NA");
        $v->AA2 ("NA");
        $v->AAPOS ("NA");
        $v->score ("");
        $v->median  ("");
        $v->seqs_rep  ("");
        $v->freq_av  ("");
        $v->freq_ceu  ("");
        $v->ensg  ("");
        $v->enst  ("");
        $v->ensp  ("");
        $v->gene_name  ("");
        $v->gene_desc ("");
        $v->ensfm  ("");
        $v->fam_desc  ("");
        $v->gene_status  ("");
        $v->fam_size  ("");
        $v->kaks_mouse  ("");
        $v->kaks_macaque  ("");
        $v->mim_status  ("");

	$v->freq_hcb("");
	$v->freq_jpt("");
	$v->freq_yri("");

	# SIMNL<2011-10-31> (allele freq, empty if absent)
	$v->avg_1k_af("");
	$v->european_1k_af("");
	$v->east_asian_1k_af("");
	$v->west_african_1k_af("");
	$v->south_asian_1k_af("");
	$v->american_1k_af("");

        return $v;
}

sub
make_variant_with_elts
{
        my ($line) = @_;

        my @elts = split (/\t/, $line);
        my $v = VariantInfo->new();

        $v->rsid (@elts[4]);
        $v->enst (@elts[5]);
        $v->ensp (@elts[6]);
        $v->region (@elts[8]);
        $v->snp_type (@elts[9]);
        $v->nt1 (@elts[10]);
        $v->nt2 (@elts[11]);
        $v->alleles_db (@elts[10] . "\/" . @elts[11]);
#print "alleles db " . $v->alleles_db . "\n";
        $v->codon1 ( @elts[14]);
        $v->codon2 ( @elts[15]);
        $v->orn_db ( @elts[3]);
        $v->AA1 ( @elts[16]);
        $v->AA2 ( @elts[17]);
        $v->AAPOS ( @elts[19]);
        $v->score ( @elts[23]);
        $v->median ( @elts[24]);
        $v->seqs_rep ( @elts[25]);
        $v->ensg (@elts[26]);
        $v->gene_name (@elts[27]);
       $v->gene_desc (@elts[28]);
        $v->ensfm (@elts[29]);
        $v->fam_desc (@elts[30]);
        $v->gene_status (@elts[31]);
        $v->fam_size (@elts[32]);
        $v->kaks_mouse (@elts[33]);
        $v->kaks_macaque (@elts[34]);
        $v->mim_status (@elts[35]);
        $v->freq_av (@elts[36]);
        $v->freq_ceu (@elts[37]);

	$v->freq_hcb(@elts[38]);
	$v->freq_jpt(@elts[39]);
	$v->freq_yri(@elts[40]);

	# SIMNL<2011-10-28> (allele freq, empty if absent)	
	$v->avg_1k_af(@elts[41]);
	$v->european_1k_af(@elts[42]);
	$v->east_asian_1k_af(@elts[43]);
	$v->west_african_1k_af(@elts[44]);
	$v->south_asian_1k_af(@elts[45]);
	$v->american_1k_af(@elts[46]);

        return $v;
}

=pod
# SIM
sub create_thousand_genome_hash() {
    open (THOUSAND_GENOME_FILE, "<$tmp/$pid\_1K_genome.txt");
    my %hash = ();
    while (my $element = <THOUSAND_GENOME_FILE>) {
	chomp $element;
	my ($key, $values) = split("\t", $element);
	$hash{$key} = $values;
    }
    close(THOUSAND_GENOME_HASH);
    return \%hash;
}
=cut

sub get1000KeyElements() {
    my ($row) = @_;
    my ($chr, $a, $b, $c, $d, $rsid, @rest) = $row;
    return ($chr, $rsid);
}


sub getAlleleFreq() {

    my ($col1, $arrayRef) = @_;
    my $col2 = $col1 + 1;
    my @array = @{$arrayRef};

    my $freq1 = $array[$col1];
    my ($allele1, $frequency1) = split(":", $freq1);
    my $allele_freq1 = "";
    if ($frequency1 =~ /(0\.\d+)/) {
	my $f = sprintf("%.3f", $frequency1);
	$allele_freq1 = $allele1 . "," . $f;
    }

    my $freq2 = $array[$col2];
    my ($allele2, $frequency2) = split(":", $freq2);
    my $allele_freq2 = "";
    if ($frequency2 =~ /(0\.\d+)/) {
	my $f = sprintf("%.3f", $frequency2);
	$allele_freq2 = $allele2 . "," . $f;
    }

    if ($allele_freq1 eq "" && $allele_freq2 eq "") {
	return "";
    }

    $result = "$allele_freq1:$allele_freq2";
    return $result;
}


# SIMNL<2011-10-31>: Now the 1000 genomes data is in the Variant_Info struct
sub return_variant_line
{
        my ($v) = @_;
        my $subst ="";
        if ($v->AA1 eq "NA") {
                $subst = "NA";
        } else {
                $subst = $v->AA1 . $v->AAPOS . $v->AA2;
        }
        my $codonchange = $v->codon1 . "-" . $v->codon2;


	my $line = "";
	if (defined($v->enst)) { $line .= $v->enst . "\t"; } else { $line .= "\t"; }
	if (defined($subst)) { $line .= $subst . "\t"; } else { $line .= "\t"; }
	if (defined($v->ensp)) { $line .= $v->ensp . "\t"; } else { $line .= "\t"; }
	if (defined($v->rsid)) { $line .= $v->rsid . "\t"; } else { $line .= "\t"; }
	if (defined($v->region)) { $line .= $v->region . "\t"; } else { $line .= "\t"; }
	if (defined($codonchange)) { $line .= $codonchange . "\t"; } else { $line .= "\t"; }
	if (defined($v->snp_type)) { $line .= $v->snp_type . "\t"; } else { $line .= "\t"; }
	if (defined($v->score)) { $line .= $v->score . "\t"; } else { $line .= "-1\t"; }
	if (defined($v->median)) { $line .= $v->median . "\t"; } else { $line .= "-1\t"; }
	if (defined($v->seqs_rep)) { $line .= $v->seqs_rep . "\t"; } else { $line .= "-1\t"; }
	if (defined($v->ensg)) { $line .= $v->ensg . "\t"; } else { $line .= "\t"; }
	if (defined($v->gene_name)) { $line .= $v->gene_name . "\t"; } else { $line .= "\t"; }
	if (defined($v->gene_desc)) { $line .= $v->gene_desc . "\t"; } else { $line .= "\t"; }
	if (defined($v->ensfm)) { $line .= $v->ensfm . "\t"; } else { $line .= "\t"; }
	if (defined($v->fam_desc)) { $line .= $v->fam_desc . "\t"; } else { $line .= "\t"; }
	if (defined($v->gene_status)) { $line .= $v->gene_status . "\t"; } else { $line .= "\t"; }
	if (defined($v->fam_size)) { $line .= $v->fam_size . "\t"; } else { $line .= "\t"; }
	if (defined($v->kaks_mouse)) { $line .= $v->kaks_mouse . "\t"; } else { $line .= "\t"; }
	if (defined($v->kaks_macaque)) { $line .= $v->kaks_macaque . "\t"; } else { $line .= "\t"; }
	if (defined($v->mim_status)) { $line .= $v->mim_status . "\t"; } else { $line .= "\t"; }

	if (defined($v->freq_av)) { $line .= $v->freq_av . "\t"; } else { $line .= "\t"; }
	if (defined($v->freq_ceu)) { $line .= $v->freq_ceu . "\t"; } else { $line .= "\t"; }
	if (defined($v->freq_hcb)) { $line .= $v->freq_hcb . "\t"; } else { $line .= "\t"; }
	if (defined($v->freq_jpt)) { $line .= $v->freq_jpt . "\t"; } else { $line .= "\t"; }
	if (defined($v->freq_yri)) { $line .= $v->freq_yri . "\t"; } else { $line .= "\t"; }

	# SIMNL<2011-10-28>: Add 1000 genomes 
	# Pushing into @fields result in 
	if ($v->avg_1k_af ne "") { $line .= $v->avg_1k_af . "\t"; } else { $line .= "\t"; }
	if ($v->european_1k_af ne "") { $line .= $v->european_1k_af . "\t"; } else { $line .= "\t"; }
	if ($v->east_asian_1k_af ne "") { $line .= $v->east_asian_1k_af . "\t"; } else { $line .= "\t"; }
	if ($v->west_african_1k_af ne "") { $line .= $v->west_african_1k_af . "\t"; } else { $line .= "\t"; }
	if ($v->south_asian_1k_af ne "") { $line .= $v->south_asian_1k_af . "\t"; } else { $line .= "\t"; }
	if ($v->american_1k_af ne "") { $line .= $v->american_1k_af . "\t"; } else { $line .= "\t"; }

=pod
        my @fields ;
         push (@fields, $v->enst, $subst, $v->ensp, $v->rsid, $v->region,
                $codonchange, $v->snp_type, $v->score, $v->median, $v->seqs_rep, $v->ensg,
                $v->gene_name, $v->gene_desc, $v->ensfm, $v->fam_desc, $v->gene_status,
                $v->fam_size, $v->kaks_mouse, $v->kaks_macaque, $v->mim_status,
                $v->freq_av, $v->freq_ceu);

	# Sim
	push (@fields, $v->freq_hcb, $v->freq_jpt, $v->freq_yri);

#	push (@fields, $v->all_1K_genome, $v->european_1K_genome, $v->east_asian_1K_genome, $v->west_african_1K_genome);
#	push (@fields, $v->south_asian_1K_genome, $v->american_1K_genome);       
	# Sim end


        my $line = join ("\t", @fields);
=cut

        return ($line);


};

sub check_ip_counts
{
	
## Check that this IP address hasn't been used too much
my $IP_address = $ENV{REMOTE_ADDR};

#       print "<HR>" . $IP_address . "<BR></HR> ";
my $remote_host = $ENV{REMOTE_HOST};

	my $ip_counts =
`cat  /home/blocks/apache/logs/access_log  | grep POST | grep $IP_address | wc -l `;
chomp($ip_counts);
if ( $ip_counts == "" ) {
        $ip_counts =
`cat /home/blocks/apache/logs/access_log  | grep POST | grep $remote_host | wc -l`;
        chomp($ip_counts);
}
#       print $ip_counts. "<BR>";
my $upper_limit = 50;
if ( $address ne "" ) {
        $upper_limit = 1000;
}
if ( $ip_counts > $upper_limit ) {
	my $content = "<p class=\"header1\">Error</p>\n";
	$content .= "<p>Your computer has exceeded its daily limit.</p>\n";
    $content .= "<p>Please download <a href=\"/\">SIFT software</a> directly to your computer or <a href=\"/sift-bin/contact.pl\">contact</a> us so that we can help you.  Thank you for using SIFT. </p>";
        &finish_script($content, -1);
}

}

  ###########################################################################
  # This function calls the template object and injects the HTML variables
  # into the template. 
  #
  # @param : the HTML content
  # @return void
  ###########################################################################
  sub finish_script {
    my($content, $code) = @_;

    my $title = "SIFT: CHR Coords Submit";
    my @stylesheets = qw(/stylesheets/main.css);
 
    my $tt = Template->new({
      ABSOLUTE => 1,
    });

    my $template_file = '/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/templates/3_column_fixed_width.tpl';
    my $vars = {
      main_content => $content,
      title => $title,
      stylesheets => \@stylesheets,
    };

    open(STDOUT);
    $tt->process($template_file, $vars) || die $tt->error();
    
    exit($code);
  } # end sub finish_script
  ###########################################################################



sub getGeneDataForMulti() {
    my ($enst) = @_;

    $sth_db_strand_tx->execute ($enst);
    @rows_strand = $sth_db_strand_tx->fetchrow_array();
    $orn_tx = @rows_strand[0];
    $tx_enst_hash{$enst} = $orn_tx;

    $sth_db_supp_geneinfo->execute($enst);
    my @rows1 = $sth_db_supp_geneinfo->fetchrow_array();
    my $ensg = @rows1[2];
    my $gene_name = @rows1[3];
    my $gene_desc = @rows1[4];
    my $ensfm = @rows1[5];
    my $fam_desc = @rows1[6];
    my $gene_status = @rows1[7];
    my $fam_size = @rows1[8];
    my $kaks_mouse = @rows1[9];
    my $kaks_macaque = @rows1[10];
    my $mim_status = @rows1[11];

    my @results = ($ensg, $gene_name, $gene_desc, $ensfm, $fam_desc, $gene_status, $fam_size, $kaks_mouse, $kaks_macaque, $min_status);
    return \@results;
} #end getGeneDataForMulti()


sub grabAlleleFreq() {
    my $freq_av1 = $_[0];
    my $favg1 = "";
    if ($freq_av1 !~ /^\s+$/) {
        my ($theAllele1_avg, $theFreq1_avg) = split(",", $freq_av1);
	# SIMNL <2011-09-30> PE bug - the if-block causes cases with 0 allele freq to disappear
	# eg. T,1.000: instead of T,1.000:C,0.000
#        if ($theFreq1_avg =~ /([0|1]\.\d+)/) {
	if ($theFreq1_avg) {
            my $f = sprintf "%0.3f", $theFreq1_avg;
            $favg1 = "$theAllele1_avg,$f";
	}
#        }
    }
    return $favg1;
}


sub getAFLine() {
    my ($rsid_formatted) = @_;
    my $af_line = "";
    if ($rsid_formatted =~ /rs/){
	$sth_db_supp_allelefreq->execute($rsid_formatted);
        @af_rows = $sth_db_supp_allelefreq->fetchrow();
        my $refseq_id = shift(@af_rows);

        while (scalar(@af_rows) != 0) {
            my $allelefreq1 = shift(@af_rows);
            my $allelefreq2 = shift(@af_rows);
            $allelefreq1 =~ s/^\s+//; $allelefreq1 =~ s/\s+$//;
            $allelefreq2 =~ s/^\s+//; $allelefreq2 =~ s/\s+$//;
            if (($allelefreq1 eq "") || $allelefreq2 eq "") {
                $af_line .= "\t";
            } else {
                my ($allele1, $freq1) = split(",", $allelefreq1);
                my ($allele2, $freq2) = split(",", $allelefreq2);
                my $formatted_freq1 = "0.000"; my $formatted_freq2 = "0.000";
                if ($freq1 != 0) { $formatted_freq1 = sprintf "%0.3f", $freq1; }
                if ($freq2 != 0) { $formatted_freq2 = sprintf "%0.3f", $freq2; }
                $af_line .= "\t$allele1,$formatted_freq1:$allele2,$formatted_freq2";
            }
        } #end while
    } #end if ($rsid_formatted =~ /rs/)
    else {
	$af_line = "\t\t\t\t\t"; # Without rsid, we do not have HapMap alllele frequency results
    }
    if ($af_line eq "") {
	$af_line = "\t\t\t\t\t";
    }
    return $af_line;
}

# SIMNL<2011-10-27>: For new database created with multiple gene annotation sources
# Ensembl, RefSeq, Known, CCDS
sub uniqueIt() {
    my $identities = $_[0];
    if ($identities ne "") { 
	my @array = split(/,/, $identities);
	my %hash = ();
	foreach my $id (@array) {
	    $hash{$id} = 0;
	}
	my @uniq = keys(%hash);
	my $uniq_identities = join(":", @uniq);

	return $uniq_identities;
    } 
    return "";
} #end uniqueIt

sub getTranscriptionProteinIDs() {
    my ($aref) = @_;
    my @annotations = @{$aref};
    # UNIQ_ID, ENST, NM, KNOWN_TRANSCRIPT, CCDS, ENSP, NP, KNOWN_PROTEIN, ENSG
    my ($uniq_id, $ensts, $nms, $knowns, $ccdses, $ensps, $nps, $known_prots, $ensgs) = @{$aref};

    my $ensembl_transcripts = &uniqueIt($ensts);
    my $refseq_transcripts =  &uniqueIt($nms);
    my $known_transcripts = &uniqueIt($knowns);
    my $ccds_transcripts = &uniqueIt($ccdses);

    my $ensembl_proteins = &uniqueIt($ensps);
    my $refseq_proteins = &uniqueIt($nps);
    my $known_proteins = &uniqueIt($known_prots);

    return ($ensembl_transcripts, $refseq_transcripts, $known_transcripts, $ccds_transcripts, 
	    $ensembl_proteins, $refseq_proteins, $known_proteins);

}#end getTranscriptionProteinIDs

sub roundDecimals() {
    my ($allele_freq) = @_;
    if (!defined($allele_freq) || $allele_freq eq "") {
	return "EMPTY";
    } 
    my ($allele, $freq) = split(":", $allele_freq);
    my $new_frequency = sprintf("%.3f", $freq);

    return "$allele:$new_frequency";

} #end roundDecimals


sub getRowData() {

    my ($row, $index, $need_ensp_exists_in_key) = @_;

    @elts = split /\t/, $row;	
    
    $chr = @elts[0];
    $coord1 = @elts[1];
    $coord2 = @elts[2];
    $orn = @elts[3];
    $nt1 = @elts[10];
    $nt2 = @elts[11];
    $enst = @elts[6];
    $rsid = @elts[4];
    if($rsid =~ /(rs\d+)\:.+?/){
	$rsid_formatted = $1;
    }
    else{
	$rsid_formatted = "";
    }	    
    $CDS = @elts[20];
    $AA1_VALID = @elts[21];
    
    
    # Sim: This is needed when using new database with 1000 Genome population columns
    $row = "";
    my @tmpArray = ();
    my $size = scalar @elts;
    for (my $i = 0; $i < 26; $i++) {
	push @tmpArray, $elts[$i];
    }
    $row = join("\t", @tmpArray);
    
    # Now use uniq_key_identity to get the "enst" or "nm"
    my $uniq_key_identity = $uniq_key_id_array[$index];
    
    
    # SIMNL<2011-10-24>: Get transcript id from Gene Annotation Map
    $db_map_stm->execute($uniq_key_identity);
    my $transcript_identities = ""; # this will replace "ENST" column
    my $protein_identities = ""; # this will replace "ENSP" column
    
    my $ens_transcripts = ""; my $ref_transcripts = ""; my $known_transcripts = ""; my $ccds_transcripts = "";
    my $ens_proteins = ""; my $ref_proteins = ""; my $known_proteins = "";
    
    while (my @annotations = $db_map_stm->fetchrow_array()){	
	($ens_transcripts, $ref_transcripts, $known_transcripts, $ccds_transcripts,
	 $ens_proteins, $ref_proteins, $known_proteins) = &getTranscriptionProteinIDs(\@annotations);	    
    }
    
    # SIMNL<2011-10-25> Now concatenate
    
    if ($ens_transcripts && $ens_transcripts ne "") { $transcript_identities .= "$ens_transcripts,"; } else { $transcript_identities .= " ,"; }
    if ($ref_transcripts && $ref_transcripts ne "") { $transcript_identities .= "$ref_transcripts,"; } else { $transcript_identities .= " ,"; }
    if ($known_transcripts && $known_transcripts ne "") { $transcript_identities .= "$known_transcripts,"; } else { $transcript_identities .= " ,"; }
    if ($ccds_transcripts && $ccds_transcripts ne "") { $transcript_identities .= "$ccds_transcripts,"; } else { $transcript_identities .= " ,"; }
    
    if ($ens_proteins && $ens_proteins ne "") { $protein_identities .= "$ens_proteins,"; } else { $protein_identities .= " ,"; }
    if ($ref_proteins && $ref_proteins ne "") { $protein_identities .= "$ref_proteins,"; } else { $protein_identities .= " ,"; }
    if ($known_proteins && $known_proteins ne "") { $protein_identities .= "$known_proteins,"; } else { $protein_identities .= " ,"; }
    
    #SIMNL <2011-10-24>: This is required downstream, to check why.
    $ensp_exists = 1;
    if ($protein_identities ne "") { $ensp_exists = 0; }
    
    # SIMNL<2011-10-24>: Now we could either have ENST or NM_, need to adjust script accordingly.
    
    my $enst = ""; my $nm = ""; # If enst exists, then we use enst. Otherwise, we use nm. If neither exist,just print empty columns.
    if ($ens_transcripts ne "") {
	my @ensts = split(/:/, $ens_transcripts);
	$enst = shift(@ensts);
	
    } 
    if ($ref_transcripts ne "") {
	my @nms = split(/:/, $ref_transcripts);
	$nm = shift(@nms);
    }
    

    my $gene_info_cols = "";
    my @rows1 = ();
    my $has_data = "NO";
    if ($enst eq "" && $nm eq "") {
	$gene_info_cols = "\t\t\t\t\t\t\t\t\t\t";
    } else { 
	$has_data = "YES";
	if ($enst && $enst ne "") {
	    $query_using_enst->execute($enst);
	    @rows1 = $query_using_enst->fetchrow_array();
	} elsif ($nm && $nm ne "") {
	    @rows1 = $query_using_nm->execute($nm);
	}
    } #end else
    
    if ($has_data eq "YES") {
	$ensg = @rows1[2];
	$gene_name = @rows1[3];
	$gene_desc = @rows1[4];
	$ensfm = @rows1[5];
	$fam_desc = @rows1[6];
	$gene_status = @rows1[7];
	$fam_size = @rows1[8];
	$kaks_mouse = @rows1[9];
	$kaks_macaque = @rows1[10];
	$mim_status = @rows1[11];	
	$gene_info_cols = "\t$ensg\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status";
    }
    
    my $af_line = &getAFLine($rsid_formatted);
    
    my @data_in_row = split(/\t/, $row);
    
    # Replace uniq_key_id with "enst"
    if ($transcript_identities =~ /,$/) { chop $transcript_identities; }
    if ($protein_identities =~ /,$/) { chop $protein_identities; }
    
    $data_in_row[5] = $transcript_identities;
    $data_in_row[6] = $protein_identities;
    
    my $modified_row = join("\t", @data_in_row);
    
    
    # SIMNL<2011-10-28>: Get 1000 genome from one_thousand_genome_hash
    my $key_for_g1k = "$chr:$coord1:$coord2:$nt1:$nt2";
    my $g1k_info = $one_thousand_genome_hash{$key_for_g1k};
    
    $row2 = $modified_row . $gene_info_cols . $af_line;
    if ($g1k_info && $g1k_info ne "") {
	$row2 .= "\t" . $g1k_info;
    } else {
	$row2 .= "\t\t\t\t\t\t";
    }

    my $key = "$chr:$coord1:$coord2";
    if ($need_ensp_exists_in_key eq "YES") {
	$key = "$chr\t$coord1\t$coord2\t$CDS\t$ensp_exists";
    }

    return ($key, $row2, $transcript_identities);

} #end getRowData

sub prepareMultiTranscriptHash() {
    my ($query_aref) = @_;
    my @query_result = @{$query_aref};
    
    my %multi_hash = (); # <k,v> = <chr:coord1:coord2, @datum>

    # SIMNL <2011-11-01>: We need to get ENST using uniq id
    # Not all uniq ids will have ENST (e.g. transcripts found only in RefSeq)    
    my @results_with_transcript_ids = ();

    if ($multi_transcripts != 0) { # We populate multi_hash only if user asks for multi.
	
	my $num_results = scalar(@query_result);
	for (my $index = 0; $index < $num_results; $index++) {  
	    my $row = $query_result[$index];
	    
	    chomp $row;

	    my ($key, $row2) = &getRowData($row, $index, "NO"); # my $key = "$chr:$coord1:$coord2"; We do not want CDS and ensp_exists in the key


	    if (!defined($multi_hash{$key})) {
		my @array = (); 
		$multi_hash{$key} = \@array;
	    }
	    my $aref = $multi_hash{$key};
	    my @arr = @{$aref};
	    push @arr, $row2;
	    $multi_hash{$key} = \@arr;
	} #end for (my $index = 0; $index < $num_results; $index++)
    } #end if (isMultiTranscripts == "YES")


    # SNL <2011-08-01> At this point in time, the %multi_hash contains all returned values from the db
    # with chr:coord1:coord2 as key
    return \%multi_hash;
}
