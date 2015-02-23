#!/usr/bin/env perl

use strict;
use warnings;
use Math::Random qw(random_normal);
use Math::Random qw(random_poisson);
use POSIX qw(ceil);
use POSIX qw(floor);
use List::Util qw(min max);
use Data::Dumper;

my $size_argv = @ARGV;
if ($size_argv != 8 || $ARGV[0] eq "-h" || $ARGV[0] eq "--help") {
	print "Script started without enough arguments:\n";
	print "Paramters are required as follows:\n";
	print "\ttranscripts\t\tfasta file containing transcript sequences on which PAR-CLIP reads are simulated\n";
	print "\tout_prefix\t\tprefix for the output files. A fastq file as well as a log file are created\n";
	print "\terror_profile\t\tfile name to the error profile with mismatch probabilities\n";
	print "\tt2c_profile_file\tfile name to T-C conversion rates, ordered by the majority of occurence within a cluster\n";
	print "\tt2c_position_file\tfile name to T-C conversion site specific probabilites. Signifies each a probability how likely a T-C conversion is at the respective cluster position\n";
	print "\tquality_file\t\tfile name to base calling qualities; PHRED64 scaling\n";
	print "\tindel_file\t\tfile name to the indel profile, giving positions-specific probabilites for insertions and deletions\n";
	print "\tbound_prob\t\tdouble value specifying the fraction of RBP-bound clusters that are simulated\n";	
	exit;
}

print "Starting PAR-CLIP read simulation script\n";

my $fasta_file = $ARGV[0];
my $output_prefix = $ARGV[1];
my $output_file = $output_prefix . ".fastq";
my $log_file = $output_prefix . ".log";
my $error_log = $output_prefix . ".err";
my $error_profile_file = $ARGV[2];
my $t2c_profile_file = $ARGV[3];
my $t2c_position_file = $ARGV[4];
my $quality_file = $ARGV[5];
my $indel_file = $ARGV[6];
# for bound_prob, standard values are 1.0 for 100%T-C, 0.5 for 50% T-C simulating real par-clip reads and 0.0 for 0% T-C as control for PAR-CLIP read structure.
my $bound_prob = $ARGV[7];


## from 0 quality to 64?: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
my @byte_to_qual = ("!", "\"", "\#", "\$", "\%", "\&", "\'", "\(", "\)", "\*", "\+", "\,", "\-", "\.", "\/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "\:", "\;", "\<", "\=", "\>", "\?", "\@", "A","B", "C","D", "E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","\[","\\","\]","\^","\_","\`","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","{","|","}","~",")");

# how many reads should be collected from a fasta sequence; position should be selected uniformly distributed given the length of the read sequence
my $select_read = 0.5;

# set to 1 to allow for indels given by indels-file
my $allow_indels = 1;

# probability for exchanging a T with a C in the read sequence
#my $prob_T2C = 0.2;
my $max_number_T2C_positions = 4;
my $max_start_pos_per_cluster = 3;
my $max_end_pos_per_cluster = 3;

# probability for exchanging a single base as sequencing error; should depend on read length
#my $prob_error = 0.0;

# read in error-profile from file:
# format: 2-dim array, jeweils ACGT und darin dann die conversion-probabilities
my (%errors, %T2C, %qual, %sd);
my $run = 0;
my @ACGT = ("A", "C", "G", "T");
open(ERROR, "<$error_profile_file");
while (my $line = <ERROR>) {
	chomp($line);
	#$errors[$run] = $line;
	my @splitted_line = split('\s+', $line);
	my $no_error = $splitted_line[($run+1)%4] + $splitted_line[($run+2)%4] + $splitted_line[($run+3)%4];
	$splitted_line[$run] = 1 - $no_error;
	$errors{$ACGT[$run]} = \@splitted_line;
	$run++;
}
close(ERROR);

$run = 0;
open(T2C, "<$t2c_profile_file");
while (my $line = <T2C>) {
	chomp($line);
	$T2C{$run} = $line;
	$run++;
}
close(T2C);

$run = 0;
open(QUAL, "<$quality_file");
while (my $line = <QUAL>) {
	chomp($line);
	my @splitted_line = split('\t', $line);
	$qual{$run} = $splitted_line[0];
	$sd{$run} = $splitted_line[1];
	$run++;
}
close(QUAL);

my (%indels);
if ($allow_indels) {
	$run = 0;
	open(INDEL, "<$indel_file");
	while (my $line = <INDEL>) {
		chomp($line);
		my @splitted_line = split('\s+', $line);
		$indels{$run} = \@splitted_line;	
		$run++;
	}
	close(INDEL);
}

my %t2c_positions;
$run = 0;
open(POS, "<$t2c_position_file");
while (my $line = <POS>) {
	chomp($line);
	$t2c_positions{$run} = $line;
	$run++;	
}
close(POS);

# read length, take normal distributed between these lengths;
my $min_length = 7;
my $max_length = 30;

# gives the probability, that a read-cluster within a transcript is "really bound by the RBP", i.e. may contain T2C mutations
# ultra is 0.6; high rate was 0.45, mid-high is 0.35, mid was 0.25, low-mid is 0.15 and low is 0.05, zero is 0.0
#my $bound_prob = 0.0;

my ($header, $sequence, $T2C_errors, $mut, $most_T2C, $most_error, $reads, $av_reads_per_cluster, $clusters, $sum_read_length, $av_read_lenght, $num_indels, $simulated_bases);

# variables that are used to calculate a statistic on created dataset;
$sum_read_length = 0;
$clusters = 0;
$av_reads_per_cluster = 0;
$most_T2C = 0;
$most_error = 0;
$T2C_errors = 0;
$mut = 0;
$reads = 0;
$num_indels = 0;
$simulated_bases = 0;

# lookups for mutations and quality scores
my @ACG = ("A", "C", "G");
my @ACT = ("A", "C", "T");
my @AGT = ("A", "G", "T");
my @CGT = ("C", "G", "T");
my %ACGT_hash = ("A"=>0, "C"=>1, "G"=>2, "T"=>3);
# quality values correspond to PHRED33: 30 - 39
my @high_qual = ("^", "_", "`", "a", "b", "c", "d", "e", "f", "g");
# quality values correspond to PHRED33: 10 - 19
my @low_qual = ("J", "K", "L", "M", "N", "O", "P", "Q", "R", "S");

# go through FASTA file of sequences, collect a complete sequence with header and run "createReads" on that sequence
my $int = 0;
open(FASTA, "<$fasta_file");
open(OUT, ">$output_file");
open(LOG, ">$log_file");
open(ERR, ">$error_log");
while (my $line = <FASTA>) {
	chomp($line);
	
	if ($line =~ m/^>/) {
		$int++;
		if ($int % 10000 == 0) {
			print "$int passed\n";
		}
		if (defined($header) && defined($sequence)) {
			&createReads($header, $sequence);
		}
		$header = $line;
		$sequence = "";
	} else {
		$sequence .= $line;
	}
}

# calculate some statics for logging-output
$av_reads_per_cluster = $reads / $clusters;
$av_read_lenght = $sum_read_length / $reads;

# maybe add average number of reads per cluster, "mutated"/"bound" clusters etc
print LOG "number reads generated: $reads\nnumber bases simulated: $simulated_bases\naverage read-length: $av_read_lenght\nnumber clusters generated: $clusters\naverage reads per cluster: $av_reads_per_cluster\nT2C mutations occured: $T2C_errors\nsequencing errors occured: $mut\nread with most T2C: $most_T2C\nread with most errors: $most_error\nnumber indels generated: $num_indels\n";

print LOG "\nSome parameters:\nselect_prob=$select_read\nread bound by RBP probability: $bound_prob\n";

close(LOG);
close(OUT);
close(FASTA);
close(ERR);


sub createReads {
	# load parameters for the method
	my $header = shift;
	my $sequence = shift;
	
	#my @header_split = split("\\|", $header);
	#my $header_start = $header_split[2];
	#my $header_end = $header_split[3];

	if (rand() < $select_read) {
		# select reads from transcript sequence only on some of the inputted transcript sequences

		# decide how many clusters and how many reads per cluster should exist
		my $num_clusters = ceil(rand() * 3);
		my @num_reads_per_cluster = map int, random_normal($num_clusters, 12, 10);
		$clusters += $num_clusters;
			
		# for each cluster, create as many reads as it is supposed by random number. read position within cluster varies, decision on whether cluster "is bound by RBP".
		my $cluster_index = 1;
		foreach (@num_reads_per_cluster) {
			my $number_reads_generate = ceil($_);
			#$reads += $number_reads_generate;
			my $cluster_pos = ceil(rand() * (length($sequence) - ($max_length)));
			next if ($cluster_pos < 10 || (length($sequence) - $cluster_pos) < $max_length);
			
			# create only some start and end positions for reads to get a more realistic "cluster-structure"
			my $num_T2C_pos = ceil(rand() * $max_number_T2C_positions);
			# key -> value: bias -> allele (0,1,2)
			my %T2C_positions_biases;
			my $num_start_pos = ceil(rand() * $max_start_pos_per_cluster);
			my $num_end_pos = ceil(rand() * $max_end_pos_per_cluster);
			
			## IMPLEMENT CHECK SO THAT END POS CANNOT BE GREATER THAN CDS SEQUENCE!!! AND CHECK CLUSTER_POS INIT SOME LINES ABOVE, WHETHER *2 IS OK?!?!??!?!?!?
			my @start_pos = map int, random_normal($num_start_pos, $cluster_pos, 1);
			my @end_pos = map int, random_normal($num_end_pos, $cluster_pos + ($max_length - $min_length), 1);
			#while (max(@start_pos) >= min(@end_pos)) {
			#	@start_pos = map int, random_normal($num_start_pos, $cluster_pos, 1);
			#	@end_pos = map int, random_normal($num_end_pos, $cluster_pos + ($max_length - $min_length), 1);
			#}
			for (my $i = 0; $i < $num_end_pos; $i++) {
				if ($end_pos[$i] >= length($sequence)) {
					print "end pos behind length of CDS\n";
					splice(@end_pos, $i, 1);
				}
			}
			
			my $cluster_bound = 0;
			if (rand() < $bound_prob) {
				$cluster_bound = 1;
				my @Tpositions;
				# calculate which Ts are mutated
				# EXCHANGE MAX AND MIN BETWEEN START AND END TO DRIVE T2C WITHIN MAIN CLUSTER OR ALSO AT OVERHANGS?!?!?!?!?
				for (my $k = max(@start_pos); $k < min(@end_pos); $k++) {
					my $base = substr($sequence, $k, 1);
					if ($base eq "T") {
						push(@Tpositions, $k);
					}
				}
				#my $numMissed = 0;
				for (my $k = 0; $k < $num_T2C_pos; $k++) {
					my $numTsInCluster = @Tpositions;
					last unless ($numTsInCluster);
					#print "numTsincluster: $numTsInCluster\n";
					
					#my $currentPositionIndex = int(rand($numTsInCluster));
					#my $currentPositionIndex;
					my $newT2CPosition;
					if ($numTsInCluster == 1) {
						$newT2CPosition = 0;
					} else {
						$newT2CPosition = &get_t2c_position(\@Tpositions, max(@start_pos));
					}
					#print "NEW: #T in cluster: $newT2CPosition and the actual position within the read: " . $Tpositions[$newT2CPosition] . "\n";
					
					#my @currentPositionIndex = map int, random_poisson(1, $numTsInCluster / 2);
					#print "index: ".$currentPositionIndex[0]."\n";
					#while ($currentPositionIndex[0] < 0 || $currentPositionIndex[0] >= $numTsInCluster) {
						#$numMissed++;
					#	@currentPositionIndex = map int, random_poisson(1, $numTsInCluster / 2);
					#}
					#print "OLD: #T in cluster: ".$currentPositionIndex[0]." and the actual position within the read: " . $Tpositions[$newT2CPosition] . "\n";
					
					#if ($currentPositionIndex == $numTsInCluster) {
					#	print "match of " . $currentPositionIndex ." and " . $numTsInCluster . "\n";
					#}
					
					#push(@T2C_positions_biases, $Tpositions[$currentPositionIndex]);
					#print "test: " . $Tpositions[$currentPositionIndex] . "\n";
					$T2C_positions_biases{$Tpositions[$newT2CPosition]} = $T2C{$k};
					
					#print "currentPosIndex: " . $currentPositionIndex . "\n";
					#splice(@Tpositions, $currentPositionIndex[0], 1);
					splice(@Tpositions, $newT2CPosition, 1);
				}
				#print $numMissed . "\n";
			}
			
			# for each read, create sequence with errors, T2C and base calling quality.
			for (my $i = 0; $i < $number_reads_generate; $i++) {
				#my $offset = ($max_length - ceil((rand() * ($max_length - $min_length)))) / 3;
				#my $start = $cluster_pos;
				#if ($start > $offset) {
				#	if (rand() > 0.5) {
				#		$offset *= (-1);
				#	}
				#	$start += $offset;
				#} else {
				#	last;
				#}
				
				# introduce errors and report them! control parameters for large/small benchmarking test set, one for paper and one for unit tests
				#my $end = $start + ($max_length - ceil((rand() * ($max_length - $min_length))));
				
				my $start = $start_pos[floor(rand() * $num_start_pos)];
				my $end = $end_pos[floor(rand() * $num_end_pos)];
				if ($end - $start > $max_length || $start >= $end) {
					# skip these start and end positions because they are outranging prerequisites, caused by the normal-distribution selection in rare cases.
					#print "start $start, end $end\n";
					next;
				}
				
				my $read_seq_wt = substr($sequence, $start, $end-$start);
				my $read_qual = "";
				my $read_seq_mut = "";
				my $num_T2C = 0;
				my $num_error = 0;
				my $indel_set = 0;
				$sum_read_length += length($read_seq_wt);
				for (my $j = 0; $j < length($read_seq_wt); $j++) {
					# additionally pre-select the up to three alleles that are t2c mutated. calculate the alleles via $CLUSTER-start-position + $i?!?!?!
					#
					#
					#
					# PREPROCESS CLUSTERSEQUENCE AND SELECT UP TO THREE ALLELES. PUT IT INTO T2C_positions_biases ARRAY AS NOTED ABOVE!!!
					#
					#
					#
					
					my $current_base = substr($read_seq_wt, $j, 1);
					# if current base is a t2c allele, check T2C mutation and CONTINUE! DO NOT LET IT CHANGE BY AN ERROR AGAIN!!!!
					
					my $test = rand();
					#my $test = 0;
					my $errors_ref = $errors{$current_base};
					my $is_T = 0;
					if ($current_base eq "T") {
						$is_T = 1;
					}
					
					#if (exists($T2C_positions_biases[$j]) && $cluster_bound) {
					if (exists($T2C_positions_biases{$start + $j}) && $cluster_bound) {
						#my $test = rand();
						#if (!exists($T2C_positions_biases{$start + $j})) {
						#	print "start $start and j $j and test $test\n";
						#}
						if ($T2C_positions_biases{$start + $j} > $test) {
							# mutate!
							#print "mean: " . $qual{$j} . "\tsd: " . $sd{$j} . "\n";
							#print "J: $j\n";
							
							$read_qual .= &get_quality($qual{$j}, $sd{$j});
							#print "byte: " . $qual_rdm[0] . "\tqual: ".$byte_to_qual[$qual_rdm[0]]. "\n";
							$read_seq_mut .= "C";
							$T2C_errors++;
							$num_T2C++;
						} else {
							# should be a T - test?!?!?!
							#$read_qual .= $high_qual[floor(rand()*9)];
							$read_qual .= &get_quality($qual{$j}, $sd{$j});
							$read_seq_mut .= substr($read_seq_wt, $j, 1);
						}
						next;
					}
					
					#my $match_prob = $$errors_ref[$ACGT_hash{$current_base}];
					#my $first_mm_prob = $$errors_ref[$ACGT_hash{$current_base}] + $$errors_ref[($ACGT_hash{$current_base}+1)%4];
					#my $second_mm_prob = $$errors_ref[$ACGT_hash{$current_base}] + $$errors_ref[($ACGT_hash{$current_base}+1)%4] + $$errors_ref[($ACGT_hash{$current_base}+2)%4];
					
					#print "check for test=$test < match=$match_prob; < first_mm_prob=$first_mm_prob and < second_mm_prob=$second_mm_prob\n";
					
					unless (defined($ACGT_hash{$current_base})) {
						print "unhandled error found. Please submit the .err error-log to the developer. Thanks!\n";
						print ERR "unrecognized base in ACGT_hash=$current_base\n";
						print ERR "Sequence_header=$header\n";
						print ERR "Sequence=$sequence\n";
						
						# select random base for $current_base as this is not covered by $ACGT_hash
						$current_base = $ACGT[floor(rand() * 4)];
					}
					
					if ($test < $$errors_ref[$ACGT_hash{$current_base}]) {
						#match
						#$read_qual .= $high_qual[floor(rand()*9)];
						$read_qual .= &get_quality($qual{$j}, $sd{$j});
						$read_seq_mut .= substr($read_seq_wt, $j, 1);
						#next;
					} elsif ($test < $$errors_ref[$ACGT_hash{$current_base}] + $$errors_ref[($ACGT_hash{$current_base}+1)%4]) {
						#mut to first base
						#$read_qual .= $low_qual[floor(rand()*9)];
						$read_qual .= &get_quality($qual{$j}, $sd{$j});
						$read_seq_mut .= $ACGT[($ACGT_hash{substr($read_seq_wt, $j, 1)}+1)%4];
						$mut++;
						$num_error++;
						next;
					} elsif ($test < $$errors_ref[$ACGT_hash{$current_base}] + $$errors_ref[($ACGT_hash{$current_base}+1)%4] + $$errors_ref[($ACGT_hash{$current_base}+2)%4]) {
						#mut to second base
						#$read_qual .= $low_qual[floor(rand()*9)];
						$read_qual .= &get_quality($qual{$j}, $sd{$j});
						$read_seq_mut .= $ACGT[($ACGT_hash{substr($read_seq_wt, $j, 1)}+2)%4];
						#unless ($is_T) {
							$mut++;
							$num_error++;
						#} else {
							#$T2C_errors++;
							$num_T2C++;
						#}
						next;
					} else {
						#mut to third base
						#$read_qual .= $low_qual[floor(rand()*9)];
						$read_qual .= &get_quality($qual{$j}, $sd{$j});
						$read_seq_mut .= $ACGT[($ACGT_hash{substr($read_seq_wt, $j, 1)}+3)%4];
						$mut++;
						$num_error++;
						next;
					}
					
					if ($allow_indels) {
						my $test_indel = rand();
						my $current_indel_pos = $indels{$j};
						#print "indel_set=$indel_set\n";
						if (!$indel_set && $test_indel <= $$current_indel_pos[0]) {
							# INSERT INSERTION AT THIS POINT OF A RANDOM 25/25/25/25 BASE.
							$read_qual .= &get_quality($qual{$j}, $sd{$j});
							$read_seq_mut .= $ACGT[floor(rand()*4)];
							$num_indels++;
							$j--;
							$indel_set = 1;
							next;
						} elsif (!$indel_set && $test_indel <= $$current_indel_pos[1]) {
							# INSERT DELETION AT THIS POINT;
							$num_indels++;
							$indel_set = 1;
							next;
						}
					}
					
					#if ($cluster_bound && substr($read_seq_wt, $j, 1) eq "T" && rand() < $prob_T2C) {
					#	# introduce T2C depending on random number per each base in read
					#	$read_qual .= $high_qual[floor(rand()*9)];
					#	$read_seq_mut .= "C";
					#	$T2C++;
					#	$num_T2C++;
					#} elsif (rand() < $prob_error) {
					#	# introduce sequencing error depending on random number per each base in read
					#	$read_qual .= $low_qual[floor(rand()*9)];
					#	$read_seq_mut .= &mutate_base(substr($read_seq_wt, $j, 1));
					#	$mut++;
					#	$num_error++;
					#} else {
					#	# copy correct base and give high quality
					#	
					#}
					
					$simulated_bases++;
					
				}
				# collect some information for statistic
				if ($num_T2C > $most_T2C) {
					$most_T2C = $num_T2C;
				}
				if ($num_error > $most_error) {
					$most_error = $num_error;
				}
				
				# print out a 4-row FASTQ read
				#my $new_header = $header_split[0] . "\|" . $header_split[1] . "\|" . ($header_start + $start) . "\|" . ($header_end + $end);
				print OUT "\@SEQ_ID:$header\|$cluster_bound-$cluster_index:$i\n$read_seq_mut\n+\n$read_qual\n";
				$reads++;
			}
			$cluster_index++;
		}
	} else {
		return;
	}
}

sub get_quality {
	my $mean_qual = shift;
	my $sd_qual = shift;
	
	my @qual_rdm = map int, random_normal(1, $mean_qual, $sd_qual);
	if ($qual_rdm[0] > 64) {
		$qual_rdm[0] = 64;
	} elsif ($qual_rdm[0] <= 2) {
		$qual_rdm[0] = 3;
	}
	
	return $byte_to_qual[$qual_rdm[0]];
}

# creates a mutation on a given base; e.g. A to C or G or T.
sub mutate_base {
	# load parameters
	my $current_base = shift;
	
	if ($current_base eq "A") {
		return $CGT[floor(rand()*3)];
	} elsif ($current_base eq "C") {
		return $AGT[floor(rand()*3)];	
	} elsif($current_base eq "G") {
		return $ACT[floor(rand()*3)];
	} elsif ($current_base eq "T") {
		return $ACG[floor(rand()*3)];	
	}	
}

# selects a position within a cluster (CAREFUL: this is the position in the T-position-array, not the actual position within the read/cluster!!!) where a new T2C allele is introduced. is based on the distribution given by the t2c-positions file from real data.
sub get_t2c_position {
	my $t_positions_ref = shift;
	my $start_shift = shift;
	my $random = rand();
	my $sum_prob = 0;
	foreach my $t_position(@$t_positions_ref) {
		$sum_prob += $t2c_positions{$t_position-$start_shift};
		#print "T2C prob for actual seq: ". $t2c_positions{$t_position-$start_shift} . "\n";
	}
	#print "random number before norm: $random\n";
	$random *= $sum_prob;
	#print "random number after norm: $random\n";
	
	my $processed = 0.0;
	my $pos = 0;
	
	foreach my $t_position(@$t_positions_ref) {
		$processed += $t2c_positions{$t_position-$start_shift};
		if ($random <= $processed) {
			return $pos;
		}
		$pos++;
	}
	
	#foreach my $dist_val (keys(%t2c_positions)) {
	#	$processed += $t2c_positions{$dist_val};
	#	if ($random <= $processed) {
	#		return $pos;
	#	}
	#	$pos++;
	#}
	
	print "ERR: get_t2c_position does not sum up to 1????\n";
	return -1;
}

# REPORT ALL CASES WHERE ERROR, T2C, CHIMERIC ARE INTRODUCED!!!!!! just count numbers as output...



































