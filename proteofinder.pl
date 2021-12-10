#!/usr/bin/perl
#Script to detect Lantibiotic proteins in bacterial proteome.
#Can be adapted to invoke motif-based detection from any proteome
# Copyright (C) 2015 by
# Deepti Vipin, Mishra Lab, CCMB ,India
# All rights reserved
# Released under MIT license (see LICENSE.txt)

use strict;
use warnings;

# Opaltoolkit should be installed in perl.
#Import the Opal libraries to perl
use OpalServices;
use OpalTypes;
# Import other libraries
use Fcntl;
use File::Basename;
use File::Spec::Functions;
use LWP::Simple;


 #Submits a FASTA file of proteome to the MEME webservice FIMO
 #tool provided by NBCR and downloads the resulting output.

my ($location, $fimoService, $req, $commandline, $fh, $motif, $fasta,
  $fastaFile, $motifFile, $fastaInputFile, $motifInputFile, $status,
  $jobid, $output_dir, $dir_count, $index_url);
# Set the location of MEME service.
# A list of NBCR services can be found at http://ws.nbcr.net/opal2/dashboard?command=serviceList
$location = "http://nbcr-222.ucsd.edu/opal2/services/FIMO_4.9.1";

# Instantiate a new service object to interact with.  Pass it the service location
$fimoService = OpalServices->new(service_url => $location);

# Instantiate a new job request
$req = JobInputType->new();

# Input name of the FASTA file of proteome.
print "Please enter the name of the file: \n" ;
$fastaFile = <STDIN>;
chop($fastaFile);

#Fasta file containing motifs in MEME format(present in directory)..
$motifFile = "input/crp0.meme.xml";

#invoke MEME suite FIMO services
$commandline = "-upseqs $fastaFile $motifFile";


$req->setArgs($commandline);


# Get contents of local files to be used as input
sysopen($fh, $fastaFile, O_RDONLY);
$fasta = do {local $/; <$fh>}; # slurp file
close $fh;
sysopen($fh, $motifFile, O_RDONLY);
$motif = do {local $/; <$fh>}; # slurp file
close $fh;

# Set name and content of remote input file.
# Two input files are used, the motif file
#and the FASTA file.
$fastaInputFile = InputFileType->new($fastaFile, $fasta);
$motifInputFile = InputFileType->new($motifFile, $motif);
$req->setInputFile($motifInputFile, $fastaInputFile);

# Launch the job and retrieve job ID
print "Launching non-blocking FIMO job\n";
$status = $fimoService->launchJob($req);
$jobid = $status->getJobID();
print "Received Job ID: ", $jobid, "\n";


# Poll for job status
print "Polling job status\n";
while (1) {
  # print current status
  print "Status:\n";
  print "\tCode: ", $status->getCode(), "\n";
  print "\tMessage: ", $status->getMessage(), "\n";
  print "\tOutput Base URL: ", $status->getBaseURL(), "\n";
# STATUS_DONE || STATUS_FAILED
  last if ($status->getCode() == 8 || $status->getCode() == 4);

  print "Waiting 30 seconds\n";
  sleep(30);

  # Query job status
  $status = $fimoService->queryStatus($jobid);
}

if ($status->getCode() == 8) { # STATUS_DONE
  # Download the output files for the job.
  print "\nDownloading Outputs:\n\n\n";
  #$output_dir = "fimo_out";
  $dir_count = 0;
  if ($status->getBaseURL() =~ m/:\/\/[^\/]+\/(.*)$/) {
    my $path = $1;
    my @path_elems = grep {length($_) > 0} split('/', $path);
    $dir_count = scalar(@path_elems);
  }
  $index_url = $status->getBaseURL()."/fimo.txt";

 getstore($index_url, 'output/Fimo.txt') or die 'Unable to get page';





#Filter FIMO output based on Amino Acid sequence
#length and the p value of the FIMO output matching lantibiotics.

my( @array,%hash,%data, @id);


#extract fimo IDs to submit to Uniprot
my $fimo='output/Fimo.txt';
open(FIMO,$fimo) or die $!;
my $fimoId = 'log/log2.txt'; # File with UniProt identifiers.
open(ID,'>'.$fimoId) or die $!;

while (<FIMO>)
{
	chomp $_;
	if ($_=~/\|(.+)\|/)

	{
		print ID "$1\n";
	}
}

close(ID);
close(FIMO);


#open files to record process
my $output='log/log3.txt';#contains data downloaded from Uniprot
open(OUT,'>'.$output) or die"Cant open $output\n";
my $info='log/log4.txt'; #extracted ID, AA length
open(INFO, '>'.$info) or die "Cant open $info";
my $protein='log/log5.txt';#contains proteins less than 80AA residues
open(PROT,'>'.$protein) or die "Cant find $protein";
my $filter='log/log6.txt';#FIMO result of proteins<80AA
open(FIL,'>'.$filter) or die $!;
my $final='output/FINAL.txt';#FINAL OUTPUT
open(FINAL, '>'.$final) or die $!;

#open file for reading
open(ID,$fimoId) or die $!;

#submit list to uniprot to get length of protein
my $base = 'http://www.uniprot.org';
my $tool = 'batch';
my $agent = LWP::UserAgent-> new;
#my ($user, $pass) = <USERNAME> <PASS>
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/",
                            [ 'file' => [$fimoId],
                              'format' => 'txt',
                            ],
                            'Content_Type' => 'form-data');

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  $response = $agent->get($response->base);
}

#get response into output file
$response->is_success ?
  print OUT $response->content :
  die 'Failed, got ' . $response->status_line .
    ' for ' . $response->request->uri . "\n";
close(OUT);

#extract data containing AA length
open(OUT, $output) or die $!;

while(<OUT>){
	chomp $_;

	if($_=~/\bID\s/)
	{
		print INFO "$_\n";
	}
}
close(OUT);

#get proteins less than 80AA
open(INFO, $info) or die $!;

while(<INFO>){
		chomp;
	@array = split /\s+/;
    push @{$hash{$array[3]}}, $_; #AA length is in the 4th column. Refer log4.txt

}
#print FIMO result for protein len<80AA
foreach my $aa (keys %hash) {
	if($aa <80){
		print PROT join "\n",@{$hash{$aa}};
		print PROT "\n";
	}
}
close (INFO);

open(PROT, $protein) or die $!;
#store the IDs of proteins<80AA to array
while (<PROT>){
	chomp;
	my @group= split/\s+/;
	push (@id,$group[1]);

}

#Extract the FIMO results of these IDs
open(FIMO,$fimo) or die $!;

chomp @id;

while (<FIMO>){
	chomp $_;
	$_=~/\|.+\|(\S+)/; #pattern tr|Q03ZY5|Q03ZY5_LEUMM
	if($1~~@id)
	{
		print FIL $_."\n";

	}
}

close (PROT);
close(FIL);

#Filter from them proteins with score> -10
open(FIL,$filter) or die $!;
while (<FIL>){
	chomp;
	my @lines= split/\s+/;
	push @{$data{$lines[5]}}, $_; #p value used as keys
}
#print FINAL "pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence\n";
foreach my $score(keys %data){
	if ($score>-10)
	{
	print FINAL join "\n", @{$data{$score}};
	print FINAL "\n";
}
}
print "Job Done!";
close(FIMO);
}
