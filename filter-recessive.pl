#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

cat my-snver-output.vcf | filter-recessive.pl 

Description:

This script filters a SNVER pooled file.  

WARNING: the recessive mutant must be the first genotype column followed by the 
reference strain.  Any other ordering will not work.  Additional columns will be ignored

Output:

    seqid                      : The contig / scaffold / chromosome name
    position                   : The position of the mutation
    reference allele frequency : The frequency of the non-reference mutation in the control pool
    mutant alllele frequency   : The frequency of the non-reference mutation in the mutant pool
    delta af                   : The allele frequency difference between the mutant and control pool
    fail - mutant in control   : We do not expect to see the mutation in the reference pool [0 = pass, 1 = fail]
    fail - mutant af too low   : Under the recessive model we expect the allele frequency to be 1.0, however 
                                 sequencing errors could result in lower allele frequencies so the cutoff is 0.75 [0 = pass, 1 = fail]
    fail - depth               : If the control or mutant pool has a depth below 5 [0 = pass, 1 = fail]
    fail flags                 : The last three columns concatenated

Contact:

   Problems or questions: zev.kronenberg\@gmail.com

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success ;


while (<STDIN>) {
    chomp;
    next if $_ =~ /^\#/;
    next if $_ =~ /NA/;
    my @l = split /\t/, $_;
    my @anc = split /:/, $l[9];
    my @mut = split /:/, $l[10];

    my $anc_af = $anc[0] / $anc[1];
    my $mut_af = $mut[0] / $mut[1];
    my $deltAf = abs($mut_af - $anc_af);
    
    my $allFlag;

    my $fail_depth = 0;
    if($anc[1] < 5 || $mut[1] < 5){
	$fail_depth = 1;
    }

    my $fail_anc_af = 0;
    if($anc_af > 0){
	$fail_anc_af = 1;
    }
    my $fail_mut_af = 0;
    if($mut_af < 0.75){
	$fail_mut_af = 1;
    }
    
    $allFlag .= $fail_depth;
    $allFlag .= $fail_anc_af;
    $allFlag .= $fail_mut_af;

    print $l[0]  , "\t";
    print $l[1]  , "\t";
    print $anc_af, "\t";
    print $mut_af, "\t";
    print $deltAf, "\t";
    print $fail_anc_af, "\t";
    print $fail_mut_af, "\t";
    print $fail_depth,  "\t";
    print $allFlag; 
    print "\n";

}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------


