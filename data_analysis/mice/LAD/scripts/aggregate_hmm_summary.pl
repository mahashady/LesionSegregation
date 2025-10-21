#!/usr/bin/perl
use strict;

my @ar1=<../results/HMM/output_HMM/*.hmm.summary>;
my @ar2=<../results/HMM/output_HMM/*.hmm.bwPlEM>;
my @ar3=<../results/HMM/output_HMM/*.hmm.bwASEM>;

my %h1;
my %h2;
my %h3;


foreach my $fi (@ar1){

    open (IN,"$fi")||die;
    while (<IN>){
	chomp;
	my @a=split/\s/;
	$h1{$a[0]}=$_;
    }
}


foreach my $fi (@ar2){
	my @na=split/\//,$fi;
	my @na1=split/\./,$na[-1];
        $na1[0].=".hmm";
	open (IN,"$fi")||die;
	while (<IN>){
	    chomp;
	    my @a=split/\s/;
	    $h2{$na1[0]}.="$a[0] ";
	}
	my $l=$h2{$na1[0]};
	substr($l,-1)='';
	$h2{$na1[0]}=$l}


foreach my $fi (@ar3){
    my @na=split/\//,$fi;
    my @na1=split/\./,$na[-1];
	$na1[0].=".hmm";
    open (IN,"$fi")||die;
    while (<IN>){
        chomp;
        my @a=split/\s/;
        $h3{$na1[0]}.="$a[0] ";
        }
    my $l=$h3{$na1[0]};
    substr($l,-1)='';
    $h3{$na1[0]}=$l}


open (OUT,">../results/HMM/HMM_summary.tsv");
foreach my $k (sort keys %h1){
	print OUT "$h1{$k} $h2{$k} $h3{$k}\n";
}