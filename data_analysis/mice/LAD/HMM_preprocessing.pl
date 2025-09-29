#!/usr/bin/perl
use strict;

my @ar=</workspace/projects/lesion_segregation/mice/data/mutations/*.nodMat>;

my %n;
my %g;

$n{'A'}=12;
$n{'C'}=13;
$n{'G'}=14;
$n{'T'}=15;
$g{12}='A';
$g{13}='C';
$g{14}='G';
$g{15}='T';

foreach my $fi (@ar){
	my %dat;
	print "$fi\n";
    open(my $IN, "<", $fi ) ||die "Can't open $fi: $!\n";
	my $samplename = $fi;
    $samplename =~ s{.*/}{};    
    $samplename =~ s/\.[^.]+$//; 
	print "$samplename\n";

    while (<$IN>){
	chomp;
	my @a=split/\,/; my $b='B';
	next unless length($a[6])==1 && length($a[7])==1;
	my $mut="$a[6]\_$a[7]"; my $cov=$a[11];
	my $mut1="$a[6]\_N";
	for my $i (12..15){if ($n{$a[7]} == $i || $n{$a[6]} == $i ){next}
			   if ($a[$i]>2){$mut.=",$g{$i}"; $b='M'; $cov.=",$a[$i]"}}
    
	${$dat{$a[0]}}{$a[1]}.="$mut $a[9] $a[10] $cov $a[4] $b $mut1\t";

	#print "$mut $a[9] $a[10] $cov $a[4]\n";
    }

	foreach my $chr (sort keys %dat){
    	open (OUT,">>/workspace/projects/lesion_segregation/mice/MRCA/results/HMM/input_HMM_test/$samplename.hmm");
    	foreach my $cor (sort {$a<=>$b} keys %{$dat{$chr}}){
			my @b=split/\t/,${$dat{$chr}}{$cor};
			for my $b (0..$#b){
	    		print OUT "$chr $cor $b[$b]\n";
			}
		}
	}
}



