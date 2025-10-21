#!/usr/bin/perl
use strict;

open (FI,"../../data/list_of_samples.txt");
#open (IN,"$fi")||die;
open (OUT,">../results/HMM/Ploidy_summary.tsv");

while (<FI>){
    chomp;
    my $na=$_;
    my $BW="../results/HMM/HMM_ploidy/$_.hmm.BW";
    my $AS="../results/HMM/HMM_ploidy/$_.hmm.CloneAS";
    my $VAF="../results/HMM/HMM_ploidy/$_.hmm.Clonesize";
    if (!(open(IN1,"$VAF") && open (IN1,"$BW") && open (IN1,"$AS"))){ print OUT "$na 1_clone 1_clone 1_clone 1_clone\n";next} #if we were not able to fit a mixture of distributions -> 1_clone
#    print "$VAF\n";
    open (IN1,"$VAF")||die;
    my %pl; 
    my $i;
    while (<IN1>){
        chomp;
	$i++;
	my @a=split/\s/;
	$pl{$i}=$_;
    }
    my $Ncl=2;
    if ($pl{'6'} <800 || $pl{'7'} <800){$Ncl=1;} # if one of the subclones contains less then 800 mutations -> 1_clone
    if ($Ncl==2 & $pl{'4'}>$pl{'5'}){$Ncl='ER'}

    open (IN1,"$BW")||die;
    my %bw1;my %bw2;
    my $i;
    my $bwp='Good';
    my $all1;
    my $all2;
    my %S1;
    my %S2;

    while (<IN1>){
        chomp;
        $i++;
        my @a=split/\s/;
	if ($i==1){
	    if($a[1]<0.67 || $a[4]<0.67 || $a[3]>0.33 || $a[6]>0.33 || $a[2]>0.6 || $a[5]>0.6 || $a[2]<0.4 || $a[5]<0.4) {$bwp='Bad';}}
	    if ($i==2){
		for my $ii (1..3){
		    my $lci=$ii+3;
		    $all1+=$a[$ii];
		    $all2+=$a[$lci];}
		
		for my $ii (1..3){
                    my $lci=$ii+3; my $c=$ii-1;
		    if ($a[$ii]<100 || $a[$lci] < 100){$bwp='Bad'}
		    $S1{$c}=$a[$ii]/$all1;
                    $S2{$c}=$a[$lci]/$all2;
		}
	    }
    }

    if ($all1<800 && $all2<800){$Ncl=1;}

    open (IN1,"$AS")||die;
    my @b; my $i;
    my $allAS; my %S12;my %exp;
    while (<IN1>){
        chomp;
        $i++;
        my @a=split/\s/;
	if ($i==1){@b=@a;
		   for my $k (0..$#a){$b[$k]=substr($a[$k],1,1); $b[$k].=substr($a[$k],4,1);}
	}

	if ($i==2){
	    for my $j (0..$#a){$allAS+=$a[$j]}
            for my $j (0..$#a){$S12{$b[$j]}=$a[$j]/$allAS}

	}}

##### rewrite adding Asymmetry;
    my $ASS='tetra';

    my %o_e; 
for my $f (0..2){
    for my $s (0..2){
#	next unless $f>$s;
#	my $sf="$s$f";
        my $fs="$f$s";
	    $exp{$fs}=($S1{$f}*$S2{$s});
	    $o_e{$fs}=0;
	    if ($exp{$fs}>0){
	        $o_e{$fs}=($S12{$fs})/($exp{$fs});}
    }}
my$o_e_20=$o_e{'20'};
my$o_e_02=$o_e{'02'};
if ($o_e{'00'}>1.5 && $o_e{'22'}>1.5 && $o_e{'20'}<0.6 && $o_e{'02'}<0.6){$ASS='Dip*'}
if ($o_e{'00'}>1.8 && $o_e{'22'}>1.8 && $o_e{'20'}<0.6 && $o_e{'02'}<0.6){$ASS='Dip'}


if ($o_e{'20'}>2 && $o_e{'02'}>2){$ASS='Symmetric'}


#    print OUT "$na $pl{6} $pl{7} $pl{4} $pl{5} $bwp $ASS A00_$S12{'00'})/($exp{'00'} A22_$S12{'22'})/($exp{'22'} A02_$S12{'02'})/($exp{'02'} A20_$S12{'20'})/($exp{'20'}\n";
print OUT "$na $bwp $o_e_20 $o_e_02 $ASS\n";

}