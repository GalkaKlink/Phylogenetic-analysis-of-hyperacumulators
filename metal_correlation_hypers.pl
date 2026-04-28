use strict;
open (IN1,"ElementHyper.GENERA.ANGIOSPERMS.tsv") or die $!;
open (OUT, ">>GenusMetalCorrelation_Hypers_100000H0.Angiosperms") or die $!;
my %me_genera;
while (<IN1>) {
    my $str = $_;
    my @all = $str =~ m/\w+/g;
    my $me = $all[0];
    my $genus = $all[1];
     push(@{$me_genera{$me}},$genus);
}


open (IN0,"AllAngiospermGenera.tsv") or die $!;
my @genera;
while (<IN0>) {
    my $str = $_;
    my @all = $str =~ m/\w+/g;
    my $genus = $all[0];
    push(@genera, $genus);
}

my @metals = keys %me_genera;
@metals = sort { $a <=> $b } @metals;
foreach my $me_n1(0..($#metals-1)) {
    my $me1 = $metals[$me_n1];
    my @genera1 = @{$me_genera{$me1}};
    foreach my $me_n2(($me_n1+1)..$#metals) {
	my $me2 = $metals[$me_n2];
	my @genera2 = @{$me_genera{$me2}};
	my $common_gen = 0;
	foreach my $gen1(@genera1) {
	    if ($gen1 ~~ @genera2) {
		$common_gen+=1;
	    }
	}

	
	my $n = 0;
	my @common_H0;
	until ($n==100000) {
	    $n +=1;
	    
	    my @new_genera1;
	    my @used1;
	    until($#new_genera1 == $#genera1) {
		my $rand = int(rand($#genera+1));
		if ($rand ~~ @used1) {
		}
		else {
		    push(@used1,$rand);
		    my $genus = $genera[$rand];
		    push(@new_genera1,$genus);
		}
	    }
	    
	    my @new_genera2;
	    my @used2;
	    until($#new_genera2 == $#genera2) {
		my $rand = int(rand($#genera+1));
		if ($rand ~~ @used2) {
		}
		else {
		    push(@used2,$rand);
		    my $genus = $genera[$rand];
		    push(@new_genera2,$genus);
		}
	    }
	    
	    my $new_common_gen = 0;
	    foreach my $gen1(@new_genera1) {
		if ($gen1 ~~ @new_genera2) {
		    $new_common_gen+=1;
		}
	    }
	    push(@common_H0,$new_common_gen);
	}
###
	my @sort_H0 = sort{$b <=> $a} @common_H0;
	
	my $H0_1 = $sort_H0[0];
	my $max_H0 = $H0_1;
	my $m=0;
	foreach my $H0(@sort_H0) {
	    if ($H0 == $H0_1) {
		$m+=1;
	    }
	    else { 
		last;
	    }
	}
	my $min_corr_pv =$m/($#common_H0+1);
      

	my $pv_uncorr = 1;
	foreach my $H0_n(0..$#sort_H0) {
	    my $H0 = $sort_H0[$H0_n];
	    if ($H0<=$common_gen) {
		$pv_uncorr = $#sort_H0-$H0_n+1;
		last;
	    }
	}


	my @sort_H0 = sort{$a <=> $b} @common_H0;
	my $H0_1 = $sort_H0[0];
	my $m=0;
	my $min_H0 = $H0_1;
	foreach my $H0(@sort_H0) {
	    if ($H0 == $H0_1) {
		$m+=1;
	    }
	    else { 
		last;
	    }
	}
	my $min_uncorr_pv =$m/($#common_H0+1);
     
	my $pv_corr = 1;
	foreach my $H0_n(0..$#sort_H0) {
	    my $H0 = $sort_H0[$H0_n];
	    if ($H0>=$common_gen) {
		$pv_corr = $#sort_H0-$H0_n+1;
		last;
	    }
	}
	my $size1 = $#genera1+1;
	my $size2 = $#genera2+1;
	my $pv_corr1 = $pv_corr/($#common_H0+1);
	my $pv_uncorr1 = $pv_uncorr/($#common_H0+1);
	print $me1."\t".$me2."\t".$pv_corr."\t".$pv_uncorr."\t".$#common_H0."\t".$pv_corr1."\t".$pv_uncorr1."\n";


	print OUT $me1."\t".$me2."\t".$size1."\t".$size2."\t".$common_gen."\t".$pv_corr1."\t".$pv_uncorr1."\t".$min_H0."\t".$max_H0."\n";
    }
} 

	
