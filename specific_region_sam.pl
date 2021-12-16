#!/usr/bin/perl -w
use strict;
my $gtf_file="Arabidopsis_thaliana.TAIR10.44.gtf";

my %gene_start;
my %gene_end;
open(IN,"$gtf_file") or die "can,t open $gtf_file";
while (<IN>) {
    chomp(my $line=$_);
    unless ($line=~m/^#/){
        my ($tag,$start,$end,$descript)=(split /\t/,$line)[2,3,4,8];
        if ($tag eq "gene") {
            $descript=~s/"//gi;
            my $gene=(split /; /,$descript)[0];
            $gene=(split / /,$gene)[1];
           # my $pos=$start."_".$end;
			$gene_start{$gene}=$start;
			$gene_end{$gene}=$end;
        } elsif(($tag eq "transcript")|($tag eq "mRNA")){
            $descript=~s/"//gi;
            my $gene=(split /; /,$descript)[0];
            $gene=(split / /,$gene)[1];
			if (exists $gene_start{$gene}){
				my $tmp=$gene_start{$gene};
				if ($start<$tmp){
					$gene_start{$gene}=$start
				}
				my $tmp2=$gene_end{$gene};
				if ($end>$tmp2){
					$gene_end{$gene}=$end;
				}	
			}else{
				$gene_start{$gene}=$start;
				$gene_end{$gene}=$end;
			}
		}
    }

}
close IN;


my $intron=$ARGV[1];
my ($gene,$ch,$intron_start,$intron_end)=split /_/,$intron;
my $start=$gene_start{$gene};
my $end=$gene_end{$gene};

open OUT, ">$ARGV[2]";
open(IN,"$ARGV[0]") or die "can't open $ARGV[0]";
while (<IN>){
	chomp(my $line=$_);
	if($line=~m/^\@/){
		print OUT "$line\n";
		next
	}
	unless(($line=~m/NH:i:1\t/)||($line=~m/NH:i:1$/)){
		next
	}
	#print TEMP "$line\n";
	my ($read_name,$ch,$read_start,$code)=(split /\t/, $line)[0,2,3,5];
	if($code=~m/D|H|S|I/){
		next
	}
	my $read_exon2_end;
	if ($code=~m/N/){
		my @n_cishu=($code=~m/([0-9]+N[0-9]+)/g);
		my $n_cishu=$#n_cishu+1;
		if ($n_cishu==1) {
				my @n_cishu=($code=~m/([0-9]+[A-Z])/gi);
                my $read_e1_length=$n_cishu[0];
                my $read_intron_length=$n_cishu[1];
                my $read_e2_length=$n_cishu[2];
                $read_e1_length=~s/[A-Z]//gi;
                $read_intron_length=~s/[A-Z]//gi;
                $read_e2_length=~s/[A-Z]//gi;
                my $read_exon1_start=$read_start;
                #my $read_intron_start=$read_start+$read_e1_length;
                #my $read_exon1_end=$read_intron_start-1;
                my $read_exon1_end=$read_start+$read_e1_length-1;
                #my $read_intron_end=$read_start+$read_e1_length+$read_intron_length-1;
                my $read_exon2_start=$read_start+$read_e1_length+$read_intron_length;
                $read_exon2_end=$read_start+$read_e1_length+$read_intron_length-1+$read_e2_length;
		}else{
                my @zuhe;
                my $zuhe_index=-1;
                my $code_temp=$code;
                for(my $i=0;$i<=10;$i++){
                    if ($code_temp=~m/N/) {
                        $code_temp=~m/([0-9]+M[0-9]+N[0-9]+M)/;
                        my $pattern=$1;
                        $zuhe_index++;
                        $zuhe[$zuhe_index]=$pattern;
                        $pattern=~m/([0-9]+M[0-9]+N)/;
                        $code_temp=~s/^$1//;
                    }else{
                        last
                    }
                }
                my $read_start1=$read_start;
                foreach my $code1(@zuhe){
                    my @n_cishu=($code1=~m/([0-9]+[A-Z])/gi);
                    my $read_e1_length=$n_cishu[0];
                    my $read_intron_length=$n_cishu[1];
                    my $read_e2_length=$n_cishu[2];
                    $read_e1_length=~s/[A-Z]//gi;
                    $read_intron_length=~s/[A-Z]//gi;
                    $read_e2_length=~s/[A-Z]//gi;
                    my $read_exon1_start=$read_start1;
                    my $read_exon1_end=$read_start1+$read_e1_length-1;
                    my $read_exon2_start=$read_start1+$read_e1_length+$read_intron_length;
                    $read_exon2_end=$read_start1+$read_e1_length+$read_intron_length-1+$read_e2_length;
					$read_start1=$read_start1+$read_e1_length+$read_intron_length;
				}
		
		}
	}else{
		unless($code=~m/M/){
			next;
		}
		$code=~m/([0-9]+)M/;
		my $read_length=$1;
		$read_exon2_end=$read_start+$read_length-1;
	}
	my $overlap=if_overlap($start,$end,$read_start,$read_exon2_end);
	if($overlap>0){
		print OUT "$line\n";
	}
}
close OUT;
sub if_overlap{
    my @pos=@_;
	my ($start,$end,$read_start,$read_end)=@pos;
	my $overlap=0;
	if (($read_start>=$start)&($read_end<=$end)){
		$overlap=1;
	}
	if (($read_start<=$start)&($read_end>=$start)){
		$overlap=1;
	}
	if (($read_start<=$end)&($read_end>=$end)){
		$overlap=1;
	}
	return $overlap
}