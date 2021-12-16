#!/usr/bin/perl -w
use strict;
my %junction_count;
my %junction_count_other;
open(IN,"$ARGV[0]") or die "can't open $ARGV[0]";
while (<IN>) {
    chomp (my $line=$_);
    my ($ch,$junction,$count,$count_other)=(split /\t/, $line)[0,1,3,4];
    $junction_count{$ch}{$junction}=$count;
    if ($count_other>0) {
        $junction_count_other{$ch}{$junction}=$count_other;
    }
    
}
close IN;
my %no_ref;
open(IN,"$ARGV[4]") or die "can't open $ARGV[4]";
while (<IN>) {
    chomp (my $line=$_);
    my ($junction,$count)=split /\t/, $line;
        my $ch=(split /_/,$junction)[0];
        $no_ref{$ch}{$junction}=$count;
    
    
}
close IN;
open OUT, ">$ARGV[3]";
open(IN,"$ARGV[1]") or die "can't open $ARGV[1]";
#1_704495_704587
while (<IN>) {
    chomp (my $line=$_);
    my ($gene,$ch,$intron,$exon1,$exon2,$psi)=(split /\t/, $line)[0,1,2,3,4,12];
    unless ($psi eq "-"){
        $psi=sprintf "%.2f",$psi;
    }
    
    my $intron2=$ch."_".$intron;
    if ($intron2 eq "$ARGV[2]") {
        my ($intron_start,$intron_end)=(split /_/,$intron);
        my ($exon1_start,$exon1_end)=split /_/,$exon1;
        my ($exon2_start,$exon2_end)=split /_/,$exon2;
        my $x_lab=$gene." ".$ch.": ".$intron_start."-".$intron_end." (PSI=$psi)";
        my $e1_i_junction=$exon1_start."_".$exon1_end."_".$intron_start."_".$intron_end;
        my $i_e2_junction=$intron_start."_".$intron_end."_".$exon2_start."_".$exon2_end;
        my $e1_e2_junction=$exon1_start."_".$exon1_end."_".$exon2_start."_".$exon2_end;
        my ($e1_i_junction_count,$i_e2_junction_count,$e1_e2_junction_count)=(0,0,0); 
        if (exists $junction_count{$ch}{$e1_i_junction}) {
            $e1_i_junction_count=$junction_count{$ch}{$e1_i_junction};
        }
        if (exists $junction_count{$ch}{$i_e2_junction}) {
            $i_e2_junction_count=$junction_count{$ch}{$i_e2_junction};
        }  
        if (exists $junction_count{$ch}{$e1_e2_junction}) {
            $e1_e2_junction_count=$junction_count{$ch}{$e1_e2_junction};
        }
        if (exists $junction_count_other{$ch}{$e1_e2_junction}) {
            foreach my $junction(keys %{$no_ref{$ch}}){
                my ($start,$end)=(split /_/,$junction)[1,2];
                if (($exon1_end==$start)&&($exon2_start!=$end)) {
                    my $no_ref_count=$no_ref{$ch}{$junction};
                    print OUT "$intron_start\t$end\t$no_ref_count\tNew junction reads\t$exon1_start,$exon1_end,$end,$exon2_start,$exon2_end\t$x_lab\n";
                }elsif(($exon1_end!=$start)&&($exon2_start==$end)){
                    my $no_ref_count=$no_ref{$ch}{$junction};
                    print OUT "$start\t$intron_end\t$no_ref_count\tNew junction reads\t$exon1_start,$start,$exon1_end,$exon2_start,$exon2_end\t$x_lab\n";
                }
                
                
            }
        }
        
        my $e1_tag=$exon1_end-1;
        my $i_tag1=$intron_start+1;
        my $i_tag2=$intron_end-1;
        my $e2_tag=$exon2_start+1;
        
        print OUT "$intron_start\t$intron_end\t$e1_e2_junction_count\tExclusion reads\t$exon1_start,$exon1_end,$exon2_start,$exon2_end\t$x_lab\n";
        print OUT "$e1_tag\t$i_tag1\t$e1_i_junction_count\tInclusion reads\t$exon1_start,$exon1_end,$exon2_start,$exon2_end\t$x_lab\n";
        print OUT "$i_tag2\t$e2_tag\t$i_e2_junction_count\tInclusion reads\t$exon1_start,$exon1_end,$exon2_start,$exon2_end\t$x_lab\n";
    }
    
}


close IN;
close OUT;
