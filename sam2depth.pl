#!/usr/bin/perl -w
use strict;
use Time::HiRes qw( time );
my $start = time();
###input parameter

my $sam_file=$ARGV[0];

    my %pos_depth;
    my @canshu=@_;
    open(IN,"$sam_file") or die "can't open $sam_file";
    while (<IN>){
        chomp(my $line=$_);
        if($line=~m/^\@/){
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
                my $read_exon1_end=$read_start+$read_e1_length-1;
                my $read_exon2_start=$read_start+$read_e1_length+$read_intron_length;
                my $read_exon2_end=$read_start+$read_e1_length+$read_intron_length-1+$read_e2_length;
                for my $pos($read_start...$read_exon1_end){
                    $pos_depth{$ch}{$pos}++;
                }
                  for my $pos($read_exon2_start...$read_exon2_end){
                    $pos_depth{$ch}{$pos}++;
                }
            }else{
                my @zuhe;
                my $zuhe_type;
                my $zuhe_type_exon;
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
                my $zuhe_num=$#zuhe+1;
                my $zuhe_order=0;
                foreach my $code1(@zuhe){
                    $zuhe_order+=1;
                    if ($zuhe_order==1) {
                        $zuhe_type="First"
                    }elsif($zuhe_order==$zuhe_num){
                        $zuhe_type="Last"
                    }else{
                        $zuhe_type="Other"
                    }
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
                    my $read_exon2_end=$read_start1+$read_e1_length+$read_intron_length-1+$read_e2_length;
                    my $read_length=$read_e1_length+$read_e2_length;
                    $read_start1=$read_start1+$read_e1_length+$read_intron_length;
                    if ($zuhe_type=~m/First|Other/){
                        for my $pos($read_exon1_start...$read_exon1_end){
                           $pos_depth{$ch}{$pos}++;
                        }
                    }else{
                        for my $pos($read_exon1_start...$read_exon1_end){
                           $pos_depth{$ch}{$pos}++;
                        }
                        for my $pos($read_exon2_start...$read_exon2_end){
                           $pos_depth{$ch}{$pos}++;
                        }
                    }
                }
            }
        }else{
            unless($code=~m/M/){
                next;
            }
            $code=~m/([0-9]+)M/;
            my $read_length=$1;
            my $read_end=$read_start+$read_length-1;
            for my $pos($read_start...$read_end){
                $pos_depth{$ch}{$pos}++;
            }
        }
    }
    open OUT,">$ARGV[1]";
   
    foreach my $ch(keys %pos_depth){
        foreach my $pos(keys %{$pos_depth{$ch}}){
                my $depth=$pos_depth{$ch}{$pos};
                print OUT "$ch\t$pos\t$depth\n";
        }
    }

    close OUT;


sub min{
    my @numbers=@_;
    my $min=$numbers[0];
    foreach my $number (@numbers){
        if ($number < $min) {
            $min=$number
        }
    }
    return $min;
}
sub max{
    my @numbers=@_;
    my $max=$numbers[0];
    foreach my $number (@numbers){
        if ($number > $max) {
            $max=$number
        }
    }
    return $max
}
sub median{
    #仅限于3个值
    my @numbers=@_;
    my @number_new=sort {$a<=>$b} @numbers;
    my $median=$number_new[1];
    return $median;
}
sub median_true{
    my @numbers=@_;
    my @number_new=sort {$a<=>$b} @numbers;
    my $num=$#number_new+1;
    my $median;
    if (($num%2)>0) {
        $median=$number_new[int($num/2)]
    }else{
        $median=($number_new[int($num/2)]+$number_new[int($num/2)-1])/2
    }
    return $median;
}
sub mean{
    my @numbers=@_;
    my $num=$#numbers+1;
    my $sum=0;
    foreach my $value(@numbers){
        $sum+=$value
    }
    my $mean_value=$sum/$num;
    return $mean_value;
}