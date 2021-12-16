#!/usr/bin/perl -w
use strict;
use Time::HiRes qw( time );
use Getopt::Long;
my $usage="$0 [options]
--sam		<str>           SAM file
--thread        <int>           the number of thread
--m             <int>           minimum length between read and exon or intron(default:8)
--outdir        <str>           the path of outdir
        ";
my($sam,$thread,$m,$outdir);
GetOptions("sam:s" =>\$sam,
        "thread:i" =>\$thread,
        "m:i" =>\$m,
        "outdir:s" =>\$outdir


);

die $usage if !defined $sam;

$m ||=8;








my $start = time();

my $sam_file=$sam;
my $th_p=$thread;
my $min_overlap=$m;
my $out_dir=$outdir;

$out_dir=~s/\/$//gi;
unless (-e $out_dir) {
    `mkdir $out_dir`;
}

###split same file
 my $sam_file_name;
if ($th_p>1) {
    use threads; 
    use Config;
    my $num=`ls -l $sam_file| awk '{ print \$5}'`;
    chomp $num;
    my $number=(split /\s+/,$num)[0];
    my $size=int($number/$th_p/1024/1024)+1;
    $size=$size."m";
    $sam_file_name=(split /\//,$sam_file)[-1];
    $sam_file_name=$out_dir."/".$sam_file_name;
    `split -C $size $sam_file -d -a 3 $sam_file_name`;
}

##no ref
sub read_count
{
    my %junction_count;
    my @canshu=@_;
    open(IN,"$canshu[0]") or die "can't open $canshu[0]";
    while (<IN>){
        chomp(my $line=$_);
        
        if($line=~m/^\@/){
            next
        }
        unless(($line=~m/NH:i:1\t/)||($line=~m/NH:i:1$/)){
            next
        }
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
                if (($read_e2_length>=$min_overlap)&&($read_e1_length>=$min_overlap)) {
                    my $junction=$ch."_".$read_exon1_end."_".$read_exon2_start;
                    $junction_count{$junction}++;
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
                    if (($read_e2_length>=$min_overlap)&&($read_e1_length>=$min_overlap)) {
                        my $junction=$ch."_".$read_exon1_end."_".$read_exon2_start;
                        $junction_count{$junction}++;
                    }
                }
            }
        }
    }
    open OUT,">$canshu[0]_junction_count_tmp.txt";
   
    foreach my $junction(keys %junction_count){ 
        my $count=$junction_count{$junction};
        print OUT "$junction\t$count\n";

    }
    close OUT;
}
# Simultaneous computation
my %junction_count_all;
if ($th_p>1) {
    my @thread_all=qw /1...$th_p/;
    for(my $i=0;$i<$th_p;$i++ ){
        my $num;
        if ($i<10) {
            $num="00".$i;
        }else{
            $num="0".$i
        }
        my $sam="$sam_file_name$num";
        $thread_all[$i] = threads->create( \&read_count,$sam);
    }    
    for(my $i=0;$i<$th_p;$i++ ){
        $thread_all[$i] ->join();
    }
    opendir (DIR, $out_dir) or die "can't open the directory!";
    my @dir = readdir DIR;
    foreach my $file (@dir){
        my $query=(split /\//, $sam_file_name)[-1];
        if ($file=~m/$query[0-9]+_junction_count_tmp\.txt/gi) {
            open(IN,"$out_dir/$file") or die "";
            while (<IN>) {
                chomp (my $line=$_);
                my ($junction,$count)=split /\t/,$line;
                $junction_count_all{$junction}+=$count;
            }
            close IN;
        }
    }
    my $file=$out_dir."/"."junction_count\.txt";
    open OUT, ">$file";
    foreach my $junction(keys %junction_count_all){
        my $count=$junction_count_all{$junction};
        print OUT "$junction\t$count\n";
    }
    close OUT;
    for (my $i=0;$i<$th_p;$i++){
        my $j;
        if ($i<10) {
            $j="00".$i;
        }elsif($i<100){
            $j="0".$i;
        }
        my $file="$sam_file_name$j";
        my $file2="$sam_file_name$j"."_junction_count_tmp.txt";
        `rm $file`;
        `rm $file2`;
    }
}else{
    
    read_count($sam_file);
    my $file_temp=$sam_file."_junction_count_tmp.txt";
    my $file=$out_dir."/"."junction_count.txt";
    `mv $file_temp $file`;
    
}
my $end = time();
printf("Execution Time: %0.02f s\n", $end - $start);

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
