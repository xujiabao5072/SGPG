use  strict;
my $chrin=shift;
my $st=shift;
my $en=shift;
my $out="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S07.Family/Reseq50/S01.Shell2/05.merge/data/$chr/$chr.$st.$en.merge";

my %samlist;
my $file="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S07.Family/Popbin/First.group.txt";
open FL,$file or die $!;
while (<FL>) {
        chomp;
        my @a=split/\s+/;
        my $sam=$a[0];
        $samlist{$sam}=0;
}
close FL;
my %snps;
my %refs;
foreach my $sam (sort keys %samlist) {
	my $vcf="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S07.Family/Reseq50/S04.Pricall/$sam/gvcf/$sam.gt.gz";
	open VCF,"zcat $vcf|" or die $!;
	<VCF>;
	while (<VCF>) {
		chomp;
		next if($_=~/\#/);
		my @a=split/\s+/;
		my $chr=$a[0];
		my $pos=$a[1];
		my $ref=$a[2];
		my $gt=$a[3];
		if($chr ne $chrin){next};
		if($pos<$st){next};
		if($pos>=$en){last};
		$snps{$chr}{$pos}{$sam}=$gt;
		$refs{$chr}{$pos}=$ref;
	}
	close VCF;
}
open OUT,">$out" or die $!;
my @sams=sort keys %samlist;
my $samstr=join "\t",@sams;
print OUT "Chr\tPos\tRef\t$samstr\n";
foreach my $chr (sort keys %snps) {
	foreach my $pos (sort {$a<=>$b} keys %{$snps{$chr}}) {
		my @vals;
		my $miss++;
		my %freq;
		foreach my $sam (@sams) {
			push @vals,$snps{$chr}{$pos}{$sam};
			if($snps{$chr}{$pos}{$sam} eq "-/-"){
				$miss++;
			}else{
				my @alle=split/\//,$snps{$chr}{$pos}{$sam};
				$freq{$alle[0]}++;
				$freq{$alle[1]}++;
			}
		}
		my $mr=$miss/scalar(@sams);
		my @alle=sort {$freq{$b}<=>$freq{$a}} keys %freq;
		if(scalar(@alle)!=2){next};
		my $maf=$freq{$alle[1]}/($freq{$alle[0]}+$freq{$alle[1]});
		next if($maf<0.05);
		my $valstr=join "\t",@vals;
		print OUT "$chr\t$pos\t$ref\t$valstr\n";
	}
}
close OUT;



