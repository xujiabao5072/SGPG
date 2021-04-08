
my %chrlen;
open CL,"/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S07.Family/Reseq50/S02.DB/Banana_A_B.format.fa.fai" or die $!;
while (<CL>) {
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $len=$a[1];
	$chrlen{$chr}=$len;
}
close CL;

foreach my $chr (sort keys %chrlen) {
	open CHR,">step08.$chr.merge.sh" or die !;
	my $len=$chrlen{$chr};
	for (my $loc=1;$loc<$len ;$loc+=1000000) {
		my $st=$loc;
		my $en=$st+1000000;
		print CHR "perl /hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S07.Family/Reseq50/S01.Shell2/05.merge/samtools_gt_merge.pl $chr $st $en\n";
	}
	close CHR;
}



