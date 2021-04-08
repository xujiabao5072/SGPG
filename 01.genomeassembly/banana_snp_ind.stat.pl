use strict;
my $sam=shift;
open OUT,">/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S04.reseq/Analysis/S10.Stat/S02.snp/S01.Data/$sam.snp.xls" or die $!;
print "/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S04.reseq/Analysis/S10.Stat/S02.snp/S01.Data/$sam.snp.xls\n";
my (%chrlist,%chrA,%chrB);
my $chrfile="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S00.DB/A+B/Banana_A_B.fa.len";
print "$chrfile\n";
open CF,$chrfile or die $!;
while (<CF>) {
	next if($_=~/^\#/);
	chomp;
	my @a=split/\s+/;
	my $chr=$a[0];
	my $len=$a[1];
	next if($chr=~/mito/ or $chr=~/random/ or $chr=~/scaffold/);
	if($chr=~/^B/){
		$chrlist{$chr}=$len;
		$chrA{$chr}+=$len;
	}else{
		$chrlist{$chr}=$len;
		$chrB{$chr}+=$len;
	}
	print "$chr\n";
}
close CF;
my %snps;
my %snpall;
foreach my $chr(sort keys %chrlist){
	my $tag="A";
	if(exists $chrB{$chr}){
		$tag="B";
	}
	my ($chrhom,$chrhet,$chrsnp);
	my $snp="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S04.reseq/Analysis/S04.Pricall/$sam/gvcf/$chr.$sam.vcf.gz";
	print "$snp\n";
	open SNP,"zcat $snp|" or die $!;
	while (<SNP>) {
		chomp;
		my @a=split/\s+/;
		my $snp=$a[1];
		my $ref=$a[3];
		my $alt=$a[4];
		my $qual=$a[5];
		my @feat=split/\:/,$a[9];
		next if ($qual<20);
		next if (length($ref)>1);
		my $flag=0;
		my @alt=split/\,/,$alt;
		foreach my $altnt (@alt) {
			if(length($altnt)>1){$flag=1};
		}
		next if($flag==1);
		my @nts=split/\//,$feat[0];
		my %alle;
		foreach my $nt (@nts) {
			$alle{$nt}=0;
		}
		my $allnum=scalar(keys %alle);
		if($allnum==1){
			$snps{$tag}{hom}++;
			$snps{$tag}{snp}++;
			$snpall{hom}++;
			$snpall{snp}++;
			$chrhom++;
			$chrsnp++;
		}else{
			$snps{$tag}{het}++;
			$snps{$tag}{snp}++;
			$snpall{het}++;
			$snpall{snp}++;
			$chrhet++;
			$chrsnp++;
		}
	}
	close SNP;
	print "Chr\t$sam\t$tag\t$chr\t$chrhom\t$chrhet\t$chrsnp\n";
	print OUT "Chr\t$sam\t$tag\t$chr\t$chrhom\t$chrhet\t$chrsnp\n";
}
print OUT "Summary\t$sam\t$snps{A}{hom}\t$snps{A}{het}\t$snps{A}{snp}\t$snps{B}{hom}\t$snps{B}{het}\t$snps{B}{snp}\t$snpall{hom}\t$snpall{het}\t$snpall{snp}\n";
close OUT;

