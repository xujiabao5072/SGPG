use strict;
my $fqlist="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S04.reseq/Lib/fq.config.list";
my %samlist;
open FQ,$fqlist or die $!;
while (<FQ>) {
	chomp;
	my @a=split/\s+/;
	my $sam=$a[0];
	$samlist{$sam}=0;
}
close FQ;
my $shelldir="/hwfssz4/BC_COM_P5/F19HTSSCWLJ4898/li201911011155549842RE0101/S04.reseq/Analysis/S10.Stat/S02.snp/S00.shell/step01.snp";
foreach my $sam (sort keys %samlist) {
	open SH,">$shelldir/step01.$sam.sh" or die $!;
	print SH "perl $shelldir/../banana_snp_ind.stat.pl $sam\n";
	close SH;
}