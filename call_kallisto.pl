use strict;
use warnings;
use Data::Dumper;

my $index = '/lustre/scratch118/infgen/team133/fr7/trichuris_projects/Mus_musculus.GRCm38.cdna.E98.all.idx';

my %samples;
while(<>){
	my @temp = split('\s+',$_);
	$temp[1] =~ s/_[12].fastq.gz//;
	$samples{$temp[0]}{$temp[1]} = '';
}

foreach my $sample (keys %samples){
	mkdir $sample;
	my $paired_lane;
	foreach my $lane (keys(%{$samples{$sample}})){
		$paired_lane .= $lane."_1.fastq.gz ".$lane."_2.fastq.gz ";		
	}
	my $kallisto = "bsub -o ".$sample."/quant.o -e ".$sample."/quant.e -R 'select[mem>=5000] rusage[mem=5000] span[hosts=1]' -M 5000 kallisto quant -i ".$index." -o ".$sample." -b 100 ".$paired_lane;
	system($kallisto);
}

#print Dumper \%samples;
