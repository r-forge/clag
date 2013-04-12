#!/usr/bin/perl -w
use POSIX;
$folder=shift;


open (IN,"$folder/inputOriginal.txt") || die "ERROR: CEBA couldn't write into  the file $folder/inputOriginal.txt";
while($line =<IN>){
	chomp $line;
	@dataline=split(" ",$line);
	$hash{$dataline[0]}{$dataline[1]}=$dataline[2];
	push(@valuesData,$dataline[2]);
}
close(IN);

@valuesDataSort=sort{$a<=>$b} @valuesData;
$valuemin=$valuesDataSort[0];
$sizeArray=(scalar @valuesDataSort);
$valuemax=$valuesDataSort[$sizeArray-1];

open (OUT,">$folder/input.txt") || die "ERROR: CEBA couldn't write into  the file $folder/input.txt";
@keyE1=keys %hash;
foreach $e1(@keyE1){
	@keyE2=keys %{$hash{$e1}};
	foreach $e2(@keyE2){
			$maxNew=1;
			$minNew=0;
			$valueN=($hash{$e1}{$e2}-$valuemin)*($maxNew-$minNew)/($valuemax-$valuemin)+$minNew;	
			print  OUT $e1." ".$e2." ".$valueN."\n";
		#	print  $hash{$e1}{$e2}." ".$valueN."\n";
	}	
}	

close (OUT);


exit;


