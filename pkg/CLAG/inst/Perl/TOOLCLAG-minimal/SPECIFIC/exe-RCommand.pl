#!/usr/bin/perl -s
use File::Find;

$path= $f;
$parameter= $k;
$delta= $d;



$exp="input.txt";
$SUBpath="";
sub recherche {
  if (-f){
    push(@acronyme,"$File::Find::dir/$_")
    if $_=~/$exp$/;
  }
}

@acronyme=();
$SUBpath="";
find(\&recherche,"$path");
foreach $acro (@acronyme)
{
	
	
	$SUBpath=substr($acro, 0,((length $acro)-10));
	print "\n\nJOB SUMMARY\n";
    print "Folder: $SUBpath\n";
	print "-----------------------------------\n\n";
	
	print "DEBUG: STEP1: Cluster";	
	system "./SPECIFIC/1-CoEvScoring.pl $SUBpath $delta";
	
	print "OK\n";

	print "DEBUG: STEP2: Cluster aggregate...";	
	system "./SPECIFIC/Aggregation.pl $SUBpath aggregation $parameter $delta";
	print "OK\n";
	
#	print "DEBUG: STEP3: CLAG displays matrix...";		
#	system "./SPECIFIC/2-MatrixNONClusterized.pl $SUBpath ";
#	$RCommandFile=$SUBpath."/"."R_COMMAND-Matrix.txt";
#	system "R --no-save < $RCommandFile > $SUBpath/DEBUG1.txt";
#	print "OK\n";
#	
#	print "DEBUG: STEP4: CLAG displays matrix Clusterized...";		
#	system "./SPECIFIC/3-MatrixClusterized.pl $SUBpath  $delta";
#	$RCommandFile=$SUBpath."/"."R_COMMAND-Matrix-Clusterized.txt";
#	system "R --no-save < $RCommandFile > $SUBpath/DEBUG2.txt";
#	print "OK\n";
#	
#	print "DEBUG: STEP5: CLAG displays matrix clusterized sym= 1, env=1...";		
#	system "./SPECIFIC/4-MatrixClusterizedParam.pl $SUBpath $delta ";
#	$RCommandFile=$SUBpath."/"."R_COMMAND-Matrix-Clusterized-Scores1.txt";
#	system "R --no-save < $RCommandFile > $SUBpath/DEBUG3.txt";
#	print "OK\n";
#	
#	
#	print "DEBUG: STEP6: Draw graph...";
#	system "./SPECIFIC/5-CallNeato.pl $SUBpath $parameter $delta";
#	print "OK\n";
#	
#
#	print "DEBUG: STEP7: CLAG displays matrix Clusterized -aggregated ...";		
#	system "./SPECIFIC/3-MatrixClusterizedAggregation.pl $SUBpath aggregation $delta";
#	$RCommandFile=$SUBpath."/"."R_COMMAND-Matrix-Aggregated.txt";
#	system "R --no-save < $RCommandFile > $SUBpath/DEBUG4.txt";
#	print "OK\n";
	
} 

exit 0;
	
	
		
