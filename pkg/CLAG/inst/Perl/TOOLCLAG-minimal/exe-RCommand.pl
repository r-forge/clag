#!/usr/bin/perl -s
 
$path= $f;
$program= $p;
if($k){
$parameter= $k;
}else{
$parameter= 0;
}
if($d){
$delta= $d;
}else{
$delta= "5/10/20/40";
}

my $res = 1;

if($program eq 1){
	print "GENERAL matrix\n";
#	system "./AnyMatrix/exe-RCommand.pl -f=$path -k=$parameter -d=$delta";
	$res = system($^X, "./AnyMatrix/exe-RCommand.pl", "-f=$path", "-k=$parameter", "-d=$delta");
}else{
	if($program eq 2){
		print "BINARY matrix\n";	
		$res = system "./Binary/exe-RCommand.pl -f=$path -k=$parameter";	
	}else{
		if($program eq 3){
			print "SPECIFIC matrix where N=M\n";
			$res = system "./SPECIFIC/exe-RCommand.pl -f=$path -k=$parameter -d=$delta";
			
		}else{
			print "ERROR";
			exit 1;
		}
	}
}
exit $res;
	

