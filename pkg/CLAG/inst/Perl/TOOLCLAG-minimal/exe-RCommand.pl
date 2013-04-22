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



if($program eq 1){
	print "GENERAL matrix\n";
	system "./AnyMatrix/exe-RCommand.pl -f=$path -k=$parameter -d=$delta";	
}else{
	if($program eq 2){
		print "BINARY matrix\n";	
		system "./Binary/exe-RCommand.pl -f=$path -k=$parameter";	
	}else{
		if($program eq 3){
			print "SPECIFIC matrix where N=M\n";
			system "./SPECIFIC/exe-RCommand.pl -f=$path -k=$parameter -d=$delta";
			
		}else{
			print "ERROR";
		}
	}
}
exit 0;
	

