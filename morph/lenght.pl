#!/usr/bin/perl

open(fh, "DeaneHV.txt") or die "Could not read file";

while($line = <fh>) {
    chomp;
    @l = split(' ', $line);
    $len=$l[2]-$l[1]+1;
    my $pdb=substr($l[0],6,6);
    $list{$pdb}=$len; 
    #print "$pdb $l[1] $l[2] $len\n";
}

close(fh);

open(fh, "stat_aF0A0S.txt") or die "Could not read file";

while($line = <fh>) {
    chomp;
    @l = split(' ', $line);
    my $pdb=substr($l[0],6,6);
    print "$pdb $list{$pdb}\n";
}

close(fh);

