#!/usr/bin/perl
use strict;
use Cwd qw();

# usage 

my $file = $ARGV[0]; # Input text file with data


my %pdbs=();


# read forward
open (MYFILE, $file);
my $if=-1;

while(<MYFILE>) {
  if ($_ =~ /MODEL/) {
    $if++;
  }
  else {
    push (@{$pdbs{$if}},$_); 
    }
 }
close(MYFILE);


# print in reverse 
for(my $i = 0; $i <= $if; $i++) {
my $size=scalar @{$pdbs{$i}};
    printf "MODEL %-7d\n", $i+1;
for(my $j = 0; $j < $size; $j++) {
    print "$pdbs{$if-$i}[$j]";
}
}

my $file = $ARGV[1]; # Input text file with data

open (MYFILE, $file);

%pdbs=();

while(<MYFILE>) {
  if ($_ =~ /MODEL/) {
    $if++;
    printf "MODEL %-7d\n", $if+1;
  }
  else {
    printf "$_";
    }
 }
close(MYFILE);







