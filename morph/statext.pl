#!/usr/bin/perl
#!/usr/local/bin/perl

#--------------------------------------------------------------
# Jussi Karlgren, jussi@sics.se
#--------------------------------------------------------------
# Modified by Mon:
# This computes some basic statistics from table text files.
#

unless(@ARGV == 6){
    print "USAGE:\n\t$0 <data> [icol] [fcol] [irow] [frow] [separator]\n";
    print "\t<data> --> Input text table.\n";
    print "\t[icol] --> Initial column\n";
    print "\t[fcol] --> Final column\n";
    print "\t[irow] --> Initial row\n";
    print "\t[frow] --> Final row. \"all\" uses all rows available.\n";
    print "\t[separator] --> Number of white spaces separating fields.\n";
    die "NOTE:\n\tAveraged range is [irow,frow] (only data displayed will be used)\n";
}

$debug = 0;
$icol=$ARGV[1]-1;
$fcol=$ARGV[2]-1;
$irow=$ARGV[3]-1;
$frow=$ARGV[4];
$sep=$ARGV[5];
$defdecpos=0; # Default decimal positions when an integer is detected

# Change this to fit whatever set of columns 
# you are studying. First column 
# is 0. Many columns ( > 2) can be specified.
@columns = ($icol .. $fcol);
my @sizecol; # This stores the maximum column size.
foreach(0 .. $fcol) { push(@sizecol,0); } # Column size array initialization
my @sizedec; # This stores the floating point decimal positions.
foreach(0 .. $fcol) { push(@sizedec,0); } # Decimal positions size array initialization


# Read data.
# $line_number is the current line number; $column the current column.
# %score will hold the value for any combination of key and column.
# %individuals will hold a "1" for every line with a non-null value. If a line
# holds some nulls but at least one non-null value the nulls will be interpreted
# as zero values.
print "READING:\n" if $debug;
$line_number = 0;
my $clen=0;
my $decpos=0;
open (FILE,"$ARGV[0]");
while (<FILE>) 
{
  chomp;
  my @line = split(/ +/, $_);
 if( $_ !~ "^#" )  # Not reading the "#" begining lines
 {
    print "$_" if $debug;
    my $ncol = 0;
    foreach my $column (@columns) 
    { 
        $clen = length($line[$column]);
	$clen -= 7 + length("$N") if $ncol == 0;
	print "$line[$column] clen= $clen\n" if $debug;
	if($clen > $sizecol[$column])
	{
	  $sizecol[$column] = $clen; # stores the maximum column size
	  $decpos=index($line[$column], '.');
	  if($decpos == -1) 
	  { $sizedec[$column] = $defdecpos; }
	  else
	  { $sizedec[$column] = $clen-$decpos-1; }
	}
	
	if ($line[$column] ne "")
	{
	    $data{$column,$line_number} = $line[$column];
#	    $individuals{$line_number} = 1;
	};

        $ncol++; # update column index
    };
    $line_number++;
 }
};
# print "sizedec= @sizedec\n";

close FILE;
print "Lines readed $line_number\n" if $debug;
foreach(@sizecol){ $_ += $sep; } # Column separator size
print "Column sizes: @sizecol\n" if $debug;
$line_number--; # this is the last available element

if ($frow eq 'all') { $frow = $line_number; }
#$frow=$frow-1;
#else { $frow = $frow - 1; }

if ($frow > $line_number)
{
  print "More rows requested ($frow) than available ($line_number), forcing maximum.";
  $frow = $line_number;
}

my $N=$frow-$irow+1;
print "irow= $irow  frow= $frow  N= $N\n" if $debug;

# Draw a separator line
print "#";
foreach(@columns)
{
  foreach(1 .. $sizecol[$_])
  {
    print "-";
  }
}
print "\n";

# Showing only the lines taken into account.
#foreach $index ($irow .. $line_number)
my $format;
foreach $index ($irow .. $frow)
{
	printf "%4d ",$index+1 if $debug; # line number begining with 1
#	print "format= $format\n";
#	exit;
#	foreach $col (@columns) { printf "%10s ",$data{$col,$index}; }
	foreach $col (@columns) 
	{ 
	  $format="%$sizecol[$col]s";
	  printf $format,$data{$col,$index}; 
	}
	printf "\n";
}

# Draw a separator line
print "#";
foreach(@columns)
{
  foreach(1 .. $sizecol[$_])
  {
    print "-";
  }
}

$sizecol[1] += $sizecol[0] - 4; # correct offset upon first col removal

my $ncol = 0;
printf("\n#Avg(%d):",$N);
# Computing Averages...
foreach my $col (@columns)
{
	$ncol++; # update column index
	my $avg = 0.0;
	foreach my $row ($irow .. $frow)
	{ $avg += $data{$col,$row}; }
#	$avg /= $frow - $irow +1;
	$avg /= $N;
	push(@avgs,$avg);
	next if($ncol == 1); # skip first column (typically text)
        $format="%$sizecol[$col].$sizedec[$col]f";
	printf $format,$avg;
}
printf("\n");
#printf(" :Averages\n");
#print "@avgs\n";
my $value;

printf("#Sig(%d):",$N);
$ncol = 0;
# Computing sigmas...
foreach my $col (@columns)
{
	$ncol++; # update column index
	my $sig=0.0;
	foreach my $row ($irow .. $frow)
	{ $sig += ($data{$col,$row} - $avgs[$col-$icol])**2; } #printf "\ndata %f\tavgs %f\n",$data{$col,$row},$avgs[$col-$icol];}
	$sig /= $N; 
	$sig = sqrt($sig);
	#push(@sigs,$sig);
	next if($ncol == 1); # skip first column (typically text)
        $format="%$sizecol[$col].$sizedec[$col]f";
	printf $format,$sig; 
}
printf("\n");
#printf(" :Sigmas\n");
#print "@sigs\n";

printf("#Med(%d):",$N);
$ncol = 0;
# Computing median...
foreach my $col (@columns)
{
	$ncol++; # update column index

	my $med;
	my @coldata=(); # New clean array 
	foreach my $row ($irow .. $frow)
	{ 
		push(@coldata,$data{$col,$row});
	}
	my @sorted = sort {$a <=> $b} @coldata;  # sort data numerically ascending
	# print "Sorted: @sorted\n";
	if($N % 2 == 1) # Odd
	{
		$med = $sorted[$N/2]; # middle element is the Median
	}
	else
	{
		$med = ( $sorted[$N/2] + $sorted[$N/2-1] ) / 2 ; # the average of the two elements in the middle is the Median
	}
	push(@meds,$med);
	next if($ncol == 1); # skip first column (typically text)
        $format="%$sizecol[$col].$sizedec[$col]f";
	printf $format,$med; 
}
printf("\n");
#printf(" :Medians\n");
#print "@meds\n";


