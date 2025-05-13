#!/usr/bin/perl -w
use strict;
#my $outfile = $ARGV[0];

my $mindist = 5;
my $mindistSmall = 3;

my @readQueue;
my $notEnd = 1;
my $stdist = 0;
#open OUT, ">$outfile";
while($notEnd) {
#    $stdist = 0; # Uhm I think.
    while($notEnd && (scalar @readQueue < 2 || $stdist < $mindist)) {
#	print "scRQ2: ", scalar @readQueue, "\n";
	$notEnd = <>;
	if($notEnd) {
	    chomp $notEnd;
	    #    my ($chr, $start, $end, $name, $qual) = split;
	    my @line = split "\t", $notEnd;
	    push @readQueue, [@line];
	    if(scalar @readQueue > 1) {
		#	    print @{$readQueue[0]};
#		print $#readQueue, "\n";
		$stdist = abs($readQueue[$#readQueue][1] - $readQueue[0][1]);
#		print "stdist:, ", $stdist, "\n";
	    }
	}
    }
#    print "RQ3: ", scalar @readQueue, "\n";
    if(scalar @readQueue == 1) { # Only happens last iteration
	next;
    }
    my $md = $readQueue[1][1] - $readQueue[0][1];
#    print "Min dist ", $md, "\n";
#    if(scalar @readQueue == 1 || $stdist > $mindist) { # No dup candidates
    if($md > $mindist) { # No dup candidates
	print join("\t", @{$readQueue[0]}), "\n";
	shift @readQueue;
    } else {
#	print "Maybe dup\n";
	print join("\t", @{$readQueue[0]}), "\n";
	my @newarr = ();
	for(my $i = 1; $i != @readQueue; ++$i) {
	    # Check if any of the reads are dups of the first one
	    my $dist1 = abs($readQueue[0][1] - $readQueue[$i][1]);
	    my $dist2 = abs($readQueue[0][2] - $readQueue[$i][2]);
#	    print "dists", $dist1, ", ", $dist2, "\n";
	    if(($dist1 < $mindist && $dist2 < $mindistSmall ) || ($dist2 < $mindist && $dist1 < $mindistSmall)) {
		print STDERR "Found dup! ", $dist1, " - ", $dist2, join ("\t", @{$readQueue[0]}), " --- ", join ("\t", @{$readQueue[$i]}), "\n";
	    } else {
		push @newarr, [@{$readQueue[$i]}];
#		print "No dup ", join ("\t", @{$readQueue[0]}), "\n-------", join ("\t", @{$newarr[$#newarr]}), "\n";
	    }
	}
	@readQueue = @newarr;
	if(scalar @readQueue > 1) {
	    $stdist = abs($readQueue[$#readQueue][1] - $readQueue[0][1]);
	}
#	print "Not end: ", $notEnd;
    }
}

while(scalar @readQueue > 1) {
    print join("\t", @{$readQueue[0]}), "\n";
    my @newarr = ();
    for(my $i = 1; $i != scalar @readQueue; ++$i) {
	# Check if any of the reads are dups of the first one
	my $dist1 = abs($readQueue[0][1] - $readQueue[$i][1]);
	my $dist2 = abs($readQueue[0][2] - $readQueue[$i][2]);
	if($dist1 > $mindist || $dist2 > $mindist) {
	    push @newarr, [@{$readQueue[$i]}];
	}
    }
    @readQueue = @newarr;
}

print join("\t", @{$readQueue[0]}), "\n";

