#!/usr/bin/perl -w
	use strict;

    use Cwd qw(getcwd);

	my ($input,$denominator,$minn,$mint) = @ARGV;

	my $dPSIcutoff = 20;
	my $pvaluecutoff = 0.01;
	my $mediangap = 1;
	#my $minn = 2;
	print "Minimal control samples = $minn\n";
	print "Minimal tumor samples = $mint\n";
	
	my %input;
	my %gn;
	my $linecount = 0;
	print "Reading $input...\n";
	open(FILE,"$input") || die "Aborting.. Can't open $input : $!\n";
    while(my $line=<FILE>){
		chomp $line;
		$linecount++;
		next if($line=~/^\#/ || $linecount == 1);
		my ($region,$gn,$exon,$type,$n,$t,$exontype,$refT,$dPSI,$p,$fdr,$nvalues,$tvalues,$dbid) = split(/\t/,$line);
		next if($gn=~/MSTRG\./);
		$gn{$dbid} = $gn;
	}
	close(FILE);
	
	my %globaldenominator;
	my %denominator;
	$linecount = 0;
	print "Reading $denominator...\n";
	open(FILE,"$denominator") || die "Aborting.. Can't open $denominator : $!\n";
    while(my $line=<FILE>){
		chomp $line;
		$linecount++;
		next if($line eq "" || $linecount < 3);
		my @columns;
		if($linecount == 3){
			@columns = split(/\t/,$line);
			next;
		}
		my @rows = split(/\t/,$line);
		my ($ID,$exonsize) = ($rows[0],$rows[1]);
		my ($gn,$dbid) = split(/\:/,$ID);
		next if(!exists $gn{$dbid});
		#15_43640497_43647293_W_ENST00000321596_1
		my ($chr,$start,$end,$type,$ENST,$num) = split(/\_/,$dbid);
		for(my $i = 2;$i < scalar @rows;$i++){
			my $value = $rows[$i];
			next if($value=~/NaN/);
			$globaldenominator{$gn{$dbid}}{$value}++;
			$denominator{$dbid} = $value;
		}
	}
	close(FILE);
	
	my %firstgap;
	foreach my $gn(sort keys %globaldenominator){
		my @values;
		foreach my $value(sort keys %{ $globaldenominator{$gn} }){
			push(@values,$value);
		}
		my $median = median(@values);
		my $first_significant_gap_start = find_first_significant_gap(@values);
		if(defined $first_significant_gap_start){
			if($mediangap == 1){
				$firstgap{$gn} = $first_significant_gap_start if($first_significant_gap_start < $median);
			}else{
				$firstgap{$gn} = $first_significant_gap_start;
			}
			
		}else{
			next;
		}
	}
	
	my $outfn = $input;
	$outfn=~s/.txt/.volcano.txt/;
	my $barfn = $input;
	$barfn=~s/.txt/.barchart.txt/;
	
	$linecount = 0;
	my %barchart;
	my @barchart;
	push(@barchart, "Exon Inclusion");
	push(@barchart, "Exon Skipping");
	push(@barchart, "Alt. 5'-splice-site");
	push(@barchart, "Alt. 3'-splice-site");
	push(@barchart, "Increased IR");
	push(@barchart, "Decreased IR");
	
	
	print "Reading $input...\n";
	open(FILE,"$input") || die "Aborting.. Can't open $input : $!\n";
	open(OUT,">$outfn") || die "Aborting.. Can't open $outfn : $!\n";
    while(my $line=<FILE>){
		chomp $line;
		$linecount++;
		print OUT "Event Region	Gene Symbol	Target Exon	Event Type	N	T	dPSI	log10(p-value)	Database ID\n" if($linecount == 1);
		next if($line eq "" || $linecount == 1);
		my ($region,$gn,$exon,$type,$n,$t,$exontype,$refT,$dPSI,$p,$fdr,$nvalues,$tvalues,$dbid) = split(/\t/,$line);
		next if($gn=~/MSTRG\./ || $dPSI == 0);
		next if($n < $minn || $t < $mint);
		my $denominator = $denominator{$dbid};
		my $log10p = sprintf("%.4f",log($p)/log(10)) * -1;
		if($type=~/SES/ || $type=~/MXS/ || $type=~/MES/){
			$type = "Exon Inclusion" if($dPSI > 0);
			$type = "Exon Skipping" if($dPSI < 0);
		}
		if($type=~/A3SS/ || $type=~/A5SS/){
			$type = "Alt. 3'-splice-site" if($type=~/A3SS/);
			$type = "Alt. 5'-splice-site" if($type=~/A5SS/);
		}
		if($type=~/IR/){
			$type = "Increased IR" if($dPSI > 0);
			$type = "Decreased IR" if($dPSI < 0);
		}
		if(exists $firstgap{$gn}){
			if($denominator <= $firstgap{$gn}){
				next;
			}else{
				print OUT "$region\t$gn\t$exon\t$type\t$n\t$t\t$dPSI\t$log10p\t$dbid\n";
			}
		}else{
			print OUT "$region\t$gn\t$exon\t$type\t$n\t$t\t$dPSI\t$log10p\t$dbid\n";
		}
		if(abs($dPSI) > $dPSIcutoff && $p < $pvaluecutoff){
			$barchart{$type}{$gn}++;
		}
	}
	close(OUT);
	close(FILE);
	
	open(OUT,">$barfn") || die "Aborting.. Can't open $barfn : $!\n";
	print OUT "File Name";
	foreach my $type(@barchart){
		print OUT "\t$type";
	}
	print OUT "\n";
	print OUT "$input";
	foreach my $type(@barchart){
		my $num = scalar keys %{ $barchart{$type} };
		print OUT "\t" . $num;
	}
	print OUT "\n";
	close(OUT);
	



# Function to find the starting point of the first significant gap
sub find_first_significant_gap {
    my @numbers = @_;

    # Sort the numbers if they are not sorted
    @numbers = sort {$a <=> $b} @numbers;

    my @gaps;
    my $total_gap = 0;

    # Calculate all gaps and the total gap size
    for(my $i = 0; $i < scalar @numbers - 1; $i++) {
        my $gap = $numbers[$i+1] - $numbers[$i];
        push @gaps, {gap => $gap, start => $numbers[$i], end => $numbers[$i+1]};
        $total_gap += $gap;
    }

	return undef if(scalar @gaps == 0);
    # Calculate the average gap size
    my $average_gap = $total_gap / scalar @gaps;

    # Find the first gap larger than the average gap size
    foreach my $gap_info (@gaps) {
        if($gap_info->{gap} > $average_gap) {
            return $gap_info->{start}; # return the starting point of the first significant gap
        }
    }

    # Return undef if no significant gap is found
    return undef;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub averagex{
        my @data = @_;
        if (not @data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@data) {
                $total += $_;
        }
        my $average = $total / @data;
        return $average;
}

sub stdev{
        my @data = @_;
        if(@data == 1){
    		return 0;
        }
        my $average = averagex(@data);
        my $sqtotal = 0;
        foreach(@data) {
            $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@data-1)) ** 0.5;
        return $std;
}

