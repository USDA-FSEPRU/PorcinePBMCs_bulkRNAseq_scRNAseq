#!/bin/env perl

############# USAGE ############
#
# perl trimming.AT.pl RNA-seq
#
# Creator : Haibo Liu, ISU
###############################
use warnings;
use strict;

my $start=0;
my $len =0;
my @lines =();

if ( scalar @ARGV != 1 )
{
   die "USAGE: Provide sequencing type as a option! $!\n";
}

my $isRNA_seq = $ARGV[0];


while (<STDIN>)
{
   chomp;
   if ( $. % 4 == 2)
   {
      $len = length($_);
      if (/A{10,}[^A]?A*$/ || /^T*[^T]?T{10,}/)
      {
         #my $full_len = length($_);
         if ( $isRNA_seq eq "RNA-seq" )
           {
               s/A+[^A]?A*$//;
               s/^T*[^T]?T{10,}//;
           }
        while ($_ =~ m/A{10,}[^A]$/ || $_ =~ m/^[^T]T{10,}/) 
        {
           if ( $isRNA_seq eq "RNA-seq" )
           {
               s/A+[^A]?A*$//;
               s/^T*[^T]?T{10,}//;
           }
        }
 
         $len = length($_);
         #my $num_A = $full_len - $len;    
         #$lines[0] = ${lines[0]}.":${num_A}A";
      }
      push @lines, $_;
   }
   elsif ($start && $. % 4 == 0)
   {
      if ($len >= 25)
      {
        my $temp = substr($_, 0, $len);
        push @lines, $temp;
        print join "\n", @lines;
        print "\n"; 
      }        
      @lines = (); 
      $len = 0;
      $start = 0;
   }
   else 
   {
      if ($. % 4 == 1)
      {
         $start = 1;
         s/([^\s]+\s+[^\s]+).+/$1/; 
      }
      push @lines, $_;
   }
}
