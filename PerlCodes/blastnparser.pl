#!/usr/bin/perl -w
#
# ====================================================================================================
# File: blastnparser.pl
# Author: Eddy
#
# This file takes the standard blastn output produced by blastn from NCBI and parses
# its contents to a tabular form for easier viewing and manipulation.
#
# Usage: perl blastnparser.pl <blastn output file>
# ====================================================================================================

use strict;
use warnings;

# Prints the usage error message of this program
sub printUsage() {
   printf("\n");
   printf("*** INCORRECT NUMBER OF ARGUMENTS OR FILE DOES NOT EXISTS OR NOT READABLE ***\n");
   printf("USAGE: perl blastnparser.pl <blastn output file>\n");
   printf("WHERE: <blastn output file>: Is the standard blastn output file specified with the '-o'\n");
   printf("                             option when running 'blastall -p blastn'\n");
   printf("\n");
} # End of sub printUsage()

# Checks if the file path passed as argument value corresponds to a plain text file with size
# non zero. Returns 0 if file is ok; otherwise return -1
sub isFileOk() {
   my ($filePath) = @_;
   my $returnCode = -1;
   
   if (-s $filePath && -T $filePath) {
      $returnCode = 0;
   } # End of if statement

   return $returnCode;
} # End of sub isFileOk()

# Removes all leading and trailing white spaces from the argument string
sub trim() {
   my ($myString) = @_;
   
   for ($myString) {
      s/^\s+//;
      s/\s+$//;
   } # End of for loop

   return $myString;
} # End of sub trim()

# Print the output headers
sub printHeaders() {
   printf("Query ID\tQuery\tQuery Length\tBlast Hit ID\tSegment #\tSubject\tSubject Length\tScore\tExpect\t");
   printf("Coverage\tIdentity\tIdent. %%\tGaps\tGaps \%%\tStrand\tQuery Start\tQuery End\tSubject Start\tSubject End\n");
} # End of sub printHeaders()

# -------------------- MAIN PROGRAM STARTS HERE ------------------- #
# Check the number of arguments passed to the program from the command file
if ($#ARGV + 1 == 1 && &isFileOk($ARGV[0]) == 0) {
   # Variables used for parsing the queries
   my $line = "";               # one line of the input
   my $queryID = 0;             # Query ID
   my $queryName = "";          # Name of the letters
   my $queryLength = 0;         # The length (bp) of the query
   my $blastHitID = 0;          # Blast hit id
   my $subjectName = "";        # The name of the subject hit
   my $subjectLength = 0;       # The length of the subject
   my $segmentID = 0;           # The segment ID count
   my $queryStart = 0;          # The query start postition
   my $queryEnd = 0;            # The query end position
   my $subjectStart = 0;        # The subject start position
   my $subjectEnd = 0;          # The subject end position
   
   # print the output headers
   &printHeaders();
  
   # open file for input and start processing
   open (INFILE, $ARGV[0]);

   while (<INFILE>) {
      $line = &trim($_);

      if (length($line) != 0) {
         if ($line =~ m/^Query= /i) {
            # Check if needs to print $queryStart, $queryEnd, $subjectStart, $subjectEnd
            if ($queryStart != 0 || $subjectStart != 0) {
               printf("%d\t%d\t%d\t%d\n", $queryStart, $queryEnd, $subjectStart, $subjectEnd);
            } # End of if statement
            
            # Reset counts and other variables
            $queryName = "";
            $queryLength = "";
            $blastHitID = 0;
            $subjectName = "";
            $subjectLength = 0;
            $segmentID = 0;
            $queryStart = 0;
            $queryEnd = 0;
            $subjectStart = 0;
            $subjectEnd = 0;

            # Increment queryID by 1
            $queryID++;
            $queryName = substr($line, 7);

            # continue reading each line from input file until (xxxx Letters)
            $line = <INFILE>;
            $line = &trim($line);
            while ($line !~ m/^\(.*letters\)$/i) {
               $queryName = join(" ", $queryName, $line);
               $line = <INFILE>;
               $line = &trim($line);
            } # End of while loop
           
            # parse the number of letters of query
            $queryLength = substr($line, 1, index($line, " "));
            $queryLength =~ s/,//g;
         } # End of if statement
        elsif ($line =~ /^\** No hits found \**$/) {
           printf("%d\t", $queryID);
           printf("%s\t", $queryName);
           printf("%d\t", $queryLength);
           printf("%s\n", $line);
        } # End of else if statement
        elsif ($line =~ /^>/) {
            # Check if needs to print $queryStart, $queryEnd, $subjectStart, $subjectEnd
            if ($queryStart != 0 || $subjectStart != 0) {
               printf("%d\t%d\t%d\t%d\n", $queryStart, $queryEnd, $subjectStart, $subjectEnd);
            } # End of if statement
         
            # Reset counts and other variables
            $segmentID = 0;
            $queryStart = 0;
            $queryEnd = 0;
            $subjectStart = 0;
            $subjectEnd = 0;
            
            $blastHitID++;
            $subjectName = substr($line, 1);

            $line = <INFILE>;
            $line = &trim($line);

            while ($line !~ m/^Length =/i) {
               $subjectName = join(" ", $subjectName, $line);
               $line = <INFILE>;
               $line = &trim($line);
            } # End of while loop

            # parse the length of the blast hit
            $subjectLength = substr($line, index($line, "=") + 1);
         } # End of else if statement
         elsif ($line =~ m/^Score =/i) {
            # Check if needs to print $queryStart, $queryEnd, $subjectStart, $subjectEnd
            if ($queryStart != 0 || $subjectStart != 0) { 
               printf("%d\t%d\t%d\t%d\n", $queryStart, $queryEnd, $subjectStart, $subjectEnd);
            } # End of if statement
                                                            
            # Reset counts and other variables
            $queryStart = 0;
            $queryEnd = 0;
            $subjectStart = 0;
            $subjectEnd = 0;

            # parsing score and expect values into @tmpArray
            my @tmpArray = split(/,/, $line);
          
            # Remove "Score =" and "Expect ="
            $tmpArray[0] =~ s/Score = //i;
            $tmpArray[1] =~ s/Expect =//i; 
            
            # trimming leading and trailing spaces
            $tmpArray[0] = &trim($tmpArray[0]);
            $tmpArray[1] = &trim($tmpArray[1]);
            $tmpArray[1] = uc($tmpArray[1]);   # Converting expected value to uppercase
            
            if ($tmpArray[1] =~ m/^E/) {
               $tmpArray[1] = "1.00" . $tmpArray[1];
            } # End of if statement
            
            # increment segmentID by 1 and print common fields
            $segmentID++;
            
            # printing common fields
            printf("%d\t", $queryID);
            printf("%s\t", $queryName);
            printf("%d\t", $queryLength);
            printf("%d\t", $blastHitID);
            printf("%d\t", $segmentID);
            printf("%s\t", $subjectName);
            printf("%d\t", $subjectLength);
            printf("%s\t", $tmpArray[0]);
            printf("%s\t", $tmpArray[1]);
         } # End of else if statement
         elsif ($line =~ m/^Identities =/i) {
            # parsing identities and gaps into @tmpArray
            my @tmpArray = split(/,/, $line);
           
            $tmpArray[0] =~ s/Identities =//i;
            $tmpArray[0] =~ s/\s//g;
            $tmpArray[0] =~ s/\)//g;
            $tmpArray[0] =~ s/\//,/g;
            $tmpArray[0] =~ s/\(/,/g;
            
            my @identities = split(/,/, &trim($tmpArray[0]));
          
            printf("%.2f%%\t", ($identities[0] / $queryLength) * 100);
            printf("%s | %s\t", $identities[0], $identities[1]);
            printf("%s\t", $identities[2]);

            # Check if it contains the gap fields
            if ($#tmpArray + 1 == 2) {
               $tmpArray[1] =~ s/Gaps =//i;
               $tmpArray[1] =~ s/\s//g;
               $tmpArray[1] =~ s/\)//g;
               $tmpArray[1] =~ s/\//,/g;
               $tmpArray[1] =~ s/\(/,/g;

               my @gaps = split(/,/, &trim($tmpArray[1]));
               printf("%s | %s\t", $gaps[0], $gaps[1]);
               printf("%s\t", $gaps[2]);
            } # End of if statement
            else {
               printf("\t\t");
            } # End of else statement
         } # End of if statement
         elsif ($line =~ m/^Strand =/i) {
            $line =~ s/Strand =//i;
            printf("%s\t", &trim($line));
         } # End of else if statement
         elsif ($line =~ m/^Query: /i) {
            my @tmpArray = split(/\s/, $line);

            if ($queryStart == 0) {
               $queryStart = $tmpArray[1];
            } # End of if statement

            $queryEnd = $tmpArray[$#tmpArray];
         } # End of else if statement
         elsif ($line =~ m/^Sbjct: /i) {
            my @tmpArray = split(/\s/, $line);

            if ($subjectStart == 0) {
               $subjectStart = $tmpArray[1];
            } # End of if statement

            $subjectEnd = $tmpArray[$#tmpArray];
         } # End of else if statement
      } # End of if statement
   } # End of while loop

   # print the $queryStart, $queryEnd, $subjectStart, $subjectEnd for the last entry
   if ($queryStart != 0 || $subjectStart != 0) { 
      printf("%d\t%d\t%d\t%d\n", $queryStart, $queryEnd, $subjectStart, $subjectEnd);
   } # End of if statement                                                    
} # End of if statement
else {
   &printUsage();
} # End of else statement

