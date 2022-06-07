#!/usr/bin/env perl

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head2 EXPORT


=head1 AUTHOR

Joseph Yesselman, jyesselm@umich.edu and Charles L. Brooks III, brookscl@umich.edu

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010 by Joseph Yesselman and Charles L. Brooks III

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
=head1 FUNCTIONS

=cut

#use 5.010000;
use strict;
use warnings;
use Carp;

#Load MATCH Libraries 
use lib $ENV{"MATCH"} . "/lib";

use MATCHBaseObject;
use MATCHFunctions ':func';
use LookUpTable;
use Storable;
use Parameters;
use MATCHParameters;
use BaseObject ':vars';

our $VERSION = '0.01';
#Exported Variables


#Non Exported Variables

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate;

$Parameters = $DefaultParameters;


# Preloaded methods go here.

srand(time ^ ($$ + ($$ << 15)));

my $FileName = $ARGV[0];

my $TopologyFilePath = $ENV{"MATCH"} . "/resources/$FileName/$FileName.rtf";
my $TopologyPatchingFile = $ENV{"MATCH"} . "/resources/$FileName/$FileName.patches";

my $Chains = SetupMoleculesFromTopology($TopologyFilePath, $TopologyPatchingFile);

open(FILE, $ENV{"MATCH"} . "/resources/Incorrect_Residue_List.txt");
my %IncorrectResidues = map { chomp($_); $_ => 1 } <FILE>;
close(FILE);

print scalar(@$Chains) . "\n";
my @CorrectChains = grep { ! exists $IncorrectResidues{$_->getMolecule(0)->getName } } @$Chains;
print scalar(@CorrectChains) . "\n";
#my $Chains = retrieve("Chains.dat");

store(\@CorrectChains, $ENV{"MATCH"} . "/t/StoredTopologys/$FileName.dat");

exit 1;

1;
__END__
