#!/usr/bin/env perl

=head1 NAME

GenerateHbondingInfo.pl

=head1 SYNOPSIS

MATCH.pl generates topology and parameter files for CHARMM force fields

=head2 EXPORT

NONE

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

use strict;
use warnings;
use Carp;

#Get MATCH libraries 
use lib $ENV{'MATCH'} . "/lib"; 

use MATCHBaseObject;
use BaseObject ':vars';

use MATCHParameters;
use MoleculeFileHandler;

use MATCHer;

#Handle CommandLine Arguments

my @Arguments;
my $ParameterFilePath;
my @Molecules;
my $IsParameter = 0;

foreach my $i (0 .. @ARGV-1) {
	
		
	if(lc($ARGV[$i]) eq "-forcefield") {
		
		croak "Invoked -forcefield but did not specify which one!" if($i+1 == @ARGV); 
		
		$ParameterFilePath = $ARGV[$i+1] . ".par"; 
		
		$IsParameter = 1;
		
	}
	
	elsif(lc($ARGV[$i]) =~ /^-h/) {
		
		
		print "Usage: \n";
		print "MATCH.pl -forcefield forcefield [parameters] filepath\n";
		print "MATCH.pl -parameterfile parameterfilepath filepath\n";
		print "MATCH.pl -parameterlist for full list of parameters\n";
		
		exit 1;
		
	}
	
	elsif(lc($ARGV[$i]) =~/^-parameterlist/) {
		
		print "Supported Parameters\n";	
	  print sprintf("%-30s", "Parameter") . "Value Type\n";
		my @Parameters = ( ['BondLengthFilePath' ,'PATH'],  ['CreatePdb', 'PATH'], ['ExitifNotInitiated', 'Binary'], ['ExitifNotTyped', 'Binary'], ['ExitifNotCharged', 'Binary'],
		['IncrementFilePath', 'PATH'], ['ImproperFilePath', 'PATH'],  ['ParameterFilePath', 'PATH'], ['RelationMatrixFilePath', 'PATH'], ['RefiningIncrementsFilePath', 'PATH'],
		['ShortenTypeFilePath', 'PATH'], ['SubstituteIncrements', 'PATH'], ['TypesFilePath', 'PATH'], ['UsingRefiningIncrements', 'PATH']
		
		
		);
		
		foreach my $Parameter (@Parameters) {
			
			print sprintf("%-30s", $Parameter->[0]) . $Parameter->[1] . "\n";
			
		}
		
		exit 1;
		
	}
	
	elsif($ARGV[$i] =~ /^-/) {
		
		croak "Invoked " . $ARGV[$i] . " but did not specify its value!" if($i+1 == @ARGV); 
		
		push @Arguments, [$ARGV[$i], $ARGV[$i+1]]; $IsParameter = 1;
		
	}
	
	else {
		
		if($IsParameter) { $IsParameter = 0; next }
		
		push @Molecules, $ARGV[$i];
		
	}
	
}

#Setup Default Parameters

my $DefaultParameters = MATCHParameters->New;

$DefaultParameters->Initiate($ParameterFilePath, @Arguments);

$Parameters = $DefaultParameters;

#Load Molecules

croak "Please specify path of File to be processed\n" unless scalar(@Molecules) > 0;

#Initiate MATCHer

my $MATCHer = MATCHer->New;

$MATCHer->Initiate;

my $PreviousTopologyName;

foreach my $MoleculeFilePath (@Molecules) {
	
	if($MoleculeFilePath =~ /.mol2/) {
		
		open(FILE, $MoleculeFilePath);
		
		my @FileContents = <FILE>;
		
		close(FILE);
		
		open(FILE, ">$MoleculeFilePath");
		
		print FILE @FileContents;
		
		print FILE "\n\n";
		
		close(FILE);		
	}

	

  my $MoleculeFileHandler = MoleculeFileHandler->New($MoleculeFilePath);

  my $Structures = $MoleculeFileHandler->BuildObjectsFromFile;

  my $TestMolecule = $Structures->[0]->getMolecule(0);

  #$TestMolecule->Initiate;

	#MATCH molecule
		
	$Parameters->setAppendingTopologyFilePath("top_" . $PreviousTopologyName . ".rtf") if defined $PreviousTopologyName;
	
	$Parameters->setAppendingParameterFilePath($PreviousTopologyName . ".prm") if defined $PreviousTopologyName;
	
	$PreviousTopologyName = undef;
	
	#$TestMolecule->setName("UNK");
	
	
	print $MATCHer->getHydrogenBondTopologyInformation($TestMolecule) . "\n";
	
	exit 1;
	
  my $Result = $MATCHer->BuildBothTopologyAndParameterFileForMolecule($TestMolecule, $PreviousTopologyName);

  if($Result == 0) {
	
	  print $TestMolecule->getFileName . " Failed!\n";
	  next;
	
	}
	
	print $TestMolecule->getFileName . " Success!\n";

  $PreviousTopologyName = $TestMolecule->getFileName unless defined $PreviousTopologyName;

  if($Parameters->getCreatePdb) {
	
	  if($TestMolecule->getName !~ /^\S+$/) {
		
		  $TestMolecule->setName("UNK");
		
	  }
			
	  $TestMolecule->setSegmentId("TEST") if ! $TestMolecule->getSegmentId;
	
	  open(FILE, ">" . $Parameters->getCreatePdb);
	
	  print FILE $TestMolecule->ToStringInPdbFormat . "\n";
	
	  print FILE "END\n";
	
	  close(FILE);
	
  }

}

