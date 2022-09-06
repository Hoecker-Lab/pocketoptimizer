package MoleculeFileHandler::PdbFileHandler;

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
use MoleculeFileHandler;

require Exporter;

our @ISA = qw(MoleculeFileHandler Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BaseObject ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all'  => [ qw( )],
											 
										 'vars' => [ qw () ],
										
										 'func' => [ qw ()]);
							

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

=head2 Initiate;

Usage: $PdbFileHandlerInstance->Initiate;

Arguments:

Synopsis:
  Setups $PdbFileHandlerInstance must be called

=cut


sub Initiate {
	
	my $Self = shift;
			
  my $EndPattern = qr/^TER/;
  my $EndComplexPattern = qr/^ENDMDL/;
	
	my $BondPattern = qr/^CONECT/;
	
	#Protein DataBank formatted
	my $Atom = qr/^(?:ATOM|HETATM)/;
  my @AtomVariables = qw(Num Name ResidueName ChainName ResidueNum X Y Z Occupancy TempFactor SegmentId);
  my @AtomSubStrs = (5, 6, 4, 1, 6, 10, 8, 8, 7, 6, -1);
	

  my $AtomLineParser =	    LineParser->New('Atom',
																			    	sub { if($_[0] !~ $Atom) { return ( ) }

																							  my $string = substr($_[0], 6); 

																							  return map { my $substr = substr($string, 0, $_);
																														 $string = substr($string, $_);
																														 $substr } @AtomSubStrs 
																					  }, 
																				  	sub { my $Hash;
																								foreach (0 .. $#AtomVariables) {
																									$_[0]->[$_] =~ s/^\s*(.*?)\s*$/$1/;
																									my $value = $_[0]->[$_];																									
																									$Hash->{$AtomVariables[$_]} = $value if defined $value;
																								}
																								return $Hash  });
  	
	my $BondLineParser = 							LineParser->New('Bond',
																							      sub { if($_[0] =~ $BondPattern) { my @spl = split /\s+/, $_[0]; shift @spl; @spl }},
																							      sub { return { Num => shift @{$_[0]}, BondedToNums => $_[0] }});
																							
  my $EndOfChainLineParser = 			  LineParser->New('EndChain',
																								    sub { $_[0] =~ $EndPattern}, 
																										sub { });
																										
	my $EndofCompelexLineParser	=		  LineParser->New('EndOfComplex',
																										sub { $_[0] =~ $EndComplexPattern},
																										sub { });
	
	$Self->{_Patterns} = [ $AtomLineParser,
	                       $BondLineParser,
	 											 $EndOfChainLineParser, 
	                       $EndofCompelexLineParser];
		
}

=head2 ReadFile;

Usage: $PdbFileHandlerInstance->ReadFile($File);

Arguments:

Synopsis:
  This is where all the data is collected from File

=cut

sub ReadFile {
	
	my ($Self, $FilePath) = @_;
	
	open(FILE, $FilePath);
	
	my @Chains;
	my $Finished = 0;
	
	while(!$Finished) {
	  push @Chains, $Self->ReadChainFromFile(*FILE); 
	  $Finished = 1 if($Chains[-1]->[0] eq 'EOF');
  }
	
	close(FILE);
	 	
	my $ComplexFlag = 0;
	
	foreach (@Chains) { 
	
		$ComplexFlag = 1 if $_->[0] eq 'EndOfComplex'; 
		
	}
	
	#Remove A and B crystalization alternatives when part of ResidueNum etc
	foreach my $Chain (@Chains) {
		
	  my $Atoms = $Chain->[1]->{'Atom'};
	
	  my %ResidueNames = map { $_->{'ResidueNum'} => 1} @$Atoms;
	
	  my %ResiduesToDelete;
	  my %SeenThisNum;
	  my %RenameResidue;

    foreach my $Atom (@$Atoms) {
	
	    #print $Atom->{'Name'} . " "  . $Atom->{'Num'} . " " . $Atom->{'ResidueName'} . " " . $Atom->{'ResidueNum'} . " " . $Atom->{'X'}. " " . $Atom->{'Y'} ." " . $Atom->{'Z'} . "\n";
	
		}
 
	
	  foreach my $ResidueName (sort { $a cmp $b } keys %ResidueNames) {
					
		  next unless(uc($ResidueName) =~ /(\d+)[A-Z]/);
		
		  my $Num = $1;
						
			if(exists $ResidueNames{$Num} || exists $SeenThisNum{$Num}) {  $ResiduesToDelete{$ResidueName} = 1; next }
			
			$RenameResidue{$ResidueName} = $Num;
									
			$SeenThisNum{$1} = 1;
			
	  }
	
	  my @GoodAtoms = grep { ! exists $ResiduesToDelete{$_->{'ResidueNum'}}} @$Atoms;
	
	  foreach my $Atom (@GoodAtoms) {
		
		  $Atom->{'ResidueNum'} = $RenameResidue{$Atom->{'ResidueNum'}} if(exists $RenameResidue{$Atom->{'ResidueNum'}});
		
	  }
	
	  $Chain->[1]->{'Atom'} = \@GoodAtoms;
		
	}
				
  if(!$ComplexFlag) {
		
	
    $Self->setChains([ map { $_->[1] } @Chains]);
    return 1;	
  }

	my @Complexes;
	my $CurrentComplex = [ ]; 
	foreach (@Chains) {
		
		if($_->[0] ne 'EndOfComplex') { push @$CurrentComplex, $_->[1] }
		
		else {
		  push @$CurrentComplex, $_->[1];
		  push @Complexes, $CurrentComplex;
		  $CurrentComplex = [ ];	
		}
		
	}

  push @Complexes, $CurrentComplex if @$CurrentComplex;
	
	if(@Complexes) {
		$Self->setComplexes(\@Complexes);
		return 1;
	}
	
	#Something bad happened!	
	confess("ReadFile Terminated Incorrectly!\n");		
}

=head2 ReadChainFromFile;

Usage: $PdbFileHandlerInstance->ReadChainFromFile($File);

Arguments:

Synopsis:
  Collects the data on a specific Chain, this allows the break up of the data

=cut

sub ReadChainFromFile {
	
		my ($Self, $FH) = @_;
		
		my $Patterns = $Self->getPatterns;

		my $Sections;
				
		my $CurrentRepresentation = "";
		
		my $DoneReading = 0;
		
		while(my $Line = <$FH>) {
			
		  if(eof($FH)) { $CurrentRepresentation = 'EOF'; last; }
												
			foreach my $LineParser (@$Patterns) {
				
				next if ! $LineParser->Parse($Line);
				
				my $DataHash = $LineParser->ReturnResults;			
								
				push @{$Sections->{$LineParser->getRepresentation}}, $DataHash if $DataHash;
											
			  $CurrentRepresentation = $LineParser->getRepresentation;
			  			
				last;
			  					
			}		
											
			last if $CurrentRepresentation eq 'EndChain' or $CurrentRepresentation eq 'EndOfComplex' ;
		
		}		
			
		return [$CurrentRepresentation, $Sections];
}



1;
__END__
