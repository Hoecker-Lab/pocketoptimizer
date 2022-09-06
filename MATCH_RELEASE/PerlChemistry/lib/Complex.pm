package Complex;

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
use BaseObject;
use Chain;

require Exporter;

our @ISA = qw(BaseObject Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

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

=head2 New

Usage: Complex->New;

Arguments:
  $Class: should be 'Complex'

Synopsis:
  creates a new Complex object which is the Highest level of datastructure organization, it contains all informatin about a system, 
all chains 

=cut

sub New {
	
	my ($Class, $StructuralHash) = @_;
	
	
	my $Self = {
		
		_Chains      => [ ], 
		_FileName		 => undef,
		_Name				 => undef,
		
	};
	
	bless $Self, $Class;
	
  $Self->BuildFromStructuralHash($StructuralHash) if defined $StructuralHash;
	
	return $Self;
	
}

=head2 BuildFromStructuralHash

Usage: $ComplexInstance->BuildFromStructuralHash;

Arguments:
  $StructuralHash: should be 'Complex'

Synopsis:
  Loads data from StructuralHash built from loading in from file

=cut

sub BuildFromStructuralHash {
	
	my ($Self, $StructuralHash) = @_;
	
	confess("This Hash is not of a Complex") if $StructuralHash->{'Structure'} ne ref($Self);
	
	my $Chains;
	
	$Chains = $StructuralHash->{'Chains'} if exists $StructuralHash->{'Chains'};
	
	my @ChainObjects;
			
	while( my ($ChainNumber, $ChainHash) = each(%$Chains)) {
				
		my $ChainObject = Chain->New($ChainHash);
		
		$ChainObject->setNum($ChainNumber);
		
		push @ChainObjects, $ChainObject;
				
	}	
	
	$Self->setChains(\@ChainObjects);
	
	
}

=head2 IdentifyChains

Usage: $ComplexInstance->IdentifyChains;

Arguments:

Synopsis:
  Tries to identify each Chain as RNA, DNA, Proteins, Water, or Ligand, Work in progress!

=cut

sub IdentifyChains {
	
	my $Self = shift;
	
	my $Chains = $Self->getChains;
	
	$_->IdentifyThisChain foreach (@$Chains);
	
}

=head2 getChainsOfType

Usage: $ComplexInstance->getChainsOfType($Type);

Arguments:
  $Type: The type of chain you wish to recover, Protein, Nucleic, Ligand, Water etc

Synopsis:
  returns the chains of specific type that are of interest

=cut

sub getChainsOfType {
	
	my ($Self, $Type) = @_;
	
	my $Chains = $Self->getChains;
	
	my @ChainsOfThisType =  grep { $_->getType && $_->getType eq $Type } @$Chains;
				
	return \@ChainsOfThisType;
	
}




1;
__END__
