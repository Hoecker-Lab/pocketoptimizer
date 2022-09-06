package LookUpTable::GridElement;

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

require Exporter;

our @ISA = qw(BaseObject);

our $VERSION = '0.01';

#Exported Variables


#Non Exported Variables


# Preloaded methods go here.

=head2 New

Usage: LookUpTable::GridElement->New;

Arguments:
  $Class: should be 'LookUpTable::GridElement'

Synopsis:
  creates a new LookUpTable::GridElement object

=cut

sub New { 
  my $Class = shift;
	
	my $Self = { 
	  _NodeReferences     => [ ],
	  _Count							=> 0,
	};
	
	bless $Self, $Class;
  return $Self;
	
};

=head2 Copy

Usage: $GridElementInstance->Copy;

Arguments:

Synopsis:
  Copies GridElement object, is only a shallow copy

=cut

sub Copy {
	
	my $Self = shift;
	
	my $NewGridElement =  LookUpTable::GridElement->New;
				
	$NewGridElement->setNodeReferences([map { $_ } @{$Self->getNodeReferences}]);
	$NewGridElement->setCount($Self->getCount);
	
	return $NewGridElement;
}

=head2 RemoveNodeReference

Usage: $GridElementInstance->RemoveNodeReference($NodeReference);

Arguments:
  $NodeReference: The Reference that needs to be removed from $GridElemenetInstance

Synopsis:
  Removes a Node reference from $GridElementInstance, this function is primarily used when an Node is being deleted from the Current lookup Table

=cut

sub RemoveNodeReference {
  my ($Self, $NodeReference) = @_;

  my $NodeReferences = $Self->getNodeReferences;

  foreach (0 .. @$NodeReferences-1) {
	  if (${$NodeReferences->[$_]} eq $$NodeReference) {
		  splice (@$NodeReferences, $_, 1); last;
	  }
  }

  $Self->setCount($Self->getCount -1);

}

=head2 AddNodeReference

Usage: $GridElementInstance->AddNodeReference($NodeReference);

Arguments:
  $NodeReference: The Reference that needs to be added from $GridElemenetInstance

Synopsis:
  Add reference to Node object to $GridElementInstance, Updates internal Count variable

=cut

sub AddNodeReference {
	my ($Self, $NodeReference) = @_;
		
	my $NodeReferences = $Self->getNodeReferences;
	
	my $Included = 0;
	foreach (@$NodeReferences) {
		
		if($$_->getAtomName eq $$NodeReference->getAtomName) { $Included = 1; last }
		
	}
	
	if(!$Included) {
		
	  push @$NodeReferences, $NodeReference;
	
	  $Self->setCount($Self->getCount + 1);
	
  }
	
}

=head2 Merge

Usage: $GridElementInstance->Merge($OtherGridElementInstance);

Arguments:
  $OtherGridElementInstance: The Other GridElementInstance that needs to be merged with the current one.

Synopsis:
  This is basically just overloading the + operator, Merge copies all references from $OtherGridElementInstance into $GridElementInstance

=cut

sub Merge {
	
	my ($Self, $OtherGridElement) = @_;
	
	confess("Cannot merge with this object, can only with other LookUpTable::GridElements!") if ref($OtherGridElement) ne 'LookUpTable::GridElement';
		
	my $OtherNodeReferences = $OtherGridElement->getNodeReferences;
			
	foreach (@$OtherNodeReferences) { $Self->AddNodeReference($_) }
			
}

=head2 Update

Usage: $GridElementInstance->Update;

Arguments:

Synopsis:
  Updates the Count internal variable

=cut

sub Update { 
	my $Self = shift;
	
	my $NodeReferences = $Self->getNodeReferences;
	
	my $Count = 0;
	
	foreach my $NodeReference (@$NodeReferences) {
		
		$Count++ if $$NodeReference->getUsed == 0;
		
	}
	
	#my $Count = scalar(grep {  $$_->getUsed == 0 } @$NodeReferences);
		
	$Self->setCount($Count);
}

1;
__END__
