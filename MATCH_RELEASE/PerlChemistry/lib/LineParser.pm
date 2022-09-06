package LineParser;

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

our @ISA = qw(BaseObject Exporter);

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

sub New {
	
	my ($Class, $Representation, $MatchingSub, $ResultsSub) = @_;

  my $Self = {
	  _Representation => $Representation,
	  _MatchingSub    => $MatchingSub,
	  _ResultsSub     => $ResultsSub,
	  _Results			  => [],
	
  };
	 
	
	bless $Self, $Class;
	
	return $Self;
	
}

sub Parse {
	
	my ($Self, $Line) = @_;
	
	confess("No Matching Sub!") if ! $Self->getMatchingSub;
	
	my $MatchingSub = $Self->getMatchingSub;
		
	my @Results;
	
	eval { @Results = &$MatchingSub($Line) };

  confess("Parse does not have a valid subroutine!") if($@); 

  if($Results[0]) {
	  $Self->setResults(\@Results);
	  return 1;
  }
	
	return 0;
	
}

sub ReturnResults {
	
	my $Self = shift;
			
	confess("No Results Sub!") if ! $Self->getResultsSub;
	
	my $ResultsSub = $Self->getResultsSub;
	
	return &$ResultsSub($Self->getResults); 
	
	
}




1;
__END__
