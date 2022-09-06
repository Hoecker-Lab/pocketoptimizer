package MATCHParameters;

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
use MATCHBaseObject;
use Parameters;

require Exporter;

our @ISA = qw(Parameters Exporter);

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
	
	my $Class = shift;
	
	my $Parameters = Parameters->New;
	
	my $Self = {
		
    _AddHydrogenBondToplogyInformation  => 0,
	  _AppendingForceField							  => undef,
		_AppendingParameterFilePath 				=> undef,
		_AppendingTopologyFilePath  				=> undef,
		_BondLengthFilePath                 => undef,
		_AtomicProtonationStates            => undef,
		_CreatePdb													=> 0,
		_DoNotCharge												=> 0,
		_ExitifNotTyped										  => 1, 
		_ExitifNotCharged                   => 1,
		_ExitifNotParameterized             => 1,
		_ForceCharging										  => 0,
		_ForceField 												=> undef,
		_IncrementTrainingFilePath				  => undef,
		_IncrementFilePath 								  => undef,
		_ImproperFilePath										=> undef,
		_ParameterFilePath 								  => undef,
	  _Parameterize					  						=> 1,
		_RefiningIncrementsFilePath			    => undef,
		_RelationMatrixFilePath 						=> undef,
		_ShortenTypeFilePath							  => undef,
		_StandAloneParameterFile            => 0, 
		_SubstituteIncrements 							=> 1,
		_TypesFilePath    							  	=> undef,
		_UsingRefiningIncrements						=> 1,
		
	};
	
	while(my ($Field, $Value) = each %$Parameters) {
		
		$Self->{$Field} = $Value;
		
	}
	
	bless $Self, $Class;
	
}

sub Initiate {
	
	my ($Self, $ParameterFilePath, @CommandLineParameters) = @_;

  my $DefaultPath = $ENV{'MATCH'};

  $ParameterFilePath = $ENV{'MATCH'} . "/resources/DefaultParameters.par" unless defined $ParameterFilePath;

  unless(-e $ParameterFilePath)  {
	
	  croak "ParameterFilePath does not exist" unless(-e $ENV{'MATCH'} . "/resources/$ParameterFilePath");
	
	  $ParameterFilePath = $ENV{'MATCH'} . "/resources/$ParameterFilePath";
	
  }

  $Self->ReadParameterFile($DefaultPath, $ParameterFilePath);

  $Self->ProcessComandLineArguments($DefaultPath, @CommandLineParameters);
  
	
}




1;
__END__
