package Parameters;

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
use BaseObject ':all';

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
	
	my $Class = shift;
	
  my $Self = {
	  _PATH 							       => undef,
	  _BondLengthFilePath        => undef,
	  _PrintOutPath              => "",
		_ExitifNotInitiated	       => 1,
		_CheckAtomProtState        => 1,
		_CheckElementBondingNumber => 1,
		_RenameAtoms               => 0,
	
	};
	
	bless $Self, $Class;
}
	

sub Initiate {
	
  my ($Self, $ParameterFilePath, @CommandLineParameters) = @_;

  my $DefaultPath = $ENV{'PerlChemistry'};

  $ParameterFilePath = $ENV{'PerlChemistry'} . "/resources/DefaultParameters.par" unless defined $ParameterFilePath;

  $Self->ReadParameterFile($DefaultPath, $ParameterFilePath);

  $Self->ProcessComandLineArguments($DefaultPath, @CommandLineParameters);

}

sub ProcessComandLineArguments {
	
	my ($Self, $DefaultPath, @CommandLineParameters) = @_;
	
	foreach my $CommandLineParameter (@CommandLineParameters) {
		
	  $Self->ProcessParameter($DefaultPath, substr($CommandLineParameter->[0],1), $CommandLineParameter->[1]);
		
	}
	
}

sub ReadParameterFile {
	
	my ($Self, $DefaultPath, $ParameterFilePath) = @_;
	
	confess("Path to ParameterFile does not exist") unless -e $ParameterFilePath;
	
	open(FILE, $ParameterFilePath);

  my %ParameterFileContents = map { my ($ParameterName, $Value) = ($_ =~ /\s*(\S+)\s*=\s*(\S+)/ );  $ParameterName => $Value } 
 														  grep { $_ =~ /^[^#]\s*(\S+)\s*=\s*(\S+)/ } <FILE>;
 
 close(FILE);

  while (my ($Parameter, $Value) = each %ParameterFileContents) {
	
    $Self->ProcessParameter($DefaultPath, $Parameter, $Value);  
	
	}
	
}

sub ProcessParameter {
	
	my ($Self, $DefaultPath, $Parameter, $Value) = @_;
	
	if($Value =~ /\// ) {
	
	  unless(-e $Value) {
		
		  confess("Parameter: " . $Parameter. " File does not exist: $Value")  unless(-e "$DefaultPath/$Value" || $Value =~ /\d{7}?.out/);
	    
	    $Value = "$DefaultPath/$Value" unless ($Value =~ /\d{7}?.out/);
#		  confess("Parameter: " . $Parameter. " File does not exist: $Value")  unless(-e "$DefaultPath/$Value");
#	    
#	    $Value = "$DefaultPath/$Value";
		
	  }
	
  }

  PrintMessage("Warning: $Parameter is an unknown parameter!", 1) unless exists $Self->{"_" . $Parameter};

  $Self->{"_" . $Parameter} = $Value;
	
}




1;
__END__
