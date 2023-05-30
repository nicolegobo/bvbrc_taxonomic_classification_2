#
# App wrapper for taxonomic classification.
# Initial version that does not internally fork and report output; instead
# is designed to be executed by p3x-app-shepherd.
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::ReadSet;
use IPC::Run;
use Cwd;
use File::Path 'make_path';
use strict;
use Data::Dumper;
use File::Basename;
use File::Temp;
use JSON::XS;
use Getopt::Long::Descriptive;

my $app = Bio::KBase::AppService::AppScript->new(\&run_classification, \&preflight);

$app->run(\@ARGV);

sub run_classification
{
    my($app, $app_def, $raw_params, $params) = @_;
    
    #
    # Set up options for tool and database.
    #
    
    my @cmd;
    my @options;

    my %db_map = (
		  bvbrc => 'bvbrc',
		  standard => 'standard');
    my $db_dir = $db_map{$params->{database}};
    if (!$db_dir)
    {
	die "Invalid database name '$params->{database}' specified. Valid values are " . join(", ", map { qq("$_") } keys %db_map);
    }
    
    if ($params->{input_type} eq 'reads')
    {
	process_read_input($app, $params, \@cmd, \@options);
    }
    else
    {
	die "Invalid input type '$params->{input_type}'";
    }

}

sub process_read_input
{
    my($app, $params, $cmd, $options) = @_;

    my @cmd = @$cmd;
    my @options = @$options;
    my $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params, 1);
    
    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($app->workspace);
    
    if (!$ok)
    {
	die "Readset failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }
    
    my $top = getcwd;
    my $staging = "$top/staging";
    my $output = "$top/output";
    make_path($staging, $output);
    $readset->localize_libraries($staging);
    $readset->stage_in($app->workspace);

    my $json_string = encode_json($params);
    # pushing the wrapper command
    print('Starting the python wrapper....');

    #### relative path from service-script dir ####
    @cmd = ("../workflow/snakefile/wrapper.py");
            push(@cmd, $json_string);

    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "wrapper command failed $?: @cmd";
    }

    save_output_files($app, $output);
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    my $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);

    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($app->workspace);

    if (!$ok)
    {
	die "Readset failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }

# this will require changes once we decide what to do about the databases
    # my $mem = "32G";
    my $mem = "160G";

    #
    # Kraken DB requires a lot more memory.
    #
    # if (lc($params->{database}) eq 'bvbrc')
    # {
	# $mem = "160G";
    # }
    
    my $time = 60 * 60 * 10;
    my $pf = {
	cpu => 6,
	memory => $mem,
	runtime => $time,
	storage => 1.1 * ($comp_size + $uncomp_size),
	policy_data => { constraint => 'kraken_db', partition => 'kraken_db' }
    };
    return $pf;
}

sub save_output_files
{
    my($app, $output) = @_;
    my %suffix_map = (
        csv => 'csv',
        err => 'txt',
        html => 'html',
        out => 'txt',
        txt => 'txt',);

    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;

    if (opendir(D, $output))
    {
    while (my $p = readdir(D))
    {
        my @cmd = ("p3-cp", "--recursive", @suffix_map, "$output/", "ws:" . $app->result_folder);
        print STDERR "saving files to workspace... @cmd\n";
        my $ok = IPC::Run::run(\@cmd);
        if (!$ok)
        {
        warn "Error $? copying output with @cmd\n";
        }
    }
    closedir(D);
    }
}