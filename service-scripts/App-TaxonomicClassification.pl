#
# App wrapper for taxonomic classification.
# Initial version that does not internally fork and report output; instead
# is designed to be executed by p3x-app-shepherd.
#

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::ReadSet;
use File::Slurp;
use IPC::Run;
use Cwd qw(abs_path getcwd);
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
    
    my %db_map = (
		  bvbrc => 'bvbrc',
		  standard => 'standard');
    my $db_dir = $db_map{$params->{database}};
    if (!$db_dir)
    {
	die "Invalid database name '$params->{database}' specified. Valid values are " . join(", ", map { qq("$_") } keys %db_map);
    }

	process_read_input($app, $params);


}

sub process_read_input
{
    my($app, $params) = @_;

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

    #
    # Modify the parameters to remove any SRA libraries and instead
    # use the SE & PE lists from the readset.
    #

    my $nparams = { single_end_libs => [], paired_end_libs => [], srr_libs => [] };;

    $readset->visit_libraries(sub { my($pe) = @_;
				    my $lib = {
					read1 => abs_path($pe->{read_path_1}),
					read2 => abs_path($pe->{read_path_2}),
					(exists($pe->{sample_id}) ? (sample_id => $pe->{sample_id}) : ())
					};
				    push(@{$nparams->{paired_end_libs}}, $lib);
				},
			      sub {
				  my($se) = @_;
				  my $lib = {
				      read => abs_path($se->{read_path}),
				      (exists($se->{sample_id}) ? (sample_id => $se->{sample_id}) : ())
				      };
				  push(@{$nparams->{single_end_libs}}, $lib);
			      },
			     );
    $params->{$_} = $nparams->{$_} foreach keys %$nparams;

    my $json_string = encode_json($params);
    # pushing the wrapper command
    print STDERR "Starting the python wrapper....\n";

    #
    # Create json config file for the execution of this workflow.
    # If we are in a production deployment, we can find the workflows
    # by looking in $KB_TOP/workflows/app-name
    # Otherwise they are in the module directory; this is indicated
    # by the value of $KB_MODULE_DIR (note this is set for both
    # deployed and dev-container builds; the deployment case
    # is determined by the existence of $KB_TOP/workflows)
    #

    my %config_vars;
    my $wf_dir = "$ENV{KB_TOP}/workflows/$ENV{KB_MODULE_DIR}";
    if (! -d $wf_dir)
    {
	$wf_dir = "$ENV{KB_TOP}/modules/$ENV{KB_MODULE_DIR}/workflow";
    }
    -d $wf_dir or die "Workflow directory $wf_dir does not exist";

    #
    # Find snakemake. We need to put this in a standard location in the runtime but for now
    # use this.
    #
    my $snakemake = "$ENV{KB_RUNTIME}/artic-ncov2019/bin/snakemake";
    $config_vars{workflow_dir} = $wf_dir;
    $config_vars{input_data_dir} = $staging;
    $config_vars{output_data_dir} = $output;
    $config_vars{snakemake} = $snakemake;
    $config_vars{params} = $params;
    $config_vars{cores} = $ENV{P3_ALLOCATED_CPU} // 2;
    write_file("$top/config.json", JSON::XS->new->pretty->canonical->encode(\%config_vars));

    my @cmd = ("python3", "$wf_dir/snakefile/wrapper.py", "$top/config.json");

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
