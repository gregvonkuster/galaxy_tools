#!/usr/bin/env perl
# Author: Eric Wafula
# Email: ekw10@psu.edu
# Institution: Penn State University, Biology Dept, Claude dePamphilis Lab
# Date: June 2018

use strict;
use warnings;
use File::Spec;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use FindBin;
use DBI;

my $home =  "$FindBin::Bin/..";

my $usage = <<__EOUSAGE__;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#                                  GENE FAMILY SCAFFOLD UPDATER
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Required Options:
#
#
#  --database_connection_string <string>    : Postgres database connection string using format
#                                             postgresql://<user>:<password>@<host>/<database name>
#
#  --proteins <string>                      : Amino acids (proteins) sequences fasta file (proteins.fasta)
#                                             This can either be an absolute path or just the file name
#
#  --coding_sequences <string>              : Corresponding coding sequences (CDS) fasta file (cds.fasta)
#
#  --scaffold <string>                      : Orthogroups or gene families proteins scaffold.  This can either be an absolute
#                                             path to the directory containing the scaffolds (e.g., /home/scaffolds/22Gv1.1)
#                                             or just the scaffold (e.g., 22Gv1.1).  If the latter, ~home/data is prepended to
#                                             the scaffold to create the absolute path.
#                                             the scaffold to create the absolute path.
#                                             If Monocots clusters (version 1.0): 12Gv1.0
#                                             If Angiosperms clusters (version 1.0): 22Gv1.0
#                                             If Angiosperms clusters (version 1.1): 22Gv1.1
#                                             If Green plants clusters (version 1.0): 31Gv1.0
#                                             If Other non PlantTribes clusters: XGvY.Z, where "X" is the number species in the scaffold,
#                                             and "Y.Z" version number such as 12Gv1.0. Please look at one of the PlantTribes scaffold
#                                               data on how data files and directories are named, formated, and organized.
#
#
#  --species_name <string>                  : Name of the species
#
#  --species_code <string>                  : Code of the species
#
#  --species_family <string>                : Family of the species
#
#  --species_order <string>                 : Order of the species
#
#  --species_group <string>                 : Group of the species
#
#  --species_clade <string>                 : Clade of the species
#
#  --rooting_order_species_code <string>    : Species code after which the new species will be placed in the rooting order config file
#
# # # # # # # # # # # # # # # # # #
#  Others Options:
#
#  --num_threads <int>                      : number of threads (CPUs) to used for HMMScan, BLASTP, and MAFFT
#                                             Default: 1
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Example Usage:
#
#  GeneFamilyScaffoldUpdater --database_connection_string postgresql://<user>:<password>@<host>/<database name>
                             --proteins proteins.fasta --coding_sequences cds.fasta --scaffold 22Gv1.1
#                            --species_name Fake genome --species_family Brassicaceae --species_order Brassicales
                             --species_group Rosids --species_clade Core Eudicots --rooting_order_species_code Phypa
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

__EOUSAGE__
    ;

# Declare and initialize variables;
my $database_connection_string;
my $proteins;
my $species_name;
my $species_code;
my $species_family;
my $species_order;
my $species_group;
my $species_clade;
my $coding_sequences;
my $scaffold;
my $rooting_order_species_code;
my $num_threads;

my $options = GetOptions ( 'database_connection_string=s' => \$database_connection_string,
                           'proteins=s' => \$proteins,
                           'species_name=s' => \$species_name,
                           'species_code=s' => \$species_code,
                           'species_family=s' => \$species_family,
                           'species_order=s' => \$species_order,
                           'species_group=s' => \$species_group,
                           'species_clade=s' => \$species_clade,
                           'coding_sequences=s' => \$coding_sequences,
                           'scaffold=s' => \$scaffold,
                           'rooting_order_species_code=s' => \$rooting_order_species_code,
                           'num_threads=i' => \$num_threads,
                         );

# # # # # # # # # # # # # # # # # # # # # # # # #  validate options and set variables  # # # # # # # # # # # # # # # # # # # # # # # # # #
# check if options are set
unless ( $options ) { die $usage; }
unless ( $database_connection_string and $proteins and $species_name and $species_code and $species_family and $species_order and $species_group and $species_clade and $coding_sequences and $scaffold and $rooting_order_species_code ) {
    print "\nOne or more required options not set\n"; die $usage; }
# get scaffold directory
my $scaffold_dir;
if (File::Spec->file_name_is_absolute($scaffold)) {
    $scaffold_dir = $scaffold;
    $scaffold = basename($scaffold);
} else {
    if ($scaffold) { $scaffold_dir = "$home/data/$scaffold"; }
    else { print "\n --scaffold option is not set\n\n"; die $usage; }
}

# validate scaffold and update type options
if ( $scaffold !~ /^\d+Gv\d+\.\d+$/) {
    print "\nOrthogroups or gene families proteins scaffold name $scaffold is not in the required format";
    print " i.e. XGvY.Z, where X is number species in the scaffold, and Y.Z version number such as 12Gv1.0.\n";
    die $usage;
}

# Find out if the received rooting_order_species_code is in
# the rooting order configuration file for the scaffold.  We
# do this before anything else since it is the least resource
# intensive and an invalid species code will force an error.
validate_rooting_order_species_code($scaffold_dir, $rooting_order_species_code);

# get scaffold clustering methods
my %methods;
my $annotation_dir = "$scaffold_dir/annot";
opendir (DIR, "$annotation_dir") or die "Can't open $annotation_dir directory\n";
while (my $filename = readdir(DIR)) {
    if ($filename =~ /(\w+)\.list$/) { $methods{$1} = $1; }
}

# set defaults
if (!$num_threads) { $num_threads = 1; }

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  sub-routine calls  # # # # # # # # # # # # # # # # # # # # # # # # # # #

print localtime()." - Starting gene family scaffold updating\n";

# Create working directory.
my $working_dir = "./geneFamilyScaffoldUpdate_dir";
if (-d $working_dir) {
    die "Exiting...!\nGene family scaffold update output directory ($working_dir) already exists!\n\n";
}
make_directory($working_dir);

# Copy original scaffold data to a working directory.
print "-- ".localtime()." - Copying original scaffold data to working directory\n\n";
my $copy_scaffold_data = system "cp -r $scaffold_dir $working_dir";
if ($copy_scaffold_data != 0) {
    print "\n-- ".localtime()." - Copying original scaffold data to working directory failed\n\n"; exit(1);
}

# Update the scaffold config files in the working directory with the new genome.
update_config_files ( $scaffold, $rooting_order_species_code, $species_name, $species_code, $species_family, $species_order, $species_group, $species_clade, $working_dir );

# Update the scaffold files in the working directory with the new genome.
foreach my $method (keys %methods) {
    sort_sequences ( $proteins, $coding_sequences, $scaffold, $method, $num_threads, $species_name, $working_dir, $scaffold_dir );
}

# Move updated scaffold data to original directory.
my $updated_scaffold_dir = "$working_dir/$scaffold";
print "-- ".localtime()." - Removing original scaffold data directory\n$scaffold_dir\n\n";
my $remove_scaffold_data = system "rm -rf $scaffold_dir";
if ($remove_scaffold_data != 0) {
    print "\n-- ".localtime()." - Removing original scaffold data directory failed\n\n"; exit(1);
}
print "-- ".localtime()." - Moving updated scaffold data from\n$updated_scaffold_dir\nto original directory\n$scaffold_dir\n\n";
my $move_scaffold_data = system "mv $updated_scaffold_dir $scaffold_dir";
if ($move_scaffold_data != 0) {
    print "\n-- ".localtime()." - Moving updated scaffold data to original directory failed\n\n"; exit(1);
}

# Update the database tables with the new genome.
update_database_tables ( $database_connection_string, $proteins, $scaffold, \%methods, $species_name, $species_family, $species_order, $species_group, $species_clade, $scaffold_dir, $working_dir );

print localtime()." - Completed gene family scaffold updating\n\n";

exit(0);

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # sub-routines # # # # # # # # # # # # # # # # # # # # # # # # # # #

sub validate_rooting_order_species_code {
    my ($scaffold_dir, $rooting_order_species_code ) = @_;
    my $rooting_order_config = "$scaffold_dir/$scaffold.rootingOrder.config";
    open (IN, $rooting_order_config) or die "Can't open $rooting_order_config\n";
    while (<IN>) {
        chomp;
        if (/^#/ || /^$/) {
            # Skip comments and section headers.
            next;
        }
        # Example line: Physcomitrella patens=Phypa
        my @F = split(/=/, $_);
        if ($F[1] eq $rooting_order_species_code) {
            return;
        }
    }
    print "\n-- ".localtime()." - Invalid rooting order species code: $rooting_order_species_code\n\n"; exit(1);
}

sub make_directory {
        my ( $new_dir ) = @_;
        if (!-d $new_dir) {
                mkdir($new_dir, 0755);
        }
}

sub update_config_files {
    my ( $scaffold, $rooting_order_species_code, $species_name, $species_code, $species_family, $species_order, $species_group, $species_clade, $working_dir ) = @_;
    print localtime()." - Updating scaffold config files in working directory\n\n";
    # Update  the rootingOrder config file.
    my $rooting_order_config = "$working_dir/$scaffold/$scaffold.rootingOrder.config";
    my $tmp_rooting_order_config = "$working_dir/$scaffold/$scaffold.rootingOrder.config.tmp";
    open (IN, $rooting_order_config) or die "Can't open $rooting_order_config\n";
    open (OUT, ">$tmp_rooting_order_config") or die "Can't open $tmp_rooting_order_config file\n";
    my $inserted = 0;
    while (<IN>) {
        # Example line: Physcomitrella patens=Phypa
        chomp;
        print OUT "$_\n";
        if (not $inserted) {
            if (not /^#/ && not /^$/) {
                my @F=split(/=/, $_);
                if ($F[1] eq $rooting_order_species_code) {
                    print OUT "$species_name=$species_code\n";
                    $inserted = 1;
                }
            }
        }
    }
    close OUT;
    close IN;
    my $update_rooting_order = system "mv $tmp_rooting_order_config $rooting_order_config >/dev/null";
    if ($update_rooting_order != 0) {
        print "\n-- ".localtime()." - Updating rooting order config in working directory failed\n\n"; exit(1);
    }
    # update  the taxaLineage config file.
    open (OUT, ">>$working_dir/$scaffold/$scaffold.taxaLineage.config") or die "Can't open $working_dir/$scaffold/$scaffold.taxaLineage.config file\n";
    print OUT "$species_name\t$species_family\t$species_order\t$species_group\t$species_clade\n";
    close OUT;

}

sub sort_sequences {
    my ( $proteins, $coding_sequences, $scaffold, $method, $num_threads, $species_name, $working_dir, $scaffold_dir ) = @_;
    my $method_dir = "$working_dir/$method";
    make_directory($method_dir);
    print localtime()." - Sorting working directory protein sequences in the $method clustering method\n";
    print "-- ".localtime()." - Running BLASTP\n\n";
    my $blastp_call = system "blastp -outfmt 6 -evalue 1e-5 -num_threads $num_threads -query $proteins -db $scaffold_dir/db/blast/$method -out $method_dir/proteins.blastp  >/dev/null";
    if ($blastp_call != 0) { print "\n-- ".localtime()." - Running BLASTP failed\n\n"; exit(1); }
    print "-- ".localtime()." - Getting best BLASTP hits\n\n";
    my $blast_results = "proteins.blastp";
    get_best_blastp_orthos ( $blast_results, $scaffold, $method, $method_dir, $scaffold_dir );
    print "-- ".localtime()." - Running HMMScan\n\n";
    my $hmmscan_call = system "hmmscan -E 1e-5 --cpu $num_threads --noali --tblout $method_dir/proteins.hmmscan -o $method_dir/hmmscan.log $scaffold_dir/db/hmm/$method $proteins >/dev/null";
    if ($hmmscan_call != 0) { print "\n-- ".localtime()." - Running HMMScan failed\n\n"; exit(1); }
    print "-- ".localtime()." - Getting best HMMScan hits\n\n";
    my $hmmscan_results = "proteins.hmmscan";
    get_best_hmmscan_orthos ( $hmmscan_results, $scaffold, $method, $method_dir, $scaffold_dir );
    get_blast_hmmscan_orthos ( $scaffold, $method, $method_dir, $scaffold_dir );
    get_orthogroup_fasta ( $proteins, $coding_sequences, $method, $method_dir, $scaffold_dir );
    update_scaffold_data( $proteins, $scaffold, $method, $num_threads, $method_dir, $species_name, $working_dir, $scaffold_dir );
}

sub  get_best_blastp_orthos {
    my ( $blast_results, $scaffold, $method, $method_dir, $scaffold_dir ) = @_;
    my (%best, %max, %list);
    open (IN, "$method_dir/$blast_results") or die "Can't open $method_dir/$blast_results file\n";
    while (<IN>) {
        chomp;
        my @F=split(/\t/, $_);
        if ($F[0] eq $F[1]) { next; }
        if (!$best{$F[0]}) {
            $best{$F[0]} = $_;
            $max{$F[0]} = $F[11];
        }
        else {
            if ($F[11] > $max{$F[0]}) {
                $best{$F[0]} = $_;
                $max{$F[0]} = $F[11];
            }
        }
    }
    close IN;
    open (IN, "$scaffold_dir/annot/$method.list") or die "Can't open $scaffold_dir/annot/$method.list file\n";
    while (<IN>) {
        chomp;
        my @F=split(/\t/, $_);
        $list{$F[1]} = $F[0];
    }
    close IN;
    open (OUT, ">$method_dir/$blast_results.bestOrthos") or die "Can't open $method_dir/$blast_results.bestOrthos file\n";
    print OUT "Gene ID\tOrthogroup ID\n";
    foreach (keys %best) {
        my @F = split(/\t/, $best{$_});
        print OUT "$F[0]\t$list{$F[1]}\n";
    }
    close OUT;
}

sub  get_best_hmmscan_orthos {
    my ( $hmmscan_results, $scaffold, $method, $method_dir, $scaffold_dir ) = @_;
    my %hits;
    open (IN, "$method_dir/$hmmscan_results") or die "Can't open $method_dir/$hmmscan_results file\n";
    while (<IN>) {
        if (/^#/){next;}
        my @F = split(/\s+/, $_);
        $hits{$F[2]}{$F[0]} = $F[5];
    }
    close IN;
    open (OUT, ">$method_dir/$hmmscan_results.bestOrthos") or die "Can't open $method_dir/$hmmscan_results.bestOrthos file\n";
    print OUT "Gene ID\tOrthogroup ID\n";
    for my $hit (keys %hits) {
        my $score = 0;
        my $best_target;
        for my $target (keys %{$hits{$hit}}) {
            if ($hits{$hit}{$target} >= $score) {
                $score = $hits{$hit}{$target};
                $best_target = $target;
            }
        }
        print OUT "$hit\t$best_target\n";
     }
    close OUT;
}

sub get_blast_hmmscan_orthos {
    my ( $scaffold, $method, $method_dir, $scaffold_dir ) = @_;
    my (%blastp, %hmmscan, %genes);
    opendir (DIR, "$method_dir") or die "Can't open $method_dir directory\n";
    while (my $filename = readdir(DIR)) {
        if ($filename =~ /^proteins\.blastp\.bestOrthos$/){
            open (IN, "$method_dir/$filename") or die "Can't open $method_dir/$filename file\n";
            while (<IN>) {
                chomp;
                if (/^Gene/) {next;}
                my @F = split(/\t/, $_);
                $blastp{$F[0]} = $F[1];
                $genes{$F[0]} = $F[0];
            }
            close IN;
        }
        if ($filename =~ /^proteins\.hmmscan\.bestOrthos$/){
            open (IN, "$method_dir/$filename") or die "Can't open $method_dir/$filename file\n";
            while (<IN>) {
                chomp;
                if (/^Gene/) {next;}
                my @F = split(/\t/, $_);
                $hmmscan{$F[0]} = $F[1];
                $genes{$F[0]} = $F[0];
            }
            close IN;
        }
    }
    closedir DIR;
    open (OUT, ">$method_dir/proteins.both.bestOrthos") or die "Can't open $method_dir/protein.both.bestOrthos file\n";
    print OUT "Gene ID\tOrthogroup ID\n";
    foreach (sort keys %genes) {
        if (!$blastp{$_} and $hmmscan{$_}) { print OUT "$_\t$hmmscan{$_}\n"; next; }
        elsif ($blastp{$_} and !$hmmscan{$_}) { print OUT "$_\t$blastp{$_}\n"; next; }
        elsif ($blastp{$_} == $hmmscan{$_}) { print OUT "$_\t$blastp{$_}\n"; next }
        else { print OUT "$_\t$hmmscan{$_}\n"; }
    }
    close OUT;
}

sub get_orthogroup_fasta {
    my ( $proteins, $coding_sequences, $method, $method_dir, $scaffold_dir ) = @_;
    print "-- ".localtime()." - Retrieving orthogroup fasta files\n\n";
    my (%orthos, %pep, %cds);
    my $orthogroup_fasta = "$method_dir/orthogroups_fasta";
    make_directory($orthogroup_fasta);
    my $orthogroup_assignment = "proteins.both.bestOrthos";
    open (IN, "$method_dir/$orthogroup_assignment") or die "Can't open $method_dir/$orthogroup_assignment file\n";
    while (<IN>) {
        chomp;
        if ($_ =~ /^Gene/) { next; }
        my @F = split(/\t/, $_);
        $orthos{$F[1]}{$F[0]} = $F[0];
    }
    close IN;
    %pep = get_sequences ($proteins);
    if ($coding_sequences) { %cds = get_sequences ($coding_sequences); }
    my ($ortho_id, $seq_id);
    foreach $ortho_id (keys %orthos) {
        open (PEP, ">$orthogroup_fasta/$ortho_id.faa") or die "Can't open $orthogroup_fasta/$ortho_id.faa file\n";
        if ($coding_sequences) { open (CDS, ">$orthogroup_fasta/$ortho_id.fna") or die "Can't open $orthogroup_fasta/$ortho_id.fna file\n"; }
        foreach $seq_id (sort keys %{$orthos{$ortho_id}}) {
            $pep{$seq_id} =~ s/.{80}(?=.)/$&\n/g;
            print PEP ">$seq_id\n$pep{$seq_id}\n";
            if ($coding_sequences) {
                $cds{$seq_id} =~ s/.{80}(?=.)/$&\n/g;
                print CDS ">$seq_id\n$cds{$seq_id}\n";
            }
        }
        close PEP;
        close CDS;
    }
    my @files;
    my $formated_fasta = "$orthogroup_fasta/formated_fasta";
    make_directory($formated_fasta);
    opendir(DIR, "$orthogroup_fasta") or die "Can't open $orthogroup_fasta directory\n";
    while(my $filename = readdir(DIR)) {
        if($filename =~ /\d+\.fna/ or $filename =~ /\d+\.faa/){
            push (@files, $filename);
        }
    }
    closedir(DIR);
    foreach my $file (@files) {
        open (IN, "$orthogroup_fasta/$file") or die "Can't open $orthogroup_fasta/$file file\n";
        open (OUT, ">$formated_fasta/$file") or die "Can't open $formated_fasta/$file file\n";
        while (<IN>) {
            chomp;
            if(/^>/){ s/\|/_/g; print OUT "$_\n"; }
            else { print OUT "$_\n"; }
        }
        close IN;
        close OUT;
    }
    my $integrated_orthogroup_fasta = "$method_dir/integrated_orthogroup_fasta";
    make_directory($integrated_orthogroup_fasta);
    integrate_orthogroup_fasta ( $formated_fasta, $method, $integrated_orthogroup_fasta, $scaffold_dir );
}

sub get_sequences {
    my ( $file ) = @_;
    my (%sequences, $id);
    open (IN, "$file") or die "Can't open $file file\n";
    while (<IN>) {
        if ($_ =~ />(\S+)/){ $id = $1; }
        else { s/\s+//g; $sequences{$id} .= $_; }
    }
    close IN;
    return %sequences;
}

sub integrate_orthogroup_fasta {
    my ( $formated_fasta, $method, $integrated_orthogroup_fasta, $scaffold_dir) = @_;
    print localtime()." - Integrating orthogroup fasta files\n\n";
    my (%pep, %cds);
    opendir (DIR, "$formated_fasta") or die "Can't open $formated_fasta directory\n";
    while ( my $filename = readdir(DIR) ) {
        if ($filename =~ /^(\d+)\.faa$/) { $pep{$1} = $1; }
        if ($filename =~ /^(\d+)\.fna$/) { $cds{$1} = $1; }
    }
    closedir DIR;
    if (keys(%cds) and (keys(%pep) != keys(%cds))) {
        die "Exiting...!\nOrthogroup classification protein and CDS fasta files not equivalent in $formated_fasta directory\n\n";
    }
    foreach my $ortho_id (keys %pep) {
        my $merging_call = system "cat $scaffold_dir/fasta/$method/$ortho_id.faa $formated_fasta/$ortho_id.faa > $integrated_orthogroup_fasta/$ortho_id.faa";
        if ($merging_call != 0) { print "\n-- ".localtime()." - Merging orthogroup $ortho_id failed\n\n"; exit(1); }
        if (keys(%cds) and $cds{$ortho_id}) {
            my $merging_call = system "cat $scaffold_dir/fasta/$method/$ortho_id.fna $formated_fasta/$ortho_id.fna > $integrated_orthogroup_fasta/$ortho_id.fna";
            if ($merging_call != 0) { print "\n-- ".localtime()." - Merging orthogroup $ortho_id failed\n\n"; exit(1); }
        }
    }
}

sub update_scaffold_data {
    my ( $proteins, $scaffold, $method, $num_threads, $method_dir, $species_name, $working_dir, $scaffold_dir ) = @_;
    print localtime()." - Updating scaffold data files\n";
    # update orthogroup annotation files
    print "-- ".localtime()." - Updating orthogroup annotation files\n\n";
    my %annot;
    open (OUT, ">>$working_dir/$scaffold/annot/$method.list") or die "Can't open $working_dir/$scaffold/annot/$method.list file\n";
    opendir (DIR, "$method_dir/orthogroups_fasta") or die "Can't open $method_dir/orthogroups_fasta directory\n";
    while ( my $filename = readdir(DIR) ) {
        my $seq_count = 0;
        my $ortho_id;
        if ($filename =~ /^(\d+)\.faa$/) {
            $ortho_id = $1;
            open (IN, "$method_dir/orthogroups_fasta/$filename") or die "Can't open $method_dir/orthogroups_fasta/$filename file\n";
            while (<IN>) {
                chomp;
                if (/^>(\S+)/) {
                    print OUT "$ortho_id\t$1\n";
                    $seq_count++;
                }
            }
            close IN;
        }
        else { next; }
        $annot{$ortho_id} =  $seq_count;
    }
    closedir DIR;
    close OUT;
    my ($fields, %avg_summary, %min_summary);
    if (File::Spec->file_name_is_absolute($proteins)) { $proteins = basename($proteins); }
    open (IN, "$working_dir/$scaffold/annot/$method.avg_evalue.summary") or die "Can't open $working_dir/$scaffold/annot/$method.avg_evalue.summary file\n";
    while (<IN>) {
        chomp;
        if (/^Orthogroup\s+ID\s+(.*)/) { $fields = "Orthogroup ID\t$species_name\t$1"; next; }
        else { /(\d+)\s+(.*)/; $avg_summary{$1} = $2; }
    }
    close IN;
    open (OUT, ">$working_dir/$scaffold/annot/$method.avg_evalue.summary") or die "Can't open $working_dir/$scaffold/annot/$method.avg_evalue.summary file\n";
    print OUT "$fields\n";
    foreach (sort {$a<=>$b} keys %avg_summary) {
        if ($annot{$_}) { print OUT "$_\t$annot{$_}\t$avg_summary{$_}\n"; }
        else { print OUT "$_\t0\t$avg_summary{$_}\n"; }
    }
    close OUT;
    open (IN, "$working_dir/$scaffold/annot/$method.min_evalue.summary") or die "Can't open $working_dir/$scaffold/annot/$method.min_evalue.summary file\n";
    while (<IN>) {
        chomp;
        if (/^Orthogroup\s+ID\s+(.*)/) { $fields = "Orthogroup ID\t$species_name\t$1"; next; }
        else { /(\d+)\s+(.*)/; $min_summary{$1} = $2; }
    }
    close IN;
    open (OUT, ">$working_dir/$scaffold/annot/$method.min_evalue.summary") or die "Can't open $working_dir/$scaffold/annot/$method.min_evalue.summary file\n";
    print OUT "$fields\n";
    foreach (sort {$a<=>$b} keys %min_summary) {
        if ($annot{$_}) { print OUT "$_\t$annot{$_}\t$min_summary{$_}\n"; }
        else { print OUT "$_\t0\t$min_summary{$_}\n"; }
    }
    close OUT;

    # update orthogroup fasta files
    print "-- ".localtime()." - Updating orthogroup fasta files\n\n";
    opendir (DIR, "$method_dir/integrated_orthogroup_fasta") or die "Can't open $method_dir/integrated_orthogroup_fasta directory\n";
    while ( my $filename = readdir(DIR) ) {
        if (($filename =~ /^(\d+)\.faa$/) or ($filename =~ /^(\d+)\.fna$/)) {
            my $update_orthogroup_fasta = system "cp $method_dir/integrated_orthogroup_fasta/$filename $working_dir/$scaffold/fasta/$method/$filename >/dev/null";
            if ($update_orthogroup_fasta != 0) { print "\n-- ".localtime()." - Updating orthogroup fasta $filename failed\n\n"; exit(1); }
        }
    }
    close IN;
    closedir DIR;

    # update orthogroup alignments
    print "-- ".localtime()." - Updating alignments\n\n";
    opendir (DIR, "$method_dir/orthogroups_fasta/formated_fasta") or die "Can't open $method_dir/orthogroups_fasta/formated_fasta directory\n";
    while ( my $filename = readdir(DIR) ) {
        if ($filename =~ /^(\d+)\.faa$/) {
            my $ortho_id = $1;
            my $align_fasta = system "mafft --thread $num_threads --add $method_dir/orthogroups_fasta/formated_fasta/$filename $working_dir/$scaffold/alns/$method/$ortho_id.aln > $method_dir/orthogroups_fasta/formated_fasta/$ortho_id.aln 2>/dev/null";
            if ($align_fasta != 0) { print "\n-- ".localtime()." - aligning orthogroup fasta file $filename failed\n\n"; exit(1); }
            my $update_orthogroup_alignment = system "mv $method_dir/orthogroups_fasta/formated_fasta/$ortho_id.aln $working_dir/$scaffold/alns/$method/$ortho_id.aln >/dev/null";
            if ($update_orthogroup_alignment != 0) { print "\n-- ".localtime()." - Updating orthogroup alignment $ortho_id.aln failed\n\n"; exit(1); }
        }
    }
    closedir DIR;

    # update orthogroup hmm profiles
    print "-- ".localtime()." - Updating hmm profiles\n\n";
    opendir (DIR, "$working_dir/$scaffold/alns/$method") or die "Can't open $working_dir/$scaffold/alns/$method directory\n";
    while ( my $filename = readdir(DIR) ) {
        if ($filename =~ /^(\d+)\.aln$/) {
            my $ortho_id = $1;
            my $update_orthogroup_hmm = system "hmmbuild -n $ortho_id --amino --cpu $num_threads $working_dir/$ortho_id.hmm $working_dir/$scaffold/alns/$method/$ortho_id.aln >/dev/null";
            if ($update_orthogroup_hmm != 0) { print "\n-- ".localtime()." - Updating orthogroup hmm profile $ortho_id.hmm failed\n\n"; exit(1); }
            my $convert_hmm_format = system "hmmconvert $working_dir/$ortho_id.hmm > $working_dir/$scaffold/hmms/$method/$ortho_id.hmm";
            if ($convert_hmm_format != 0) { print "\n-- ".localtime()." - Converting orthogroup hmm profile $ortho_id.hmm format failed\n\n"; exit(1); }
            my $remove_tmp_file = system "rm  $working_dir/$ortho_id.hmm >/dev/null";
            if ($remove_tmp_file != 0) { print "\n-- ".localtime()." - Could not remove temporary hmm profile $ortho_id.hmm\n\n"; exit(1); }
        }
    }

    # update orthogroup blast and hmm databases
    print "-- ".localtime()." - Updating blast and hmm databases\n\n";
    my $update_blast_database = system "cat $method_dir/orthogroups_fasta/*.faa >> $working_dir/$scaffold/db/blast/$method";
    if ($update_blast_database != 0) { print "\n-- ".localtime()." - Updating blast database failed\n\n"; exit(1); }
    my $index_blast_database = system "makeblastdb -in $working_dir/$scaffold/db/blast/$method -dbtype prot >/dev/null";
    if ($index_blast_database != 0) { print "\n-- ".localtime()." - Indexing blast database failed\n\n"; exit(1); }
    my $update_hmm_database = system "cat $working_dir/$scaffold/hmms/$method/*.hmm > $working_dir/$scaffold/db/hmm/$method";
    if ($update_hmm_database != 0) { print "\n-- ".localtime()." - Updating hmm database failed\n\n"; exit(1); }
    my $index_hmm_database = system "hmmpress -f $working_dir/$scaffold/db/hmm/$method >/dev/null";
    if ($index_hmm_database != 0) { print "\n-- ".localtime()." - Indexing hmm database failed\n\n"; exit(1); }
}

sub update_database_tables {
    my ( $database_connection_string, $proteins, $scaffold, $methods, $species_name, $species_family, $species_order, $species_group, $species_clade, $scaffold_dir, $working_dir ) = @_;
    print localtime()." - Updating for database tables \n";
    my ( $species_code, %method_genes, %gene_sequences, %dna, %aa );
    # database connection and variables, the format of database_connection_string is
    # postgresql://<user>:<password>@<host>/<database name>
    my @conn_part = split(/:\/\//, $database_connection_string);
    my $conn_part_str = $conn_part[1];
    my $driver = "Pg";
    my @conn_part2 = split(/\//, $conn_part_str);
    my $database = $conn_part2[1];
    my $dsn = "DBI:$driver:dbname = $database;host = 127.0.0.1;port = 5432";
    @conn_part2 = split(/:/, $conn_part_str);
    my $userid = $conn_part2[0];
    @conn_part2 = split(/@/, $conn_part_str);
    my $conn_part2_str = $conn_part2[0];
    my @conn_part3 = split(/:/, $conn_part2_str);
    my $password = $conn_part3[1];
    my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) or die "$DBI::errstr\nError : Unable to open database galaxy_plant_tribes\n";
    print "Successfully connected to database $database\n\n";
    my $gsot_association_prep_file = "$working_dir/gene_scaffold_orthogroup_taxon_association.tsv";

    # Output a prep file that stores information for updating
    # the gene_scaffold_orthogroup_taxon_association table.
    open (ASSOC, ">$gsot_association_prep_file") or die "Can't open $gsot_association_prep_file file\n";
    print ASSOC "gene_id\tscaffold_id\tclustering_method\torthogroup_id\tspecies_name\n";
    # get new species name and code
    if (File::Spec->file_name_is_absolute($proteins)) { $proteins = basename($proteins); }
    $species_name =~ s/\_/ /g;
    open(IN, "$scaffold_dir/$scaffold.rootingOrder.config") or die "Can't open $scaffold_dir/$scaffold.rootingOrder.config file\n";
    while (<IN>){
        chomp;
        if (/^\#/ or /^\s+/ or /^\[/){ next; }
        if (/(\w+\s+\w+)\=(\w+)/) { if ($species_name eq $1) { $species_code = $2;} }
    }
    close IN;
    foreach my $clustering_method (keys %$methods) {
        # Updating orthogroup database table
        print "-- ".localtime()." - Updating $clustering_method records for the plant_tribes_orthogroup database table\n\n";
        my ( $stmt, $sth, $rv, $scaffold_id );
        $stmt = qq(SELECT id FROM plant_tribes_scaffold WHERE scaffold_id = '$scaffold' AND clustering_method = '$clustering_method';);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute() or die $DBI::errstr;
        if ($rv < 0) { print $DBI::errstr; }
        while (my @row = $sth->fetchrow_array()) {
            $scaffold_id = $row[0];
        }
        my $scaffold_annotation_dir = "$scaffold_dir/annot";
        opendir (DIR, $scaffold_annotation_dir) or die "Can't open $scaffold_annotation_dir directory\n";
        while ( my $filename = readdir(DIR) ) {
            if ($filename =~ /$clustering_method.min_evalue\.summary/) {
                open (IN, "$scaffold_annotation_dir/$filename") or die "Can't open $scaffold_annotation_dir/$filename file\n";
                while (<IN>){
                    chomp;
                    if (/^Orthogroup/){ next; }
                    my @fields = split(/\t/, $_);
                    my $num_species = 0;
                    my $num_genes = 0;
                    $scaffold =~ /(\d+)Gv\d+\.\d+/; # 22Gv1.1
                    my $genomes = $1 + 1;
                    for (1..$genomes){
                        if ($fields[$_] > 0){ $num_species++; }
                        $num_genes += $fields[$_];
                    }
                    $stmt = qq(UPDATE plant_tribes_orthogroup SET num_species = $num_species, num_genes = $num_genes WHERE orthogroup_id = $fields[0] AND scaffold_id = $scaffold_id;);
                    $rv = $dbh->do($stmt) or die $DBI::errstr;
                    if($rv < 0) {
                        print $DBI::errstr;
                    }
                }
                close IN;

            }
            if ($filename =~ /$clustering_method\.list/) {
                open (IN, "$scaffold_annotation_dir/$filename") or die "Can't open $scaffold_annotation_dir/$filename file\n";
                while (<IN>){
                    chomp;
                     my @fields = split(/\t/, $_);
                     my @gene_id = split(/\|/, $fields[1]);
                     if ($gene_id[1] =~ /$species_code/) {
                         $fields[1] =~ s/\|/_/g;
                         $method_genes{$clustering_method}{$fields[0]}{$fields[1]} = $fields[1];
                     }
                }
                close IN;
            }
        }
        close DIR;
        print "$scaffold $clustering_method records successfully updated in the plant_tribes_orthogroup table\n\n";
        # Updating taxon database table
        print "-- ".localtime()." - Inserting $clustering_method records into the plant_tribes_taxon database table\n\n";
        my $num_genes = 0;
        foreach my $ortho_id (keys %{$method_genes{$clustering_method}}){
            foreach (keys %{$method_genes{$clustering_method}{$ortho_id}}){
                $num_genes++;
            }
        }
        open(IN, "$scaffold_dir/$scaffold.taxaLineage.config") or die "Can't open $scaffold_dir/$scaffold.taxaLineage.config file\n";
        while (<IN>){
            chomp;
            my @fields = split(/\t/, $_);
            if ($fields[0] ne $species_name) { next; }
            $stmt = qq(INSERT INTO plant_tribes_taxon (species_name, scaffold_id, num_genes, species_family, species_order, species_group, species_clade) VALUES ('$fields[0]', $scaffold_id, $num_genes, '$fields[1]', '$fields[2]', '$fields[3]', '$fields[4]'));
            $rv = $dbh->do($stmt) or die $DBI::errstr;
        }
        close IN;
        print "$species_name $scaffold $clustering_method record successfully inserted into the plant_tribes_taxon table\n\n";
        my ($dna_id, $aa_id);
        my $orthogroups_fasta_dir = "$working_dir/$clustering_method/orthogroups_fasta/formated_fasta";
        opendir (DIR, $orthogroups_fasta_dir) or die "Can't open $orthogroups_fasta_dir directory\n";
        while ( my $filename = readdir(DIR) ) {
            if ($filename =~ /^(\d+)\.fna$/) {
                my $ortho_id = $1;
                open(IN, "$orthogroups_fasta_dir/$filename") or die "Can't open $orthogroups_fasta_dir/$filename file\n";
                while(<IN>){
                    chomp;
                    if (/^>(\S+)/){
                        $dna_id = $1;
                        print ASSOC "$dna_id\t$scaffold\t$clustering_method\t$ortho_id\t$species_name\n";
                        next;
                    }
                    else { s/\s+//g; $dna{$dna_id} .= $_; }
                }
                close IN;
            }
            if ($filename =~ /^(\d+)\.faa$/) {
                open(IN, "$orthogroups_fasta_dir/$filename") or die "Can't open $orthogroups_fasta_dir/$filename file\n";
                while(<IN>){
                    chomp;
                    if (/^>(\S+)/){ $aa_id = $1; next; }
                    else { s/\s+//g; $aa{$aa_id} .= $_; }
                }
                close IN;
            }
        }
        close DIR;
    }
    close ASSOC;
    # Updating gene database table
    print "-- ".localtime()." - Inserting records into the plant_tribes_gene database table\n\n";
    foreach my $gene_id (sort keys %dna) {
        my $stmt = qq(INSERT INTO plant_tribes_gene (gene_id, dna_sequence, aa_sequence) VALUES ('$gene_id', '$dna{$gene_id}', '$aa{$gene_id}'));
        my $rv = $dbh->do($stmt) or die $DBI::errstr;
    }
    print "$species_name $scaffold records successfully inserted into the plantr_tribes_gene table\n\n";
    # Updaing gene-scaffold-orthogroup-taxon-association database table
    print "-- ".localtime()." - Inserting  records into the gene_scaffold_orthogroup_taxon_association database table\n\n";
    open(IN, "$gsot_association_prep_file") or die "Can't open $gsot_association_prep_file file\n";
    while(<IN>){
        chomp;
        if (/^gene_id/){ next; }
        my @fields = split(/\t/, $_);
        my ( $stmt, $sth, $rv, $scaffold_id, $orthogroup_id, $taxon_id, $gene_id );
        $stmt = qq(SELECT id FROM plant_tribes_scaffold WHERE scaffold_id = '$fields[1]' AND clustering_method = '$fields[2]';);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute() or die $DBI::errstr;
        if ($rv < 0) { print $DBI::errstr; }
        while (my @row = $sth->fetchrow_array()) {
            $scaffold_id = $row[0];
        }
        $stmt = qq(SELECT id FROM plant_tribes_orthogroup WHERE orthogroup_id = $fields[3] AND scaffold_id = '$scaffold_id';);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute() or die $DBI::errstr;
        if ($rv < 0) { print $DBI::errstr; }
        while (my @row = $sth->fetchrow_array()) {
            $orthogroup_id = $row[0];
        }
        $stmt = qq(SELECT id FROM plant_tribes_taxon WHERE species_name = '$species_name' AND scaffold_id = '$scaffold_id';);
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute() or die $DBI::errstr;
        if ($rv < 0) { print $DBI::errstr; }
        while (my @row = $sth->fetchrow_array()) {
            $taxon_id = $row[0];
        }
        $stmt = qq(SELECT id FROM plant_tribes_gene WHERE gene_id = '$fields[0]' );
        $sth = $dbh->prepare( $stmt );
        $rv = $sth->execute() or die $DBI::errstr;
        if ($rv < 0) { print $DBI::errstr; }
        while (my @row = $sth->fetchrow_array()) {
            $gene_id = $row[0];
        }
        $stmt = qq(INSERT INTO gene_scaffold_orthogroup_taxon_association (gene_id, scaffold_id, orthogroup_id, taxon_id) VALUES ($gene_id, $scaffold_id, $orthogroup_id, $taxon_id));
        $rv = $dbh->do($stmt) or die $DBI::errstr;

    }
    close IN;
    print "$species_name $scaffold records successfully inserted into the gene_scaffold_orthogroup_taxon_association table\n\n";
    $dbh->disconnect();
}
