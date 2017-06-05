#!/usr/bin/perl -w

#
# This script runs the Sfold package. Do not call this script
# directly. Call the shell script ./sfold instead.
#

# set auto flush
$| = 1;

use strict;
use Getopt::Std;
use vars qw/ %opt /;
use Cwd qw/ cwd /;

# set additional arguments for the Sfold core program;
# the default is an empty string
my($sfold_extra_arg) = "";
my($sclass_extra_arg) = "";

my($cmdenv) = "";
my($seqfile) = "";
my($outdir) = "";
my($basedir) = &cwd();
my($seqlen);
my($n_sample);

&checkenv();

my($sfoldhdr) = `$cmdenv $ENV{'SFOLDEXE'} | $ENV{'PATH_GREP'} "Sfold Executable Code"`;
chomp($sfoldhdr);

&usage() if (scalar(@ARGV) == 0);
&initarg();

print "\n$sfoldhdr\n";

#
# Prepare to run Sfold partition function and sampling
#
print " Running partition function and statistical sampling module... ";
&run_sfold();
print "done\n";
print " (general output in $outdir/sfold.out)\n";

#
# Prepare to run Sfold clustering unless user chooses to skip it
#
if (defined($opt{a}) && $opt{a}) {
  print " Running structure clustering module... ";
  &run_sclass();
  print "done\n";
  print " (general output in $outdir/sclass.out)\n";
}

print "Process completed. See output files in '$outdir'.\n\n";

exit;


#
# check environment settings to make sure we
# have got everything we need to run Sfold
#
sub checkenv {
  # check hardware architecture and operating system
  die " Error: environment variable 'OSTYPE' not defined"
    if (!defined($ENV{'OSTYPE'}));
  die " Error: environment variable 'MACH' not defined"
    if (!defined($ENV{'MACH'}));
  die " Error: unsupported hardware/OS found"
    if ($ENV{'OSTYPE'} eq "UNKNOWN" || $ENV{'MACH'} eq "UNKNOWN");

  # check Sfold directories
  die " Error: invalid environment variable 'SFOLDDIR'"
    if (!defined($ENV{'SFOLDDIR'}) || !-e $ENV{'SFOLDDIR'});
  die " Error: invalid environment variable 'SFOLDBIN'"
    if (!defined($ENV{'SFOLDBIN'}) || !-e $ENV{'SFOLDBIN'});
  die " Error: invalid environment variable 'SFOLDLIB'"
    if (!defined($ENV{'SFOLDLIB'}) || !-e $ENV{'SFOLDLIB'});
  die " Error: invalid environment variable 'SFOLDPAR'"
    if (!defined($ENV{'SFOLDPAR'}) || !-e $ENV{'SFOLDPAR'});

  # check Sfold binaries
  die " Error: Sfold does not run on your platform"
    if (!defined($ENV{'SFOLDEXE'}) || !-e $ENV{'SFOLDEXE'});
  die " Error: invalid environment variable 'FINDFEEXE'"
    if (!defined($ENV{'FINDFEEXE'}) || !-e $ENV{'FINDFEEXE'});
  die " Error: invalid environment variable 'GETU0EXE'"
    if (!defined($ENV{'GETU0EXE'}) || !-e $ENV{'GETU0EXE'});
  die " Error: invalid environment variable 'SCLASSEXE'"
    if (!defined($ENV{'SCLASSEXE'}) || !-e $ENV{'SCLASSEXE'});
  die " Error: invalid environment variable 'BPROBEXE'"
    if (!defined($ENV{'BPROBEXE'}) || !-e $ENV{'BPROBEXE'});
  die " Error: invalid environment variable 'FINDDISTEXE'"
    if (!defined($ENV{'FINDDISTEXE'}) || !-e $ENV{'FINDDISTEXE'});

  # check external programs needed by Sfold
  die " Error: external program 'awk' not found"
    if (!defined($ENV{'PATH_AWK'}) || !-e $ENV{'PATH_AWK'});
  die " Error: external program 'grep' not found"
    if (!defined($ENV{'PATH_GREP'}) || !-e $ENV{'PATH_GREP'});
  die " Error: external program 'perl' not found"
    if (!defined($ENV{'PATH_PERL'}) || !-e $ENV{'PATH_PERL'});
  die " Error: external program 'R' not found"
    if (!defined($ENV{'PATH_R'}) || !-e $ENV{'PATH_R'});

  return 0;
}


#
# initialize command line arguments processing
#
sub initarg {
  my($str) = 'a:f:hel:m:o:r:w:i:';
  getopts("$str", \%opt) || &usage();
  &usage() if ($opt{h});

  # do checking on arguments that we use in this script;
  # for all other arguments, we just relay them to the
  # Sfold core program for doing the checking
  $opt{a} = 1 if (!defined($opt{a})); # set default value if not defined
  die " Error: invalid '-a' option"
    if (defined($opt{a}) && $opt{a} !~ /^[01]$/);

  $opt{i} = 0 if (!defined($opt{i})); # set default value if not defined
  die " Error: invalid '-i' option"
    if (defined($opt{i}) && $opt{i} !~ /^[0123]$/);


  if (!defined($opt{o})) {
    $outdir = join("", $basedir, "/output");
  } else {
    $outdir = $opt{o};
  }
  if (defined($opt{r})) {
    $sfold_extra_arg .= " -r $opt{r}" if ($opt{r} =~ /^\d+$/);
  }

  die " Error: a file named '$outdir' already exists; \
        unable to create output directory"
    if (-e $outdir && !-d $outdir);

  die " Error: unable to locate MFE structure file '$opt{m}'"
    if (defined($opt{m}) && !-e $opt{m});

  die " Error: sequence file not specified"
    if (!defined($ARGV[0]));
  die " Error: unable to locate sequence file '$ARGV[0]'"
    if (!-e $ARGV[0]);
  $seqfile = $ARGV[0];
}


#
# display help information for this program and then exit
#
sub usage {
  print <<EndOfOutput;

$sfoldhdr

Usage: $0 [options]... sequence_file
Options:
  -a <0 or 1>         Run clustering on the sampled ensemble [default=1]
  -f <string>         Name of file containing folding constraints. Constraint
                      syntax follows what is used in mfold 3.1
                      [default=no constraint]
  -h                  Display this information
  -l <+ve integer>    Maximum distance between paired bases [default=no limit]
  -m <string>         Name of file containing the MFE structure in GCG connect
                      format. If provided, Sfold clustering module will
                      determine the cluster to which this structure belongs.
  -o <string>         Name of directory to which output files are written.
                      Directory will be created if it does not already exist.
                      Existing files will be overwritten.
                      [default=$basedir/output]
  -w <+ve integer>    Length of antisense oligos [default=20]
  -e                  Do not obliterate sample.out to save space [default=do]
  -i <0,1,2,3>        1=do Sirna, 2=do Soligo, 3=both, 0=neither [default=0]

IMPORTANT: Use of this program is restricted to non-commercial internal
           research use only. Please read the enclosed license agreement
           before using this program. For more info, please visit us at:
           http://sfold.wadsworth.org

EndOfOutput

  exit;
}


#
# Run Sfold core program - partition function and sampling
#
sub run_sfold {
  my($sfold_arg) = $sfold_extra_arg;
  my($sfoldout) = "";
  my($scale);

  # set option to specify constraint file
  $sfold_arg .= " -f $opt{f}" if (defined($opt{f}));

  # set local folding length
  $sfold_arg .= " -l $opt{l}" if (defined($opt{l}));

  # set output directory
  $sfold_arg .= " -o $outdir";

  # set parameter directory
  $sfold_arg .= " -p $ENV{'SFOLDPAR'}";

  # determine scaling factor
  if (defined($opt{l})) {
    $scale = `$ENV{'GETU0EXE'} $seqfile $opt{l}`;
  } else {
    $scale = `$ENV{'GETU0EXE'} $seqfile`;
  }
  die " Error: $ENV{'GETU0EXE'} $seqfile $opt{l} failed" if ($?);
  chomp($scale);
  $sfold_arg .= " -s $scale";

  # set antisense oligo length
  $sfold_arg .= " -w $opt{w}" if (defined($opt{w}));
  
  # set verbose option
  $sfold_arg .= " -e" if (defined($opt{e}));

  # set option to enable siRNA module
  $sfold_arg .= " -i $opt{i}" if (defined($opt{i}));

  # create output directory if not already exist
  if (!-e $outdir) {
    mkdir($outdir, 0755) ||
      die " Error: unable to create output directory '$outdir'";
  }
  $sfoldout = join("", $outdir, "/sfold.out");

  system("$cmdenv $ENV{'SFOLDEXE'} $sfold_arg $seqfile > $sfoldout 2>&1");

  my($endofout) = `$ENV{'PATH_GREP'} "Total cpu time" $sfoldout` || "";
  chomp($endofout);

  die "\n Error: $ENV{'SFOLDEXE'} failed. See $sfoldout for more info."
    if ($? || $endofout eq "");

  $seqlen = `$ENV{'PATH_GREP'} "Sequence length n" $sfoldout | $ENV{'PATH_AWK'} '{print \$5}'` || "";
  chomp($seqlen);

  $n_sample = `$ENV{'PATH_GREP'} "Number of structures in a sample" $sfoldout | $ENV{'PATH_AWK'} '{print \$8}'` || "";
  chomp($n_sample);

  return 0;
}


#
# Run Sfold clustering program
#
sub run_sclass {
  my($sclass_arg) = $sclass_extra_arg;
  my($sclassout) = join("", $outdir, "/sclass.out");
  my($stdout) = join("", $outdir, "/sclass.stdout");
  my($clustdir) = join("", $outdir, "/clusters");
  my($chfile) = join("", $clustdir, "/ch.index.out");

  # make sure output directory already exists
  die " Error: output directory '$outdir' does not exist" if (!-e $outdir);

  # prepare directory for storing clustering output
  if (-e $clustdir) {
    system("rm -rf $clustdir > /dev/null 2>&1");
  }
  mkdir($clustdir, 0755) ||
    die " Error: unable to create clusters directory '$clustdir'";

  $sclass_arg .= " -b $outdir/bp.out";
  $sclass_arg .= " -c $clustdir";
  $sclass_arg .= " -e $outdir";
  $sclass_arg .= " -f $outdir/fe.out";
  $sclass_arg .= " -p $outdir/bprob.out";
  $sclass_arg .= " -s $seqfile";
  $sclass_arg .= " -w 1";
  $sclass_arg .= " -m $opt{m}" if (defined($opt{m}));

  system("$cmdenv $ENV{'SCLASSEXE'} $sclass_arg > $stdout 2>&1");
  die " Error: $ENV{'SCLASSEXE'} failed" if ($?);

  # check whether the clustering output is complete...
  my($tmpout) = `$ENV{'PATH_GREP'} "Optimal number of clusters" $stdout` || "";
  chomp($tmpout);
  die " Error: incomplete output from $ENV{'SCLASSEXE'}" if ($tmpout eq "");

  # check whether sclass returns only a single cluster (only when all input
  # structures are identical)
  $tmpout = `$ENV{'PATH_GREP'} "CH index =" $stdout | $ENV{'PATH_AWK'} '{print \$4}'`
              || "";
  chomp($tmpout);
  if ($tmpout !~ /^[\d\.]+[\n\r]*$/) {
    # identical sampled structures
    system("mv -f $stdout $sclassout > /dev/null 2>&1");

    return 0;
  }

  # multiple clusters returned, collect statistics

  # write file for CH index plot
  open(SCLASS_OUT, "$stdout") || die " Error: unable to open '$stdout'";
  open(CHINDEX_OUT, ">$chfile") || die " Error: unable to create '$chfile'";
  while (<SCLASS_OUT>) {
    next if (!/Determining optimal number of clusters.../);

    # started determining the number of clusters
    while (<SCLASS_OUT>) {
      chomp;
      if (/^\s*\d+\s+[\d\.]+\s*$/) {
        print CHINDEX_OUT "$_\n";
      } elsif ($_ ne "") {
        last;
      }
    }

    last;
  }
  close(CHINDEX_OUT);
  close(SCLASS_OUT);

  # load cluster sizes
  my($mfe_clust) = -1;
  if (defined($opt{m})) {
    $tmpout = `$ENV{'PATH_GREP'} 'mfold optimal structure appears in cluster' $stdout` || "";
    $tmpout =~ /mfold optimal structure appears in cluster (\d+)/;
    $mfe_clust = $1;
  }

  my(@clust_sizes) = ();
  open(GREP, "$ENV{'PATH_GREP'} 'Probability' $stdout |")
    || die " Error: unable to grep probabilities";
  $tmpout = 1;
  while (<GREP>) {
    /Probability\s*:\s*\d+\/\d+\s*=\s*([\d\.]+)/;
    my($prob) = $1;

    if (defined($opt{m}) && $tmpout == $mfe_clust) {
      push @clust_sizes, "$prob*";
    } else {
      push @clust_sizes, "$prob ";
    }

    $tmpout++;
  }
  close(GREP);

  # set the optimal number of clusters
  my($nclust) = scalar(@clust_sizes);

  # if the mfe should be by itself in its own cluster, then
  # we add the special entry 0.000 at the end of the array
  # after we have determined the total number of clusters,
  # i.e. nclust does not include this special class of the mfe.
  if (defined($opt{m}) && $mfe_clust == 0) {
    push @clust_sizes, "0.000*";
  }


  # extract blocks of cluster-level info & distance matrices
  # from sclass standard output
  my(@clustinfo) = ();
  my($clustinfoindex) = 0;
  my(@min_dmat, @max_dmat, @avg_dmat, @cc_dmat);
  @min_dmat = @max_dmat = @avg_dmat = @cc_dmat = ();

  my($sclasshdr);
  open(SCLASS_OUT, "$stdout") || die " Error: unable to open '$stdout'";
  while (<SCLASS_OUT>) {
    chomp;

    if (/Sfold Executable Code/) {
      $sclasshdr = <SCLASS_OUT>;
      chomp($sclasshdr);
    } elsif (/^Cluster\s*:\s*\d+$/) {
      $clustinfo[$clustinfoindex] = "$_\n";
      while (<SCLASS_OUT>) {
        chomp;
        last if (/^\s*$/);
        $clustinfo[$clustinfoindex] .= "$_\n";
      }
      $clustinfoindex++;
    } elsif (/Minimum Distance Matrix between Clusters/i) {
      while (<SCLASS_OUT>) {
        chomp;
        last if (/^\s*$/);
        my(@e) = split;
        foreach $tmpout (@e) {
          push @min_dmat, $tmpout;
        }
      }
    } elsif (/Maximum Distance Matrix between Clusters/i) {
      while (<SCLASS_OUT>) {
        chomp;
        last if (/^\s*$/);
        my(@e) = split;
        foreach $tmpout (@e) {
          push @max_dmat, $tmpout;
        }
      }
    } elsif (/Average Distance Matrix between Clusters/i) {
      while (<SCLASS_OUT>) {
        chomp;
        last if (/^\s*$/);
        my(@e) = split;
        foreach $tmpout (@e) {
          push @avg_dmat, $tmpout;
        }
      }
    } elsif (/Distance matrix between cluster centroids/i) {
      while (<SCLASS_OUT>) {
        chomp;
        last if (/^\s*$/);
        my(@e) = split;
        foreach $tmpout (@e) {
          push @cc_dmat, $tmpout;
        }
      }
    }
  }
  close(SCLASS_OUT);

  # compute average distance between sampled structures and representative structure
  $tmpout = $n_sample + 1;
  my($avgdist_mfe, $avgdist_ec);

  if (defined($opt{m})) {
    # compute average distance between sampled structures and MFE structure
    system("$ENV{'FINDDISTEXE'} -b $outdir/bp.out -f $outdir/fe.out -p $outdir/bprob.out -i $tmpout -m $opt{m} > $outdir/bp.dist.from.optimal.out");
    die " Error: $ENV{'FINDDISTEXE'} failed" if ($?);
    $avgdist_mfe = &find_avg_dist("$outdir/bp.dist.from.optimal.out");
  }

  # compute average distance between sampled structures and ensemble centroid
  system("$ENV{'FINDDISTEXE'} -b $outdir/bp.out -f $outdir/fe.out -p $outdir/bprob.out -i $tmpout -m $outdir/ecentroid.ct > $outdir/bp.dist.from.ecentroid.out");
  die " Error: $ENV{'FINDDISTEXE'} failed" if ($?);
  $avgdist_ec = &find_avg_dist("$outdir/bp.dist.from.ecentroid.out");


  # prepare statistics for each cluster
  for (my $i=1; $i<=$nclust; $i++) {
    my($maxclust) = 20;
    my($jid) = join("", "0"x(length($maxclust) - length($i)), $i);

    my($clist) = join("", $clustdir, "/c", $jid, ".list");
    my($cbp) = join("", $clustdir, "/c", $jid, ".bp.out");
    my($cfe) = join("", $clustdir, "/c", $jid, ".fe.out");
    my($c2dhist) = join("", $clustdir, "/c", $jid, ".2dhist.out");
    my($ccenbp) = join("", $clustdir, "/c", $jid, ".ccentroid.bp");
    my($ccenfe) = join("", $clustdir, "/c", $jid, ".ccentroid.fe");
    my($ccenct) = join("", $clustdir, "/c", $jid, ".ccentroid.ct");
    my($ccenbpdist) = join("", $clustdir, "/c", $jid, ".bp.dist.from.ccentroid.out");
    my($ccenavgdist) = join("", $clustdir, "/c", $jid, ".avg.dist.from.ccentroid.out");
    my($cprob) = $clust_sizes[$i-1];
    $cprob =~ s/[^\d\.]//g;

    if ($cprob != 0) {
      my($bid) = $cprob*$n_sample + 1;
      system("$ENV{'SFOLDBIN'}/bp_count.pl $cbp > $c2dhist");

      system("$ENV{'FINDDISTEXE'} -b $cbp -f $cfe -p $outdir/bprob.out -i $bid -m $ccenct > $ccenbpdist");
      die " Error: $ENV{'FINDDISTEXE'} failed" if ($?);
#       = &find_avg_dist("$ccenbpdist");
#      system("$ENV{'SFOLDBIN'}/find.avg.dist.pl $ccenbpdist > $ccenavgdist");
    }
  }


  # finally, output results to sclass.out
  open(OUTFILE, ">$sclassout") || die " Error: unable to create '$sclassout'";

  # print Sfold version header
  print OUTFILE "\n$sfoldhdr\n";
  print OUTFILE "$sclasshdr\n\n";

  print OUTFILE "Sequence length            = $seqlen\n";

  $tmpout = `$ENV{'PATH_GREP'} 'Optimal number of clusters' $stdout`;
  print OUTFILE "$tmpout";

  print OUTFILE "Cluster sizes              = @clust_sizes\n";
  print OUTFILE "(* marks the cluster in which the MFE structure is located)\n"
    if (defined($opt{m}));
  print OUTFILE "\n";

  if (defined($opt{m})) {
    $tmpout = `$ENV{'PATH_GREP'} initially $opt{m}`;
    $tmpout =~ /\[initially\s+([\-\d\.]+)\]/;

    print OUTFILE "Free energy of the MFE structure     = ",
      sprintf("%.2f", $1), "\n";
  }

  $tmpout = `$ENV{'PATH_GREP'} initially $outdir/ecentroid.ct`;
  $tmpout =~ /\[initially\s+([\-\d\.]+)\]/;
  my($ec_energy) = $1;
  print OUTFILE "Free energy of the ensemble centroid = ",
    sprintf("%.2f", $ec_energy), "\n";

  if (defined($opt{m})) {
    $tmpout = `$ENV{'PATH_GREP'} 'Boltzmann probability' $stdout`;
    $tmpout =~ /Boltzmann probability of mfold optimal structure\s*=\s*([^\n]+)/;
    print OUTFILE "Boltzmann probability of the MFE structure     = $1\n";
  }

  my($ec_bprob) = `$ENV{'BPROBEXE'} $outdir/bprob.out $ec_energy` || "";
  chomp($ec_bprob);
  print OUTFILE "Boltzmann probability of the ensemble centroid = $ec_bprob\n";
  print OUTFILE "\n";

  print OUTFILE "Average base pair distance between sample and\n";
  print OUTFILE ".....MFE structure     = $avgdist_mfe\n" if (defined($opt{m}));
  print OUTFILE ".....ensemble centroid = $avgdist_ec\n";

  if (defined($opt{m})) {
    my($improve_p);
    if ($avgdist_mfe != 0) {
      $improve_p = 0 - ($avgdist_ec - $avgdist_mfe) * 100.0 / $avgdist_mfe;
    } else {
      $improve_p = 0;
    }
    print OUTFILE "Percentage improvement = ",
      sprintf("%.4f", $improve_p), "\n";
  }
  print OUTFILE "\n";
  print OUTFILE "\n";

  # print cluster-level info to file
  for (my $i=1; $i<=$nclust; $i++) {
    print OUTFILE "$clustinfo[$i-1]\n";
  }
  print OUTFILE "\n";

  print OUTFILE "Minimum distance matrix between clusters\n";
  my($spacing) = join("", "%", length($seqlen)+2, "s");
  for (my $i=1; $i<=$nclust; $i++) {
    for (my $j=1; $j<=$nclust; $j++) {
      my($val) = $min_dmat[($i-1)*$nclust + $j - 1];
      $val = sprintf("$spacing", $val);
      print OUTFILE "$val";
    }

    print OUTFILE "\n";
  }
  print OUTFILE "\n";


  print OUTFILE "Maximum distance matrix between clusters\n";
  $spacing = join("", "%", length($seqlen)+2, "s");
  for (my $i=1; $i<=$nclust; $i++) {
    for (my $j=1; $j<=$nclust; $j++) {
      my($val) = $max_dmat[($i-1)*$nclust + $j - 1];
      $val = sprintf("$spacing", $val);
      print OUTFILE "$val";
    }

    print OUTFILE "\n";
  }
  print OUTFILE "\n";

  print OUTFILE "Average distance matrix between clusters\n";
  $spacing = join("", "%", length($seqlen)+5, "s");
  for (my $i=1; $i<=$nclust; $i++) {
    for (my $j=1; $j<=$nclust; $j++) {
      my($val) = $avg_dmat[($i-1)*$nclust + $j - 1];
      $val = sprintf("$spacing", $val);
      print OUTFILE "$val";
    }

    print OUTFILE "\n";
  }
  print OUTFILE "\n";


  # print out distances among cluster centroids
  print OUTFILE "Distance matrix between cluster centroids\n";
  $spacing = join("", "%", length($seqlen)+2, "s");
  for (my $i=1; $i<=$nclust; $i++) {
    for (my $j=1; $j<=$nclust; $j++) {
      my($val) = $cc_dmat[($i-1)*$nclust + $j - 1];
      $val = sprintf("$spacing", $val);
      print OUTFILE "$val";
    }

    print OUTFILE "\n";
  }
  print OUTFILE "\n";

  close(OUTFILE);

  # remove intermediate sclass output
  unlink "$stdout";

  return 0;
}


# compute average base-pair distance from find.dist output
sub find_avg_dist {
  my($fname) = @_;
  my($sum) = 0;
  my($count) = 0;

  die " Error: unable to locate file '$fname'" if (!-e $fname);
  open(FINDDISTOUT, "$fname") || die " Error: unable to open '$fname'";
  while (<FINDDISTOUT>) {
    next if (!/^\s*\d+\s+\d+\s*$/);
    my($sid, $d) = split;

    $sum += $d;
    $count++;
  }
  close(FINDDISTOUT);

  die " Error: '$fname' not in valid format" if ($count < 1);
  if ($count == 1) {
    return 0;
  } else {
    # do not count the last line from find.dist output
    return ($sum / ($count-1));
  }
}
