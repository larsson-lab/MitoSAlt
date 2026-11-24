#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min max);
use File::Path qw(make_path);
use FindBin qw($Bin);
use Cwd qw(getcwd);
use Sys::Hostname qw(hostname);

#VERSION FILE NEXT TO THE SCRIPT; BUILDS MAY APPEND THE COMMIT AS SEMVER BUILD METADATA (E.G. 1.1.2+abc1234)
my $VERSION = 'unknown';
if(open(my $vf,'<',"$Bin/VERSION")){
  my $line = <$vf>;
  close($vf) or undef $line;
  if(defined $line){chomp($line); $VERSION = $line if length $line;}
}

#RUN SUMMARY log/<study name>.summary.tsv, WRITTEN ON SUCCESS AND ON FAILURE
my ($summary_file,%summary);
my @summary_keys = qw(version sample mt_reads mt_alignments mt_mean_depth split_reads clusters events_reported warnings);
#WRITTEN TO A TEMPORARY FILE AND RENAMED, SO A FAILED WRITE NEVER LEAVES A TRUNCATED SUMMARY; RETURNS FALSE ON FAILURE
sub write_summary{
  my ($status,$code,$msg) = @_;
  return 1 unless $summary_file;
  my $tmp = "$summary_file.tmp";
  open(my $fh,'>',$tmp) or return 0;
  print $fh "status\t$status\n";
  print $fh "exit_code\t$code\n" if defined $code;
  if(defined $msg){$msg =~ s/:\s*\n\s*/: /; $msg =~ s/\s*\n\s*/; /g; print $fh "message\t$msg\n";}
  print $fh "$_\t$summary{$_}\n" for grep {defined $summary{$_}} @summary_keys;
  close($fh) or return 0;
  return rename($tmp,$summary_file);
}

#EXIT CODES: 1 = INVALID INPUT OR CONFIGURATION, 2 = A PIPELINE STEP FAILED,
#3 = NO MT READS OR ALIGNMENTS LEFT AFTER A STEP (BAD INPUT OR REFERENCE, NOT A NEGATIVE RESULT)
sub abort{
  my ($code,$msg) = @_;
  #SUMMARY FIRST: IF THE LOG PIPE IS GONE (E.G. TEE KILLED WITH THE JOB) THE PRINT CAN FAIL
  write_summary('failed',$code,$msg);
  local $SIG{PIPE} = 'IGNORE';
  print STDERR scalar(localtime).": ERROR: $msg\n";
  exit $code;
}

#CONFIG TOOL ENTRIES CAN BE COMMANDS WITH ARGUMENTS, E.G. "bedtools sort"
sub find_executable{
  my ($cmd) = split(' ',shift);
  return -f $cmd && -x $cmd if $cmd =~ m{/};
  return grep {-f "$_/$cmd" && -x "$_/$cmd"} split(/:/,$ENV{PATH});
}

#PIPEFAIL SO THAT A FAILING COMMAND INSIDE A PIPE IS NOT MASKED BY THE LAST ONE
sub run{
  my ($cmd) = @_;
  #HOST SHELL STARTUP FILES (E.G. LMOD VIA BASH_ENV UNDER SINGULARITY) MUST NOT RUN FOR EACH COMMAND
  delete local @ENV{qw(BASH_ENV ENV)};
  my $rc = system('bash','-o','pipefail','-c',$cmd);
  abort(2,"Command failed (".cmd_status($rc)."): $cmd") if $rc;
  return;
}

#LIKE run() BUT RETURNS THE COMMAND OUTPUT
sub capture{
  my ($cmd) = @_;
  delete local @ENV{qw(BASH_ENV ENV)};
  open(my $fh,'-|','bash','-o','pipefail','-c',$cmd) or abort(2,"Command failed (could not start: $!): $cmd");
  my $out = do {local $/; <$fh>};
  close($fh) or abort(2,"Command failed (".cmd_status($?)."): $cmd");
  return defined $out ? $out : '';
}

sub cmd_status{
  my ($rc) = @_;
  return $rc == -1 ? "could not start: $!" : ($rc & 127) ? "killed by signal ".($rc & 127) : "exit code ".($rc >> 8);
}

sub count_lines{
  my ($file) = @_;
  open(my $fh,'<',$file) or abort(2,"Cannot open $file: $!");
  my $n = 0;
  $n++ while <$fh>;
  close($fh) or abort(2,"Cannot read $file: $!");
  return $n;
}

sub fai_contig_length{
  my ($fai,$contig) = @_;
  open(my $fh,'<',$fai) or return 0;
  while(<$fh>){my @f = split /\t/; return $f[1] if $f[0] eq $contig;}
  return 0;
}

if(@ARGV && $ARGV[0] eq '--version'){print "MitoSAlt $VERSION\n"; exit 0;}

#INPUT
my $config_file = $ARGV[0];
my $p1 = $ARGV[1];
my $tag = $ARGV[2];

my $usage = "usage: perl $0 <config_file> <fastq file> <study name>";
abort(1,"Configuration file missing, $usage") unless $config_file;
abort(1,"Read file missing, $usage") unless $p1;
abort(1,"Study name not given, $usage") unless $tag;
#THE STUDY NAME IS USED UNQUOTED IN SHELL COMMANDS AND FILE NAMES
abort(1,"Study name '$tag' may only contain letters, digits, '.', '_' and '-'") unless $tag =~ /^[\w.-]+$/;

make_path(qw(bam bw indel log plot tab));
$summary_file = "log/$tag.summary.tsv";
%summary = (version => $VERSION, sample => $tag, warnings => 0);
#A SUMMARY LEFT AT "running" MEANS THE PROCESS WAS KILLED (E.G. SIGKILL BY THE OOM KILLER)
write_summary('running') or abort(1,"Cannot write $summary_file: $!");

#ONE LOG PER RUN: OUR MESSAGES AND THE STDOUT/STDERR OF EVERY TOOL GO TO log/<study name>.log AND THE TERMINAL
open(STDOUT,'|-','tee','-i',"log/$tag.log") or abort(1,"Cannot open log/$tag.log: $!");
open(STDERR,'>&',\*STDOUT) or abort(1,"Cannot redirect STDERR to log: $!");
$| = 1;
print scalar(localtime).": MitoSAlt $VERSION\n";
print "Command: $0 @ARGV\n";
print "Host: ".hostname()."\n";
print "Working directory: ".getcwd()."\n";
$SIG{__WARN__} = sub {$summary{warnings}++; print STDERR scalar(localtime).": WARNING: $_[0]";};

#PERL ERRORS AND TERMINATION SIGNALS ARE REPORTED AS FAILED STEPS
$SIG{__DIE__} = sub {return if $^S; my $m = shift; chomp($m); abort(2,"Perl error: $m");};
$SIG{$_} = sub {abort(2,"Terminated by SIG$_[0]");} for qw(INT TERM HUP);

abort(1,"Configuration file not found: $config_file") unless -f $config_file;
abort(1,"FASTQ file not found or empty: $p1") unless -s $p1;

#LOAD CONFIGURATION FILE
open(CONFIG, "<$config_file") or abort(1,"Cannot read configuration file $config_file: $!");
my $User_Preferences;
while (<CONFIG>) {
    chomp;                  # no newline
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences->{$var} = $value;
    my $line = $_;
    print $line."\n";
}
close(CONFIG) or abort(1,"Cannot read configuration file $config_file: $!");

my @step_keys = qw(dna enriched nu_mt rmtmp o_mt i_del cn_mt);
my @param_keys = qw(threads refchr msize exclude orihs orihe orils orile score_threshold evalue_threshold split_length paired_distance deletion_threshold_min deletion_threshold_max breakthreshold cluster_threshold breakspan sizelimit hplimit flank split_distance_threshold);
my @missing = grep {!defined $User_Preferences->{$_} || $User_Preferences->{$_} eq ''} (@step_keys,@param_keys);
abort(1,"Missing keys in $config_file: @missing") if @missing;
for my $key (@step_keys){
  abort(1,"Key '$key' in $config_file must be 'yes' or 'no', got '$User_Preferences->{$key}'") unless $User_Preferences->{$key} =~ /^(yes|no)$/;
}
#NON-NUMERIC VALUES WOULD BE USED AS 0 OR A NUMERIC PREFIX WITHOUT AN ERROR
my @not_numeric = grep {$_ ne 'refchr' && $User_Preferences->{$_} !~ /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/} @param_keys;
abort(1,"Keys in $config_file must be numbers: ".join(", ",map {"$_ = '$User_Preferences->{$_}'"} @not_numeric)) if @not_numeric;

#PROGRAMS
my $hisat2 = $User_Preferences->{hisat2};
my $lastal = $User_Preferences->{lastal};
my $lastsp = $User_Preferences->{lastsp};
my $mfcv = $User_Preferences->{mfcv};
my $b2fq = $User_Preferences->{b2fq};
my $reformat = $User_Preferences->{reformat};
my $samtools = $User_Preferences->{samtools};
my $sambamba = $User_Preferences->{sambamba};
my $gcov = $User_Preferences->{gcov};
my $intersectBed = $User_Preferences->{intersectBed};
my $sortBed = $User_Preferences->{sortBed};
my $clusterBed = $User_Preferences->{clusterBed};
my $randomBed = $User_Preferences->{randomBed};
my $groupBy = $User_Preferences->{groupBy};
my $bg2bw = $User_Preferences->{bg2bw};

#DATABASES
my $hsindex = $User_Preferences->{hsindex};
my $faindex = $User_Preferences->{faindex};
my $lastindex = $User_Preferences->{lastindex};
my $mtfaindex = $User_Preferences->{mtfaindex};
my $gsize = $User_Preferences->{gsize};
my $MT_fasta = $User_Preferences->{MT_fasta};

#COMPUTATION
my $threads = $User_Preferences->{threads};

#MITOCHONDRIA FEATURES
my $refchr = $User_Preferences->{refchr};
my $msize = $User_Preferences->{msize};
my $orihs = $User_Preferences->{orihs};
my $orihe = $User_Preferences->{orihe};
my $orils = $User_Preferences->{orils};
my $orile = $User_Preferences->{orile};
my $exclude = $User_Preferences->{exclude};
my $dloop1 = $exclude;
my $dloop2 = $msize - $exclude;
my $hash;

#SCORING AND FILTERING FEATURES
my $score_threshold = $User_Preferences->{score_threshold};
my $evalue_threshold = $User_Preferences->{evalue_threshold};
my $split_length = $User_Preferences->{split_length};
my $paired_distance = $User_Preferences->{paired_distance};
my $deletion_threshold_min = $User_Preferences->{deletion_threshold_min};
my $deletion_threshold_max = $User_Preferences->{deletion_threshold_max};
my $breakthreshold = $User_Preferences->{breakthreshold};
my $cluster_threshold =  $User_Preferences->{cluster_threshold};
my $breakspan = $User_Preferences->{breakspan};
my $sizelimit = $User_Preferences->{sizelimit};
my $hplimit = $User_Preferences->{hplimit};
my $flank = $User_Preferences->{flank};
my $split_distance_threshold = $User_Preferences->{split_distance_threshold};

#STEPS
my $dna = $User_Preferences->{dna}; #USE SPECIFIC HISAT PARAMETERS FOR DNA AND RNA SEQUENCING
my $enriched = $User_Preferences->{enriched}; #IF THE SEQUENCING IS MITOCHONDRIAL DNA ENRICHED THEN SKIP THE INITIAL NUCLEAR GENOME ALIGNMENT STEP
my $nu_mt = $User_Preferences->{nu_mt}; #MAPPING TO NUCLEAR AND MITOCHONDRIAL GENOME WITH BOWTIE2
my $rmtmp = $User_Preferences->{rmtmp}; #REMOVE TEMPORARY FILES
my $o_mt = $User_Preferences->{o_mt}; #MITOCHONDRIAL READ EXTRACTION AND REMAPPING TO ONLY MITOCHONDRIAL GENOME WITH LASTAL
my $i_del = $User_Preferences->{i_del}; #IDENTIFICATION AND CLUSTERING OF DELETIONS/DUPLICATIONS
my $cn_mt = $User_Preferences->{cn_mt}; #ESTIMATION OF MT DNA COPY NUMBER

#EXIT IF THE CONFIGURATION STEPS DO NOT MATCH
if(($nu_mt eq 'yes' && $enriched eq 'yes')||($nu_mt eq 'no' && $o_mt eq 'no' && $i_del eq 'no')||($cn_mt eq 'yes' && $nu_mt eq 'no')||($cn_mt eq 'yes' && $o_mt eq 'no')||($i_del eq 'yes' && $o_mt eq 'no' && $nu_mt eq 'yes')){
  abort(1,"Invalid STEPS combination in $config_file (enriched=$enriched nu_mt=$nu_mt o_mt=$o_mt i_del=$i_del cn_mt=$cn_mt), see README");
}

#CHECK TOOLS AND REFERENCE FILES NEEDED BY THE ENABLED STEPS
my (@tools,@files,@problems);
if($nu_mt eq 'yes'){
  push @tools, qw(hisat2 samtools);
  push @files, qw(faindex);
  push @problems, "hsindex: no HISAT2 index found at $hsindex" unless (-s "$hsindex.1.ht2" || -s "$hsindex.1.ht2l");
}
if($o_mt eq 'yes'){
  push @tools, qw(lastal lastsp mfcv samtools gcov bg2bw);
  push @files, qw(mtfaindex);
  push @problems, "lastindex: no LAST index found at $lastindex" unless -s "$lastindex.prj";
}
if($cn_mt eq 'yes'){
  push @tools, qw(sambamba randomBed groupBy);
  push @files, qw(gsize);
}
if($i_del eq 'yes'){
  push @tools, qw(sortBed clusterBed intersectBed);
  push @files, qw(MT_fasta);
  push @problems, "Rscript not found in PATH" unless find_executable('Rscript');
  push @problems, "delplot.R not found in $Bin" unless -f "$Bin/delplot.R";
  push @problems, "tab/$tag.tab.gz not found, run with o_mt = yes first" if $o_mt eq 'no' && !-s "tab/$tag.tab.gz";
}
my %seen;
for my $key (grep {!$seen{$_}++} @tools){
  my $cmd = $User_Preferences->{$key};
  if(!defined $cmd || $cmd eq ''){push @problems, "$key: missing from configuration"; next;}
  push @problems, "$key: '$cmd' not found or not executable" unless find_executable($cmd);
}
for my $key (@files){
  my $file = $User_Preferences->{$key};
  if(!defined $file || $file eq ''){push @problems, "$key: missing from configuration"; next;}
  push @problems, "$key: file not found or empty: $file" unless -s $file;
}
#A WRONG CONTIG NAME (E.G. MT VS chrM) OTHERWISE GIVES EMPTY RESULTS WITHOUT ANY ERROR
push @problems, "refchr '$refchr' is not a contig in $faindex" if $nu_mt eq 'yes' && -s $faindex && !fai_contig_length($faindex,$refchr);
push @problems, "refchr '$refchr' is not a contig in $mtfaindex" if $o_mt eq 'yes' && -s $mtfaindex && !fai_contig_length($mtfaindex,$refchr);
abort(1,"Configuration check failed:\n  ".join("\n  ",@problems)) if @problems;

#MAP TO NU+MT GENOME
if($nu_mt eq 'yes' && $enriched eq 'no'){
  print scalar(localtime).": Map to NU+MT genome\n";
  if($dna eq 'yes'){run("$hisat2 -p $threads -x $hsindex --no-temp-splicesite --no-spliced-alignment --max-intronlen 5000 -U $p1 -S tmp_$tag.sam");}
  if($dna eq 'no'){run("$hisat2 -p $threads -x $hsindex -U $p1 -S tmp_$tag.sam");}

  #EXTRACT READS
  print scalar(localtime).": Extract reads\n";
  run("$samtools view -@ $threads -bt $faindex tmp_$tag.sam|$samtools sort -@ $threads -o tmp_$tag.bam -"); 
  run("$samtools index tmp_$tag.bam");
  run("$samtools idxstats tmp_$tag.bam > indel/$tag.count.txt");
  open(my $idx,'<',"indel/$tag.count.txt") or abort(2,"Cannot open indel/$tag.count.txt: $!");
  while(<$idx>){my @f = split /\t/; $summary{mt_reads} = $f[2] if $f[0] eq $refchr;}
  close($idx) or abort(2,"Cannot read indel/$tag.count.txt: $!");
  abort(3,"No reads mapped to $refchr in the nuclear+MT alignment, check the input FASTQ file") unless $summary{mt_reads};
  run("$samtools view -u -@ $threads -f 4 tmp_$tag.bam > tmp1_$tag.bam"); #flag to extract unmapped reads
  run("$samtools view -u -@ $threads tmp_$tag.bam $refchr> tmp2_$tag.bam");
  run("$samtools merge -@ $threads tmp_final_$tag.bam tmp1_$tag.bam tmp2_$tag.bam");

  #EXTRACT FASTQ
  print scalar(localtime).": Create FASTQ\n";
  run("$samtools fastq tmp_final_$tag.bam > tmp_$tag.fq"); #ADDS PAIR /1, /2 TO THE READ NAMES HENCE BBMAP REFORMAT NOT REQUIRED
  abort(3,"No reads extracted for MT remapping, tmp_$tag.fq is empty") unless -s "tmp_$tag.fq";
  #system("$reformat in=tmp.fq out=tmp1.fq overwrite=true -uniquenames -Xmx100g &>> log/$tag.log");
}

if($nu_mt eq 'yes' && $o_mt eq 'yes' && $enriched eq 'no'){
  #REMAP ON MT GENOME
  print scalar(localtime).": Map to MT genome\n";
  run("$lastal -Q1 -e80 -P$threads $lastindex tmp_$tag.fq|$lastsp > tmp_$tag.maf");
  run("$mfcv sam -d tmp_$tag.maf|$samtools view -@ $threads -bt $mtfaindex -|$samtools sort -@ $threads -o bam/$tag.bam -");
  run("$samtools index bam/$tag.bam");
  check_mt_bam();
  run("$mfcv tab tmp_$tag.maf > tab/$tag.tab");

  #COMPRESS TAB
  run("gzip -f tab/$tag.tab");

  #GENERATE BIGWIG
  print scalar(localtime).": Generate Bigwig\n";
  run("$gcov -split -ibam bam/$tag.bam -bg|sort -k1,1 -k2,2n > tmp_$tag.bg");
  run("$bg2bw tmp_$tag.bg $mtfaindex bw/$tag.bw");
  mt_mean_depth("tmp_$tag.bg");
  unlink("tmp_$tag.bg");
}

if($nu_mt eq 'no' && $o_mt eq 'yes' && $enriched eq 'yes'){
  #REMAP ON MT GENOME
  print scalar(localtime).": Map to MT genome\n";
  #system("$reformat in=$p1 in2=$p2 out=tmp_$tag.fq overwrite=true addslash=t trimreaddescription=t spaceslash=f -Xmx100g &>> log/$tag.log");
  run("$lastal -Q1 -e80 -P$threads $lastindex $p1|$lastsp > tmp_$tag.maf");
  run("$mfcv sam -d tmp_$tag.maf|$samtools view -@ $threads -bt $mtfaindex -|$samtools sort -@ $threads -o bam/$tag.bam -");
  run("$samtools index bam/$tag.bam");
  check_mt_bam();
  run("$mfcv tab tmp_$tag.maf > tab/$tag.tab");

  #COMPRESS TAB
  run("gzip -f tab/$tag.tab");

  #GENERATE BIGWIG
  print scalar(localtime).": Generate Bigwig\n";
  run("$gcov -split -ibam bam/$tag.bam -bg|sort -k1,1 -k2,2n > tmp_$tag.bg");
  run("$bg2bw tmp_$tag.bg $mtfaindex bw/$tag.bw");
  mt_mean_depth("tmp_$tag.bg");
  unlink("tmp_$tag.bg");
}

if($cn_mt eq 'yes' && $enriched eq 'no' && $o_mt eq 'yes' && $nu_mt eq 'yes'){
  run("$sambamba depth window -w $msize bam/$tag.bam|cut -f 1,5 > indel/$tag.cnmt.txt");
  abort(2,"$sambamba depth window wrote no coverage rows to indel/$tag.cnmt.txt") unless count_lines("indel/$tag.cnmt.txt") > 1;
  run("$randomBed -l 1000 -n 3000 -g $gsize > tmp_random_$tag.bed");
  abort(2,"$randomBed produced no regions from $gsize, check that it lists nuclear contigs and their lengths") unless -s "tmp_random_$tag.bed";
  run("$sambamba depth region -L tmp_random_$tag.bed -t $threads -o tmp_random_$tag.cov tmp_$tag.bam");
  #BEDTOOLS GROUPBY EXITS 0 ON AN EMPTY INPUT, SO CHECK FOR COVERAGE ROWS BELOW THE HEADER
  abort(3,"No nuclear reads found in the random regions used for copy number estimation: the input may contain only MT reads (then set enriched = yes, nu_mt = no and cn_mt = no), or the contig names in $gsize do not match $faindex") unless count_lines("tmp_random_$tag.cov") > 1;
  run("$groupBy -g 5 -c 8 -o mean -i tmp_random_$tag.cov >> indel/$tag.cnmt.txt");
}

#IDENTIFY DELETIONS/DUPLICATIONS
if($i_del eq 'yes'){
  my $infile = "tab/$tag.tab.gz";
  my $bedfile = "indel/$tag.bed";
  my $breakpointfile = "indel/$tag.breakpoint";
  my $clusterfile = "indel/$tag.cluster";
  
  my $check_paired = &check_paired($infile);
  print scalar(localtime).": Build split read hash\n";
    $hash = &build_hash($infile,$check_paired);
  abort(3,"No alignments in $infile passed the score and e-value thresholds") unless $hash && %$hash;
  &remove_duplicates();

  print scalar(localtime).": Generate non-split read BED\n";
  my $nosplitbed_file = &print_bed($hash,$tag,$refchr);

  print scalar(localtime).": Process hash to get best deletion/duplication candidates\n";
  my $delhash = &process_hash($hash,$bedfile,$breakpointfile,$msize,$tag);
  $summary{split_reads} = count_lines($breakpointfile);

  print scalar(localtime).": Build split read clusters\n";
  my $clusterhash = &get_cluster($delhash,$breakthreshold,$tag);

  print scalar(localtime).": Generate and print results\n";
  my $clustercount = &build_results($hash,$clusterhash,$infile,$refchr,$cluster_threshold,$tag);
  &print_result($clustercount,$clusterfile);
  $summary{clusters} = count_lines($clusterfile);

  #PLOT DELETIONS
  print scalar(localtime).":Plot deletions/duplications\n";
  my $tsvfile = "indel/$tag.tsv";
  #A PLOT FROM AN EARLIER RUN WOULD OTHERWISE REMAIN WHEN THIS RUN FINDS NO EVENTS
  unlink $tsvfile,"plot/$tag.pdf";
  run("Rscript '$Bin/delplot.R' $msize $orihs $orihe $orils $orile $sizelimit indel/$tag.cluster indel/$tag.breakpoint $tag $hplimit $MT_fasta $flank");
  abort(2,"delplot.R did not write $tsvfile") unless -s $tsvfile;
  $summary{events_reported} = count_lines($tsvfile) - 1;
}

#EXACT NAMES: A GLOB SUCH AS tmp_S1* WOULD ALSO REMOVE THE FILES OF SAMPLE S10 IN THE SAME DIRECTORY
if($rmtmp eq 'yes'){
  my @tmp_sample = ((map {"tmp_$tag.$_"} qw(sam bam bam.bai fq maf bg)),(map {"tmp_${tag}$_"} qw(_nosplit.bed _bps.bed _bpe.bed _bps.cls _bpe.cls .split.start.bed .split.end.bed .intersect.start.bed .intersect.end.bed)));
  if($nu_mt eq 'yes'){unlink(@tmp_sample,"tmp1_$tag.bam","tmp2_$tag.bam","tmp_final_$tag.bam");}
  if($o_mt eq 'yes' && $enriched eq 'yes'){unlink @tmp_sample;}
  if($i_del eq 'yes'){unlink @tmp_sample;}
  if($cn_mt eq 'yes'){unlink("tmp_random_$tag.bed","tmp_random_$tag.cov");}
}

write_summary('ok',0) or abort(2,"Cannot write $summary_file: $!");
print scalar(localtime).":Finished\n";
###########################################################SUBROUTINES################################################################
sub check_mt_bam{
  my $n = capture("$samtools view -c -F 4 bam/$tag.bam");
  chomp($n);
  $summary{mt_alignments} = $n;
  abort(3,"No reads aligned to the MT genome in bam/$tag.bam") unless $n;
  return;
}

#MEAN DEPTH OVER THE MT CONTIG FROM THE BEDGRAPH, TO JUDGE WHETHER A RUN WITHOUT EVENTS HAD ENOUGH COVERAGE
sub mt_mean_depth{
  my ($bg) = @_;
  my $len = fai_contig_length($mtfaindex,$refchr) or return;
  open(my $fh,'<',$bg) or abort(2,"Cannot open $bg: $!");
  my $sum = 0;
  while(<$fh>){my @f = split /\t/; $sum += ($f[2]-$f[1])*$f[3] if $f[0] eq $refchr;}
  close($fh) or abort(2,"Cannot read $bg: $!");
  $summary{mt_mean_depth} = sprintf("%.1f",$sum/$len);
  return;
}

sub check_paired{
  my $check_paired = 0;
  my $path = shift;
  my $line = capture("zcat $path|tail -n 1");
  my @elements = split(/\t/, $line);
  #PAIRED (/1, /2) OR UNPAIRED READ NAMES ARE DETECTED FROM THE LAST LINE ONLY, SO IT MUST BE AN ALIGNMENT ROW
  if($line =~ /^#/ || $line !~ /\S/){
    my $rows = capture("zcat $path|awk '!/^#/{n++} END{print n+0}'");
    abort(3,"No alignments in $path") unless $rows > 0;
  }
  abort(2,"Cannot detect paired or unpaired read names: the last line of $path is not an alignment row") if $line =~ /^#/ || @elements < 12;

  $check_paired = 1 if $elements[6]=~m/\/.$/;
  return $check_paired;
}

sub build_hash{
  my $file = shift;
  my $check_paired = shift;

  open(IN, "gunzip -c $file |") or abort(2,"Cannot read $file: $!");
  my $hash;
  my $check_hash;
  while(<IN>){
    my $line = $_;
    chomp($line);
    next if $line=~m/^\#/;

    my @elements = split(" ",$line);
    my $score = $elements[0];
    my $chr = $elements[1];

    my $start = $elements[2];
    my $length = $elements[3];
    my $end = $start + $length;

    my $read_start = $elements[7];
    my $id;
    my $mate = 1;

    if($check_paired==0){
      $id = $elements[6];
      if($id=~m/\_.$/){
        $mate = $id;
        $mate =~s/.*\_(.)$/$1/;
        $id =~s/\_.$//;
      }
    }
    elsif($check_paired == 1){
      my @name = split(/\//,$elements[6]);       $id = $name[0];
      $mate = $name[1];  
    }

    my $strand = $elements[9];
    my $evalue = $elements[12];
    $evalue=~s/.*\=//;
    next if exists $check_hash->{$id}->{$start.$length};
    $check_hash->{$id}->{$start.$length}++;
    $mate = 2 if $mate > 1;
    #print "$id\t$mate\t$start\t$end\t$strand\t$len\t$score\t$eval\n"; 

    #FILTER BY SCORE AND MAP EVALUE
    next if $score < $score_threshold;
    next if $evalue > $evalue_threshold;
    
    push(@{$hash->{$id}->{$mate}->{starts}},$start);
    push(@{$hash->{$id}->{$mate}->{ends}},$end);
    push(@{$hash->{$id}->{$mate}->{rstarts}},$read_start);
    push(@{$hash->{$id}->{$mate}->{strands}},$strand);
    push(@{$hash->{$id}->{$mate}->{scores}},$score);
    push(@{$hash->{$id}->{$mate}->{evalue}},$evalue);
    push(@{$hash->{$id}->{$mate}->{lengths}},$length);
  }
  undef $check_hash;
  close(IN) or abort(2,"Failed to read $file");
  return $hash;
}

sub remove_duplicates{
  my $check_duplicates;

  for my $id (keys %{$hash}){
    my $starts1='NULL';
    my $ends1='NULL';
    my $lengths1='NULL';
    my $starts2='NULL';
    my $ends2='NULL';
    my $lengths2='NULL';
    my $count1 = 0;
    my $count2 = 0;

    if(exists $hash->{$id}->{1}){
      $starts1 = join("",@{$hash->{$id}->{1}->{starts}});
      $ends1 = join("",@{$hash->{$id}->{1}->{ends}});
      $lengths1 = join("",@{$hash->{$id}->{1}->{lengths}});
      $count1 = @{$hash->{$id}->{1}->{starts}};

    }
    if(exists $hash->{$id}->{2}){
      $starts2 = join("",@{$hash->{$id}->{2}->{starts}});
      $ends2 = join("",@{$hash->{$id}->{2}->{ends}});
      $lengths2 = join("",@{$hash->{$id}->{2}->{lengths}});
      $count2 = @{$hash->{$id}->{2}->{starts}};
    }

    my $signature1 = $starts1.":".$ends1.":".$lengths1;
    my $signature2 = $starts2.":".$ends2.":".$lengths2;

    my @signatures = ($signature1,$signature2);
    @signatures = sort @signatures;
    my $signature = join("",@signatures);    


    if(exists $check_duplicates->{$signature} && $count1>0 && $count2>0){
      delete $hash->{$id};
    }
    else{
      $check_duplicates->{$signature}++;
    }
  }
}

#PRINT NON-SPLIT BED
sub print_bed{
  my $hash = shift;
  my $tag = shift;
  my $refchr = shift;
  my $filename  = "tmp_".$tag."_nosplit.bed";
  open(ABED,">$filename") or abort(2,"Cannot open $filename: $!");
  for my $id (keys %{$hash}){
     #print $id."\n";
 
    for my $read(keys %{$hash->{$id}}){
      my $name = $id."_".$read;

      #CHECK IF READ IS SPLIT ALIGNED
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next if $count > 1;
      my $start = @{$hash->{$id}->{$read}->{starts}}[0];    
      my $end = @{$hash->{$id}->{$read}->{ends}}[0];
      print ABED "$refchr\t$start\t$end\t$name\n"; 
    }        
  }
  close(ABED) or abort(2,"Cannot write $filename: $!");
  run("sort -k2,2n $filename -o $filename");
  return $filename;
}

sub process_hash{
  my $hash = shift;
  my $bedfile = shift;
  my $breakpointfile = shift;
  my $size = shift;
  my $tag = shift;
  my $bps = "tmp_".$tag."_bps.bed";
  my $bpe = "tmp_".$tag."_bpe.bed";
  my $delhash=&process_hash1($hash,$bedfile,$breakpointfile,$msize,$bps,$bpe);

  for my $id (keys %{$hash}){
    #print $id."\n";
    for my $read(keys %{$hash->{$id}}){
      #CHECK IF READ IS SPLIT ALIGNED
      my $name = $id."_".$read;
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next unless $count == 2;
      #next unless $name eq 'C1LUFACXX130228:4:1312:19820:41898_1';
      #CHECK IF THE SPLIT READ FRAGMENTS ARE ABOVE LENGTH THRESHOLD
      my $min_len = min @{$hash->{$id}->{$read}->{lengths}};
      next unless $min_len >= $split_length;

      #GET DISTANCE BETWEEN THE FRAGMENTS AND THE BREAKPOINTS
      my @read_starts = @{$hash->{$id}->{$read}->{starts}};    
      my @read_ends = @{$hash->{$id}->{$read}->{ends}};
      my @read_strands = @{$hash->{$id}->{$read}->{strands}};
      my @read_lengths = @{$hash->{$id}->{$read}->{lengths}};
      my @read_local_starts = @{$hash->{$id}->{$read}->{rstarts}};

      #LOOK FOR SPLIT READS WHICH MAP INVERSELY SPANNING THE D'LOOP
      my $read_check = 'no';
      if($read_starts[0]<$read_starts[1] && $read_local_starts[0]>$read_local_starts[1]){$read_check = 'yes';}
      elsif($read_starts[0]>$read_starts[1] && $read_local_starts[0]<$read_local_starts[1]){$read_check = 'yes';}
      elsif($read_starts[0]<=$read_starts[1] && $read_local_starts[0]<$read_local_starts[1] && $read_starts[1]<$read_ends[0] ){$read_check = 'yes';}
      elsif($read_starts[0]>=$read_starts[1] && $read_local_starts[0]>$read_local_starts[1] && $read_ends[1]>$read_starts[0] ){$read_check = 'yes';}
      
      #CHECK IF THE SPLIT READ FRAGMENTS SPAN THE DLOOP
      my $min_start = min @{$hash->{$id}->{$read}->{starts}};    
      my $max_end = max @{$hash->{$id}->{$read}->{ends}};
      next if $min_start <= $dloop1 && $max_end >= $dloop2;   

      my ($size,$start,$end) = &get_frag_distance(\@read_starts,\@read_ends,$read_check,$msize);
      next if $size < $deletion_threshold_min;
      next if $size > $deletion_threshold_max;
      next if $read_strands[0] ne $read_strands[1];#SPLIT READS IN OPPOSITE STRANDS MAY INDICATE INVERSION
      #print "$name\t$start\t$end\t$size\n";

      $start = 1 if $start == 0;
      $end = 1 if $end == 0;   

      #CHECK IF PAIRED READ EXISTS AND ITS ALIGNMENT POSITION SUPPORTS THE SPLIT READ
      my $pair;
      $pair = 1 if $read == 2;
      $pair = 2 if $read == 1;
      my $paired_support = 'no';
      my $distance_paired_support = &paired_support($hash,$id,$read,$pair);
      $paired_support = 'yes' if $distance_paired_support <= $paired_distance;

      #LOOK FOR SPLIT READS WITH LARGE UNMAPPED AREA
      my ($split_distance) = &get_split_distance(\@read_local_starts,\@read_lengths);
      next if $split_distance > $split_distance_threshold;

      #PRINT SPLIT READ IN BED FORMAT
      my @read_scores = @{$hash->{$id}->{$read}->{scores}};
      my @len = &generate_bed($id,$read,\@read_starts,\@read_ends,\@read_lengths,\@read_scores,$read_strands[0]);

      my $readid = $id."_".$read;
      my $clusterid = "cluster_".$count;
      $delhash->{$readid}->{breakstart}=$len[2];
      $delhash->{$readid}->{breakend}=$len[3];
      $delhash->{$readid}->{breaksize}=$len[4];
      $delhash->{$readid}->{readcheck}=$read_check;
      $delhash->{$readid}->{lenstart}=$len[0];
      $delhash->{$readid}->{lenend}=$len[1];
      
      $start = $len[2];
      $end = $len[3];
      $size = $len[4];

      print BP "$refchr\t$name\t$size\t$start\t$end\t$read_lengths[0]\t$read_lengths[1]\t$paired_support\t$distance_paired_support\t$read_check\n";
      print BPS "$refchr\t$start\t$start\t$name\t0\t+\n" if $read_check eq 'no';
      print BPS "$refchr\t$start\t$start\t$name\t0\t-\n" if $read_check eq 'yes';
      print BPE "$refchr\t$end\t$end\t$name\t0\t+\n" if $read_check eq 'no';
      print BPE "$refchr\t$end\t$end\t$name\t0\t-\n" if $read_check eq 'yes';
    }
  }
  close(BP) or abort(2,"Cannot write $breakpointfile: $!");
  close(BED) or abort(2,"Cannot write $bedfile: $!");
  close(BPS) or abort(2,"Cannot write $bps: $!");
  close(BPE) or abort(2,"Cannot write $bpe: $!");
  return $delhash;
}

sub process_hash1{
  my $hash = shift;
  my $bedfile = shift;
  my $breakpointfile = shift;
  my $size = shift;
  my $bps = shift;
  my $bpe = shift;

  my $delhash;
 
  open(BP,">$breakpointfile") or abort(2,"Cannot open $breakpointfile: $!");
  open(BED,">$bedfile") or abort(2,"Cannot open $bedfile: $!");
  open(BPS,">$bps") or abort(2,"Cannot open $bps: $!");
  open(BPE,">$bpe") or abort(2,"Cannot open $bpe: $!");

  for my $id (keys %{$hash}){
    #print $id."\n";
    for my $read(keys %{$hash->{$id}}){
      #CHECK IF READ IS SPLIT ALIGNED
      my $name = $id."_".$read;
      my $count = @{$hash->{$id}->{$read}->{starts}};
      next unless $count == 3;

      #GET DISTANCE BETWEEN THE FRAGMENTS AND THE BREAKPOINTS
      my @read_starts = @{$hash->{$id}->{$read}->{starts}};    
      my @read_ends = @{$hash->{$id}->{$read}->{ends}};
      my @read_strands = @{$hash->{$id}->{$read}->{strands}};
      my @read_lengths = @{$hash->{$id}->{$read}->{lengths}};
      my @read_scores = @{$hash->{$id}->{$read}->{scores}};
      my @read_local_starts = @{$hash->{$id}->{$read}->{rstarts}};


      my %strands=map {$_ => 1} @read_strands;
      my $count_read_strands = keys(%strands);
      next unless $count_read_strands == 1;

      my $as=join(';',@read_starts);my $ae=join(';',@read_ends);my $asr=join(';',@read_strands);my $als=join(';',@read_local_starts);
      my @pos = (1000,1000,1000,1000);
      #CHECK FOR FALSE SPLITS AT ORIGIN WITH A REAL DELETION
      $pos[0] = 0 if grep {$_ == 0} @read_starts;
      $pos[1] = $msize if grep {$_ == $msize} @read_starts;
      $pos[2] = 0 if grep {$_ == 0} @read_ends;
      $pos[3] = $msize if grep {$_ == $msize} @read_ends;

      #print "$name\t$as\t$ae\t$asr\t$als\t$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\n";
      if(($pos[0] != 1000||$pos[1] != 1000)&&($pos[2] != 1000||$pos[3] != 1000)){
        my $sort_hash;
        for(my $i=0;$i<3;$i++){
          $sort_hash->{$read_local_starts[$i]}->{start}=$read_starts[$i];
          $sort_hash->{$read_local_starts[$i]}->{end}=$read_ends[$i];
          $sort_hash->{$read_local_starts[$i]}->{strand}=$read_strands[$i];
          $sort_hash->{$read_local_starts[$i]}->{length}=$read_lengths[$i];
          $sort_hash->{$read_local_starts[$i]}->{score}=$read_scores[$i];
        }
        my $order_count=1;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          my $start_check=$sort_hash->{$local_start}->{start};
          my $end_check=$sort_hash->{$local_start}->{end};
          my $pos_check='no';
          $pos_check='yes' unless($start_check==0||$start_check==$msize||$end_check==0||$end_check==$msize);
          if($pos_check eq 'yes'){last;}
          $order_count++;
        }

        my $delete_count=3 if $order_count==1;
        $delete_count=1 if $order_count==3;
        $order_count=1;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          delete $sort_hash->{$local_start} if($order_count==$delete_count);
          $order_count++;
        }
        my @read_startsN;
        my @read_endsN;
        my @read_strandsN;
        my @read_local_startsN;
        my @read_scoresN;
        my @read_lengthsN;
        foreach my $local_start (sort { $a <=> $b } keys %{$sort_hash}){
          my $start =   $sort_hash->{$local_start}->{start};
          my $end =   $sort_hash->{$local_start}->{end};
          my $score =   $sort_hash->{$local_start}->{score};
          my $length =   $sort_hash->{$local_start}->{length};
          my $strand =   $sort_hash->{$local_start}->{strand};
          push @read_startsN,$start;
          push @read_endsN,$end;
          push @read_local_startsN,$local_start;
          push @read_strandsN,$strand;
          push @read_scoresN,$score;
          push @read_lengthsN,$length;
        }

        #CHECK IF THE SPLIT READ FRAGMENTS ARE ABOVE LENGTH THRESHOLD
        my $min_len = min @read_lengthsN;
        next unless $min_len >= $split_length;
  

        #LOOK FOR SPLIT READS WHICH MAP INVERSELY SPANNING THE D'LOOP
        my $read_check = 'no';
        if($read_startsN[0]<$read_startsN[1] && $read_local_startsN[0]>$read_local_startsN[1]){$read_check = 'yes';}
        elsif($read_startsN[0]>$read_startsN[1] && $read_local_startsN[0]<$read_local_startsN[1]){$read_check = 'yes';}
        my ($size,$start,$end) = &get_frag_distance(\@read_startsN,\@read_endsN,$read_check,$msize);
        next if $size < $deletion_threshold_min;
        next if $size > $deletion_threshold_max;

        #CHECK IF PAIRED READ EXISTS AND ITS ALIGNMENT POSITION SUPPORTS THE SPLIT READ
        my $pair;
        $pair = 1 if $read == 2;
        $pair = 2 if $read == 1;
        my $paired_support = 'no';
        my $distance_paired_support = &paired_support($hash,$id,$read,$pair);
        $paired_support = 'yes' if $distance_paired_support <= $paired_distance;

        $start = 1 if $start == 0;
        $end = 1 if $end == 0;   

        #LOOK FOR SPLIT READS WITH LARGE UNMAPPED AREA
        my ($split_distance) = &get_split_distance(\@read_local_startsN,\@read_lengthsN);
        next if $split_distance > 5;

        #PRINT SPLIT READ IN BED FORMAT
        my @len = &generate_bed1($id,$read,\@read_startsN,\@read_endsN,\@read_lengthsN,\@read_scoresN,$read_strandsN[0]);
        my $readid = $id."_".$read;
        my $clusterid = "cluster_".$count;
        $delhash->{$readid}->{breakstart}=$start;
        $delhash->{$readid}->{breakend}=$end;
        $delhash->{$readid}->{breaksize}=$size;
        $delhash->{$readid}->{readcheck}=$read_check;
        $delhash->{$readid}->{lenstart}=$len[0];
        $delhash->{$readid}->{lenend}=$len[1];
        
        #$start = $len[2];
        #$end = $len[3];
        #$size = $len[4];

        print BP "$refchr\t$name\t$size\t$start\t$end\t$read_lengths[0]\t$read_lengths[1]\t$paired_support\t$distance_paired_support\t$read_check\n";
        print BPS "$refchr\t$start\t$start\t$name\t0\t+\n" if $read_check eq 'no';
        print BPS "$refchr\t$start\t$start\t$name\t0\t-\n" if $read_check eq 'yes';
        print BPE "$refchr\t$end\t$end\t$name\t0\t+\n" if $read_check eq 'no';
        print BPE "$refchr\t$end\t$end\t$name\t0\t-\n" if $read_check eq 'yes';

      }
    }
  }
  return $delhash;
}

#CHECK IF THE PAIR OF A SPLIT READ LIES WITHIN A GIVEN THRESHOLD DISTANCE TO A FRAGMENT OF THE SPLIT READ
sub paired_support{
   my $hash = shift;
   my $id = shift;
   my $read = shift;
   my $pair = shift;
   my $distance = 1000;

   if(exists $hash->{$id}->{$pair}){
     my $distance1 = 0;
     my $distance2=0;

     my @read_starts = @{$hash->{$id}->{$read}->{starts}};    
     my @read_ends = @{$hash->{$id}->{$read}->{ends}};
    
     my $pair_start = @{$hash->{$id}->{$pair}->{starts}}[0];    
     my $pair_end = @{$hash->{$id}->{$pair}->{ends}}[0];

     $distance1 = $read_starts[0] - $pair_end if $read_starts[0] > $pair_end;
     $distance1 = $pair_start - $read_ends[0] if $pair_start > $read_ends[0];

     $distance2 = $read_starts[1] - $pair_end if $read_starts[1] > $pair_end;
     $distance2 = $pair_start - $read_ends[1] if $pair_start > $read_ends[1];
     $distance = min ($distance1,$distance2);
   }
   return $distance;
}

#GET THE SIZE OF PUTATIVE DELETION AND THE BREAKPOINTS
sub get_frag_distance{
  my $starts = shift;
  my $ends = shift;
  my $read_check = shift;
  my $msize = shift;
  my @res;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];
  
  if($read_check eq 'no'){
    if($frag1_start > $frag2_end){
      my $size = $frag1_start - $frag2_end; 
      @res = ($size,$frag2_end,$frag1_start);
    }
    elsif($frag2_start > $frag1_end){
      my $size = $frag2_start - $frag1_end;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag1_start < $frag2_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag2_start < $frag1_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag1_start <= $frag2_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_start + $frag1_end ;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag2_start <= $frag1_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag1_end + $frag2_start ;
      @res = ($size,$frag2_start,$frag1_end);
    }
  }

  if($read_check eq 'yes'){
    if($frag1_start > $frag2_end){
      my $size = ($msize -$frag1_end) + $frag2_start; 
      @res = ($size,$frag2_start,$frag1_end);
    }
    elsif($frag2_start > $frag1_end){
      my $size = ($msize-$frag2_end) + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag1_start <= $frag2_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_start + $frag1_end ;
      @res = ($size,$frag1_end,$frag2_start);
    }
    elsif($frag2_start <= $frag1_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag1_end + $frag2_start ;
      @res = ($size,$frag2_start,$frag1_end);
    }
    elsif($frag1_start < $frag2_start && $frag2_end <= $frag1_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
    elsif($frag2_start < $frag1_start && $frag1_end <= $frag2_end){
      my $size = $msize - $frag2_end + $frag1_start;
      @res = ($size,$frag1_start,$frag2_end);
    }
  }
  return @res;
}

#GET THE DISTANCE BETWEEN SPLIT READS
sub get_split_distance{
  my $starts = shift;
  my $lengths = shift;
  my @res;

  my $local1_start = @{$starts}[0];
  my $local2_start = @{$starts}[1];

  my $local1_length = @{$lengths}[0];
  my $local2_length = @{$lengths}[1];
  
  if($local1_start>$local2_start){
    my $distance = $local1_start - $local2_start - $local2_length;
    @res = ($distance);
  }
  if($local2_start>$local1_start){
    my $distance = $local2_start - $local1_start - $local1_length;
    @res = ($distance);
  }  
  return @res;
}

#GENERATE A BED FILE FOR THE SPLIT READS FOR IGV VISUALIZATION
sub generate_bed{
  my $id = shift;
  my $read = shift;
  my $starts = shift;
  my $ends = shift;
  my $lengths = shift;
  my $scores = shift;
  my $strand = shift;
  
  my $size = shift;
  my $bstart = shift;
  my $bend = shift;
  my @len;

  my $start;
  my $end;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];

  my $frag1_length = @{$lengths}[0];
  my $frag2_length = @{$lengths}[1];

  my $frag1_score = @{$scores}[0];
  my $frag2_score = @{$scores}[1];
  my $score = int(($frag1_score+$frag2_score)/2);

  if($frag1_start > $frag2_end && $strand eq '-'){
    $start = $frag2_start;
    $end = $frag1_end;
    $bstart = $frag2_end;
    $bend = $frag1_start;
    $size = $frag1_start - $frag2_end;
    my $block_start=$frag1_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";  
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }
  elsif($frag1_start > $frag2_end && $strand eq '+'){
    $start = $frag2_start;
    $end = $frag1_end;
    $bstart = $frag2_start;
    $bend = $frag1_end;
    $size = 1 + $msize- $frag1_end + $frag2_start;
    my $block_start=$frag1_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";  
    @len = ($frag1_length,$frag2_length,$bstart,$bend,$size);
  }
  elsif($frag2_start > $frag1_end && $strand eq '+'){
    $start = $frag1_start;
    $end = $frag2_end; 
    $bstart = $frag1_end;
    $bend = $frag2_start;
    $size = $frag2_start - $frag1_end;
    my $block_start=$frag2_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";  
    @len = ($frag1_length,$frag2_length,$bstart,$bend,$size);
  }
  elsif($frag2_start > $frag1_end && $strand eq '-'){
    $start = $frag1_start;
    $end = $frag2_end; 
    $bstart = $frag1_start;
    $bend = $frag2_end;
    $size = 1 + $msize - $frag2_end + $frag1_start;
    my $block_start=$frag2_start-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";  
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }
  elsif($frag2_start < $frag1_end && $strand eq '+'){
    $start = $frag1_start;
    $end = $frag2_end; 
    $bstart = $frag2_start;
    $bend = $frag1_end; 
    $size = $msize - $frag1_end + $frag2_start;
    my $block_start=$frag1_end-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";  
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }  
  elsif($frag1_start < $frag2_end && $strand eq '-'){
    $start = $frag2_start;
    $end = $frag1_end; 
    $bstart = $frag1_start;
    $bend = $frag2_end; 
    $size = $msize - $frag2_end + $frag1_start;
    my $block_start=$frag2_end-$start;
    $start = 1 if $start == 0;$end = 1 if $end == 0;$bstart = 1 if $bstart == 0;$bend = 1 if $bend == 0;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";  
    @len = ($frag2_length,$frag1_length,$bstart,$bend,$size);
  }  
  
  return @len;
}

#GENERATE A BED FILE FOR THE SPLIT READS FOR IGV VISUALIZATION (reads which split twice, once due to deletion/duplication and once due to genome circularity)
sub generate_bed1{
  my $id = shift;
  my $read = shift;
  my $starts = shift;
  my $ends = shift;
  my $lengths = shift;
  my $scores = shift;
  my $strand = shift;
  my @len;

  my $start;
  my $end;

  my $frag1_start = @{$starts}[0];
  my $frag2_start = @{$starts}[1];

  my $frag1_end = @{$ends}[0];
  my $frag2_end = @{$ends}[1];

  my $frag1_length = @{$lengths}[0];
  my $frag2_length = @{$lengths}[1];

  my $frag1_score = @{$scores}[0];
  my $frag2_score = @{$scores}[1];
  my $score = int(($frag1_score+$frag2_score)/2);

  if($frag1_start > $frag2_end){
    $start = $frag2_start;
    $end = $frag1_end;
    my $block_start=$frag1_start-$start;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag2_length\,$frag1_length\t0\,$block_start\n";  
    @len = ($frag2_length,$frag1_length);
  }
  elsif($frag2_start > $frag1_end){
    $start = $frag1_start;
    $end = $frag2_end; 
    my $block_start=$frag2_start-$start;
    print BED "$refchr\t$start\t$end\t$id\_$read\t$score\t$strand\t$start\t$end\t0\t2\t$frag1_length\,$frag2_length\t0\,$block_start\n";  
    @len = ($frag1_length,$frag2_length);
  }
  return @len;
}

sub get_cluster{
  my $delhash = shift;
  my $breakthreshold = shift;
  my $tag = shift;
  my $bps = "tmp_".$tag."_bps.bed";
  my $bpe = "tmp_".$tag."_bpe.bed";
  my $bpsc = "tmp_".$tag."_bps.cls";
  my $bpec = "tmp_".$tag."_bpe.cls";
  my $clusterhash;

  run("$sortBed -i $bps|$clusterBed -s -d $breakthreshold -i stdin > $bpsc");
  run("$sortBed -i $bpe|$clusterBed -s -d $breakthreshold -i stdin > $bpec");

  open(BPS,"<$bpsc") or abort(2,"Cannot open $bpsc: $!");
  my @bps=<BPS>;
  my $bps_hash;
  foreach my $line(@bps){
    chomp $line;
    my @fields=split(/\t/,$line);
    my $readid=$fields[3];
    my $clusterid=$fields[6];
    $bps_hash->{$readid} = $clusterid;
  }
  close(BPS) or abort(2,"Cannot read $bpsc: $!");

  open(BPE,"<$bpec") or abort(2,"Cannot open $bpec: $!");
  my @bpe=<BPE>;
  my $bpe_hash;
  foreach my $line(@bpe){
    chomp $line;
    my @fields=split(/\t/,$line);
    my $readid=$fields[3];
    my $clusterid=$fields[6];
    $bpe_hash->{$readid} = $clusterid;
  }
  close(BPE) or abort(2,"Cannot read $bpec: $!");
 
  for my $readid(keys %{$delhash}){
    next unless exists $bps_hash->{$readid};
    next unless exists $bpe_hash->{$readid};
    my $clusterids = $bps_hash->{$readid};
    my $clusteride = $bpe_hash->{$readid};
    my $clusterid="cluster".$clusterids.$clusteride;
    $clusterhash->{$clusterid}->{$readid} = $delhash->{$readid};
  }
  return $clusterhash;
}

sub build_results{
  my $hash = shift;
  my $clusterhash = shift;
  my $infile = shift;
  my $refchr = shift;
  my $cluster_read_count_threshold = shift;
  my $tag = shift;
  my $clustercount;  

  my $nosplitbed_file = "tmp_".$tag."_nosplit.bed";;
  for my $clusterid (keys %{$clusterhash}){
    my $cluster_read_count = scalar keys %{$clusterhash->{$clusterid}}; 
    next if $cluster_read_count < $cluster_read_count_threshold; 

    print "Check $clusterid\n";
    my $clustercheck;
    $clustercount->{$clusterid}->{wt} = 0;
    $clustercount->{$clusterid}->{mt} = scalar keys %{$clusterhash->{$clusterid}};    

    my $splitbed_fileS = "tmp_".$tag.".split.start.bed";
    my $splitbed_fileE = "tmp_".$tag.".split.end.bed";

    open(SBEDS,">$splitbed_fileS") or abort(2,"Cannot open $splitbed_fileS: $!");
    open(SBEDE,">$splitbed_fileE") or abort(2,"Cannot open $splitbed_fileE: $!");

    for my $readid (keys %{$clusterhash->{$clusterid}}){
      my $breakstart = $clusterhash->{$clusterid}->{$readid}->{breakstart};
      my $breakend = $clusterhash->{$clusterid}->{$readid}->{breakend};
      my $lenstart = $clusterhash->{$clusterid}->{$readid}->{lenstart};
      my $lenend = $clusterhash->{$clusterid}->{$readid}->{lenend};

      push(@{$clustercount->{$clusterid}->{starts}},$breakstart);
      push(@{$clustercount->{$clusterid}->{ends}},$breakend);
      push(@{$clustercount->{$clusterid}->{names}},$readid);
      push(@{$clustercount->{$clusterid}->{lenstarts}},$lenstart);
      push(@{$clustercount->{$clusterid}->{lenends}},$lenend);

      print SBEDS "$refchr\t$breakstart\t$breakstart\n";
      print SBEDE "$refchr\t$breakend\t$breakend\n";
    }
    close(SBEDS) or abort(2,"Cannot write $splitbed_fileS: $!");
    close(SBEDE) or abort(2,"Cannot write $splitbed_fileE: $!");

    run("sort -u -k2,2n -k3 $splitbed_fileS -o $splitbed_fileS");
    run("sort -u -k2,2n -k3 $splitbed_fileE -o $splitbed_fileE");

    my $intersectbed_fileS = "tmp_".$tag.".intersect.start.bed";
    my $intersectbed_fileE = "tmp_".$tag.".intersect.end.bed";


    run("$intersectBed -wo -sorted -a $nosplitbed_file -b $splitbed_fileS > $intersectbed_fileS");
    open(IS,"<$intersectbed_fileS") or abort(2,"Cannot open $intersectbed_fileS: $!");
    while(<IS>){
      my $line = $_;
      chomp($line);
      my @crd = split(/\t/,$line);
      my $start = $crd[1];
      my $end = $crd[2];
      my $name = $crd[3];
      my $breakstart = $crd[5];
      next if exists $clustercheck->{$name};
      my $diffstart = $breakstart - $start;
      my $diffend = $end - $breakstart;
      next if $diffstart <= $breakspan || $diffend <= $breakspan;
      $clustercount->{$clusterid}->{wt1}++;
      $clustercheck->{$name}++;
    }
    close(IS) or abort(2,"Cannot read $intersectbed_fileS: $!");
    unlink $splitbed_fileS,$intersectbed_fileS;

    run("$intersectBed -wo -sorted -a $nosplitbed_file -b $splitbed_fileE > $intersectbed_fileE");
    open(IE,"<$intersectbed_fileE") or abort(2,"Cannot open $intersectbed_fileE: $!");
    while(<IE>){
      my $line = $_;
      chomp($line);
      my @crd = split(/\t/,$line);
      my $start = $crd[1];
      my $end = $crd[2];
      my $name = $crd[3];
      my $breakend = $crd[5];
      next if exists $clustercheck->{$name};
      my $diffstart = $breakend - $start;
      my $diffend = $end - $breakend;
      next if $diffstart <= $breakspan || $diffend <= $breakspan;
      $clustercount->{$clusterid}->{wt2}++;
      $clustercheck->{$name}++;
    }
    close(IE) or abort(2,"Cannot read $intersectbed_fileE: $!");
    unlink $splitbed_fileE,$intersectbed_fileE;
  }
  return $clustercount;
}

sub print_result{
  my $clustercount = shift;
  my $clusterfile = shift;
  open(CF,">$clusterfile") or abort(2,"Cannot open $clusterfile: $!");
  for my $clusterid(keys %{$clustercount}){
    my $names = $clustercount->{$clusterid}->{names};
    my $starts = $clustercount->{$clusterid}->{starts};
    my $ends = $clustercount->{$clusterid}->{ends};
    my $lenstarts = $clustercount->{$clusterid}->{lenstarts};
    my $lenends = $clustercount->{$clusterid}->{lenends};
    my $mt = $clustercount->{$clusterid}->{mt};
    my $wt1 = 0;
    $wt1 = $clustercount->{$clusterid}->{wt1} if exists $clustercount->{$clusterid}->{wt1};
    my $wt2 = 0;
    $wt2 = $clustercount->{$clusterid}->{wt2} if exists $clustercount->{$clusterid}->{wt2};
    my $sum_wt = $wt1 + $wt2;
    my $mean_wt = sprintf "%.0f",$sum_wt/2;

    my $names_print = join("\,",@$names);
    my $starts_print = join("\,",@$starts);
    my $ends_print = join("\,",@$ends);

    #if fragments overlap, then lenstart and and lenend are committed to NA
    my $lenstarts_print = "NA";
    my $lenends_print = "NA";
    if(defined(@$lenstarts[0])){$lenstarts_print = join("\,",@$lenstarts);}
    if(defined(@$lenends[0])){$lenends_print = join("\,",@$lenends);}

    my $perc_hp = $mt*100/($mt+$mean_wt);

    print CF "$clusterid\t$names_print\t$starts_print\t$ends_print\t$lenstarts_print\t$lenends_print\t$mt\t$mean_wt\t$perc_hp\n";
  }
  close(CF) or abort(2,"Cannot write $clusterfile: $!");
}
