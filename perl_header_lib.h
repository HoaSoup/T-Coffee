char *PerlScriptName[]={"rec_sum.pl","count.pl","p\
rocess_list.pl","make_license.pl","CCsed.script","\
msa2bootstrap.pl","tc_generic_method.pl","mmseqs2p\
rf.pl","dynamic.pl","rnapdb2protpdb.pl","generic_m\
ethod.tc_method","clustalw_method.tc_method","extr\
act_from_pdb","install.pl","clean_cache.pl","natur\
e_protocol.pl","mocca","dalilite.pl","wublast.pl",\
"blastpgp.pl","ncbiblast_lwp.pl","wublast_lwp.pl",\
"RNAplfold2tclib.pl","fasta_seq2RNAplfold_template\
file.pl","fasta_seq2hmmtop_fasta.pl","fasta_seq2co\
nsan_aln.pl","clustalw_aln2fasta_aln.pl","seq2name\
_seq.pl","seq2intersection.pl","msf_aln2fasta_aln.\
pl","msa.pl","upp.pl","clustalo.pl","dca.pl","blas\
t_aln2fasta_aln.pl","blast_xml2fasta_aln.pl","fast\
a_aln2fasta_aln_unique_name.pl","newick2name_list.\
pl","excel2fasta.pl","nameseq2fasta.pl","any_file2\
unix_file.pl","EndList"};char *PerlScriptFile[]={"\
use File::Copy;\nuse Env qw(HOST);\nuse Env qw(HOM\
E);\nuse Env qw(USER);\n$x_field=0;\n$y_field=1;\n\
$y_field_set=1;\n$nyf=1;\n\n$interval=0;\n$file=\"\
stdin\";\n\n$print_avg=1;\n$print_sd=0;\n$print_su\
m=0;\n$print_n=0;\nforeach $value ( @ARGV)\n    {\\
n	if ($value ne $ARGV[$np]) \n	    {\n	    ;\n	   \
 }\n	elsif($value eq \"-s\")\n	     {\n	       $st\
ep=$ARGV[++$np];\n	       $np++;\n	     }\n	elsif(\
$value eq \"-print_all\")\n	    {\n	    $print_sd=\
$print_avg=$print_n=$print_sum=1;\n	    $np++;\n	 \
   }\n	elsif($value eq \"-print_sum\")\n	    {\n	 \
   $print_sum=1;\n	    $print_avg=0;\n	    $np++;\\
n	    }\n	elsif($value eq \"-print_n\")\n	    {\n	\
    $print_n=1;\n	    $print_avg=0;\n	    $np++;\n\
	    }\n	elsif($value eq \"-print_avg\")\n	    {\n\
	    $print_avg=1;\n	    $print_avg=0;\n	    $np++\
;\n	    }\n	elsif($value eq \"-sd\")\n	    {\n	   \
 $print_sd=1;\n	    $print_avg=0;\n	    $np++;\n	 \
   }\n	elsif($value eq \"-h\")\n	    {\n	    $head\
er=1;\n	    $np++;\n	    }\n	elsif ($value eq \"-i\
\")\n	    {\n	    $interval= $ARGV[++$np];\n	    $\
np++;\n    	    }\n	elsif ($value eq \"-r2\")\n	  \
  {\n	      $r2=1;\n	      \n	      $np  = $ARGV[+\
+$np];\n	      $nsim= $ARGV[++$np];\n	      $np++;\
\n    	    }\n	elsif ($value eq \"-r\")\n	    {\n	\
    $min= $ARGV[++$np];\n	    $max= $ARGV[++$np];\\
n	    $np++;\n    	    }\n	\n	elsif ($value eq \"-\
x\")\n	    {\n	    $x_field= $ARGV[++$np]-1;\n	   \
 $np++;\n    	    }\n	elsif ($value eq \"-y\")\n	 \
   {\n	    $nyf=0;  \n	    while ($ARGV[$np+1] && \
!($ARGV[$np+1]=~/\\-/))\n	      {\n		$y_field[$nyf\
++]=$ARGV[++$np]-1;\n		$y_field_set=1;\n	      }\n\
\n	    $np++;\n    	    }\n	elsif ($value eq \"-fi\
le\")\n	    {\n	    $file= $ARGV[++$np];\n	    $fi\
le_set=1;\n	    $np++;\n    	    }       \n	elsif \
( $value eq \"h\" ||  $value eq \"-h\" || $value e\
q \"-H\" || $value eq \"-help\" || $value eq \"hel\
p\")\n	  {\n	    print STDOUT \"data_analyse: Anal\
yse and discretization of data\\n\";\n	    print S\
TDOUT \"       -file:    <file containing the data\
 to analyze>,.<def=STDIN>\\n\";\n	    print STDOUT\
 \"       -x: <field containing the X>,...........\
....<Def=0>\\n\";\n	    print STDOUT \"       -y: \
<field containing the Y>,...............<Def=1>\\n\
\";\n	    print STDOUT \"       -i:<Interval size \
on the X>,...............<Def=0>\\n\";\n	    print\
 STDOUT \"       -i:<0:only one interval>\\n\";\n	\
    print STDOUT \"       -r:<Range of the X>\\n\"\
;\n	    print STDOUT \"       -s:<Step on the  X, \
0 means non sliding bins>\\n\";\n	    print STDOUT\
 \"       -sd: print standard deviation on the Y\"\
;\n	    print STDOUT \"       -h  : print column h\
eader \\n\";\n	    exit (0);\n	  }\n	elsif ($value\
=~/-/)\n	  {\n	    print \"$value is not a valid F\
LAG[FATAL]\\n\";\n	    exit (0);\n	   } \n	elsif (\
$list eq \"\") \n	    {\n	    $file=$ARGV[$np];\n	\
    $np++;\n	    }\n	\n	\n      }\n\n\n\n\n\nif ($\
file eq \"stdin\")\n	{\n	$remove_file=1;\n	$file=\\
"tmp$$\";\n	open (F, \">$file\");\n	while (<STDIN>\
)\n		{\n		print F $_;\n		}\n	close (F);\n	 \n	;}\n\
\n\n\nif ($interval && $step)\n  {\n    my $nl;\n \
   open(F,$file);\n    while (<F>)\n      {\n	$lin\
e=$_;\n	\n	if (!/\\S/){next;}\n	@list=($line=~/(\\\
S+)/g);\n	$val{$nl}{x}=$list[$x_field];\n	$val{$nl\
}{y}=$list[$y_field[0]];\n	$nl++\n      }\n    clo\
se (F);\n    \n    for (my $a=$min; $a<($max+$inte\
rval); $a+=$step)\n      {\n	my ($avgx, $avgy, $cn\
);\n	\n	my $rmin=$a-$interval;\n	my $rmax=$a;\n	$c\
n=0;\n	for (my $b=0; $b<$nl; $b++)\n	  {\n	    my \
$x=$val{$b}{x};\n	    my $y=$val{$b}{y};\n	    if \
($x<=$rmax && $x>=$rmin)\n	      {\n		$avgx+=$x;\n\
		$avgy+=$y;\n		$cn++;\n		$tcn++;\n		$val{$b}{used\
}=1;\n	      }\n	  }\n	if ($cn)\n	  {\n	    $avgx/\
=$cn;\n	    $avgy/=$cn;\n	  }\n	printf \"%.3f %.3f\
 %.3f\\n\", $avgx, $avgy, $avgx-$avgy;\n      }\n \
   for (my $a=0; $a<$nl; $a++)\n      {\n	if ( !$v\
al{$a}{used})\n	  {\n	    print \"---$val{$a}{x}; \
$val{$a}{y}\\n\";\n	  }\n      }\n  }\nelse\n  {\n\
    if ($interval && $max)\n      {\n	$interval_si\
ze=($max-$min)/$interval;\n      }\n    elsif ($in\
terval)\n      {\n	open(F,$file);  \n	my $set_max=\
0;\n	my $set_min=0;\n	while (<F>)\n	  {\n	    my $\
v=$_;\n	    chomp($v);\n	    print \"--$v--\";\n	 \
   \n	    if ($v<$min ||!$set_min){$set_min=1;$min\
=$v;}\n	    if ($v>$max ||!$set_max){$set_max=1;$m\
ax=$v;}\n	  }\n	close (F);\n	print \"$min $max uuu\
u\";\n	$interval_size=($max-$min)/$interval;\n    \
  }\n    open(F,$file);  \n    while (<F>)\n      \
{\n	$line=$_;\n	if (!/\\S/){next;}\n	@list=($line=\
~/(\\S+)/g);\n	\n	if ($interval==0){$bin=0;}\n	els\
e{$bin=int (($list[$x_field]-$min)/($interval_size\
));}\n	\n	\n	if ($bin && $bin==$interval){$bin--;}\
\n	for ( $a=0; $a<$nyf; $a++)\n	  {\n	    $sum{$a}\
{$bin}+=$list[$y_field[$a]];\n	    $sum2{$a}{$bin}\
+=$list[$y_field[$a]]*$list[$y_field[$a]];\n	    $\
n{$a}{$bin}++;\n	  }\n      }\n    \n    if (!$int\
erval){$interval=1;}\n    for ( $a=0; $a<$interval\
; $a++)\n      {\n	printf ( \"%4d %4d \", $interva\
l_size*$a, $interval_size*($a+1));\n	for ( $b=0; $\
b<$nyf; $b++)	\n	  {\n	    $i=$interval*$a;\n	    \
if ( $n{$b}{$a}==0)\n	      {\n		$avg=0;\n		$sd=0;\
\n	      }\n	    else\n	      {\n		$avg=$sum{$b}{$\
a}/$n{$b}{$a};\n		$sd=sqrt($sum2{$b}{$a}*$n{$b}{$a\
}-$sum{$b}{$a}*$sum{$b}{$a})/($n{$b}{$a}*$n{$b}{$a\
});\n	      }\n	    if ($print_n) {printf \"%15.4f\
 \", $n{$b}{$a};}\n	    if ($print_sum){printf \"%\
15.4f \", $sum{$b}{$a};}\n	    if ($print_avg){pri\
ntf \"%15.4f \", $avg}\n	    if ($print_sd) {print\
f \"%15.4f \", $sd;}\n	  }\n	printf (\"\\n\");\n  \
    }\n  }\n\nif ( $remove_file){unlink $file;}\n"\
,"use File::Copy;\nuse Env qw(HOST);\nuse Env qw(H\
OME);\nuse Env qw(USER);\n\nforeach $v (@ARGV){$cl\
.=$v;}\n\n\nif ( $cl=~/-k(\\d+)/){$k=$1;}\nelse {$\
k=1;}\nif ( $cl=~/-w(\\d+)/){$w=$1;}\nelse {$w=-1;\
}\nif ( $cl=~/-p(\\d+)/){$p=$1;}\nelse {$p=-1;}\n\\
nwhile (<STDIN>)\n  {\n    @l=($_=~/(\\S+)/g);\n  \
  $v=$l[$k-1];\n    if ( !$h{$v}){@ll=($v, @ll);}\\
n    \n    if ( $w==-1)\n      {$h{$v}++;}\n    el\
se\n      {$h{$v}+=$l[$w-1];}\n\n    if ($p!=-1){$\
print{$v}=$l[$p-1];}\n\n  }\nforeach $v (@ll)\n  {\
\n    print \"$v $print{$v} $h{$v}\\n\";\n  }\n","\
\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw\
(USER);\n$random_tag=int (rand 10000)+1;\n$unique_\
prefix=\"$$.$HOST.$random_tag\";\n$queue=\"distill\
ery.and.mid\";\n$monitor=0;\n$stderr_file=\"/dev/n\
ull\";\n$stdio_file=\"/dev/null\";\n$log_file=\"/d\
ev/null\";\n$pause_time=0;\n$max_sub_jobs=60;\n$mi\
n_sub_jobs=30;\n$output_all=0;\n$var='\\$';\n\nfor\
each $value ( @ARGV)\n    {\n	if ($value ne $ARGV[\
$np]) \n	    {\n	    ;\n	    }\n	elsif ($value eq \
\"-max_sub_jobs\")\n	    {\n	    $max_sub_jobs= $A\
RGV[++$np];\n	    $np++;\n    	    }	\n	elsif ($va\
lue eq \"-min_sub_jobs\" )\n	    {\n	    $min_sub_\
jobs= $ARGV[++$np];\n	    $np++;\n    	    }\n	els\
if ($value eq \"-para\")\n	    {\n	    $para=1;\n	\
    $monitor=1;\n	    $np++;\n    	    }\n	elsif (\
$value eq \"-monitor\") \n	    {\n	    $monitor=1;\
\n	    $np++;\n	    }\n	elsif ($value eq \"-no_mon\
itor\") \n	    {\n	    $monitor=0;\n	    $np++;\n	\
    }\n	elsif ($value eq \"-queue\")\n	    {\n	   \
 $queue=$ARGV[++$np];\n	    $np++;\n	    }	\n	elsi\
f ($value eq \"-stderr_file\")\n	    {\n	    $stde\
rr_file=$ARGV[++$np];\n	    $np++;\n	    }\n	elsif\
 ($value eq \"-stdio_file\")\n	    {\n	    $stdio_\
file=$ARGV[++$np];\n	    $np++;\n	    }\n	elsif ($\
value eq \"-output_all\")\n	    {\n	    $output_al\
l=1;\n	    $np++;\n	    }\n	elsif ($value eq \"-pa\
use\") \n	    {\n	    $pause_time=$ARGV[++$np];\n	\
    $np++;\n	    }\n	elsif ($value eq \"-log\")\n	\
      {\n	       $log=1;\n	       \n	       if ($A\
RGV[$np+1]=~/\\-\\S+/) \n	          {\n		  $log_fi\
le=\"stderr\";\n	          }\n	       else \n	    \
      {\n		  $log_file=$ARGV[++$np]; \n		  ++$np;\\
n		 \n	          }\n	      }\n	elsif ( $value eq \\
"-com\")\n	    {\n		\n		if (!$ARGV[$np+1]=~/^\\'/)\
 { $com=$ARGV[++$np];}\n		else {$com=$ARGV[++$np];\
}\n\n	     $np++;\n	    }\n	elsif ( $value eq \"-c\
heck\")\n	  {\n	    \n	    if (!$ARGV[$np+1]=~/^\\\
'/) { $check=$ARGV[++$np];}\n	    else {$check=$AR\
GV[++$np];}\n	    $np++;\n	  }\n	elsif ($com eq \"\
\") \n	    {\n	    $com_set=1;\n	    $com=$ARGV[$n\
p];\n	    \n	    $np++;\n	    }\n	elsif ($list eq \
\"\") \n	    {\n	    $list_set=1;\n	    $list=$ARG\
V[$np];\n	    $np++;\n	    }\n	elsif ( $var_set eq\
 \"\")\n	    {\n	    $var_set=1;\n	    $var=$ARGV[\
$np];\n	    $np++;\n	    }\n	}\n\n\n\n\nif ( $com \
eq \"\"){print \"You Need to Provide a Command [FA\
TAL]\\n\";\n	      die;\n	     }\n\n\n\nif ($list_\
set==0) \n    {\n    $x= int (rand 100000)+1;\n   \
 $tmp_file_name=\"tmp_file_$x\";\n    open ( TMP, \
\">$tmp_file_name\");\n    while (<STDIN>)\n      \
{\n	print TMP $_;\n      }\n    close (TMP);\n    \
open (F, $tmp_file_name);\n    }\nelse \n    {\n  \
  open (F, $list);\n    }\n\nif ($para==0) \n    {\
\n\n     @tc_list= <F>;\n     close (F); \n     \n\
     foreach $val(@tc_list) \n	    {\n	      \n	  \
    \n	      \n	      $loc_com=$com;\n	      if ($\
check){$loc_check=$check;}\n	      \n	      @i_val\
=($val=~/([^\\s]+)/g);\n	      \n	      if ( $#i_v\
al==0)\n		{\n		  if ($check){$loc_check=~s/$var/$i\
_val[0]/g;}\n		  $loc_com=~s/$var/$i_val[0]/g;\n		\
}\n	      else\n		{\n		  for ($n=1; $n<=$#i_val+1;\
$n++ )\n		    {\n		      \n		      $sub=\"$var$n\"\
;\n		      \n		      $loc_com=~s/$sub/$i_val[$n-1]\
/g;\n		      if ($check){$loc_check=~s/$var/$i_val\
[0]/g;}\n		    }\n		}\n	      if ( $check && -e $l\
oc_check)\n		{\n		  print STDERR \"skipping $loc_c\
om...\\n\";\n		  }\n	      else\n		{\n		  system \\
"$loc_com\";\n		}\n	    }\n    exit;\n    }\n\nels\
if ($para==1) \n    {\n    print STDERR \"do paral\
lel execution of: \\\"$com $list\\\"\\n\";\n    \n\
    if ($log==1) \n	{\n	if ($log_file eq \"stdout\\
" || $log_file eq \"stderr\" ) \n		{\n		$log_file=\
\"\";\n	        }\n\n        else \n		{\n		system \
\"echo LOG FILE> $log_file\";\n		\n	        }\n	}\\
n    else	\n	{\n	open ( OUT, \">/dev/null\");\n	}\\
n	\n    \n    $id=0;\n    $n_sub=0;\n    while ($v\
al=<F>) \n	    {	    	    \n	    $job_log[$id]=\"$\
HOME/tmp/$unique_prefix.$id.log_file\";\n	    \n	 \
   $job=$unique_prefix.\"_$id\";\n	    open (JOB, \
\">$job\");\n	    \n	    $loc_com=$com;\n	    chop\
 $val;\n\n	    $loc_com=~s/\\$/$val/g;\n	 \n	    p\
rint JOB \"#!/bin/csh\\n\";\n	    print JOB \"#\\$\
 -cwd\\n\";\n	    print JOB \"#\\$ -N $unique_pref\
ix\\n\";\n	    if ($queue && !($queue eq \" \")) {\
print JOB \"#\\$ -l $queue\\n\";}\n	    print JOB \
\"#\\n\";	    \n            print JOB \"$loc_com\\\
n\";\n	    print JOB \"echo FINISHED  >> $job_log[\
$id]\\n\";\n	    print JOB \"pwd\\n\";\n	    \n	  \
  close (JOB);\n	    if ( $output_all==1)\n		{\n		\
system \"qsub $job >  $unique_prefix\";		\n	      \
  }\n	    else\n		{system \"qsub $job -e $stderr_f\
ile -o $stdio_file >$unique_prefix\";	        \n	 \
       } \n\n\n\n	    print STDERR \"$id: $output_\
all\\n\";\n	    $n_sub++;\n	    if ( $max_sub_jobs\
 && $n_sub==$max_sub_jobs) \n		{\n		$n_sub=monitor\
_process($min_sub_jobs,@job_log); 		 \n		\n	      \
  }	\n	   \n            unlink $unique_prefix;\n	 \
   sleep $pause_time;\n	    $id++;\n	    }\n\n    \
close (OUT);\n    close (F);\n\n    print STDERR \\
"Your $id Jobs Have Been Submited (NAME=$unique_pr\
efix)\\n\";\n    monitor_process (0, @job_log);\n \
   foreach $file(@job_log) {if (-e $file) {unlink(\
$file);}}\n    \n    }\n\nsub monitor_process ( @j\
ob_list)\n    {\n    my (@job_list)=@_;\n    my $m\
in_sub_jobs=shift (@job_list);\n    my $n_sub_jobs\
;\n    my $finished;\n    my $n=0;\n\n    $n_sub_j\
obs=-1;\n    $finished=0;\n    print STDERR \"\\nM\
onitor Batch: [$min_sub_jobs]\";\n       \n    whi\
le (!$finished && (($n_sub_jobs>$min_sub_jobs)|| $\
n_sub_jobs==-1) ) \n	{\n	$finished=1;\n	$n_sub_job\
s=0;\n	$n=0;\n	foreach $file (@job_list)\n	       \
 {\n	\n		if (-e $file){;}\n		else \n		    {\n		   \
 $finished=0; $n_sub_jobs++;\n	            }\n	   \
     }\n	system \"sleep 1\";\n        }\n    \n   \
 return $n_sub_jobs;\n    }\n    \n    \nif ($tmp_\
file_name){unlink($tmp_file_name);}\n","\n\nforeac\
h ($np=0; $np<=$#ARGV; $np++)\n    {\n    $value=$\
ARGV[$np];\n\n    if ($value eq \"-file\")\n      \
{\n      $file= $ARGV[++$np];\n      }\n    elsif \
($value eq \"-type\")\n      {\n        $type= $AR\
GV[++$np];\n      }\n    elsif ($value eq \"-insti\
tute\")\n      {\n        $institute= $ARGV[++$np]\
;\n      }\n    elsif ($value eq \"-author\")\n   \
   {\n        $author= $ARGV[++$np];\n      }\n   \
 elsif ($value eq \"-date\")\n      {\n        $da\
te= $ARGV[++$np];\n      }\n     elsif ($value eq \
\"-program\")\n      {\n        $program= $ARGV[++\
$np];\n      }\n    elsif ($value eq \"-email\")\n\
      {\n        $email= $ARGV[++$np];\n      }\n \
   else\n      {\n	print \"$value is an unkown arg\
ument[FATAL]\\n\";\n	exit (1);\n      }\n  }\n\n\n\
\nopen F, $file || die;\nprint $INSTITUTE;\nif ( $\
type eq \"c\"){print \"/**************************\
****COPYRIGHT NOTICE******************************\
*/\\n\";}\nif ( $type eq \"perl\"){print \"#######\
#######################COPYRIGHT NOTICE###########\
###################/\\n\";}\nif ( $type eq \"txt\"\
){print \"-------------------------------COPYRIGHT\
 NOTICE------------------------------/\\n\";}\n\n\\
nwhile (<F>)\n  {\n  s/\\$INSTITUTE/$institute/g;\\
n  s/\\$AUTHOR/$author/g;\n  s/\\$DATE/$date/g;\n \
 s/\\$PROGRAM/$program/g;  \n  s/\\$EMAIL/$email/g\
;  \n  if ( $type eq \"txt\"){print $_;}\n  elsif \
($type eq \"c\"){chop $_; print \"\\/*$_*\\/\\n\";\
}\n  elsif ($type eq \"perl\"){print \"\\#$_\";}\n\
}\nclose (F);\nif ( $type eq \"c\"){print \"/*****\
*************************COPYRIGHT NOTICE*********\
**********************/\\n\";}\nif ( $type eq \"pe\
rl\"){print \"##############################COPYRI\
GHT NOTICE##############################/\\n\";}\n\
if ( $type eq \"txt\"){print \"-------------------\
------------COPYRIGHT NOTICE----------------------\
--------/\\n\";}\n\n","\nwhile (<>)	\n	{\n	s/\\=cc\
/123456789/g;\n	s/\\bcc/\\$\\(CC\\)/g;\n	s/1234567\
89/\\=cc/g;\n	print $_;\n	}\n\n","$version=\"1.00\\
";\n$rseed= int(rand(100000))+1;\n\n\nif ( $#ARGV=\
=-1)\n  {\n    print \"msa2bootstrap -i <input_fil\
e> -input <seq|msa|matrix|tree> -n <N-Boostrap> -o\
 <outtree> -tmode <nj|upgma|parsimony|ml> -dmode <\
kimura> -alignpg <t_coffee | muscle | clustalw> -r\
tree <file> -stype <prot|cdna|dna> -recompute -sys\
tem <cygwin|unix>\";\n    print \"\\n\\t-i: input \
file, can be sequneces, msa, matrix, trees, type i\
s specified via -input\";\n    print \"\\n\\t-inpu\
t: Type of input data\";\n    print \"\\n\\t\\tmsa\
: msa in fasta format\";\n    print \"\\n\\t\\tseq\
: compute an msa with -alignpg\";\n    print \"\\n\
\\t\\tmatrix: phylipp distance matrix fed directly\
 to method -tmode [caveat: tmode=nj or upgma]\";\n\
    print \"\\n\\t\\ttree: list of newick trees di\
rectly fed to consence in order to generate a boot\
straped tree\";\n    \n    print \"\\n\\t-n: numbe\
r of bootstrap replicates\";\n    print \"\\n\\t-o\
: name of the output tree. Files are not overwritt\
en. Use -recompute to overwrite existing file\";\n\
    print \"\\n\\t-tmode: tree mode: nj|upgma|pars\
imony|ml\";\n    print \"\\n\\t-dmode: distance mo\
de\";\n    print \"\\n\\t-alignpg: program for ali\
gning sequences (t_coffee=default)\";\n    print \\
"\\n\\t-rtree: replicate tree file (default: no fi\
le)\";\n    print \"\\n\\t-rmsa: replicate msa fil\
e (default: no file)\";\n    print \"\\n\\t-rmat: \
replicate matrix file (default: no file)\";\n    p\
rint \"\\n\\t-stype: sequence type: protein, dna o\
r cdna\";\n    print \"\\n\\t-recompute: force fil\
es to be overwritten\";\n    print \"\\n\\t-system\
: cygwin|unix\";\n      \n\n    \n    &my_exit (EX\
IT_FAILURE);\n  }\nforeach $arg (@ARGV){$command.=\
\"$arg \";}\n\nprint \"CLINE: $command\\n\";\n$thr\
eshold=100;\n$trim_msa=0;\n$stype=\"prot\";\nprint\
 \"msa2bootstrap \";\n\n$system=\"cygwin\";\nif(($\
command=~/\\-system (\\S+)/))\n  {\n    $system=$1\
;\n    if ( $system eq \"cygwin\")\n      {\n	$exe\
c_extension=\".exe\";\n      }\n    elsif ( $syste\
m eq \"unix\")\n      {\n	$exec_extension=\"\";\n	\
print \"system=Unix\";die;\n      }\n    else\n   \
   {\n	print \"msa2boostrap: -system=$system is an\
 unknown mode [FATAL]\\n\"; die;\n      }\n    \n \
   print \"-system $system \";\n  }\nif(($command=\
~/\\-stype (\\S+)/))\n  {\n    $stype=$1;\n  }\npr\
int \"-stype=$stype \";\n\n\n\nif(($command=~/\\-i\
 (\\S+)/))\n  {\n    $msa=$1;\n    print \"-i $msa\
 \";\n  }\n\nif(($command=~/\\-rtree (\\S+)/))\n  \
{\n    $rtree=$1;\n    print \"-rtree=$rtree \";\n\
  }\n\nif(($command=~/\\-rmsa (\\S+)/))\n  {\n    \
$rmsa=$1;\n  }\nif(($command=~/\\-rmat (\\S+)/))\n\
  {\n    $rmat=$1;\n  }\n$input=\"seq\";\nif(($com\
mand=~/\\-input (\\S+)/))\n  {\n    $input=$1;\n  \
}\nprint \"-input=$input \";\n\n$dmode=\"kimura\";\
\nif(($command=~/\\-dmode (\\S+)/))\n  {\n    $dmo\
de=$1;\n  }\nprint \"-dmode=$dmode \";\n$alignpg=\\
"muscle\";\nif(($command=~/\\-alignpg (\\S+)/))\n \
 {\n    $alignpg=$1;\n  }\nprint \"-alignpg=$dmode\
 \";\n\n$tmode=\"nj\";\nif(($command=~/\\-tmode (\\
\S+)/))\n  {\n    $tmode=$1;\n  }\nprint \"-tmode=\
$tmode \";\n$recompute=0;\nif(($command=~/\\-recom\
pute/))\n  {\n    $recompute=1;\n    print \"-reco\
mpute \";\n  }\n\n$out=$msa;\n$out=~s/\\..*//;\n$o\
ut.=\".bph\";\nif(($command=~/\\-o (\\S+)/))\n  {\\
n    $out=$1;\n    \n  }\nprint \"-out=$out \";\ni\
f (-e $out && !$recompute)\n  {\n    print \"\\nNo\
 Computation Required $out already exists\\n\";\n \
   &my_exit (EXIT_SUCCESS);\n    \n  }\n\n$n=100;\\
nif(($command=~/\\-n (\\d+)/))\n  {\n    $n=$1;\n \
 }\nprint \"-n=$n \";\n$seed=3;\nif(($command=~/\\\
-s (\\d+)/))\n  {\n    $seed=$1;\n  }\nprint \"-s=\
$seed \";\n\nif(($command=~/\\-run_name (\\d+)/))\\
n  {\n    $suffix=$1;\n  }\nelse\n  {\n    $msa=~/\
([^.]+)/;\n    $suffix=$1;\n  }\nprint \"-run_name\
=$suffix\\n\";\n\n\nif ( $input eq \"seq\")\n  {\n\
    $seq=$msa;\n    $msa=\"$suffix.prot_msa\";\n  \
  \n    if ($stype eq \"cdna\")\n      {\n	$cdna_s\
eq=$seq;\n	$clean_cdna_seq=&vtmpnam();\n	$seq=&vtm\
pnam();\n	`t_coffee -other_pg seq_reformat -in $cd\
na_seq -action +clean_cdna >$clean_cdna_seq`;\n	`t\
_coffee -other_pg seq_reformat -in $clean_cdna_seq\
 -action +translate >$seq`;\n	\n      }\n\n    if \
(!-e $msa || $recompute)\n      {\n	print \"\\n###\
##   Compute an MSA With $alignpg\\n\";\n	\n	if ( \
$alignpg eq \"t_coffee\")\n	  {`$alignpg $seq -out\
file=$msa >/dev/null 2>/dev/null`;}\n	elsif ( $ali\
gnpg eq \"muscle\")\n	  {\n	    `$alignpg -in $seq\
 > $msa 2>/dev/null`;\n	  }\n	elsif ( $alignpg eq \
\"clustalw\")\n	  {\n	    `$alignpg -infile=$seq -\
outfile=$msa -quicktree >/dev/null 2>/dev/null`;\n\
	  }\n	elsif ( $align eq \"mafft\")\n	  {\n	    `$\
alignpg $seq > $msa >/dev/null 2>/dev/null`;\n	  }\
\n	else\n	  {\n	    `$alignpg -in=$seq -outfile=$m\
sa`;\n	  }\n      }\n    if (!-e $msa)\n      {\n	\
print \"\\nError: $alignpg Could Not produce the M\
SA $msa [FATAL]\\n\";\n      }\n\n    if ($stype e\
q \"cdna\")\n      {\n	$msa2=\"$suffix.cdna_msa\";\
\n	`t_coffee -other_pg seq_reformat -in $clean_cdn\
a_seq -in2 $msa -action +thread_dna_on_prot_aln -o\
utput fasta_aln  >$msa2`;\n	$msa=$msa2;\n      }\n\
    \n    $input=\"msa\";\n  }\n\n\n\n$seqboot_o=&\
vtmpnam();\n$seqboot_c=&vtmpnam();\n\n$protdist_o=\
&vtmpnam();\n$protdist_c=&vtmpnam();\nif ( $input \
eq \"msa\")\n  {\n    if ($tmode eq \"nj\" || $tmo\
de eq \"upgma\"){$input=\"matrix\";}\n    \n    $l\
msa= &vtmpnam ();\n    `t_coffee -other_pg seq_ref\
ormat -in $msa -output phylip_aln > $lmsa`;\n    \\
n    if ( -e \"outfile\"){unlink (\"outfile\");}\n\
    # run seqboot\n  \n    if ( $n>1)\n      {\n	p\
rint \"Run SeqBoot .....\";\n	open (F, \">$seqboot\
_c\");\n	print F \"$lmsa\\nR\\n$n\\nY\\n$seed\\n\"\
;\n	close (F);\n	`seqboot$exec_extension  < $seqbo\
ot_c`;\n	if ( -e \"outfile\"){ print \"[OK]\\n\";}\
\n	else { print \"[FAILED]\\n\";&my_exit (EXIT_FAI\
LURE);}\n	`mv outfile $seqboot_o`;\n      }\n    e\
lse\n      {\n	`cp $lmsa $seqboot_o`;\n      }\n\n\
    if ($rmsa){`cp $seqboot_o $rmsa`;}\n    \n    \
if ($tmode eq \"nj\" || $tmode eq \"upgma\")\n    \
  {\n	if ( $stype eq \"prot\")\n	  {\n	    # run p\
rotdist\n	    print \"Run Protdist [dmode=$dmode]\\
";\n	    if ($dmode eq \"kimura\")\n	      {\n		$d\
mode=\"P\\nP\\nP\";\n	      }\n	    else\n	      {\
\n		print \"\\n$dmode is an unknown mode for Protd\
ist [FATAL:msa2bootstrap.pl]\\n\";\n		&my_exit (EX\
IT_FAILURE);\n	      }\n	    open (F, \">$protdist\
_c\");\n	    if ($n>1){print F \"$seqboot_o\\n$dmo\
de\\nM\\nD\\n$n\\nY\\n\";}\n	    else {printf F \"\
$seqboot_o\\n$dmode\\nY\\n\";}\n	    close (F);\n	\
    `protdist$exec_extension  < $protdist_c`;\n	  \
  if ( -e \"outfile\"){ print \"[OK]\\n\";}\n	    \
else { print \"[FAILED]\\n\";&my_exit (EXIT_FAILUR\
E);}\n	    `mv outfile $protdist_o`;\n	 \n	  }\n	e\
lsif ( $stype eq \"cdna\" || $stype eq \"dna\")\n	\
  {\n	    print \"Run dnadist [dmode=default\";\n	\
    open (F, \">$protdist_c\");\n	    if ($n>1){pr\
int F \"$seqboot_o\\nM\\nD\\n$n\\nY\\n\";}\n	    e\
lse {printf F \"$seqboot_o\\nY\\n\";}\n	    close \
(F);\n	    `protdist$exec_extension  < $protdist_c\
`;\n	    if ( -e \"outfile\"){ print \"[OK]\\n\";}\
\n	    else { print \"[FAILED]\\n\";&my_exit (EXIT\
_FAILURE);}\n	    `mv outfile $protdist_o`;\n	  }\\
n      }\n  }\nelsif ( $input eq \"matrix\")\n  {\\
n    $protdist_o=&vtmpnam();\n    print \"MSA: $ms\
a\\n\";\n    `cp $msa $protdist_o`;\n    $n=1;\n  \
}\n\n\n\n\n\n$nb_o=&vtmpnam();\n$nb_c=&vtmpnam();\\
nif ($input eq \"matrix\" && $tmode ne \"parsimony\
\" && $tmode ne \"ml\")\n  {\n    print \"Run neig\
hbor [tmode=$tmode]\";\n\n    if ($tmode eq \"nj\"\
)\n      {\n	$tmode=\"\\nN\\nN\";\n      }\n    el\
sif ( $tmode eq \"upgma\")\n      {\n	$tmode = \"\\
\nN\";\n      }\n    else\n      {\n	print \"\\n E\
RROR: $tmode is an unknown tree computation mode\\\
n\";\n	&my_exit (EXIT_FAILURE);\n      }\n\n    op\
en (F, \">$nb_c\");\n    if ($n>1){print F \"$prot\
dist_o$tmode\\nM\\n$n\\n$seed\\nY\\n\";}\n    else\
 {print F \"$protdist_o$tmode\\nY\\n\";}\n    clos\
e (F);\n\n    `neighbor$exec_extension  < $nb_c`;\\
n    if ( -e \"outtree\"){ print \"[Neighbor OK]\\\
n\";}\n    else { print \"[FAILED]\\n\";&my_exit (\
EXIT_FAILURE);}\n    `mv outtree $nb_o`;\n    unli\
nk (\"outfile\");\n  }\nelsif ($input eq \"msa\" &\
& $tmode eq \"parsimony\")\n  {\n    if ( -e \"out\
file\"){unlink (\"outfile\");}\n    if ( -e \"outt\
ree\"){unlink (\"outtree\");}\n    \n    if ($styp\
e eq \"prot\")\n      {\n	print \"Run protpars [tm\
ode=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1)\
{print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\n\
Y\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n\
	close (F);\n	`protpars$exec_extension  < $nb_c`;\\
n      }\n    elsif ( $stype eq \"dna\" || $stype \
eq \"cdna\")\n      {\n	print \"Run dnapars [tmode\
=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){pr\
int F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\\
n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	cl\
ose (F);\n	`dnapars$exec_extension  < $nb_c`;\n   \
   }\n    if ( -e \"outtree\"){ print \"[OK]\\n\";\
}\n    else { print \"[FAILED]\\n\";&my_exit (EXIT\
_FAILURE);}\n    `mv outtree $nb_o`;\n   unlink (\\
"outfile\");\n  }\nelsif ($input eq \"msa\" && $tm\
ode eq \"ml\")\n  {\n    if ( -e \"outfile\"){unli\
nk (\"outfile\");}\n    if ( -e \"outtree\"){unlin\
k (\"outtree\");}\n    \n    if ($stype eq \"prot\\
")\n      {\n	print \"Error: ML impossible with Pr\
otein Sequences [ERROR]\";\n	&my_exit (EXIT_FAILUR\
E);\n      }\n    elsif ( $stype eq \"dna\" || $st\
ype eq \"cdna\")\n      {\n	print \"Run dnaml [tmo\
de=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){\
print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\
\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	\
close (F);\n	`dnaml$exec_extension  < $nb_c`;\n   \
   }\n    if ( -e \"outtree\"){ print \"[OK]\\n\";\
}\n    else { print \"[FAILED]\\n\";&my_exit (EXIT\
_FAILURE);}\n    `mv outtree $nb_o`;\n   unlink (\\
"outfile\");\n  }\n\n\nelse\n  {\n    `cp $msa $nb\
_o`;\n    $n=2;\n  }\n\nif ($rmsa && -e $seqboot_o\
){print \"\\nOutput List of $n Replicate MSA: $rms\
a\\n\";`cp $seqboot_o $rmsa`;}\nif ($rmat && -e $p\
rotdist_o){print \"\\nOutput List of $n Replicate \
MATRICES: $rmat\\n\";`cp $protdist_o $rmat`;}\nif \
($rtree && -e $nb_o){print \"\\nOutput List of $n \
Replicate TREES: $rtree\\n\";`cp $nb_o $rtree`;}\n\
\n\n\n$con_o=&vtmpnam();\n$con_c=&vtmpnam();\nif (\
$n >1)\n  {\n    print \"Run Consense.....\";\n   \
 open (F, \">$con_c\");\n    print F \"$nb_o\\nY\\\
n\";\n    close (F);\n    `consense$exec_extension\
  < $con_c`;\n    if ( -s \"outtree\"  > 0) { prin\
t \"[OK]\\n\";}\n    else { print \"[FAILED]\\n\";\
&my_exit (EXIT_FAILURE);}\n    `mv outtree $con_o`\
;\n    unlink (\"outfile\");\n  }\nelse\n  {\n    \
`cp $nb_o $con_o`;\n  }\n\n\n`cp $con_o $out`;\nif\
 ( !-e $out)\n  {\n    print \"Tree Computation fa\
iled [FAILED]\\n\";\n    &my_exit (EXIT_FAILURE);\\
n  }\nelsif ($n>1)\n  {\n    print \"\\nOutput Boo\
tstrapped Tree: $out\\n\";\n    $avg=`t_coffee -ot\
her_pg seq_reformat -in $out -action +avg_bootstra\
p`;\n    $avg=~s/\\n//g;\n    print \"$avg\\n\";\n\
  }\nelse\n  {\n    print \"\\nOutput Tree: $out\\\
n\";\n  }\n\nopen (F, \"$out\");\nwhile (<F>)\n  {\
\n    \n    $tree.=$_;\n  }\nclose (F);\n$tree=~s/\
\\n//g;\nprint \"BPH: $tree\\n\";\n\n\n&my_exit (E\
XIT_SUCCESS);\n\nsub my_exit \n  {\n    my $m=@_[0\
];\n    &clean_vtmpnam();\n    exit ($m);\n  }\nsu\
b vtmpnam \n  {\n    my $file;\n\n\n    $ntmp++;\n\
    $file=\"tmp4msa2bootstrap.$rseed.$$.$ntmp\";\n\
    \n    push (@tmpfile, $file);\n    return $fil\
e;\n  }\nsub clean_vtmpnam \n  {\n    my $t;\n    \
foreach $t (@tmpfile)\n      {\n	if ( -e $t){unlin\
k ($t)};\n      }\n  }\n","use Env;\nuse FileHandl\
e;\nuse Cwd;\nuse File::Path;\nuse Sys::Hostname;\\
n\n\nour $PIDCHILD;\nour $ERROR_DONE;\nour @TMPFIL\
E_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=0\
;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCESS=0;\nou\
r $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_generic_me\
thod.pl\";\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT_S\
TARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\"};\nou\
r $debug_generic_method=$ENV{\"DEBUG_GENERIC_METHO\
D\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\n\
if (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=\
$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV\
{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (i\
sshellpid(getppid())){lock4tc(getppid(), \"LLOCK\"\
, \"LSET\", \"$$\\n\");}\nour %RECODE;\nour $RECOD\
E_N;\n\n\n\n\nour $BLAST_MAX_NRUNS=2;\nour $COMMAN\
D;\nour $PIDCHILD;\n\n$REF_EMAIL=\"\";\n$tmp_dir=\\
"\";\n$init_dir=\"\";\n\n\n$test=0;\nif ($test==1)\
\n  {\n    $SERVER=\"NCBI\";\n    $query=$ARGV[0];\
\n    $hitf=$ARGV[1];\n    %s=read_fasta_seq($quer\
y);\n    @sl=keys(%s);\n    &blast_xml2profile (\"\
xx\", $s{$sl[0]}{seq},$maxid,$minid,$mincov, $hitf\
);\n    myexit ($EXIT_FAILURE);\n  }\n\nforeach $v\
(@ARGV){$cl.=\"$v \";}\n$COMMAND=$cl;\n($mode)=&my\
_get_opt ( $cl, \"-mode=\",1,0);\n\n($A)=(&my_get_\
opt ( $cl, \"-name1=\",0,0));\n($B)=(&my_get_opt (\
 $cl, \"-name2=\",0,0));\n($TMPDIR)=(&my_get_opt (\
 $cl, \"-tmpdir=\",0,0));\n($CACHE)=(&my_get_opt (\
 $cl, \"-cache=\",0,0));\n($SERVER)=((&my_get_opt \
( $cl, \"-server=\",0,0)));\n($EMAIL)=((&my_get_op\
t ( $cl, \"-email=\",0,0)));\n\nif (!$A){$A=\"A\";\
}\nif (!$B){$B=\"B\";}\n\n\nif (!$TMPDIR)\n  {\n  \
  $HOME=$ENV{HOME};\n    if ($ENV{TMP_4_TCOFFEE}){\
$TMPDIR=$ENV{TMP_4_TCOFFEE};}\n    else{$TMPDIR=\"\
$HOME/.t_coffee/tmp/\";}\n  }\nif ( ! -d $TMPDIR)\\
n  {\n    mkdir $TMPDIR;\n  }\nif ( ! -d $TMPDIR)\\
n  {\n    print \"ERROR: Could not create temporar\
y dir: $TMPDIR\\n\";\n    myexit ($EXIT_FAILURE);\\
n  }\n\n$EMAIL=~s/XEMAILX/\\@/g;\nif (!$EMAIL)\n  \
{\n    if ($ENV{EMAIL_4_TCOFFEE}){$EMAIL=$ENV{EMAI\
L_4_TCOFFEE};}\n    elsif ($ENV{EMAIL}){$EMAIL=$EN\
V{EMAIL};}\n    else {$EMAIL=$REF_EMAIL;}\n  }\n\n\
($maxid,$minid,$mincov,$trim)=(&my_get_opt ( $cl, \
\"-maxid=\",0,0, \"-minid=\",0,0,\"-mincov=\",0,0,\
 \"-trim=\",0,0));\n\nif (!$cl=~/\\-maxid\\=/){$ma\
xid=95;}\nif (!$cl=~/\\-minid\\=/){$minid=35;}\nif\
 (!$cl=~/\\-mincov\\=/){$mincov=80;}\nif (!$cl=~/\\
\-trim\\=/){$trim;}\n\n\n\n\nif ($mode eq \"seq_ms\
a\")\n  {\n    &seq2msa($mode,&my_get_opt ( $cl, \\
"-infile=\",1,1, \"-method=\",1,2, \"-param=\",0,0\
,\"-outfile=\",1,0, \"-database=\",0,0));\n  }\nel\
sif ($mode eq \"blast2prf\")\n  {\n\n    blast2prf\
 (&my_get_opt ( $cl, \"-infile=\",0,0,\"-seqfile=\\
",0,0,\"-outfile=\",0,0));\n  }\nelsif ( $mode eq \
\"tblastx_msa\")\n  {\n    &seq2tblastx_lib ($mode\
,&my_get_opt ( $cl, \"-infile=\",1,1, \"-outfile=\\
",1,0));\n  }\nelsif ( $mode eq \"tblastpx_msa\")\\
n  {\n    &seq2tblastpx_lib ($mode,&my_get_opt ( $\
cl, \"-infile=\",1,1, \"-outfile=\",1,0));\n  }\ne\
lsif ( $mode eq \"thread_pair\")\n  {\n    &seq2th\
read_pair($mode,&my_get_opt ( $cl, \"-infile=\",1,\
1, \"-pdbfile1=\",1,1, \"-method=\",1,2,\"-param=\\
",0,0, \"-outfile=\",1,0, ));\n  }\nelsif ( $mode \
eq \"pdbid_pair\")\n  {\n    &seq2pdbid_pair($mode\
,&my_get_opt ( $cl, \"-pdbfile1=\",1,0, \"-pdbfile\
2=\",1,0, \"-method=\",1,2,\"-param=\",0,0, \"-out\
file=\",1,0, ));\n  }\nelsif ( $mode eq \"pdb_pair\
\")\n  {\n    &seq2pdb_pair($mode,&my_get_opt ( $c\
l, \"-pdbfile1=\",1,1, \"-pdbfile2=\",1,1, \"-meth\
od=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0, ));\\
n  }\nelsif ( $mode eq \"rnapdb_pair\")\n{\n    &s\
eq2rnapdb_pair($mode,&my_get_opt ( $cl, \"-pdbfile\
1=\",1,1, \"-pdbfile2=\",1,1, \"-method=\",1,2,\"-\
param=\",0,0, \"-outfile=\",1,0, ));\n}\nelsif ( $\
mode eq \"profile_pair\")\n  {\n     &seq2profile_\
pair($mode,&my_get_opt ( $cl, \"-profile1=\",1,1, \
\"-profile2=\",1,1, \"-method=\",1,2,\"-param=\",0\
,0, \"-outfile=\",1,0 ));\n  }\nelsif ($mode eq \"\
pdb_template_test\")\n  {\n    &blast2pdb_template\
_test ($mode,&my_get_opt ( $cl, \"-infile=\",1,1))\
;\n\n  }\nelsif ($mode eq \"psi_template_test\")\n\
  {\n    &psiblast2profile_template_test ($mode,&m\
y_get_opt ( $cl, \"-seq=\",1,1,\"-blast=\",1,1));\\
n\n  }\n\nelsif ( $mode eq \"pdb_template\")\n  {\\
n    &blast2pdb_template ($mode,&my_get_opt ( $cl,\
 \"-infile=\",1,1, \"-database=\",1,0, \"-method=\\
",1,0, \"-outfile=\",1,0,\"-pdb_type=\",1,0));\n  \
}\n\nelsif ( $mode eq \"profile_template\")\n  {\n\
\n    &seq2profile_template ($mode,&my_get_opt ( $\
cl, \"-infile=\",1,1, \"-database=\",1,0, \"-metho\
d=\",1,0, \"-outfile=\",1,0));\n  }\nelsif ( $mode\
 eq \"psiprofile_template\")\n  {\n    &seq2profil\
e_template ($mode,&my_get_opt ( $cl, \"-infile=\",\
1,1, \"-database=\",1,0, \"-method=\",1,0, \"-outf\
ile=\",1,0));\n  }\nelsif ( $mode eq \"RNA_templat\
e\")\n  {\n    &seq2RNA_template ($mode,&my_get_op\
t ( $cl, \"-infile=\",1,1,\"-pdbfile=\",1,1,\"-out\
file=\",1,0));\n  }\nelsif ( $mode eq \"tm_templat\
e\")\n  {\n    &seq2tm_template ($mode,&my_get_opt\
 ( $cl, \"-infile=\",1,1,\"-arch=\",1,1,\"-psv=\",\
1,1, \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \\
"psitm_template\")\n  {\n    &seq2tm_template ($mo\
de,&my_get_opt ( $cl, \"-infile=\",1,1, \"-arch=\"\
,1,1,\"-psv=\",1,1, \"-outfile=\",1,0,\"-database=\
\",1,0));\n  }\nelsif ( $mode eq \"ssp_template\")\
\n  {\n    &seq2ssp_template ($mode,&my_get_opt ( \
$cl, \"-infile=\",1,1,\"-seq=\",1,1,\"-obs=\",1,1,\
 \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \"psi\
ssp_template\")\n  {\n    &seq2ssp_template ($mode\
,&my_get_opt ( $cl, \"-infile=\",1,1,\"-seq=\",1,1\
,\"-obs=\",1,1, \"-outfile=\",1,0));\n  }\n\n\n\ne\
lse\n  {\n    myexit(flush_error( \"$mode is an un\
known mode of tc_generic_method.pl\"));\n  }\nmyex\
it ($EXIT_SUCCESS);\n\n\nsub seq2ssp_template\n  {\
\n  my ($mode, $infile,$gor_seq,$gor_obs,$outfile)\
=@_;\n  my %s, %h;\n  my $result;\n  my (@profiles\
);\n  &set_temporary_dir (\"set\",$infile,\"seq.pe\
p\");\n  %s=read_fasta_seq (\"seq.pep\");\n\n\n  o\
pen (R, \">result.aln\");\n\n  #print stdout \"\\n\
\";\n  foreach $seq (keys(%s))\n    {\n\n      ope\
n (F, \">seqfile\");\n      $s{$seq}{seq}=uc$s{$se\
q}{seq};\n      print (F \">$s{$seq}{name}\\n$s{$s\
eq}{seq}\\n\");\n      close (F);\n      $lib_name\
=\"$s{$seq}{name}.ssp\";\n      $lib_name=&clean_f\
ile_name ($lib_name);\n\n      if ($mode eq \"ssp_\
template\"){&seq2gor_prediction ($s{$seq}{name},$s\
{$seq}{seq}, \"seqfile\", $lib_name,$gor_seq, $gor\
_obs);}\n      elsif ($mode eq \"psissp_template\"\
)\n	{\n	  &seq2msa_gor_prediction ($s{$seq}{name},\
$s{$seq}{seq},\"seqfile\", $lib_name,$gor_seq, $go\
r_obs);\n	}\n\n      if ( !-e $lib_name)\n	{\n	  m\
yexit(flush_error(\"GORIV failed to compute the se\
condary structure of $s{$seq}{name}\"));\n	  myexi\
t ($EXIT_FAILURE);\n	}\n      else\n	{\n	  print s\
tdout \"!\\tProcess: >$s{$seq}{name} _E_ $lib_name\
 \\n\";\n	  print R \">$s{$seq}{name} _E_ $lib_nam\
e\\n\";\n	}\n      unshift (@profiles, $lib_name);\
\n    }\n  close (R);\n  &set_temporary_dir (\"uns\
et\",$mode, $method,\"result.aln\",$outfile, @prof\
iles);\n}\n\nsub seq2tm_template\n  {\n  my ($mode\
,$infile,$arch,$psv,$outfile,$db)=@_;\n  my %s, %h\
;\n  my $result;\n  my (@profiles);\n  &set_tempor\
ary_dir (\"set\",$infile,\"seq.pep\");\n  %s=read_\
fasta_seq (\"seq.pep\");\n\n\n  open (R, \">result\
.aln\");\n\n  #print stdout \"\\n\";\n  foreach $s\
eq (keys(%s))\n    {\n      open (F, \">seqfile\")\
;\n      print (F \">$s{$seq}{name}\\n$s{$seq}{seq\
}\\n\");\n      close (F);\n      $lib_name=\"$s{$\
seq}{name}.tmp\";\n      $lib_name=&clean_file_nam\
e ($lib_name);\n\n      if ($mode eq \"tm_template\
\")\n	{\n	  &safe_system (\"t_coffee -other_pg fas\
ta_seq2hmmtop_fasta.pl -in=seqfile -out=$lib_name \
-arch=$arch -psv=$psv\");\n	}\n      elsif ( $mode\
 eq \"psitm_template\")\n	{\n	  &seq2msa_tm_predic\
tion ($s{$seq}{name},$s{$seq}{seq}, $db, \"seqfile\
\", $lib_name,$arch, $psv);\n	}\n      if ( !-e $l\
ib_name)\n	{\n	  myexit(flush_error(\"hmmtop faile\
d to compute the secondary structure of $s{$seq}{n\
ame}\"));\n	  myexit ($EXIT_FAILURE);\n	}\n      e\
lse\n	{\n	  print stdout \"!\\tProcess: >$s{$seq}{\
name} _T_ $lib_name\\n\";\n	  print R \">$s{$seq}{\
name} _T_ $lib_name\\n\";\n	}\n      unshift (@pro\
files, $lib_name);\n    }\n  close (R);\n  &set_te\
mporary_dir (\"unset\",$mode, $method,\"result.aln\
\",$outfile, @profiles);\n}\n\n\n\nsub seq2RNA_tem\
plate\n  {\n    \n    my ($mode, $infile, $pdbfile\
, $outfile)=@_;\n    my %s, %h ;\n    my $result;\\
n    my (@profiles);\n    my ($seq_mode, $pdb_mode\
, $pwd);\n    \n    #use $seq_mode to estimate the\
 template of sequences WITHOUT a PDB\n    #use $pd\
b_mode to estimate the template of sequences WITH \
   a PDB\n\n    $seq_mode=$ENV{\"SEQ2TEMPLATE4_F_\\
"};\n    $pdb_mode=$ENV{\"PDB2TEMPLATE4_F_\"};\n  \
  \n    if (!$pdb_mode){$pdb_mode=\"find_pair-p\";\
}\n    if (!$seq_mode){$seq_mode=\"RNAplfold\";}\n\
    \n    my $cwd = cwd();\n    &set_temporary_dir\
 (\"set\",$infile,\"seq.pep\");\n    %s=read_fasta\
_seq (\"seq.pep\");\n    %pdb_template_h = &read_t\
emplate_file($pdbfile);\n    my $pdb_chain;\n    \\
n       \n    open (R, \">result.aln\");\n    #pri\
nt stdout \"\\n\";\n    foreach $seq (keys(%s))\n \
     {\n	\n	open (F, \">seqfile\");\n	print (F \">\
$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n	close (F);\
\n	$pdb_chain = $pdb_template_h{$seq};\n	$lib_name\
=\"$s{$seq}{name}.rfold\";\n	$lib_name=&clean_file\
_name ($lib_name);\n	if ($pdb_template_h{$seq} eq \
\"\")\n	  {\n	    if    ($seq_mode eq \"RNAplfold\\
"){RNAplfold2lib (\"seqfile\", \"$lib_name\");}\n	\
    elsif ($seq_mode eq \"no\"){$lib_name=0;}\n	  \
  else\n	      {\n		myexit(add_error (EXIT_FAILURE\
,$$,$$,getppid(), \"seq2RNA_template failure::meth\
od $seq_mode not available for sequences without P\
DB structures\"));\n	      }\n	  }\n	elsif ($pdb_t\
emplate_h{$seq} ne \"\")\n	  {\n	    my $pdbf;\n	 \
   if    ( -e \"$cwd/$pdb_chain\"   ){$pdbf=\"$cwd\
/$pdb_chain\"; }\n	    elsif ( -e  $pdb_chain     \
    ){$pdbf=\"$pdb_chain\";      }\n	    elsif ( -\
e  \"$CACHE$pdb_chain\" ){$pdbf=\"$CACHE$pdb_chain\
\";}\n	    elsif ( -e  \"$CACHE/$pdb_chain\"){$pdb\
f=\"$CACHE/$pdb_chain\";}\n	    else\n	      {\n		\
myexit(flush_error(\"Could not read $pdb_chain \")\
);\n	      }\n\n	    if($pdb_mode eq \"x3dna-ssr\"\
)\n	      {\n		x3dnassr2lib (\"seqfile\", \"$pdbf\\
", \"$lib_name\");\n	      }\n	    elsif ($pdb_mod\
e eq \"find_pair-p\")\n	      {\n		x3dna_find_pair\
2lib (\"seqfile\", \"$pdbf\", \"$lib_name\", \"fin\
d_pair -p\");\n	      }\n	    elsif ($pdb_mode eq \
\"find_pair\")\n	      {\n		x3dna_find_pair2lib (\\
"seqfile\", \"$pdbf\", \"$lib_name\", \"find_pair\\
");\n	      }\n	    elsif ($pdb_mode eq \"RNAplfol\
d\")\n	      {\n		RNAplfold2lib (\"seqfile\", \"$l\
ib_name\");\n	      }\n	    elsif ($pdb_mode eq \"\
no\"){$lib_name=0;}\n	    else\n	      {\n		myexit\
(add_error (EXIT_FAILURE,$$,$$,getppid(), \"seq2RN\
A_template failure::Could not find method $pdb_mod\
e\"));\n	      }\n	  }\n	if ($lib_name)\n	  {\n	  \
  print stdout \"!\\tProcess: >$s{$seq}{name} _F_ \
$lib_name\\n\";\n	    print R \">$s{$seq}{name} _F\
_ $lib_name\\n\";\n	    unshift (@profiles, $lib_n\
ame);\n	  }\n      }\n    close (R);\n    &set_tem\
porary_dir (\"unset\",$mode, $method,\"result.aln\\
",$outfile, @profiles);\n  }\n\n\n\nsub psiblast2p\
rofile_template_test\n  {\n  my ($mode, $seq,$blas\
t)=@_;\n  my %s, %h, ;\n  my ($result,$psiblast_ou\
tput,$profile_name,@profiles);\n  my $trim=0;\n  m\
y $maxid=100;\n  my $minid=0;\n  my $mincov=0;\n  \
my $maxcov=100;\n\n  %s=read_fasta_seq ($seq);\n  \
open (R, \">result.aln\");\n\n  #print stdout \"\\\
n\";\n  foreach $seq (keys(%s))\n    {\n\n      op\
en (F, \">seqfile\");\n      print (F \">$A\\n$s{$\
seq}{seq}\\n\");\n      close (F);\n      $psiblas\
t_output=$blast;\n      if ( -e $psiblast_output)\\
n	{\n	  %profile=blast_xml2profile($s{$seq}{name},\
 $s{$seq}{seq},$maxid, $minid,$mincov,$psiblast_ou\
tput);\n\n\n\n	  $profile_name=\"$s{$seq}{name}.pr\
f\";\n	  $profile_name=&clean_file_name ($profile_\
name);\n	  unshift (@profiles, $profile_name);\n	 \
 output_profile ($profile_name, \\%profile, $trim)\
;\n	  print stdout \"!\\tProcess: >$s{$seq}{name} \
_R_ $profile_name [$profile{n} Seq.] [$SERVER/blas\
t/$db][$CACHE_STATUS]\\n\";\n	  print R \">$s{$seq\
}{name} _R_ $profile_name\\n\";\n	}\n    }\n  clos\
e (R);\n\n  die;\n}\nsub seq2profile_template\n   \
 {\n      my ($mode, $infile, $db, $method, $outfi\
le)=@_;\n      if    ($method eq \"psiblast\"){ret\
urn psiblast2profile_template ($mode, $infile, $db\
, $method, $outfile);}\n      elsif ($method eq \"\
blastp\")   {return psiblast2profile_template ($mo\
de, $infile, $db, $method, $outfile);}\n      elsi\
f ($method eq \"hh\")      {return hh2profile_temp\
late ($mode, $infile, $db, $method, $outfile);}\n \
   }\n\nsub psiblast2profile_template\n  {\n  my (\
$mode, $infile, $db, $method, $outfile)=@_;\n  my \
%s, %h, ;\n  my ($result,$psiblast_output,$profile\
_name,@profiles);\n  &set_temporary_dir (\"set\",$\
infile,\"seq.pep\");\n  %s=read_fasta_seq (\"seq.p\
ep\");\n  open (R, \">result.aln\");\n\n  #print s\
tdout \"\\n\";\n  foreach $seq (keys(%s))\n    {\n\
      open (F, \">seqfile\");\n      print (F \">$\
A\\n$s{$seq}{seq}\\n\");\n      close (F);\n      \
$psiblast_output=&run_blast ($s{$seq}{name},$metho\
d, $db, \"seqfile\",\"outfile\");\n\n      if ( -e\
 $psiblast_output)\n	{\n	  my %profile=blast_xml2p\
rofile($s{$seq}{name}, $s{$seq}{seq},$maxid, $mini\
d,$mincov,$psiblast_output);\n	  unlink ($psiblast\
_output);\n	  \n	  $profile_name=\"$s{$seq}{name}.\
prf\";\n	  $profile_name=&clean_file_name ($profil\
e_name);\n	  unshift (@profiles, $profile_name);\n\
	  output_profile ($profile_name, \\%profile, $tri\
m);\n	  \n	  print stdout \"!\\tProcess: >$s{$seq}\
{name} _R_ $profile_name [$profile{n} Seq.] [$SERV\
ER/blast/$db][$CACHE_STATUS]\\n\";\n	  print R \">\
$s{$seq}{name} _R_ $profile_name\\n\";\n	  \n	  \n\
	}\n      \n    }\n  close (R);\n  \n  \n\n  &set_\
temporary_dir (\"unset\",$mode, $method,\"result.a\
ln\",$outfile, @profiles);\n}\n\nsub hh2profile_te\
mplate\n  {\n\n  #for each sequence, build a profi\
le, in FASTA, with ungapped querry on top  \n  my \
($mode, $infile, $db, $method, $outfile)=@_;\n  my\
 %s, %h, ;\n  my ($result,$psiblast_output,$profil\
e_name,@profiles);\n  &set_temporary_dir (\"set\",\
$infile,\"seq.pep\");\n  %s=read_fasta_seq (\"seq.\
pep\");\n  open (R, \">result.aln\");\n  \n  my $h\
h=$ENV{\"HHSEARCH_4_TCOFFEE\"};\n  if (!$hh)\n    \
{\n      print \"ERROR: HHSEARCH_4_TCOFFEE is not \
set\\n\";\n      myexit ($EXIT_FAILURE);\n    }\n \
 \n  #print stdout \"\\n\";\n  foreach $seq (keys(\
%s))\n    {\n      my ($profile_name, $nseq);\n   \
   open (F, \">seqfile\");\n      print (F \">$A\\\
n$s{$seq}{seq}\\n\");\n      close (F);\n      \n \
     #This function should input a querry and a da\
tabase and return as output a fasta MSA with quesr\
y on top\n      $profile_name=\"$s{$seq}{name}.prf\
\";\n      $profile_name=&clean_file_name ($profil\
e_name);\n      unshift (@profiles, $profile_name)\
;\n      \n      \n      safe_system  (\"$hh -name\
=$s{$seq}{name} -method=search -db=$db -seq=seqfil\
e -outfile=$profile_name\");\n      if (-e $profil\
e_name){$nseq=fasta2nseq($profile_name);}\n      \\
n      print stdout \"!\\tProcess: >$s{$seq}{name}\
 _R_ $profile_name [$nseq Seq.] [$method/$db][$CAC\
HE_STATUS]\\n\";\n      print R \">$s{$seq}{name} \
_R_ $profile_name\\n\";\n    }\n  close (R);\n  &s\
et_temporary_dir (\"unset\",$mode, $method,\"resul\
t.aln\",$outfile, @profiles);\n}\n\nsub blast2pdb_\
template_test\n    {\n      my ($mode,$infile)=@_;\
\n      my ($maxid,$minid,$mincov);\n      $maxid=\
100;\n      $minid=0;\n      $mincov=0;\n\n      p\
rint \"$infile\\n\";\n\n      %p=blast_xml2profile\
($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$min\
cov,$infile);\n      $c=1;\n      print stdout \"!\
\\tProcess: >$s{$seq}{name} [$SERVER/blast/$db][$C\
ACHE_STATUS]\\n\";\n      while (!$found && $c<$p{\
n})\n	{\n	  $pdbid=&id2pdbid($p{$c}{identifyer});\\
n	  if ( length ($pdbid)>5){$pdbid=id2pdbid($p{$c}\
{definition});}\n\n	  if ( length ($pdbid)>5)\n	  \
  {\n	      myexit(add_error (EXIT_FAILURE,$$,$$,g\
etppid(), \"BLAST_FAILURE::Could Not Parse PDBID (\
$p{$c}{identifyer},$p{$c}{definition})\"));\n	    \
}\n\n\n	  if (!&pdb_is_released($pdbid))\n	    {\n\
	      print stdout \"\\t\\t**$pdbid [WARNIG: PDB \
NOT RELEASED or WITHDRAWN]\\n\";\n	      $c++;\n	 \
   }\n	  elsif (!&pdb_has_right_type ($pdbid,$type\
))\n	    {\n	      my $ptype=&pdb2type ($pdbid);\n\
	      my $etype=&type2etype($type);\n\n	      pri\
nt stdout \"\\t\\t**$pdbid [$ptype cannot be used \
(expected: $etype)]\\n\";\n	      $c++;\n	    }\n	\
  else\n	    {\n	      $found=1;\n	    }\n	}\n\n  \
    if ($found)\n	{\n	  print stdout \"\\t\\t >$s{\
$seq}{name} _P_ $pdbid\\n\";\n	}\n      else\n	{\n\
	  print stdout \"\\t\\t >$s{$seq}{name} No Templa\
te Selected\\n\";\n	}\n      die;\n    }\nsub blas\
t2pdb_template\n  {\n  my ($mode, $infile, $db, $m\
ethod, $outfile,$type)=@_;\n  my %s, %h, ;\n  my (\
$result,$blast_output);\n  &set_temporary_dir (\"s\
et\",$infile,\"seq.pep\");\n  %s=read_fasta_seq (\\
"seq.pep\");\n  open (R, \">result.aln\");\n\n\n  \
#print stdout \"\\n\";\n  foreach $seq (keys(%s))\\
n    {\n      my $c;\n      my $found;\n\n      op\
en (F, \">seqfile\");\n      print (F \">$A\\n$s{$\
seq}{seq}\\n\");\n      close (F);\n\n      $blast\
_output=&run_blast ($s{$seq}{name},$method, $db, \\
"seqfile\",\"outfile\");\n\n      %p=blast_xml2pro\
file($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,\
$mincov,$blast_output);\n      unlink ($blast_outp\
ut);\n\n      $c=1;\n      print stdout \"!\\tProc\
ess: >$s{$seq}{name} [$SERVER/blast/$db][$CACHE_ST\
ATUS]\\n\";\n      while (!$found && $c<$p{n})\n	{\
\n	  $pdbid=&id2pdbid($p{$c}{identifyer});\n	  if \
( length ($pdbid)>5){$pdbid=id2pdbid($p{$c}{defini\
tion});}\n\n	  if ( length ($pdbid)>5)\n	    {\n	 \
     myexit(add_error (EXIT_FAILURE,$$,$$,getppid(\
), \"BLAST_FAILURE::Could Not Parse PDBID ($p{$c}{\
identifyer},$p{$c}{definition})\"));\n	    }\n\n\n\
	  if (!&pdb_is_released($pdbid))\n	    {\n	      \
print stdout \"\\t\\t**$pdbid [PDB NOT RELEASED or\
 WITHDRAWN]\\n\";\n	      $c++;\n	    }\n	  elsif \
(!&pdb_has_right_type ($pdbid,$type))\n	    {\n	  \
    my $ptype=&pdb2type ($pdbid);\n	      my $etyp\
e=&type2etype($type);\n\n	      print stdout \"\\t\
\\t**$pdbid [$ptype cannot be used (expected: $ety\
pe)]\\n\";\n	      $c++;\n	    }\n	  else\n	    {\\
n	      $found=1;\n	    }\n	}\n\n      if ($found)\
\n	{\n	  print R \">$s{$seq}{name} _P_ $pdbid\\n\"\
;\n	  print stdout \"\\t\\t >$s{$seq}{name} _P_ $p\
dbid\\n\";\n	}\n      else\n	{\n	  print R \">$s{$\
seq}{name}\\n\";\n	  print stdout \"\\t\\t >$s{$se\
q}{name} No Template Selected\\n\";\n	}\n    }\n  \
close (R);\n  &set_temporary_dir (\"unset\",$mode,\
 $method,\"result.aln\",$outfile);\n}\nsub type2et\
ype\n  {\n    my $type=shift;\n    my $etype;\n\n \
   if ( $type=~/n/){$etype.=\"NMR \";}\n    if ( $\
type=~/d/){$etype.=\"diffraction \";}\n	if ( $type\
=~/e/){$etype.=\"EM \";}\n    if ( $type=~/m/){$et\
ype.=\"model \";}\n    return $etype;\n  }\nsub pd\
b2type\n  {\n     my $pdb=shift;\n     my $f=vtmpn\
am();\n\n     my $value= &safe_system (\"t_coffee \
-other_pg extract_from_pdb -model_type $pdb > $f\"\
);\n     my $r=&file2string ($f);\n     chomp($r);\
\n     return $r;\n   }\nsub pdb_has_right_type\n \
 {\n    my $pdb=shift;\n    my $type=shift;\n\n   \
 my $f=vtmpnam();\n\n    my $value= &safe_system (\
\"t_coffee -other_pg extract_from_pdb -model_type \
$pdb > $f\");\n    my $r=&file2string ($f);\n    c\
homp($r);\n\n\n    if ( $r eq \"NMR\" && $type=~/n\
/){return 1;}\n    elsif ( $r eq \"diffraction\" &\
& $type=~/d/){return 1;}\n	elsif ( $r eq \"EM\" &&\
 $type=~/e/){return 1;}\n    elsif ( $r eq \"model\
\" && $type=~/m/){return 1;}\n    else {return 0;}\
\n  }\nsub pdb_is_released\n  {\n    my $pdb=shift\
;\n    my $f=vtmpnam();\n\n    $value= &safe_syste\
m (\"t_coffee -other_pg extract_from_pdb -is_relea\
sed_pdb_name $pdb > $f\");\n    my $r=&file2string\
 ($f);\n    chomp($r);\n    return $r;\n  }\nsub b\
last_msa\n  {\n    my ($blast,$infile,$db,$outfile\
)=@_;\n    my ($a, %s1, %s, %qs, %qs1);\n    my $s\
eqfile;\n    my $SEQ=new FileHandle;\n    my $seqf\
ile=\"seqfile\";\n    my @txt;\n\n\n    %s1=&read_\
fasta_seq ($db);\n    %s=&fasta_hash2index_hash(%s\
1);\n    %qs1=&read_fasta_seq ($infile);\n    %qs=\
&fasta_hash2index_hash(%qs1);\n\n\n    #&safe_syst\
em (\"formatdb -i $db\");\n    if ($blast eq \"bla\
stp\"){&safe_system  (\"blastall -i $infile -d $db\
 -m7 -p blastp -o io\");}\n    elsif ($blast eq \"\
blastn\"){&safe_system  (\"blastn -query $infile -\
db $db -outfmt 5 -word_size 4 -out io\");}\n\n    \
&set_blast_type (\"io\");\n\n\n    my %FB=&xml2tag\
_list (\"io\", \"Iteration\");\n    open (F, \">$o\
utfile\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\
\n    print F \"$s{n}\\n\";\n    for ( my $a=0; $a\
<$s{n}; $a++)\n      {\n	print F \"$s{$a}{name} $s\
{$a}{len} $s{$a}{seq}\\n\";\n      }\n\n\n    for \
( my $a=0; $a<$FB{n}; $a++)\n      {\n	my %p=blast\
_xml2profile ($qs{$a}{name}, $qs{$a}{seq},100, 0, \
0, $FB{$a}{body});\n	my $query=$p{0}{name};\n	my $\
i= $s1{$query}{order}+1;\n	for (my $b=1; $b<$p{n};\
 $b++)\n	  {\n	    my $l=length ($p{$b}{Qseq});\n	\
    my $hit=$p{$b}{definition};\n	    my $Qstart=$\
p{$b}{Qstart};\n	    my $Hstart=$p{$b}{Hstart};\n	\
    my $identity=$p{$b}{identity};\n	    my @lrQ=s\
plit (//,$p{$b}{Qseq});\n	    my @lrH=split (//,$p\
{$b}{Hseq});\n\n	    my $j= $s1{$hit}{order}+1;\n	\
    #if ( $j==$i){next;}\n	    printf F \"# %d %d\\
\n\", $i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} \
($Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($\
c=0; $c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n	\
	my $rH=$lrH[$c];\n		my $n=0;\n\n		if ($rQ ne \"-\\
"){$n++, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $H\
start++;}\n\n		if ( $n==2)\n		  {\n		    printf F \
\"\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,$identity\
;\n		  }\n	      }\n	  }\n      }\n    print F \"!\
 SEQ_1_TO_N\\n\";\n    close (F);\n    return $out\
put;\n  }\n\nsub blast_msa_old\n  {\n    my ($infi\
le,$outfile)=@_;\n    my ($a, %seq);\n    %s1=&rea\
d_fasta_seq ($infile);\n    foreach $s (keys (%s1)\
)\n      {\n	$i=$s1{$s}{order};\n	$s{$i}{name}=$s;\
\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=length(\
 $s{$i}{seq});\n	$s{n}++;\n      }\n    &safe_syst\
em (\"formatdb -i $infile\");\n    &safe_system (\\
"blastall -i $infile -d $infile -m7 -o io\");\n   \
 &set_blast_type (\"io\");\n\n    %FB=&xml2tag_lis\
t (\"io\", \"Iteration\");\n\n    open (F, \">$out\
file\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n\
    print F \"$s{n}\\n\";\n    for ( $a=0; $a<$s{n\
}; $a++)\n      {\n	print F \"$s{$a}{name} $s{$a}{\
len} $s{$a}{seq}\\n\";\n      }\n    for ( $a=0; $\
a<$FB{n}; $a++)\n      {\n	%p=blast_xml2profile ($\
s{$a}{name}, $s{$a}{seq},100, 0, 0, $FB{$a}{body})\
;\n	for ($b=1; $b<$p{n}; $b++)\n	  {\n	    my $l=l\
ength ($p{$b}{Qseq});\n	    my $hit=$p{$b}{definit\
ion};\n	    my $Qstart=$p{$b}{Qstart};\n	    my $H\
start=$p{$b}{Hstart};\n	    my $identity=$p{$b}{id\
entity};\n	    my @lrQ=split (//,$p{$b}{Qseq});\n	\
    my @lrH=split (//,$p{$b}{Hseq});\n	    my $i= \
$s1{$s{$a}{name}}{order}+1;\n	    my $j= $s1{$hit}\
{order}+1;\n	    #if ( $j==$i){next;}\n	    printf\
 F \"# %d %d\\n\", $i, $j;\n	    #  print  F \"\\n\
$p{$b}{Qseq} ($Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\
\n	    for ($c=0; $c<$l; $c++)\n	      {\n		my $rQ\
=$lrQ[$c];\n		my $rH=$lrH[$c];\n		my $n=0;\n\n		if\
 ($rQ ne \"-\"){$n++, $Qstart++;}\n		if ($rH ne \"\
-\"){$n++; $Hstart++;}\n\n		if ( $n==2)\n		  {\n		\
    printf F \"\\t%d %d %d\\n\", $Qstart-1, $Hstar\
t-1,$identity;\n		  }\n	      }\n	  }\n      }\n  \
  print F \"! SEQ_1_TO_N\\n\";\n    close (F);\n  \
  return $output;\n\n  }\n\nsub seq2msa\n  {\n    \
my ($mode, $infile, $method, $param, $outfile,$dat\
abase)=@_;\n    &set_temporary_dir (\"set\",$infil\
e,\"seq.pep\", $database, \"db.pep\");\n    $param\
.=\" >/dev/null 2>&1 \";\n\n\n    #make sure test.\
pep is in FASTA\n    &safe_system (\"t_coffee -oth\
er_pg seq_reformat -in seq.pep -output fasta_seq >\
 x\");\n    `mv x seq.pep`;\n\n    if ( $method eq\
 \"blastp\")\n      {\n	&blast_msa (\"blastp\",\"s\
eq.pep\",$database,\"result.aln\");\n      }\n    \
elsif ( $method eq \"blastn\")\n      {\n	&blast_m\
sa (\"blastn\",\"seq.pep\",$database,\"result.aln\\
");\n      }\n\n    elsif ( $method eq \"muscle\")\
\n      {\n	`muscle -in seq.pep -out result.aln $p\
aram`;\n      }\n    elsif ( $method eq \"probcons\
\")\n      {\n	`probcons seq.pep >result.aln 2>/de\
v/null`;\n      }\n    elsif ( $method eq \"mafft\\
")\n      {\n	`mafft --quiet --localpair --maxiter\
ate 1000 seq.pep> result.aln  2>/dev/null`\n      \
}\n    elsif ( $method=~/prank/)\n      {\n	`$meth\
od -d=seq.pep -o=result.aln -quiet 2>/dev/null`;\n\
	`mv result.aln.1.fas result.aln`;\n      }\n    e\
lsif ($method eq \"clustalo\")\n      {\n	`clustal\
o -i seq.pep > result.aln`;\n      }\n\n    else\n\
      {\n	`$method -infile=seq.pep -outfile=result\
.aln`;\n      }\n\n    &set_temporary_dir (\"unset\
\",$mode, $method,\"result.aln\",$outfile);\n    m\
yexit ($EXIT_SUCCESS);\n  }\n\nsub seq2thread_pair\
\n  {\n    my ($mode, $infile, $pdbfile1, $method,\
 $param, $outfile)=@_;\n    &set_temporary_dir (\"\
set\",$infile,\"seq.pep\",$pdbfile1,\"struc.pdb\")\
;\n    if ($method eq \"fugueali\")\n      {\n	#En\
v Variable that need to be defined for Fugue\n	if \
(!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIST}=\"DUM\
MY\";}\n	if (!$ENV{HOMSTRAD_PATH})  {$ENV{HOMSTRAD\
_PATH}=\"DUMMY\";}\n	if (!$ENV{HOMS_PATH}){$ENV{HO\
MS_PATH}=\"DUMMY\";}\n\n	`joy struc.pdb >x 2>x`;\n\
	&check_file(\"struc.tem\", \"Joy failed [FATAL:$P\
ROGRAM/$method]\");\n	`melody -t struc.tem >x 2>x`\
;\n	&check_file(\"struc.tem\", \"Melody failed [FA\
TAL:$PROGRAM/$method]\");\n	`fugueali -seq seq.pep\
 -prf struc.fug -print > tmp_result.aln`;\n\n	&che\
ck_file(\"tmp_result.aln\", \"Fugue failed [FATAL:\
$PROGRAM/$method]\");\n	&safe_system (\"t_coffee -\
other_pg seq_reformat -in tmp_result.aln -output f\
asta_aln >result.aln\");\n      }\n    elsif ( $me\
thod eq \"t_coffee\")\n      {\n	&safe_system (\"t\
_coffee -in Pstruc.pdb Sseq.pep Mslow_pair -outfil\
e result.aln -quiet\");\n      }\n    else\n      \
{\n	&safe_system (\"$method -infile=seq.pep -pdbfi\
le1=struc.pdb -outfile=result.aln $param>x 2>x\");\
\n      }\n    &set_temporary_dir (\"unset\",$mode\
,$method,\"result.aln\",$outfile);\n    myexit ($E\
XIT_SUCCESS);\n  }\nsub seq2pdbid_pair\n  {\n    m\
y ($mode, $pdbfile1, $pdbfile2, $method, $param, $\
outfile)=@_;\n    my ($name);\n\n\n    &set_tempor\
ary_dir (\"set\");\n    $name=$pdbfile1.\" \".$pdb\
file2;\n\n    if (    &cache_file(\"GET\",\"\",\"$\
name\",\"$method\",\"dali\",$outfile,\"EBI\"))\n  \
    {return $outfile;}\n    else\n      {\n	if ($m\
ethod eq \"daliweb\")\n	  {\n	    $pdbfile1=~/(...\
.)(.)/;\n	    $id1=$1; $c1=$2;\n\n	    $pdbfile2=~\
/(....)(.)/;\n	    $id2=$1; $c2=$2;\n\n	    $comma\
nd=\"t_coffee -other_pg dalilite.pl --pdb1 $id1 --\
chainid1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$E\
MAIL  >dali_stderr 2>dali_stderr\";\n	    $dali=`$\
command`;\n\n	    open (F, \"dali_stderr\");\n	   \
 while (<F>)\n	      {\n		if ( /JobId: dalilite-(\\
\S+)/)\n		{\n		  $jobid=$1;\n		}\n	      }\n	    c\
lose (F);\n	    unlink (\"dali_stderr\");\n\n	    \
$output1=\"dalilite-$jobid.txt\";\n	    if ( -e $o\
utput1)\n	      {\n		unlink ($output1);\n		&url2fi\
le (\"http://www.ebi.ac.uk/Tools/es/cgi-bin/jobres\
ults.cgi/dalilite/dalilite-$jobid/aln.html\", \"ou\
tput2\");\n\n		if ( -e \"output2\")\n		  {\n		    \
my ($seq1, $seq2);\n		    $seq1=$seq2=\"\";\n\n		 \
   open (F, \"output2\");\n		    while (<F>)\n		  \
    {\n			$l=$_;\n			if ( $l=~/Query\\s+(\\S+)/)\n\
			  {\n			    $seq1.=$1;\n			  }\n			elsif ( $l=~\
/Sbjct\\s+(\\S+)/)\n			  {\n			    $seq2.=$1;\n			\
  }\n		      }\n		    close (F);\n		    unlink (\"\
output2\");\n		    if ($seq1 ne \"\" && $seq2 ne \\
"\")\n		      {\n			$output3=\">$A\\n$seq1\\n>$B\\\
n$seq2\\n\";\n			$output3=~s/\\./-/g;\n			open (F,\
 \">result.aln\");\n			print F \"$output3\";\n			c\
lose (F);\n		      }\n		  }\n	      }\n	  }\n     \
 }\n    &cache_file(\"SET\",\"\",\"$name\",\"$meth\
od\",\"dali\",\"result.aln\",\"EBI\");\n    &set_t\
emporary_dir (\"unset\",$mode, $method, \"result.a\
ln\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\\
nsub seq2pdb_pair\n  {\n    my ($mode, $pdbfile1, \
$pdbfile2, $method, $param, $outfile)=@_;\n\n    &\
set_temporary_dir (\"set\",$pdbfile1,\"pdb1.pdb\",\
$pdbfile2,\"pdb2.pdb\");\n    if ($method eq \"t_c\
offee\")\n      {\n	&safe_system (\"t_coffee -in P\
pdb1.pdb Ppdb2.pdb -quiet -outfile=result.aln\");\\
n      }\n    elsif ( $method eq \"DaliLite\")\n  \
    {\n	if ( &safe_system (\"DaliLite -pairwise pd\
b1.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\n	 \
    my ($seq1, $seq2);\n	     $seq1=$seq2=\"\";\n\\
n	     open (F, \"tmp1\");\n	     while (<F>)\n	  \
     {\n		 $l=$_;\n		 if ( $l=~/Query\\s+(\\S+)/)\\
n		   {\n		     $seq1.=$1;\n		   }\n		 elsif ( $l=\
~/Sbjct\\s+(\\S+)/)\n		   {\n		     $seq2.=$1;\n		\
   }\n	       }\n	     close (F);\n	     unlink (\\
"tmp1\");\n	     if ($seq1 ne \"\" && $seq2 ne \"\\
")\n	       {\n		 my $output3=\">$A\\n$seq1\\n>$B\\
\n$seq2\\n\";\n		 $output3=~s/\\./-/g;\n		 open (F\
, \">result.aln\");\n		 print F \"$output3\";\n		 \
close (F);\n	       }\n	   }\n	else\n	  {\n	    pr\
int \"ERROR: DalLite failed to align the considere\
d structures[tc_generic_method.pl]\\n\";\n	  }\n  \
    }\n    elsif ( $method eq \"TMalign\")\n      \
{\n	if ( &safe_system (\"TMalign pdb1.pdb pdb2.pdb\
 >tmp1\")==$EXIT_SUCCESS)\n	  {\n	    `tail -4 tmp\
1 > tmp2`;\n\n	    open (F, \"tmp2\");\n	    while\
 (<F>)\n	      {\n		unshift(@l, $_);\n	      }\n	 \
   close (F);\n	    open (F, \">result.aln\");\n	 \
   $l[3]=~s/[^a-zA-Z0-9-]/\\-/g;\n	    $l[1]=~s/[^\
a-zA-Z0-9-]/\\-/g;\n	    print F \">$A\\n$l[3]\\n>\
$B\\n$l[1]\\n\";\n	    close (F);\n	  }\n	else\n	 \
 {\n	    print \"ERROR: TMalign failed to align th\
e considered structures[tc_generic_method.pl]\\n\"\
;\n	    `rm result.aln >/dev/null 2>/dev/null`;\n	\
  }\n      }\n    elsif ( $method eq \"mustang\")\\
n      {\n	if ( &safe_system (\"mustang -i pdb1.pd\
b pdb2.pdb -F fasta >/dev/null 2>/dev/null\")==$EX\
IT_SUCCESS)\n	  {\n	    `mv results.afasta result.\
aln`;\n	  }\n	else\n	  {\n	    print \"ERROR: must\
ang failed to align the considered structures[tc_g\
eneric_method.pl]\\n\";\n	    `rm result.aln >/dev\
/null 2>/dev/null`;\n	  }\n      }\n    else\n    \
  {\n	if ( &safe_system (\"$method -pdbfile1=pdb1.\
pep -pdbfile2=pdb2.pdb -outfile=result.aln $param>\
x 2>x\")==$EXIT_SUCCESS)\n	  {\n	    `mv results.a\
fasta result.aln`;\n	  }\n	else\n	  {\n	    print \
\"ERROR: $method failed to align the considered st\
ructures[tc_generic_method.pl]\\n\";\n	    `rm res\
ult.aln >/dev/null 2>/dev/null`;\n	  }\n      }\n \
   &set_temporary_dir (\"unset\",$mode, $method, \\
"result.aln\",$outfile);\n    myexit ($EXIT_SUCCES\
S);\n  }\n\nsub seq2rnapdb_pair\n  {\n    my ($mod\
e, $pdbfile1, $pdbfile2, $method, $param, $outfile\
)=@_;\n    \n    if ($method eq \"runsara.py\")\n \
     {\n	my $path=$ENV{PATH};\n	\n	if ($ENV{X3DNA_\
4_SARA}){$ENV{PATH}=\"$ENV{X3DNA_4_SARA}:$path\";}\
\n	\n	open(TMP,\"<$pdbfile1\");\n	my $count = 0;\n\
	my $line;\n	while (<TMP>)\n	  {\n	    $line = $_;\
\n	    if ($count ==1)\n	      {\n		last;\n	      \
}\n	    $count += 1;\n	  }\n	\n	\n	$chain1 = subst\
r($line,length($line)-3,1);\n	\n	close TMP;\n	open\
(TMP,\"<$pdbfile2\");\n	my $count = 0;\n	while (<T\
MP>)\n	  {\n	    $line = $_;\n	    if ($count ==1)\
\n	      {\n		last;\n	      }\n	    $count += 1;\n\
	  }\n	$chain2 = substr($line,length($line)-3,1);\\
n	close TMP;\n	\n	$tmp_file=&vtmpnam();\n	\n	safe_\
system(\"runsara.py $pdbfile1 $chain1 $pdbfile2 $c\
hain2 -s -o $tmp_file --limitation 5000 > /dev/nul\
l 2> /dev/null\");\n	if ($ENV{X3DNA_4_SARA}){$ENV{\
PATH}=$path;}\n	\n	open(TMP,\"<$tmp_file\") or die\
 \"cannot open the sara tmp file:$!\\n\";\n	open(O\
UT,\">$outfile\") or die \"cannot open the $outfil\
e file:$!\\n\";\n	\n	my $switch = 0;\n	my $seqNum \
= 0;\n	foreach my $line (<TMP>)\n	  {\n	    next u\
nless ($line=~/SARAALI/);\n	    if ($line=~/>/)\n	\
      {\n		$switch =0;\n		print OUT \">seq$seqNum\\
\n\";\n		$seqNum++;\n	      }\n	    if ($switch < \
2){\n	      $switch++;\n	      next;\n	    }\n	   \
 \n	    if ($line =~/REMARK\\s+SARAALI\\s+([^\\*]+\
)\\*/)\n	      {\n		my $string = $1;\n		print OUT \
\"$string\\n\";\n	      }\n	  }\n	close TMP;\n	clo\
se OUT;\n	unlink($tmp_file);\n      }\n  }\nsub se\
q2profile_pair\n  {\n    my ($mode, $profile1, $pr\
ofile2, $method, $param, $outfile)=@_;\n    \n    \
\n    if ($method eq \"clustalw\")\n      {\n	`clu\
stalw -profile1=$profile1 -profile2=$profile2 -out\
file=$outfile`;\n      }\n    elsif ( $method eq \\
"clustalo\")\n      {\n	\n	`clustalo --p1 $profile\
1 --p2 $profile2 -o $outfile --force`;\n      }\n \
   elsif ( $method eq \"hhalign\")\n      {\n	hhal\
ign ( $profile1,$profile2,$outfile,$param);\n     \
 }\n    else\n      {\n	`$method -profile1=$profil\
e1 -profile2=$profile2 -outfile=$outfile $param> /\
dev/null 2>/dev/null`;\n      }\n    myexit ($EXIT\
_SUCCESS);\n  }\n\nsub pg_is_installed\n  {\n    m\
y @ml=@_;\n    my ($r, $p, $m);\n    my $supported\
=0;\n\n    my $p=shift (@ml);\n    if ($p=~/::/)\n\
      {\n	if (safe_system (\"perl -M$p -e 1\")==$E\
XIT_SUCCESS){return 1;}\n	else {return 0;}\n      \
}\n    else\n      {\n	my $cwhich=vtmpnam();\n	$r=\
`which $p >$cwhich 2>/dev/null`;\n	$r=file2string \
($cwhich);\n	if ($r=~/^\\//){$r=1;}\n	else {$r=0;}\
\n	\n	if ($r==0 && is_blast_package ($p)){return p\
g_is_installed (\"legacy_blast.pl\");}\n	else {ret\
urn $r;}\n      }\n  }\n\nsub is_blast_package\n  \
{\n    my $p=shift;\n    if ( $p=~/blastp/){return\
 1;}\n    elsif ($p=~/blastall/){return 1;}\n    e\
lsif ($p=~/blastn/){return 1;}\n    elsif ($p=~/bl\
astx/){return 1;}\n    elsif ($p=~/formatdb/){retu\
rn 1;}\n    else {return 0;}\n  }\n\nsub check_int\
ernet_connection\n  {\n    my $internet;\n    my $\
tmp;\n    &check_configuration ( \"wget\");\n\n   \
 $tmp=&vtmpnam ();\n\n    if     (&pg_is_installed\
    (\"wget\")){`wget www.google.com -O$tmp >/dev/\
null 2>/dev/null`;}\n    elsif  (&pg_is_installed \
   (\"curl\")){`curl www.google.com -o$tmp >/dev/n\
ull 2>/dev/null`;}\n\n    if ( !-e $tmp || -s $tmp\
 < 10){$internet=0;}\n    else {$internet=1;}\n   \
 if (-e $tmp){unlink $tmp;}\n\n    return $interne\
t;\n  }\nsub check_pg_is_installed\n  {\n    my @m\
l=@_;\n    my $r=&pg_is_installed (@ml);\n    if (\
!$r && $p=~/::/)\n      {\n	print STDERR \"\\nYou \
Must Install the perl package $p on your system.\\\
nRUN:\\n\\tsudo perl -MCPAN -e 'install $pg'\\n\";\
\n      }\n    elsif (!$r)\n      {\n	myexit(flush\
_error(\"\\nProgram $p Supported but Not Installed\
 on your system\"));\n      }\n    else\n      {\n\
	return 1;\n      }\n  }\nsub set_temporary_dir\n \
 {\n    my @list=@_;\n    my $dir_mode, $a, $mode,\
 $method;\n\n    $dir_mode=shift (@list);\n\n\n   \
 if ( $dir_mode eq \"set\")\n      {\n	$initial_di\
r=cwd();\n	if ( !$tmp_dir)\n	  {\n	    $rand=rand \
(100000);\n	    $tmp_dir=\"$TMPDIR/tmp4tcoffee_pro\
file_pair_dir_$$\\_P_$rand\";\n	  }\n	if ( !-d $tm\
p_dir)\n	  {\n	    push (@TMPDIR_LIST, $tmp_dir);\\
n	    `mkdir $tmp_dir`;\n	  }\n\n	for ( $a=0; $a<=\
$#list; $a+=2)\n	      {\n		if (-e $list[$a]){ `cp\
 $list[$a] $tmp_dir/$list[$a+1]`;}\n	      }\n	chd\
ir $tmp_dir;\n      }\n    elsif ( $dir_mode eq \"\
unset\")\n      {\n	$mode=shift (@list);\n	$method\
=shift (@list);\n\n	if (!-e $list[0])\n	  {\n	   m\
yexit(flush_error(\"Program $method failed to prod\
uce $list[1]\" ));\n	    myexit ($EXIT_FAILURE);\n\
	  }\n	else\n	  {\n	    chdir $initial_dir;\n	    \
# `t_coffee -other_pg seq_reformat -in $tmp_dir/$l\
ist[0] -output fasta_aln -out $tmp_dir/result2.aln\
`;\n	    `cp $tmp_dir/$list[0] $tmp_dir/result2.al\
n`;\n	    if ( $list[1] eq \"stdout\")\n	      {\n\
		open (F, \"$tmp_dir/result2.aln\");\n		while (<F\
>){print $_;}close(F);\n	      }\n	    else\n	    \
  {\n		`mv $tmp_dir/result2.aln $list[1]`;\n	     \
 }\n	    shift (@list); shift (@list);\n	    forea\
ch $f (@list)\n	      {\n		if (-e (\"$tmp_dir/$f\"\
)){`mv $tmp_dir/$f .`;}\n	      }\n	  }\n      }\n\
  }\n\n\n\n\nsub my_get_opt\n  {\n    my @list=@_;\
\n    my ($cl, $a, $argv, @argl);\n\n    \n    @ar\
gl=();\n    $cl=shift @list;\n    for ( my $a=0; $\
a<=$#list; $a+=3)\n      {\n	my $option=$list[$a];\
\n	my $optional=$list[$a+1];\n	my $status=$list[$a\
+2];\n	my $argv=\"\";\n	if ($cl=~/$option(\\S+)/){\
$argv=$1;}\n	@argl=(@argl,$argv);\n\n\n	#$optional\
:0=>optional\n	#$optional:1=>must be set\n	#$statu\
s: 0=>no requirement\n	#$status: 1=>must be an exi\
sting file\n	#$status: 2=>must be an installed pac\
kage\n	\n\n	if ($optional==0){;}\n	elsif ( $option\
al==1 && $argv eq \"\")\n	  {\n	    myexit(flush_e\
rror( \"ERROR: Option $option must be set\"));\n	 \
   myexit ($EXIT_FAILURE);\n	  }\n	if ($status==0)\
{;}\n	elsif ($status ==1 && $argv ne \"\" && !-e $\
argv)\n	  {\n	    myexit(flush_error( \"File [$arg\
v] must exist\"));\n	    myexit ($EXIT_FAILURE);\n\
	  }\n	elsif ( $status==2 && $argv ne \"\" && &che\
ck_pg_is_installed ($argv)==0)\n	  {\n	    myexit(\
flush_error( \" $argv is not installed\"));\n	    \
myexit ($EXIT_FAILURE);\n	  }\n      }\n    return\
 @argl;\n    }\n\nsub check_file\n  {\n    my ($fi\
le, $msg)=@_;\n\n    if ( !-e $file)\n      {\n	my\
exit(flush_error(\"$msg\"));\n      }\n    }\nsub \
hhalign\n  {\n    my ($aln1, $aln2, $outfile, $par\
am)=@_;\n    my $hh=$ENV{\"HHALIGN_4_TCOFFEE\"};\n\
    \n    \n    if ($hh)\n      {\n	\n	#external_h\
halign\n	# set via HHALIGN_4_TCOFFEE\n	#<pg> -prof\
ile1 <fasta_prf with seq1 top> -profile2 <fasta pr\
ofile with seq2 top> -outfile < fasta alignmentof \
seq1 and 2 | tc_lib of seq 1 and 2>\n	\n	safe_syst\
em (\"$hh -method=align -profile1=$aln1 -profile2=\
$aln2 -outfile=$outfile\");\n      }\n    else\n  \
    {\n	&local_hhalign ($aln1, $aln2, $outfile, $p\
aram);\n      }\n  }\n\n    \n    \nsub local_hhal\
ign\n  {\n    my ($aln1, $aln2, $outfile, $param)=\
@_;\n    my $h1, $h2;\n\n    $h{0}{index}=0;\n    \
$h{1}{index}=1;\n\n    $h{0}{aln}=$aln1;\n    $h{1\
}{aln}=$aln2;\n\n\n\n    %{$h{0}}=aln2psi_profile \
(%{$h{0}});\n    %{$h{1}}=aln2psi_profile (%{$h{1}\
});\n\n    $param=~s/#S/ /g;\n    $param=~s/#M/\\-\
/g;\n    $param=~s/#E/\\=/g;\n\n\n\n    $command=\\
"hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.\
tmp -rank 1 -mapt 0 $param\";\n    `$command`;\n\n\
  #  `hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $out\
file.tmp -rank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;\n\n\
\n    # To run global use the following\n\n    ope\
n (I, \"$outfile.tmp\");\n    open (O, \">$outfile\
\");\n    $h{0}{cons}=s/\\./x/g;\n    $h{1}{cons}=\
s/\\./x/g;\n\n    print O \"! TC_LIB_FORMAT_01\\n2\
\\n$h{0}{name} $h{0}{len} $h{0}{seq}\\n$h{1}{name}\
 $h{1}{len} $h{1}{seq}\\n#1 2\\n\";\n\n    while (\
<I>)\n      {\n	if (/(\\d+)\\s+(\\d+)\\s+(\\d+)/)\\
n	  {\n	    print O \"\\t$h{0}{$1}\\t$h{1}{$2}\\t$\
3\\n\";\n	  }\n      }\n    print O \"! SEQ_1_TO_N\
\\n\";\n\n    close (O);\n    close (I);\n  }\n\ns\
ub aln2psi_profile\n  {\n    my (%h)=@_;\n    my (\
$aln,$i,$hv, $a, @c, $n);\n\n\n    $i=$h{index};\n\
    $aln=$h{aln};\n\n    `cp $aln $$.hhh_aln`;\n  \
  $command=\"t_coffee -other_pg seq_reformat -in $\
aln -output hasch\";\n    $hv=`$command`;chomp ($h\
v);\n\n    $h{a2m}=\"$tmp/$hv.tmp4hhpred.a2m\";\n \
   $h{a3m}=\"$tmp/$hv.tmp4hhpred.a3m\";\n    if ( \
-e $h{a3m}){;}\n    else\n      {\n	$x=`which hhco\
nsensus`;\n	`hhconsensus  -M 50 -i $h{aln} -oa2m $\
h{a2m}`;\n	if (!-e $h{a2m})\n	  {\n	    print STDE\
RR \"Program tc_generic_method.pl FAILED to run:\\\
n\\thhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}\";\
\n	    myexit ($EXIT_FAILURE);\n	  }\n\n	`hhconsen\
sus  -M 50 -i $h{aln} -oa3m $h{a3m}`;\n	if (!-e $h\
{a3m})\n	  {\n	    print STDERR \"Program tc_gener\
ic_method.pl FAILED to run:\\n\\thhconsensus  -M 5\
0 -i $h{aln} -oa3m $h{a3m}\";\n	    myexit ($EXIT_\
FAILURE);\n	  }\n       `buildali.pl $h{a3m} -n 1`\
;\n      }\n\n\n    $h{a2m_seq}=`head -n 2 $h{a2m}\
 | grep -v \">\"`;chomp ($h{a2m_seq});\n    $h{a3m\
_seq}=`head -n 2 $h{a3m} | grep -v \">\"`;chomp ($\
h{a3m_seq});\n    $h{cons}=$h{a2m_seq};\n    $h{se\
q}=`head -n 2 $h{aln} | grep -v \">\"`;chomp ($h{s\
eq});\n\n\n\n    @c=split (//, $h{cons});\n    $h{\
len}=$#c+1;\n    for ($n=0,$a=0, $b=0; $a<$h{len};\
$a++)\n      {\n	if ( $c[$a]=~/[A-Z]/)\n	  {\n	   \
 $h{++$n}=++$b;\n\n	  }\n	elsif ( $c[$a]=~/[a-z\\.\
]/)\n	  {\n	    ++$b;\n	  }\n      }\n\n    $name=\
`head -n 2 $h{aln} | grep \">\"`;\n    $name=~/\\>\
(\\S+)/;\n    $h{name}=$1;\n\n    `cp $h{a2m} $i.a\
2m`;\n    `cp $h{a3m} $i.a3m`;\n    `cp $h{aln} $i\
.hh_aln`;\n\n    return %h;\n  }\nsub read_fasta_s\
eq_index\n  {\n    my $f=@_[0];\n    my %hseq;\n  \
  my (@seq, @com, @name);\n    my ($a, $s,$nseq);\\
n\n    open (F, $f);\n    while (<F>)\n      {\n	$\
s.=$_;\n      }\n    close (F);\n\n\n    @name=($s\
=~/>(\\S*).*\\n[^>]*/g);\n\n    @seq =($s=~/>.*.*\\
\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)\
/g);\n\n\n    $nseq=$#name+1;\n\n    for ($a=0; $a\
<$nseq; $a++)\n      {\n	my $s;\n	my $n=$name[$a];\
\n	$hseq{$a}{name}=$n;\n	$seq[$a]=~s/[^A-Za-z]//g;\
\n	$hseq{$a}{order}=$a;\n	$hseq{$a}{seq}=$seq[$a];\
\n	$hseq{$a}{com}=$com[$a];\n\n      }\n    return\
 %hseq;\n  }\nsub read_fasta_seq\n  {\n    my $f=@\
_[0];\n    my %hseq;\n    my (@seq, @com, @name);\\
n    my ($a, $s,$nseq);\n\n    open (F, $f);\n    \
while (<F>)\n      {\n	$s.=$_;\n      }\n    close\
 (F);\n\n\n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n\
\n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =(\
$s=~/>\\S*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name\
+1;\n\n    for ($a=0; $a<$nseq; $a++)\n      {\n	m\
y $s;\n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$\
seq[$a]=~s/[^A-Za-z]//g;\n	$hseq{$n}{order}=$a;\n	\
$hseq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a]\
;\n\n      }\n    return %hseq;\n  }\n\n\nsub read\
_fasta_aln\n  {\n    my $f=@_[0];\n    my %hseq;\n\
    my (@seq, @com, @name);\n    my ($a, $s,$nseq)\
;\n\n    open (F, $f);\n    while (<F>)\n      {\n\
	$s.=$_;\n      }\n    close (F);\n\n\n    @name=(\
$s=~/>(\\S*).*\\n[^>]*/g);\n\n    @seq =($s=~/>.*.\
*\\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]\
*)/g);\n\n\n    $nseq=$#name+1;\n\n    for ($a=0; \
$a<$nseq; $a++)\n      {\n	my $s;\n	my $n=$name[$a\
];\n	$hseq{$n}{name}=$n;\n	$seq[$a]=~s/[^A-Za-z-.(\
)[\\]]//g;\n	$hseq{$n}{order}=$a;\n	$hseq{$n}{seq}\
=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n\n      }\n\
    return %hseq;\n  }\n\nsub recode_name2\n{\n	my\
 ($in)=shift;\n	my $mode=shift;\n\n	my %seq;\n	my \
$new_name;\n\n	if (! -e $in){return;}\n\n	#needed \
by ClustalOmega to avoid very long names\n	open (I\
NFILE, \"+<$in\");\n\n	my $line;\n\n	if ($mode eq \
\"code\")\n	{\n		chomp($line = <INFILE>);\n		my $l\
ine_length = length($line);\n		$new_name=++$RECODE\
_N;\n		$new_name=\">$new_name\";\n		my $new_length\
 = length($new_name);\n		$RECODE {$new_name}=$line\
;\n		for ($count = $new_length; $count < $line_len\
gth; $count++)\n		{\n			$new_name .= \" \";\n		}\n\
		$new_name=\"$new_name\\n\";\n		seek INFILE, 0, 0\
\n			or die \"could not seek: $!\";\n		print INFIL\
E \"$new_name\";\n	}\n	else\n	{\n		my $n_found = 0\
;\n		my $file_pos=0;\n		$file_pos=tell INFILE;\n		\
while (<INFILE>)\n		{\n			$line=$_;\n			$line =~ s\
/\\s*$//;\n\n			$old_name= $RECODE{$line};\n			if \
($old_name ne \"\")\n			{\n				seek INFILE, $file_\
pos, 0\n					or die \"could not seek: $!\";\n				p\
rint INFILE \"$old_name\\n\";\n				$file_pos++;\n	\
			if ($file_pos == 2)\n				{\n					print \"stop\\\
n\";\n					break;\n				}\n			}\n			$file_pos=tell \
INFILE;\n		}\n\n	}\n\n\n	close INFILE;\n}\n\n\nsub\
 recode_name\n{\n	my ($in)=shift;\n	my $mode=shift\
;\n	my $f=new FileHandle;\n	my %seq;\n	my $new_nam\
e;\n\n	if (! -e $in){return;}\n\n	#needed by Clust\
alOmega to avoid very long names\n	%seq=read_fasta\
_aln ($in);\n\n	open ($f, \">$in\");\n	foreach my \
$s (keys(%seq))\n	{\n		if ($mode eq \"code\")\n		{\
\n			$new_name=++$RECODE_N;\n			$RECODE {$new_name\
}=$seq{$s}{name};\n		}\n		else\n		{\n			$new_name=\
$RECODE{$seq{$s}{name}};\n		}\n		print $f \">$new_\
name\\n$seq{$s}{seq}\\n\";\n	}\n	close $f;\n}\n\ns\
ub fasta_hash2index_hash\n  {\n    my %s1=@_;\n   \
 my %s;\n    foreach my $s (keys (%s1))\n      {\n\
	my $i=$s1{$s}{order};\n	$s{$i}{name}=$s;\n	$s{$i}\
{seq}=$s1{$s}{seq};\n	$s{$i}{len}=length( $s{$i}{s\
eq});\n	$s{n}++;\n      }\n    return %s;\n  }\nsu\
b file_contains\n  {\n    my ($file, $tag, $max)=(\
@_);\n    my ($n);\n    $n=0;\n\n    if ( !-e $fil\
e && ($file =~/$tag/)) {return 1;}\n    elsif ( !-\
e $file){return 0;}\n    else\n      {\n	open (FC,\
 \"$file\");\n	while ( <FC>)\n	  {\n	    if ( ($_=\
~/$tag/))\n	      {\n		close (FC);\n		return 1;\n	\
      }\n	    elsif ($max && $n>$max)\n	      {\n	\
	close (FC);\n		return 0;\n	      }\n	    $n++;\n	\
  }\n      }\n    close (FC);\n    return 0;\n  }\\
n\n\nsub file2string\n  {\n    my $f=@_[0];\n    m\
y $string, $l;\n    open (F,\"$f\");\n    while (<\
F>)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=\
$l;\n      }\n    close (F);\n    $string=~s/\\r\\\
n//g;\n    $string=~s/\\n//g;\n    return $string;\
\n  }\n\n\nsub tag2value\n  {\n\n    my $tag=(@_[0\
]);\n    my $word=(@_[1]);\n    my $return;\n\n   \
 $tag=~/$word=\"([^\"]+)\"/;\n    $return=$1;\n   \
 return $return;\n  }\n\nsub hit_tag2pdbid\n  {\n \
   my $tag=(@_[0]);\n    my $pdbid;\n\n    $tag=~/\
id=\"(\\S+)\"/;\n    $pdbid=$1;\n    $pdbid=~s/_//\
;\n    return $pdbid;\n  }\nsub id2pdbid\n  {\n   \
 my $in=@_[0];\n    my $id;\n\n    $in=~/(\\S+)/;\\
n    $id=$in;\n    $id=~s/PDB/pdb/g;\n\n    if ($i\
d =~/pdb(.*)/){$id=$1;}\n    elsif ( $id=~/(\\S+)\\
\s+mol:protein/){$id=$1;}\n    $id=~s/[:|��_]/\
/g;\n    return $id;\n  }\nsub set_blast_type\n  {\
\n    my $file =@_[0];\n    if (&file_contains ($f\
ile,\"EBIApplicationResult\",100)){$BLAST_TYPE=\"E\
BI\";}\n    elsif (&file_contains ($file,\"NCBI_Bl\
astOutput\",100)) {$BLAST_TYPE=\"NCBI\";}\n    els\
e\n      {\n	$BLAST_TYPE=\"\";\n      }\n    retur\
n $BLAST_TYPE;\n  }\nsub is_valid_blast_xml\n    {\
\n      my $file=shift;\n      my $line;\n\n\n    \
  if ( !-e $file) {return 0;}\n      $line=&file2t\
ail ($file,100);\n\n      if ( $line=~/<\\/EBIAppl\
icationResult/ || $line=~/<\\/NCBI_BlastOutput/ ||\
 $line=~/<\\/BlastOutput/ ){return 1;}\n      retu\
rn 0;\n    }\nsub file2blast_flavor\n      {\n	my \
$file=shift; \n	if (&file_contains ($file,\"EBIApp\
licationResult\",100)){return \"EBI\";}\n	elsif (&\
file_contains ($file,\"NCBI_BlastOutput\",100)){re\
turn \"NCBI\";}\n	else {return \"UNKNOWN\";}\n    \
  }\nsub blast2prf\n	{\n	  my ($blastF, $seqF,$out\
file)=@_;\n	  my (%s, %profile);\n	  my ($result,$\
psiblast_output,$profile_name,@profiles);\n	  %s=r\
ead_fasta_seq_index ($seqF);\n	  my ($z1,$z1m)=unc\
ompress($blastF);\n	  %profile=blast_xml2profile($\
s{0}{name}, $s{0}{seq},$maxid, $minid,$mincov,$bla\
stF);\n	  output_profile ($outfile, \\%profile, $t\
rim);\n	  compress($z1,$z1m);\n	}\nsub blast_xml2p\
rofile\n  {\n    my ($name,$seq,$maxid, $minid, $m\
incov, $file)=(@_);\n    my (%p, $a, $string, $n);\
\n\n\n\n    if ($BLAST_TYPE eq \"EBI\" || &file_co\
ntains ($file,\"EBIApplicationResult\",100)){%p=eb\
i_blast_xml2profile(@_);}\n    elsif ($BLAST_TYPE \
eq \"NCBI\" || &file_contains ($file,\"NCBI_BlastO\
utput\",100)){%p=ncbi_blast_xml2profile(@_);}\n   \
 else\n      {\n	myexit(add_error ( $$,$$,getppid(\
), \"BLAST_FAILURE::unkown XML\",$CL));\n      }\n\
    for ($a=0; $a<$p{n}; $a++)\n      {\n	my $name\
=$p{$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n	$p{\
$name}{index}=$a;\n      }\n    return %p;\n  }\ns\
ub ncbi_tblastx_xml2lib_file\n  {\n    my  ($outli\
b,$string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$i,$n\
hits,@identifyerL);\n    my (%ITERATION);\n\n    o\
pen (F, \">>$outlib\");\n\n    $seq=~s/[^a-zA-Z]//\
g;\n    $L=length ($seq);\n\n    %ITERATION=xml2ta\
g_list ($string, \"Iteration\");\n    for ($i=0; $\
i<$ITERATION{n};$i++)\n      {\n	my ($qindex, $qle\
n, %hit, $string);\n	$string=$ITERATION{$i}{body};\
\n\n	$qindex=xmltag2value($string,\"Iteration_iter\
-num\");\n	$qlen  =xmltag2value($string,\"Iteratio\
n_query-len\");\n	%hit=&xml2tag_list  ($string, \"\
Hit\");\n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	 \
   my ($string);\n	    $string=$hit{$a}{body};\n\n\
	    $hindex=xmltag2value($string,\"Hit_accession\\
")+1;\n	    if ($hindex<=$qindex){next;}\n	    els\
e  {print F  \"# $qindex $hindex\\n\";}\n\n\n	    \
$hlen=xmltag2value  ($string,\"Hit_len\");\n	    %\
HSP=&xml2tag_list  ($string, \"Hsp\");\n\n	    for\
 ($b=0; $b<$HSP{n}; $b++)\n	      {\n		my ($string\
, $qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e);\n		$string=\
$HSP{$b}{body};\n\n		$qs=xmltag2value  ($string,\"\
Hsp_query-from\");\n		$qe=xmltag2value  ($string,\\
"Hsp_query-to\");\n		$qf=xmltag2value  ($string,\"\
Hsp_query-frame\");\n\n		$hs=xmltag2value  ($strin\
g,\"Hsp_hit-from\");\n		$he=xmltag2value  ($string\
,\"Hsp_hit-to\");\n		$hf=xmltag2value  ($string,\"\
Hsp_hit-frame\");\n\n		$s=xmltag2value  ($string,\\
"Hsp_identity\");\n		$l=xmltag2value  ($string,\"H\
sp_align-len\");\n		$s=int(($s*100)/$l);\n\n		if (\
$qf>0)\n		  {$rqs=$qs; $rqe=$qe;}\n		else\n		  {\n\
		    $rqe=($qlen-$qs)+1;\n		    $rqs=($qlen-$qe)+\
1;\n		  }\n\n		if ($hf>0)\n		  {$rhs=$hs; $rhe=$he\
;}\n		else\n		  {\n		    $rhe=($hlen-$hs)+1;\n		  \
  $rhs=($hlen-$he)+1;\n		  }\n		for ($d=0,$e=$rqs;\
 $e<$rqe; $e++,$d++)\n		  {\n		    my ($r1,$r2);\n\
		    $r1=$e;\n		    $r2=$rhs+$d;\n		    print F \\
" $r1 $r2 $s 0\\n\";\n		  }\n	      }\n	  }\n     \
 }\n    print F \"! SEQ_1_TO_N\\n\";\n\n    close \
(F);\n    return %lib;\n  }\n\nsub ncbi_tblastpx_x\
ml2lib_file\n  {\n    my  ($outlib,$string,%s)=(@_\
);\n    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identify\
erL);\n    my (%ITERATION,%hdes, %qdes);\n\n    op\
en (F, \">>$outlib\");\n\n    $seq=~s/[^a-zA-Z]//g\
;\n    $L=length ($seq);\n\n    %ITERATION=xml2tag\
_list ($string, \"Iteration\");\n    for ($i=0; $i\
<$ITERATION{n};$i++)\n      {\n	my ($qindex, $qlen\
, %hit, $string);\n	$string=$ITERATION{$i}{body};\\
n\n	$qdef=xmltag2value($string,\"Iteration_query-d\
ef\");\n	%qdes=&tblastpx_name2description($qdef,%s\
);\n	$qlen  =xmltag2value($string,\"Iteration_quer\
y-len\");\n	%hit=&xml2tag_list  ($string, \"Hit\")\
;\n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	    my \
($string);\n	    $string=$hit{$a}{body};\n	    $hd\
ef=xmltag2value($string,\"Hit_def\");\n	    %hdes=\
&tblastpx_name2description($hdef,%s);\n	    if ($h\
des{index}<=$qdes{index}){next;}\n	    else  {prin\
t F  \"# $qdes{index} $hdes{index}\\n\";}\n\n\n	  \
  $hlen=xmltag2value  ($string,\"Hit_len\");\n	   \
 %HSP=&xml2tag_list  ($string, \"Hsp\");\n\n	    f\
or ($b=0; $b<$HSP{n}; $b++)\n	      {\n		my ($stri\
ng, $l,$qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e, @s1, @s\
2);\n		$string=$HSP{$b}{body};\n\n		$qs=xmltag2val\
ue  ($string,\"Hsp_query-from\");\n		$qe=xmltag2va\
lue  ($string,\"Hsp_query-to\");\n		$qf=$qdes{fram\
e};\n		$qseq=xmltag2value  ($string,\"Hsp_qseq\");\
\n\n		$hs=xmltag2value  ($string,\"Hsp_hit-from\")\
;\n		$he=xmltag2value  ($string,\"Hsp_hit-to\");\n\
		$hf=$hdes{frame};\n		$hseq=xmltag2value  ($strin\
g,\"Hsp_hseq\");\n\n		$s=xmltag2value  ($string,\"\
Hsp_identity\");\n		$l=xmltag2value  ($string,\"Hs\
p_align-len\");\n		$s=int(($s*100)/$l);\n		@s1=tbl\
astpx_hsp2coordinates($qseq,$qs,$qe,%qdes);\n		@s2\
=tblastpx_hsp2coordinates($hseq,$hs,$he,%hdes);\n\\
n\n		for ($f=0; $f<=$#s1; $f++)\n		  {\n		    if (\
$s1[$f]==-1 || $s2[$f]==-1){next;}\n		    else\n		\
      {\n			print F \" $s1[$f] $s2[$f] $s 0\\n\";\\
n		      }\n		  }\n	      }\n	  }\n      }\n    pr\
int F \"! SEQ_1_TO_N\\n\";\n\n    close (F);\n    \
return %lib;\n  }\nsub tblastpx_hsp2coordinates\n \
 {\n    my ($seq, $s, $e, %des)=@_;\n    my @list;\
\n    my @sa;\n    my @gap=(-1,-1,-1);\n\n    $s=$\
des{start}+3*($s-1);\n\n    if ($des{strand} eq \"\
d\"){;}\n    else {$s=($des{length}-$s)+1;}\n\n   \
 foreach $c (split (//,$seq))\n      {\n	if ( $c e\
q '-'){push (@list,@gap);}\n	elsif ($des{strand} e\
q \"d\")\n	  {\n	    push(@list,$s++,$s++,$s++);\n\
	  }\n	else\n	  {\n	    push(@list, $s--,$s--,$s--\
);\n	  }\n      }\n    return @list;\n  }\n\nsub t\
blastpx_name2description\n  {\n    my ($name, %s)=\
@_;\n    my @at=split(\"__\", $name);\n    my %des\
;\n\n    $des{name}=$at[0];\n    $des{strand}=$at[\
1];\n\n    $des{start}=$at[2];\n    $des{end}=$at[\
3];\n    $des{length}=$at[4];\n    $des{index}=$s{\
$at[0]}{order}+1;\n    return %des;\n  }\nsub ncbi\
_blast_xml2profile\n  {\n    my ($name,$seq,$maxid\
, $minid, $mincov, $string)=(@_);\n    my ($L,$l, \
$a,$b,$c,$d,$nhits,@identifyerL);\n\n    $seq=~s/[\
^a-zA-Z]//g;\n    $L=length ($seq);\n    \n    #Th\
is is causing the NCBI parser to fail when Iterati\
on_query-def is missing\n    #%query=&xml2tag_list\
 ($string, \"Iteration_query-def\");\n    #$name=$\
query{0}{body};\n\n    %hit=&xml2tag_list ($string\
, \"Hit\");\n\n    \n    for ($nhits=0,$a=0; $a<$h\
it{n}; $a++)\n      {\n	my ($ldb,$id, $identity, $\
expectation, $start, $end, $coverage, $r);\n	my (%\
ID,%DE,%HSP);\n\n	$ldb=\"\";\n\n	%ID=&xml2tag_list\
 ($hit{$a}{body}, \"Hit_id\");\n	$identifyer=$ID{0\
}{body};\n\n	%DE=&xml2tag_list ($hit{$a}{body}, \"\
Hit_def\");\n	$definition=$DE{0}{body};\n\n	%HSP=&\
xml2tag_list ($hit{$a}{body}, \"Hsp\");\n	for ($b=\
0; $b<$HSP{n}; $b++)\n	  {\n	    my (%START,%END,%\
E,%I,%Q,%M);\n\n\n	    %START=&xml2tag_list ($HSP{\
$b}{body}, \"Hsp_query-from\");\n	    %HSTART=&xml\
2tag_list ($HSP{$b}{body}, \"Hsp_hit-from\");\n\n	\
    %LEN=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_al\
ign-len\");\n	    %END=  &xml2tag_list ($HSP{$b}{b\
ody}, \"Hsp_query-to\");\n	    %HEND=  &xml2tag_li\
st ($HSP{$b}{body}, \"Hsp_hit-to\");\n	    %E=&xml\
2tag_list     ($HSP{$b}{body}, \"Hsp_evalue\");\n	\
    %I=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_id\
entity\");\n	    %Q=&xml2tag_list     ($HSP{$b}{bo\
dy}, \"Hsp_qseq\");\n	    %M=&xml2tag_list     ($H\
SP{$b}{body}, \"Hsp_hseq\");\n\n	    for ($e=0; $e\
<$Q{n}; $e++)\n\n	      {\n		$qs=$Q{$e}{body};\n		\
$ms=$M{$e}{body};\n\n		$expectation=$E{$e}{body};\\
n		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$L\
EN{$e}{body}*100;\n		$start=$START{$e}{body};\n		$\
end=$END{$e}{body};\n		$Hstart=$HSTART{$e}{body};\\
n		$Hend=$HEND{$e}{body};\n\n		$coverage=($L)?(($e\
nd-$start)*100)/$L:0;\n		if ($identity>$maxid || $\
identity<$minid || $coverage<$mincov)\n		  {\n		  \
  next;\n		  }\n		@lr1=(split (//,$qs));\n		@lr2=(\
split (//,$ms));\n		$l=$#lr1+1;\n		for ($c=0;$c<$L\
;$c++){$p[$nhits][$c]=\"-\";}\n		for ($d=0,$c=0; $\
c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\n		    if (\
 $r=~/[A-Za-z]/)\n		      {\n\n			$p[$nhits][$d + \
$start-1]=$lr2[$c];\n			$d++;\n		      }\n		  }\n	\
	$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\n		$Qsta\
rtL[$nhits]=$start;\n		$HstartL[$nhits]=$Hstart;\n\
		$identityL[$nhits]=$identity;\n		$endL[$nhits]=$\
end;\n		$definitionL[$nhits]=$definition;\n		$iden\
tifyerL[$nhits]=$identifyer;\n		$comment[$nhits]=\\
"$ldb|$identifyer [Eval=$expectation][id=$identity\
%][start=$Hstart end=$Hend]\";\n		$nhits++;\n	    \
  }\n	  }\n      }\n\n\n    $profile{n}=0;\n    $p\
rofile{$profile{n}}{name}=$name;\n    $profile{$pr\
ofile{n}}{seq}=$seq;\n    $profile {n}++;\n\n    f\
or ($a=0; $a<$nhits; $a++)\n      {\n	$n=$a+1;\n\n\
	$profile{$n}{name}=\"$name\\_$a\";\n	$profile{$n}\
{seq}=\"\";\n	$profile{$n}{Qseq}=$Qseq[$a];\n	$pro\
file{$n}{Hseq}=$Hseq[$a];\n	$profile{$n}{Qstart}=$\
QstartL[$a];\n	$profile{$n}{Hstart}=$HstartL[$a];\\
n	$profile{$n}{identity}=$identityL[$a];\n	$profil\
e{$n}{definition}=$definitionL[$a];\n	$profile{$n}\
{identifyer}=$identifyerL[$a];\n	$profile{$n}{comm\
ent}=$comment[$a];\n\n	for ($b=0; $b<$L; $b++)\n	 \
 {\n	    if ($p[$a][$b])\n	      {\n		$profile{$n}\
{seq}.=$p[$a][$b];\n	      }\n	    else\n	      {\\
n		$profile{$n}{seq}.=\"-\";\n	      }\n	  }\n    \
  }\n\n    $profile{n}=$nhits+1;\n    return %prof\
ile;\n  }\nsub ebi_blast_xml2profile\n  {\n    my \
($name,$seq,$maxid, $minid, $mincov, $string)=(@_)\
;\n    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,\
$identifyer);\n\n\n\n    $seq=~s/[^a-zA-Z]//g;\n  \
  $L=length ($seq);\n    %hit=&xml2tag_list ($stri\
ng, \"hit\");\n\n    for ($nhits=0,$a=0; $a<$hit{n\
}; $a++)\n      {\n	my ($ldb,$id, $identity, $expe\
ctation, $start, $end, $coverage, $r);\n	my (%Q,%M\
,%E,%I);\n\n	$ldb=&tag2value ($hit{$a}{open}, \"da\
tabase\");\n	$identifyer=&tag2value ($hit{$a}{open\
}, \"id\");\n\n	$description=&tag2value ($hit{$a}{\
open}, \"description\");\n\n	%Q=&xml2tag_list ($hi\
t{$a}{body}, \"querySeq\");\n	%M=&xml2tag_list ($h\
it{$a}{body}, \"matchSeq\");\n	%E=&xml2tag_list ($\
hit{$a}{body}, \"expectation\");\n	%I=&xml2tag_lis\
t ($hit{$a}{body}, \"identity\");\n\n\n	for ($b=0;\
 $b<$Q{n}; $b++)\n	  {\n\n	    $qs=$Q{$b}{body};\n\
	    $ms=$M{$b}{body};\n\n	    $expectation=$E{$b}\
{body};\n	    $identity=$I{$b}{body};\n\n\n	    $s\
tart=&tag2value ($Q{$b}{open}, \"start\");\n	    $\
end=&tag2value ($Q{$b}{open}, \"end\");\n	    $sta\
rtM=&tag2value ($M{$b}{open}, \"start\");\n	    $e\
ndM=&tag2value ($M{$b}{open}, \"end\");\n	    $cov\
erage=(($end-$start)*100)/$L;\n\n	   # print \"$id\
: ID: $identity COV: $coverage [$start $end]\\n\";\
\n\n	    if ($identity>$maxid || $identity<$minid \
|| $coverage<$mincov){next;}\n	    # print \"KEEP\\
\n\";\n\n\n	    @lr1=(split (//,$qs));\n	    @lr2=\
(split (//,$ms));\n	    $l=$#lr1+1;\n	    for ($c=\
0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n	    for ($d\
=0,$c=0; $c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n	\
	if ( $r=~/[A-Za-z]/)\n		  {\n\n		    $p[$nhits][$\
d + $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	    \
  }\n\n	    $Qseq[$nhits]=$qs;\n	    $Hseq[$nhits]\
=$ms;\n	    $QstartL[$nhits]=$start;\n	    $Hstart\
L[$nhits]=$Hstart;\n	    $identityL[$nhits]=$ident\
ity;\n	    $endL[$nhits]=$end;\n	    $definitionL[\
$nhits]=$definition;\n	    $identifyerL[$nhits]=$i\
dentifyer;\n	    $comment[$nhits]=\"$ldb|$identify\
er [Eval=$expectation][id=$identity%][start=$start\
M end=$endM]\";\n	    $nhits++;\n	  }\n      }\n\n\
    $profile{n}=0;\n    $profile{$profile{n}}{name\
}=$name;\n    $profile{$profile{n}}{seq}=$seq;\n  \
  $profile {n}++;\n\n    for ($a=0; $a<$nhits; $a+\
+)\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$nam\
e\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}\
{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\\
n	$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n\
}{Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=$\
identityL[$a];\n	$profile{$n}{definition}=$definit\
ionL[$a];\n	$profile{$n}{identifyer}=$identifyerL[\
$a];\n	$profile{$n}{comment}=$comment[$a];\n\n	for\
 ($b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n\
	      {\n		$profile{$n}{seq}.=$p[$a][$b];\n	     \
 }\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\\
";\n	      }\n	  }\n      }\n    $profile{n}=$nhit\
s+1;\n\n    return %profile;\n  }\nsub output_prof\
ile\n  {\n    my ($outfile,$profileR, $trim)=(@_);\
\n    my ($a);\n    my %profile=%$profileR;\n    m\
y $P= new FileHandle;\n    my $tmp=vtmpnam();\n\n \
   open ($P, \">$tmp\");\n    for ($a=0; $a<$profi\
le{n}; $a++)\n      {\n	print $P \">$profile{$a}{n\
ame} $profile{$a}{comment}\\n$profile{$a}{seq}\\n\\
";\n      }\n    close ($P);\n\n    if ( $trim)\n \
     {\n	if ($ENV{psitrim_mode_4_TCOFFEE} eq \"tri\
m\")# old trimming\n	  {\n	    if ($trim>0)\n	    \
  {\n		&safe_system (\"t_coffee -other_pg seq_refo\
rmat -in $tmp -action +trim _aln_n$trim\\_K1 -outp\
ut fasta_aln -out $outfile\");\n	      }\n	    els\
e\n	      {\n		&safe_system (\"t_coffee -other_pg \
seq_reformat -in $tmp -action +trim _aln_%%$trim\\\
_K1 -output fasta_aln -out $outfile\");\n	      }\\
n	  }\n	else # newtrimming\n	  {\n	    my $tm;\n	 \
   if ($ENV{psitrim_tree_4_TCOFFEE}){$tm=\"-treemo\
de=\".$ENV{psitrim_tree_4_TCOFFEE};}\n	    else {$\
tm=\"-treemode=codnd\";}\n	    if ($trim>0)\n	    \
  {\n		&safe_system (\"t_coffee -other_pg seq_refo\
rmat -in $tmp $tm -keep 1 -action +regtrim $trim -\
output fasta_aln -out $outfile\");\n	      }\n	   \
 else\n	      {\n		&safe_system (\"t_coffee -other\
_pg seq_reformat -in $tmp $tm -keep 1 -action +reg\
trim $trim\\% -output fasta_aln -out $outfile\");\\
n	      }\n	   \n	  }\n	\n      }\n    else\n     \
 {\n	&safe_system (\"mv $tmp $outfile\");\n      }\
\n    return;\n  }\nsub blast_xml2hit_list\n  {\n \
   my $string=(@_[0]);\n    return &xml2tag_list (\
$string, \"hit\");\n  }\nsub xmltag2value\n  {\n  \
  my ($string_in, $tag)=@_;\n    my %TAG;\n    %TA\
G=xml2tag_list ($string_in, $tag);\n    return $TA\
G{0}{body};\n  }\n\nsub xml2tag_list\n  {\n    my \
($string_in,$tag)=@_;\n    my $tag_in, $tag_out;\n\
    my %tag;\n\n    if (-e $string_in)\n      {\n	\
$string=&file2string ($string_in);\n      }\n    e\
lse\n      {\n	$string=$string_in;\n      }\n    $\
tag_in1=\"<$tag \";\n    $tag_in2=\"<$tag>\";\n   \
 $tag_out=\"/$tag>\";\n    $string=~s/>/>##1/g;\n \
   $string=~s/</##2</g;\n    $string=~s/##1/<#/g;\\
n    $string=~s/##2/#>/g;\n    @l=($string=~/(\\<[\
^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$n=-1;\n\n\
\n\n    foreach $t (@l)\n      {\n\n	$t=~s/<#//;\n\
	$t=~s/#>//;\n\n	if ( $t=~/$tag_in1/ || $t=~/$tag_\
in2/)\n	  {\n\n	    $in=1;\n	    $tag{$tag{n}}{ope\
n}=$t;\n	    $n++;\n\n	  }\n	elsif ($t=~/$tag_out/\
)\n	  {\n\n\n	    $tag{$tag{n}}{close}=$t;\n	    $\
tag{n}++;\n	    $in=0;\n	  }\n	elsif ($in)\n	  {\n\
\n	    $tag{$tag{n}}{body}.=$t;\n	  }\n      }\n\n\
    return %tag;\n  }\n\n\nsub seq2gor_prediction\\
n  {\n    my ($name, $seq,$infile, $outfile, $gor_\
seq, $gor_obs)=(@_);\n    my ($l);\n\n    `gorIV -\
prd $infile -seq $gor_seq -obs $gor_obs > gor_tmp`\
;\n    open (GR, \">$outfile\");\n    open (OG, \"\
gor_tmp\");\n\n    while (<OG>)\n      {\n\n	$l=$_\
;\n	if ($l=~/\\>/){print GR \"$l\";}\n	elsif ( $l=\
~/Predicted Sec. Struct./)\n	  {\n	    $l=~s/Predi\
cted Sec. Struct\\.//;\n	    print GR \"$l\";\n	  \
}\n      }\n    close (GR);\n    close (OG);\n    \
return;\n  }\nsub seq2msa_tm_prediction\n  {\n    \
my ($name, $seq, $db, $infile, $outfile, $arch, $p\
sv)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s,\
 $l);\n    my $R2=new FileHandle;\n    my $method=\
\"psitm\";\n    my $SERVER=\"EBI\";\n\n    $blast_\
output=&run_blast ($name,\"blastp\", $db, $infile,\
 \"outfile\");\n\n    if (&cache_file(\"GET\",$inf\
ile,$name,$method,$db,$outfile,$SERVER))\n      {\\
n	print \"\\tPSITM: USE Cache\\n\";\n	return $outf\
ile;\n      }\n    else\n      {\n	$CACHE_STATUS=\\
"COMPUTE CACHE\";\n	%p=blast_xml2profile($name,$se\
q,$maxid, $minid,$mincov,$blast_output);\n\n\n	ope\
n (F, \">tm_input\");\n	for (my $a=0; $a<$p{n}; $a\
++)\n	  {\n	    my $s;\n\n	    $s=$p{$a}{seq};\n	 \
   $s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\\
n\";\n	    #print stdout \">$p{$a}{name}\\n$s\\n\"\
;\n	  }\n	close (F);\n	print \"\\tPSITM: kept  $p{\
n} Homologues for Sequence $p{0}{name}\\n\";\n	&sa\
fe_system (\"t_coffee -other_pg fasta_seq2hmmtop_f\
asta.pl -in=tm_input -out=$outfile -output=cons -c\
ov=70 -trim=95 -arch=$arch -psv=$psv\");\n	unlink \
(\"tm_input\");\n	&cache_file(\"SET\",$infile,$nam\
e,$method,$db,$outfile,$SERVER);\n	return;\n      \
}\n  }\n\n\nsub seq2msa_gor_prediction\n  {\n    m\
y ($name, $seq,$infile, $outfile, $gor_seq, $gor_o\
bs)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s,\
 $l);\n    my $R2=new FileHandle;\n    my $db=\"un\
iprot\";\n    my $method=\"psigor\";\n    my $SERV\
ER=\"EBI\";\n\n    $blast_output=&run_blast ($name\
,\"blastp\", \"uniprot\", $infile, \"outfile\");\n\
\n    if (&cache_file(\"GET\",$infile,$name,$metho\
d,$db,$outfile,$SERVER, $psiJ))\n      {\n	print \\
"\\tPSIGOR: USE Cache\\n\";\n	return $outfile;\n  \
    }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE\
 CACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid\
, $minid,$mincov,$blast_output);\n\n\n	open (F, \"\
>gor_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	 \
 {\n	    my $s;\n\n	    $s=$p{$a}{seq};\n	    $s=u\
c($s);\n	    print F \">$p{$a}{name}\\n$s\\n\";\n	\
    #print stdout \">$p{$a}{name}\\n$s\\n\";\n	  }\
\n	close (F);\n	print \"\\tGORTM: kept  $p{n} Homo\
logues for Sequence $p{0}{name}\\n\";\n	&safe_syst\
em (\"t_coffee -other_pg fasta_seq2hmmtop_fasta.pl\
 -in=gor_input -out=$outfile -output=cons -cov=70 \
-trim=95 -gor_seq=$gor_seq -gor_obs=$gor_obs -mode\
=gor\");\n	unlink (\"tm_input\");\n	&cache_file(\"\
SET\",$infile,$name,$method,$db,$outfile,$SERVER);\
\n	return;\n      }\n  }\n\n\n\nsub run_blast\n  {\
\n    my ($name, $method, $db, $infile, $outfile, \
$run)=(@_);\n    if (!$run){$run=1;}\n    my $erro\
r_log=vtmpnam();\n    my $cl_db;\n    my $psiJ=($E\
NV{psiJ_4_TCOFFEE})?$ENV{psiJ_4_TCOFFEE}:1;\n    \\
n    my $psiJFlag=\"-j$psiJ\";\n   \n    if (&cach\
e_file(\"GET\",$infile,$name,$method,$db,$outfile,\
$SERVER, $psiJ) && is_valid_blast_xml ($outfile))\\
n      {return $outfile;}\n    else\n      {\n	$CA\
CHE_STATUS=\"COMPUTE CACHE\";\n	if ( $SERVER eq \"\
EBI_SOAP\")\n	  {\n	    &check_configuration (\"EM\
AIL\",\"SOAP::Light\",\"INTERNET\");\n\n	    $cl_m\
ethod=$method;\n	    if ($cl_method =~/wu/)\n	    \
  {\n		if ( $method eq \"psiblast\" || $psiJ>1){my\
exit(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BL\
AST_FAILURE::$SERVER does not Support psiblast mod\
e ($psiJFlag)\",$CL));}\n		$cl_method=~s/wu//;\n		\
if ( $cl_method eq \"psiblast\" || $psiJ>1)\n		  {\
\n		    add_warning($$,$$,\"PSI BLAST cannot be us\
ed with the wuBLAST Client. Use server=EBI Or serv\
er=LOCAL. blastp will be used instead\");\n		    $\
cl_method=\"blastp\";\n		  }\n\n		$command=\"t_cof\
fee -other_pg wublast.pl --email $EMAIL $infile -D\
 $db -p $cl_method --outfile $outfile -o xml>/dev/\
null 2>$error_log\";\n		&safe_system ( $command);\\
n		if (-e \"$outfile.xml\") {`mv $outfile.xml $out\
file`;}\n	      }\n	    else\n	      {\n		if ($cl_\
method eq \"psiblast\"){$cl_method =\"blastp $psiJ\
Flag\";}\n\n		$command=\"t_coffee -other_pg blastp\
gp.pl --email $EMAIL $infile -d $db --outfile $out\
file -p $cl_method --mode PSI-Blast>/dev/null 2>$e\
rror_log\";\n		&safe_system ( $command);\n\n		if (\
-e \"$outfile.xml\") {`mv $outfile.xml $outfile`;}\
\n	      }\n	  }\n	elsif ($SERVER eq \"EBI_REST\" \
|| $SERVER eq \"EBI\")\n	  {\n	    $cl_method=$met\
hod;\n	    &check_configuration(\"EMAIL\",\"XML::S\
imple\", \"INTERNET\");\n	    if ($db eq \"uniprot\
\"){$db1=\"uniprotkb\";}\n	    else {$db1=$db;}\n\\
n	    \n	    if ($cl_method =~/wu/)\n	      {\n		$\
cl_method=~s/wu//;\n\n		\n		if ( $cl_method eq \"p\
siblast\"){$cl_method=\"blastp\";}\n\n		$command=\\
"t_coffee -other_pg wublast_lwp.pl --email $EMAIL \
-D $db1 -p $cl_method --outfile $outfile --align 5\
 --stype protein $infile>/dev/null 2>error_log\";\\
n	      }\n	    else\n	      {\n		#if ( $cl_method\
 =~/psiblast/){$cl_method =\"blastp $psiJFlag\";}\\
n		if ( $cl_method =~/psiblast/){$cl_method =\"bla\
stp\";}\n		$command=\"t_coffee -other_pg ncbiblast\
_lwp.pl --email $EMAIL --database $db1 --program $\
cl_method --outfile $outfile --alignments 5 --styp\
e protein $infile>/dev/null 2>$error_log\";\n		#DE\
BUG\n		#$command=\"t_coffee -other_pg ncbiblast_lw\
p.pl --email $EMAIL -D $db1 -p $cl_method --outfil\
e $outfile --align 5 --stype protein $infile\";\n	\
	\n		my $maxrun=5;#number of crashes accepetd\n		m\
y $nrun;\n		my $keep_going=1;\n		while ($keep_goin\
g)\n		  {\n		    \n		    #print \"----- $command [\
$nrun]\\n\";\n		    $nrun++;\n		    $keep_going=0;\
\n		    &safe_system ( $command,5);\n		    \n		   \
 my $success=0;\n		    $success =$success || -e \"\
$outfile.out.xml\";\n		    $success =$success || -\
e \"$outfile.xml.xml\";\n		    $success =$success \
|| -e \"$outfile.out..xml\";\n		    $success =$suc\
cess || -e \"$outfile.xml..xml\";\n		    \n		    i\
f (!$success && ($nrun<$maxrun || -e \"$outfile.ou\
t.txt\"))\n		      {\n			$keep_going=1;\n			add_wa\
rning($$,$$,\"[ncbiblast_lwp.pl] [$command] failed\
 to produce xml output -- will ne tried again [$nr\
un]\");\n		      }\n		  }\n		\n		if (-e \"$outfile\
.out.xml\") {`mv $outfile.out.xml $outfile`;}\n		e\
lsif (-e \"$outfile.xml.xml\"){`mv $outfile.xml.xm\
l $outfile`;}\n		elsif (-e \"$outfile.out..xml\") \
{`mv $outfile.out..xml $outfile`;}\n		elsif (-e \"\
$outfile.xml..xml\"){`mv $outfile.xml..xml $outfil\
e`;}\n		else\n		  {\n		    add_warning($$,$$,\"[nc\
biblast_lwp.pl] [$command] failed to produce xml o\
utput\");\n		  }\n	      }\n	  }\n	elsif ($SERVER \
eq \"NCBI\")\n	  {\n	    &check_configuration (\"I\
NTERNET\");\n	    if ($db eq \"uniprot\"){$cl_db=\\
"swissprot\";}\n	    else {$cl_db=$db;}\n\n	    if\
 ( $method eq \"psiblast\" || $psiJ>1){myexit(add_\
error (EXIT_FAILURE,$$,$$,getppid(), \"BLAST_FAILU\
RE::$SERVER does not Support psiblast mode ($psiJF\
lag)\",$CL));}\n	    my $cl_method=$method;\n	    \
\n	    &check_configuration ($cl_method);  \n	    \
$command=\"$cl_method -db $cl_db -query $infile -o\
ut $outfile -outfmt 5 -remote\";\n	    &safe_syste\
m ($command);\n	  }\n	elsif ($SERVER =~/CLIENT_(.*\
)/)\n	  {\n	    my $client=$1;\n	    if ( $method \
eq \"psiblast\" || $psiJ>1){myexit(add_error (EXIT\
_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::$SERVER\
 does not Support psiblast mode ($psiJFlag)\",$CL)\
);}\n	    $command=\"$client -p $method -d $db -i \
$infile -o $outfile -m 7\";\n	    &safe_system ($c\
ommand);\n	  }\n	elsif ( $SERVER eq \"LOCAL_blasta\
ll\")\n	  {\n	    &check_configuration (\"blastall\
\");\n	    if ( $method eq \"psiblast\" || $psiJ>1\
){myexit(add_error (EXIT_FAILURE,$$,$$,getppid(), \
\"BLAST_FAILURE::$SERVER does not Support psiblast\
 mode ($psiJFlag)\",$CL));}\n	    $command=\"blast\
all -d $db -i $infile -o $outfile -m7 -p blastp\";\
\n	    &safe_system ($command);\n	  }\n	elsif ( $S\
ERVER eq \"LOCAL\")\n	  {\n	    my $legacy=0;\n	  \
  if ($ENV{\"BLAST_DB_DIR\"}) \n	      {\n	    	$x\
=$ENV{\"BLAST_DB_DIR\"};\n		$cl_db=\"$x/$db\";\n	 \
     }\n	    else\n	      {\n		$cl_db=$db;\n	     \
 }\n	    \n	    ##\n	    ## BLAST+ provide differe\
nt binaries names and CLI options\n	    ## Use the\
 'legacy_blast.pl' to keep compatibility with old \
blast commands\n	    ##\n	    $path=`which legacy_\
blast.pl 2>/dev/null`;  \n	    $path=`dirname $pat\
h`; \n	    chomp($path);\n	    \n	    if    (!$leg\
acy && ($method eq \"blastp\" || $method eq \"psib\
last\"))\n	      {\n		\n		&check_configuration(\"p\
siblast\");\n		$command=\"psiblast -db $cl_db -que\
ry $infile -num_iterations $psiJ -out $outfile -ou\
tfmt 5\";\n	      }\n	    elsif ($legacy && $metho\
d eq \"blastp\")\n	     {\n	       &check_configur\
ation(\"legacy_blast.pl\");\n	       $command=\"le\
gacy_blast.pl blastpgp --path $path -d $cl_db -i $\
infile -o $outfile -m7 $psiJFlag\";		\n	     }\n	 \
   elsif ($legacy && $method eq \"psiblast\")\n	  \
    {\n		&check_configuration(\"legacy_blast.pl\")\
;\n		$command=\"legacy_blast.pl blastpgp --path $p\
ath -d $cl_db -i $infile -o $outfile -m7 $psiJFlag\
\";\n	      }\n	    elsif ($method eq \"blastn\")\\
n	      {\n		&check_configuration(\"legacy_blast.p\
l\");\n		$command=\"legacy_blast.pl blastall --pat\
h $path -p blastn -d $cl_db -i $infile -o $outfile\
 -m7 -W6\";\n	      }\n	    &safe_system ($command\
);\n	  }\n	else\n	  {\n\n	    myexit(add_error (EX\
IT_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::Unkno\
wnServer\",$CL));\n	  }\n\n\n	#Check that everythi\
ng went well\n\n	if ( !-e $outfile || !&is_valid_b\
last_xml($outfile))\n	  {\n\n	    if ( -e $outfile\
)\n	      {\n		add_warning ($$,$$,\"Corrupted Blas\
t Output (Run $run)\");\n		unlink($outfile);\n	   \
   }\n	    if ( -e $error_log)\n	      {\n\n		my $\
error_msg=file2string ($error_log);\n\n		if ( $err\
or_msg =~/enter a valid email/)\n		  {\n		    myex\
it(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAS\
T_FAILURE::Invalid_or_rejected_email::$EMAIL\", \"\
$command\"));\n		  }\n	      }\n	    if ( $run==$B\
LAST_MAX_NRUNS)\n	      {\n\n		myexit(add_error (E\
XIT_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::Unkn\
ownReason\", \"$command\"));\n	      }\n	    else\\
n	      {\n		my $out;\n		if ($SERVER eq \"NCBI\") \
{$SERVER=\"EBI\"; }\n		elsif ($SERVER eq \"EBI\"){\
$SERVER=\"NCBI\";}\n		add_warning ($$,$$,\"Blast f\
or $name failed (Run: $run out of $BLAST_MAX_NRUNS\
. Use $SERVER)\");\n		$out=&run_blast ($name, $met\
hod, $db,$infile, $outfile, $run+1);\n		if ($SERVE\
R eq \"NCBI\") {$SERVER=\"EBI\"; }\n		elsif ($SERV\
ER eq \"EBI\"){$SERVER=\"NCBI\";}\n		return $out;\\
n	      }\n	  }\n\n	&cache_file(\"SET\",$infile,$n\
ame,$method,$db,$outfile,$SERVER, $psiJ);\n	return\
 $outfile;\n      }\n  }\n\nsub cache_file\n  {\n \
   my ($cache_mode,$infile,$name,$method,$db, $out\
file,$server,$it)=(@_);\n    my $cache_file;\n    \
#Protect names so that they can be turned into leg\
al filenames\n    $name=&clean_file_name ($name);\\
n    if (!$it){$it=1;}\n    if ($db=~/\\//)\n     \
 {\n	$db=~/([^\\/]+)$/;\n	$db=$1;\n      }\n\n    \
$cache_file=\"$CACHE/$name.$method.$db.$server.$it\
.tmp\";\n    #print \"Look for $cache_file [$cache\
_mode][$CACHE] \\n\";\n    if ($infile ne \"\"){$c\
ache_file_infile=\"$CACHE/$name.$method.$db.$serve\
r.$it.infile.tmp\";}\n\n    if ($cache_mode eq \"G\
ET\")\n      {\n	if ($CACHE eq \"\" || $CACHE eq \\
"no\" || $CACHE eq \"ignore\"  || $CACHE eq \"loca\
l\" || $CACHE eq \"update\"){return 0;}\n	elsif ( \
!-d $CACHE)\n	  {\n	    print STDERR \"ERROR: Cach\
e Dir: $CACHE Does not Exist\";\n	    return 0;\n	\
  }\n	else\n	  {\n	    my ($z1,$z1m)=uncompress($c\
ache_file_infile);\n	    my ($z2,$z2m)=uncompress(\
$cache_file);\n	    \n	    if ( -e $cache_file && \
&fasta_file1_eq_fasta_file2($infile,$cache_file_in\
file)==1)\n	      {\n		`cp $cache_file $outfile`;\\
n		$CACHE_STATUS=\"READ CACHE\";\n		compress($z1,$\
z1m);\n		compress($z2,$z2m);\n		\n		return 1;\n	  \
    }\n	  }\n      }\n    elsif ($cache_mode eq \"\
SET\")\n      {\n	if ($CACHE eq \"\" || $CACHE eq \
\"no\" || $CACHE eq \"ignore\"  || $CACHE eq \"loc\
al\" || $CACHE eq \"update\"){return 0;}\n	elsif (\
 !-d $CACHE)\n	  {\n	    print STDERR \"ERROR: Cac\
he Dir: $CACHE Does not Exist\";\n	    return 0;\n\
	  }\n	elsif (-e $outfile)\n	  {\n	    `cp $outfil\
e $cache_file`;\n	    if ($cache_file_infile ne \"\
\"){ `cp $infile $cache_file_infile`;}\n	    retur\
n 1;\n	  }\n      }\n    $CACHE_STATUS=\"COMPUTE C\
ACHE\";\n    return 0;\n  }\nsub file1_eq_file2\n \
 {\n    my ($f1, $f2)=@_;\n    if ( $f1 eq \"\"){r\
eturn 1;}\n    elsif ( $f2 eq \"\"){return 1;}\n  \
  elsif ( !-e $f1){return 0;}\n    elsif ( !-e $f2\
){return 0;}\n    elsif ($f1 eq \"\" || $f2 eq \"\\
" || `diff $f1 $f2` eq \"\"){return 1;}\n\n    ret\
urn 0;\n  }\nsub clean_file_name\n  {\n    my $nam\
e=@_[0];\n\n    $name=~s/[^A-Za-z1-9.-]/_/g;\n    \
return $name;\n  }\nsub url2file\n  {\n    my ($ad\
dress, $out)=(@_);\n\n    if (&pg_is_installed (\"\
wget\"))\n	{\n	  return &safe_system (\"wget $addr\
ess -O$out >/dev/null 2>/dev/null\");\n	}\n    els\
if (&pg_is_installed (\"curl\"))\n      {\n	return\
 &safe_system (\"curl $address -o$out >/dev/null 2\
>/dev/null\");\n      }\n    else\n      {\n	myexi\
t(flus_error(\"neither curl nor wget are installed\
. Imnpossible to fectch remote file\"));\n	exit ($\
EXIT_FAILURE);\n      }\n  }\nsub fasta_file1_eq_f\
asta_file2\n  {\n    my ($f1, $f2)=@_;\n    my (%s\
1, %s2);\n    my @names;\n    %s1=read_fasta_seq (\
$f1);\n    %s2=read_fasta_seq ($f2);\n\n    @names\
=(keys (%s1));\n\n    foreach $n (keys(%s1))\n    \
  {\n	my $ss1=lc($s1{$n}{seq});\n	my $ss2=lc($s2{$\
n}{seq});\n	if ($ss1 ne $ss2){return 0;}\n      }\\
n    foreach $n (keys(%s2))\n      {\n	my $ss1=lc(\
$s1{$n}{seq});\n	my $ss2=lc($s2{$n}{seq});\n	if ($\
ss1 ne $ss2){return 0;}\n      }\n    \n    return\
 1;\n  }\n\n\n\nsub read_template_file\n  {\n    m\
y $pdb_templates = @_[0];\n    my $tmp=new FileHan\
dle;\n    open ($tmp, \"<$pdb_templates\");\n    m\
y %temp_h;\n    my ($skip,$seq, $temp);\n\n    #su\
pports both a simple [seq pdb] format and the regu\
lar fasta like template format\n    while (<$tmp>)\
\n      {\n	\n	$line = $_;\n	if ($line=~/\\>(\\S+)\
\\s_._\\s(\\S+)/){$temp_h{$1}= $2;}\n	elsif ($line\
 =~/(\\S+)\\s(\\S+)/){$temp_h{$1}= $2;}\n      }\n\
    close($tmp);\n    return %temp_h;\n  }\n\n\n\n\
\n\n\nsub seq2tblastx_lib\n  {\n    my ($mode, $in\
file, $outfile)=@_;\n    my (%s, $method,$nseq);\n\
\n    $method=$mode;\n    &set_temporary_dir (\"se\
t\",$infile,\"infile\");\n    %s=read_fasta_seq(\"\
infile\");\n\n\n    foreach $seq (keys(%s))\n     \
 {\n	$slist[$s{$seq}{order}]=$s{$seq}{seq};\n	$sna\
me[$s{$seq}{order}]=$s{$seq}{name};\n	$slen[$s{$se\
q}{order}]=length ($s{$seq}{seq});\n      }\n    $\
nseq=$#sname+1;\n    open (F, \">outfile\");\n    \
print F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"\
$nseq\\n\";\n    for ($a=0; $a<$nseq;$a++)\n      \
{\n	print F \"$sname[$a] $slen[$a]  $slist[$a]\\n\\
"\n      }\n    close (F);\n    &safe_system (\"fo\
rmatdb -i infile -p F\");\n    &safe_system (\"bla\
stall -p tblastx -i infile -d infile -m 7 -S1>blas\
t.output\");\n\n    ncbi_tblastx_xml2lib_file (\"o\
utfile\", file2string (\"blast.output\"));\n    &s\
et_temporary_dir (\"unset\",$mode, $method, \"outf\
ile\",$outfile);\n    myexit ($EXIT_SUCCESS);\n   \
 }\nsub seq2tblastpx_lib\n  {\n    my ($mode, $inf\
ile, $outfile)=@_;\n    my (%s, $method,$nseq);\n \
   $method=$mode;\n    &set_temporary_dir (\"set\"\
,$infile,\"infile\");\n    %s=read_fasta_seq(\"inf\
ile\");\n\n    foreach $seq (keys(%s))\n      {\n	\
$slist[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s\
{$seq}{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{or\
der}]=length ($s{$seq}{seq});\n      }\n    $nseq=\
$#sname+1;\n    open (F, \">outfile\");\n    print\
 F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\
\\n\";\n    for ($a=0; $a<$nseq;$a++)\n      {\n	p\
rint F \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n  \
    }\n    close (F);\n    &safe_system(\"t_coffee\
 -other_pg seq_reformat -in infile -output tblastx\
_db1 > tblastxdb\");\n    &safe_system (\"formatdb\
 -i tblastxdb -p T\");\n    &safe_system (\"blasta\
ll -p blastp -i tblastxdb -d tblastxdb -m7 >blast.\
output\");\n    ncbi_tblastpx_xml2lib_file (\"outf\
ile\", file2string (\"blast.output\"), %s);\n    &\
set_temporary_dir (\"unset\",$mode, $method, \"out\
file\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  \
  }\n\nsub x3dna_find_pair2lib\n      {\n      my \
($seq, $pdb, $lib, $pg)=@_;\n      my $outfile1=\"\
dssr-2ndstrs.dbn\";\n      my $outfile2=\"simple.o\
utput\";\n      my $f= new FileHandle;\n      my (\
$rnaSS,$pdbSS);\n      my $command;\n      my %s_p\
db;\n      my %s_seq;\n      my $pdbseq=vtmpnam(NU\
LL);\n      \n      #$pg: \"find_pair\" OR \"find_\
pair -p\"\n      \n      if (!pg_is_installed (\"f\
ind_pair\"))\n	{\n	  add_warning ($$,$$, \"x3dna/f\
ind_pair could not be used to extract RNA secondar\
y structures. Secondary structures will be extract\
ed by x3dna-ssr instead -- Install the find-pair m\
odule of x3dna  [http://x3dna.org/]\");\n	  return\
 x3dnassr2lib ($seq, $pdb, $lib);\n	}\n      \n\n \
     #get PDB sequence\n      safe_system (\"t_cof\
fee -other_pg extract_from_pdb $pdb -seq >$pdbseq\\
");\n      \n      #get find_pair contacts\n      \
$command=\"$pg $pdb simple.output > /dev/null 2>/d\
ev/null\";\n      safe_system ($command);\n\n     \
 if (($command=~/find_pair -p/)){$outfile2=\"allpa\
irs.ana\";}\n      else {$outfile2=\"simple.output\
\";}\n      \n      if ( !-e $outfile2)\n	{\n	  my\
exit(flush_error(\"x3dna failed to compute the sec\
ondary structure file $outfile2 for $pdb\"));\n	  \
myexit ($EXIT_FAILURE);\n	}\n      \n\n      #Hand\
le situations when the pdb sequence differs from t\
he RNA sequence\n      #my @out=file2array($outfil\
e1);\n      %s_pdb=read_fasta_seq_index ($pdbseq);\
\n      %s_seq=read_fasta_seq_index ($seq);\n     \
 my $rnaS=uc($s_seq{0}{seq});\n      my $pdbS=uc($\
s_pdb{0}{seq});\n      my $vienna;\n      my @lu;\\
n    \n      if ($rnaS ne $pdbS)\n	{\n	  \n	  my (\
$rna,$pdb);\n	  $rnaSS=$rnaS;\n	  $pdbSS=$pdbS;\n	\
  $rnaSS=~s/T/U/g;\n	  $pdbSS=~s/T/U/g;\n	  ($rnaS\
S,$pdbSS)=nw ($rnaS, $pdbS);\n	  \n	  my @rnaA =sp\
lit (//,$rnaSS);\n	  my @pdbA=split (//,$pdbSS);\n\
	  my $l=@rnaA;\n	  \n	  #print \"\\n--- $s_seq{0}\
{name} $rnaSS\\n--- $s_seq{0}{name} $pdbSS\\n\\n\"\
;\n	  \n	  for (my $b=0,my $a=0; $a<$l; $a++)\n	  \
  {\n	      if   ($rnaA[$a] ne '-' && $pdbA[$a] ne\
 '-'){$lu[++$pdb]=++$rna;}\n	      elsif($rnaA[$a]\
 eq '-'){$lu[++$pdb]=-1;}\n	      elsif($pdbA[$a] \
eq '-'){++$rna;}\n	    }\n	}\n      else\n	{\n	  f\
or (my $a=0; $a<=length ($rnaS); $a++)\n	    {\n	 \
     $lu[$a]=$a;\n	    }\n	}\n      my $l=length (\
$rnaS);\n      open ($f, \">$lib\");\n      print \
$f \"! TC_LIB_FORMAT_01\\n\";\n      print $f \"1\\
\n\";\n      print $f \"$s_seq{0}{name} $l $rnaS\\\
n\";\n      print $f \"!CMT: [SOURCE] >$s_seq{0}{n\
ame} 3D contact library Generated by $pg (x3dna)\\\
n\";\n      print $f \"#1 1\\n\";\n      \n      m\
y $ne;\n      my @array=file2array($outfile2);\n  \
    for (my $a=0; $a<5; $a++){shift (@array);}\n  \
    while (!($array[0]=~/####/))\n	{\n	  my $line=\
 shift (@array);\n	  my @l=($line=~/(\\d+)/g);\n	 \
 \n	 \n	  my $f1=$lu[$l[0]];\n	  my $s1=$lu[$l[1]]\
;\n\n	  #print \"\\n$line\\n$l[0] --> $f1\\n$l[1] \
--> $s1\\n\\n\"; \n	  \n	  if (!$f1 || !$s1)\n	   \
 {\n	      print \"\\n1---- $rnaSS\\n2---- $pdbSS\\
\n$line\\n[$l[0] --- $l[1]]<---->[$f1 --- $s1]\\n\\
";\n	      myexit(flush_error(\"Error while parsin\
g s3dna::find_pair output\"));\n	    }\n	  elsif (\
$f1==-1 || $s1==-1){;}\n	  elsif ($f1<$s1){print $\
f \"$f1 $s1 100\\n\";}\n	  else {print $f \"$s1 $f\
1 100\\n\";$ne++;}\n	}\n      print $f \"! SEQ_1_T\
O_N\\n\";\n      close ($f);\n      return;\n    }\
\nsub RNAplfold2lib\n  {\n    my ($seq, $lib)=@_;\\
n    my $f= new FileHandle;\n    \n    &safe_syste\
m (\"t_coffee -other_pg RNAplfold2tclib.pl -in=$se\
q -out=$lib\");\n    \n    if ( !-e $lib)\n	{\n	 m\
yexit(flush_error(\"RNAplfold failed to compute th\
e secondary structure of $s{$seq}{name}\"));\n	 my\
exit ($EXIT_FAILURE);\n       }\n    open ($f, \">\
>$lib\");\n    print $f \"!CMT: [SOURCE] 2D contac\
t library Generated by RNAPlfold (Vienna Package)\\
\n\";\n    close $f;\n    return;\n  }\nsub x3dnas\
sr2lib\n    {\n      my ($seq, $pdb, $lib)=@_;\n  \
    my $outfile=\"dssr-2ndstrs.dbn\";\n      my $f\
= new FileHandle;\n      \n\n      if (!pg_is_inst\
alled (\"x3dna-ssr\"))\n	{\n	  add_warning ($$,$$,\
 \"x3dna-ssr could not be used to extract RNA seco\
ndary structures. Secondary structures will be pre\
dicted ab-initio instead with RNAPlfold -- Install\
 s3dna [http://x3dna.org/] \");\n	  return RNAplfo\
ld2lib ($seq,$lib);\n	}\n      \n      safe_system\
 (\"x3dna-ssr -i=$pdb >/dev/null 2>/dev/null\");\n\
      if ( !-e $outfile)\n	{\n	  myexit(flush_erro\
r(\"x3dna-ssr failed to compute the secondary stru\
cture file \"));\n	  myexit ($EXIT_FAILURE);\n	}\n\
\n      #Handle situations when the pdb sequence d\
iffers from the RNA sequence\n      @out=file2arra\
y($outfile);\n      my %s=read_fasta_seq ($seq);\n\
      my @names=keys (%s);\n      my $rnaS=uc($s{$\
names[0]}{seq});\n      my $pdbS=uc($out[1]);\n   \
   my $vienna;\n      \n      #x3dna returns non l\
egitimate nucleotides\n       $pdbS=~s/[^AGCTU]//g\
;\n      \n      if ($rnaS ne $pdbS)\n	{\n	  my ($\
rna,$pdb);\n	  my $rnaSS=$rnaS;\n	  my $pdbSS=$pdb\
S;\n	  $rnaSS=~s/T/U/g;\n	  $pdbSS=~s/T/U/g;\n	  (\
$rnaSS,$pdbSS)=nw ($rnaSS, $pdbSS);\n	  my @rnaA =\
split (//,$rnaSS);\n	  my @pdbA=split (//,$pdbSS);\
\n	  my @SS=split (//, $out[2]);\n	  \n	  my $l=@r\
naA;\n	  for (my $b=0,my $a=0; $a<$l; $a++)\n	    \
{\n	      if   ($rnaA[$a] ne '-' && $pdbA[$a] ne '\
-'){$vienna.=$SS[$b++];}\n	      elsif($rnaA[$a] e\
q '-'){$b++;}\n	      elsif($pdbA[$a] eq '-'){$vie\
nna.='.';}\n	    }\n	}\n      else\n	{\n	  $vienna\
=$out[2];\n	}\n    \n\n      open ($f, \">seq\");\\
n      print $f \">$names[0]\\n$rnaS\\n\";\n      \
close $f;\n      \n      open ($f, \">str\");\n   \
   print $f \">$names[0]\\n$vienna\\n\";\n      cl\
ose $f;\n      \n      safe_system (\"t_coffee -ot\
her_pg seq_reformat -in seq -in2 str -output vienn\
a2tc_lib >$lib\");\n      if ( !-e $lib)\n	    {\n\
	      myexit(flush_error(\"seq_reformat failed to\
 convert your secondary structure\"));\n	      mye\
xit ($EXIT_FAILURE);\n	    }\n      \n      open (\
$f, \">>$lib\");\n      print $f \"!CMT: [SOURCE] \
>$names[0] 2D Contact library generated by x3dna-s\
sr\\n\";\n      #print $f \"! Vienna_Format: >$nam\
es[0]\\n\";\n      #print $f \"! Vienna_Format: $v\
ienna\\n\";\n      \n      close $f;\n      return\
;\n    }\n\n\nsub file2head\n      {\n	my $file = \
shift;\n	my $size = shift;\n	my $f= new FileHandle\
;\n	my $line;\n	open ($f,$file);\n	read ($f,$line,\
 $size);\n	close ($f);\n	return $line;\n      }\ns\
ub file2tail\n      {\n	my $file = shift;\n	my $si\
ze = shift;\n	my $f= new FileHandle;\n	my $line;\n\
\n	open ($f,$file);\n	seek ($f,$size*-1, 2);\n	rea\
d ($f,$line, $size);\n	close ($f);\n	return $line;\
\n      }\n\n\nsub vtmpnam\n      {\n	my $r=rand(1\
00000);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	\
  {\n	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIS\
T, $f);\n	return $f;\n      }\n\nsub myexit\n  {\n\
    my $code=@_[0];\n    if ($CLEAN_EXIT_STARTED==\
1){return;}\n    else {$CLEAN_EXIT_STARTED=1;}\n  \
  ### ONLY BARE EXIT\n    exit ($code);\n  }\nsub \
set_error_lock\n    {\n      my $name = shift;\n  \
    my $pid=$$;\n\n\n      &lock4tc ($$,\"LERROR\"\
, \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\n\
      return;\n    }\nsub set_lock\n  {\n    my $p\
id=shift;\n    my $msg= shift;\n    my $p=getppid(\
);\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$m\
sg\\n\");\n  }\nsub unset_lock\n   {\n\n    my $pi\
d=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRELEASE\\
",\"\");\n  }\nsub shift_lock\n  {\n    my $from=s\
hift;\n    my $to=shift;\n    my $from_type=shift;\
\n    my $to_type=shift;\n    my $action=shift;\n \
   my $msg;\n\n    if (!&lock4tc($from, $from_type\
, \"LCHECK\", \"\")){return 0;}\n    $msg=&lock4tc\
 ($from, $from_type, \"LREAD\", \"\");\n    &lock4\
tc ($from, $from_type,\"LRELEASE\", $msg);\n    &l\
ock4tc ($to, $to_type, $action, $msg);\n    return\
;\n  }\nsub isshellpid\n  {\n    my $p=shift;\n   \
 if (!lock4tc ($p, \"LLOCK\", \"LCHECK\")){return \
0;}\n    else\n      {\n	my $c=lock4tc($p, \"LLOCK\
\", \"LREAD\");\n	if ( $c=~/-SHELL-/){return 1;}\n\
      }\n    return 0;\n  }\nsub isrootpid\n  {\n \
   if(lock4tc (getppid(), \"LLOCK\", \"LCHECK\")){\
return 0;}\n    else {return 1;}\n  }\nsub lock4tc\
\n	{\n	  my ($pid,$type,$action,$value)=@_;\n	  my\
 $fname;\n	  my $host=hostname;\n\n	  if ($type eq\
 \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$host.lock4tco\
ffee\";}\n	  elsif ( $type eq \"LERROR\"){ $fname=\
\"$LOCKDIR/.$pid.$host.error4tcoffee\";}\n	  elsif\
 ( $type eq \"LWARNING\"){ $fname=\"$LOCKDIR/.$pid\
.$host.warning4tcoffee\";}\n\n	  if ($debug_lock)\\
n	    {\n	      print STDERR \"\\n\\t---lock4tc(tc\
g): $action => $fname =>$value (RD: $LOCKDIR)\\n\"\
;\n	    }\n\n	  if    ($action eq \"LCHECK\") {ret\
urn -e $fname;}\n	  elsif ($action eq \"LREAD\"){r\
eturn file2string($fname);}\n	  elsif ($action eq \
\"LSET\") {return string2file ($value, $fname, \">\
>\");}\n	  elsif ($action eq \"LRESET\") {return s\
tring2file ($value, $fname, \">\");}\n	  elsif ($a\
ction eq \"LRELEASE\")\n	    {\n	      if ( $debug\
_lock)\n		{\n		  my $g=new FileHandle;\n		  open (\
$g, \">>$fname\");\n		  print $g \"\\nDestroyed by\
 $$\\n\";\n		  close ($g);\n		  safe_system (\"mv \
$fname $fname.old\");\n		}\n	      else\n		{\n		  \
unlink ($fname);\n		}\n	    }\n	  return \"\";\n	}\
\n\nsub file2string\n	{\n	  my $file=@_[0];\n	  my\
 $f=new FileHandle;\n	  my $r;\n	  open ($f, \"$fi\
le\");\n	  while (<$f>){$r.=$_;}\n	  close ($f);\n\
	  return $r;\n	}\nsub file2array\n	{\n	  my $file\
=@_[0];\n	  my $f=new FileHandle;\n	  my @r;\n	  o\
pen ($f, \"$file\");\n	  while (<$f>){@r=(@r,$_);}\
\n	  close ($f);\n	  return @r;\n	}\nsub string2fi\
le\n    {\n    my ($s,$file,$mode)=@_;\n    my $f=\
new FileHandle;\n\n    open ($f, \"$mode$file\");\\
n    print $f  \"$s\";\n    close ($f);\n  }\n\nBE\
GIN\n    {\n      srand;\n\n      $SIG{'SIGUP'}='s\
ignal_cleanup';\n      $SIG{'SIGINT'}='signal_clea\
nup';\n      $SIG{'SIGQUIT'}='signal_cleanup';\n  \
    $SIG{'SIGILL'}='signal_cleanup';\n      $SIG{'\
SIGTRAP'}='signal_cleanup';\n      $SIG{'SIGABRT'}\
='signal_cleanup';\n      $SIG{'SIGEMT'}='signal_c\
leanup';\n      $SIG{'SIGFPE'}='signal_cleanup';\n\
\n      $SIG{'SIGKILL'}='signal_cleanup';\n      $\
SIG{'SIGPIPE'}='signal_cleanup';\n      $SIG{'SIGS\
TOP'}='signal_cleanup';\n      $SIG{'SIGTTIN'}='si\
gnal_cleanup';\n      $SIG{'SIGXFSZ'}='signal_clea\
nup';\n      $SIG{'SIGINFO'}='signal_cleanup';\n\n\
      $SIG{'SIGBUS'}='signal_cleanup';\n      $SIG\
{'SIGALRM'}='signal_cleanup';\n      $SIG{'SIGTSTP\
'}='signal_cleanup';\n      $SIG{'SIGTTOU'}='signa\
l_cleanup';\n      $SIG{'SIGVTALRM'}='signal_clean\
up';\n      $SIG{'SIGUSR1'}='signal_cleanup';\n\n\\
n      $SIG{'SIGSEGV'}='signal_cleanup';\n      $S\
IG{'SIGTERM'}='signal_cleanup';\n      $SIG{'SIGCO\
NT'}='signal_cleanup';\n      $SIG{'SIGIO'}='signa\
l_cleanup';\n      $SIG{'SIGPROF'}='signal_cleanup\
';\n      $SIG{'SIGUSR2'}='signal_cleanup';\n\n   \
   $SIG{'SIGSYS'}='signal_cleanup';\n      $SIG{'S\
IGURG'}='signal_cleanup';\n      $SIG{'SIGCHLD'}='\
signal_cleanup';\n      $SIG{'SIGXCPU'}='signal_cl\
eanup';\n      $SIG{'SIGWINCH'}='signal_cleanup';\\
n\n      $SIG{'INT'}='signal_cleanup';\n      $SIG\
{'TERM'}='signal_cleanup';\n      $SIG{'KILL'}='si\
gnal_cleanup';\n      $SIG{'QUIT'}='signal_cleanup\
';\n\n      our $debug_lock=$ENV{\"DEBUG_LOCK\"};\\
n\n\n\n\n      foreach my $a (@ARGV){$CL.=\" $a\";\
}\n      if ( $debug_lock ){print STDERR \"\\n\\n\\
\n********** START PG: $PROGRAM *************\\n\"\
;}\n      if ( $debug_lock ){print STDERR \"\\n\\n\
\\n**********(tcg) LOCKDIR: $LOCKDIR $$ **********\
***\\n\";}\n      if ( $debug_lock ){print STDERR \
\"\\n --- $$ -- $CL\\n\";}\n\n\n\n\n    }\nsub flu\
sh_error\n  {\n    my $msg=shift;\n    $msg.=\" [t\
c_generic_method.pl/FATAL]\";\n    return add_erro\
r ($EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);\n  \
}\nsub add_error\n  {\n    my $code=shift;\n    my\
 $rpid=shift;\n    my $pid=shift;\n    my $ppid=sh\
ift;\n    my $type=shift;\n    my $com=shift;\n\n \
   $ERROR_DONE=1;\n    lock4tc ($rpid, \"LERROR\",\
\"LSET\",\"$pid -- ERROR: $type\\n\");\n    lock4t\
c ($$, \"LERROR\",\"LSET\", \"$pid -- COM: $com\\n\
\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pid\
 -- STACK: $ppid -> $pid\\n\");\n\n    return $cod\
e;\n  }\nsub add_warning\n  {\n    my $rpid=shift;\
\n    my $pid =shift;\n    my $command=shift;\n   \
 my $msg=\"$$ -- WARNING: $command\\n\";\n    prin\
t STDERR \"$msg\";\n    lock4tc ($$, \"LWARNING\",\
 \"LSET\", $msg);\n  }\n\nsub signal_cleanup\n  {\\
n    print dtderr \"\\n**** $$ (tcg) was killed\\n\
\";\n    &cleanup;\n    exit ($EXIT_FAILURE);\n  }\
\nsub clean_dir\n  {\n    my $dir=@_[0];\n    if (\
 !-d $dir){return ;}\n    elsif (!($dir=~/tmp/)){r\
eturn ;}#safety check 1\n    elsif (($dir=~/\\*/))\
{return ;}#safety check 2\n    else\n      {\n	`rm\
 -rf $dir`;\n      }\n    return;\n  }\nsub cleanu\
p\n  {\n    #print stderr \"\\n----tc: $$ Kills $P\
IDCHILD\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n   \
 my $p=getppid();\n    $CLEAN_EXIT_STARTED=1;\n\n\\
n\n    if (&lock4tc($$,\"LERROR\", \"LCHECK\", \"\\
"))\n      {\n	my $ppid=getppid();\n	if (!$ERROR_D\
ONE)\n	  {\n	    &lock4tc($$,\"LERROR\", \"LSET\",\
 \"$$ -- STACK: $p -> $$\\n\");\n	    &lock4tc($$,\
\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\");\n	  \
}\n      }\n    my $warning=&lock4tc($$, \"LWARNIN\
G\", \"LREAD\", \"\");\n    my $error=&lock4tc($$,\
  \"LERROR\", \"LREAD\", \"\");\n    #release erro\
r and warning lock if root\n\n    if (isrootpid() \
&& ($warning || $error) )\n      {\n\n	print STDER\
R \"**************** Summary *************\\n$erro\
r\\n$warning\\n\";\n\n	&lock4tc($$,\"LERROR\",\"RE\
LEASE\",\"\");\n	&lock4tc($$,\"LWARNING\",\"RELEAS\
E\",\"\");\n      }\n\n\n    foreach my $f (@TMPFI\
LE_LIST)\n      {\n	if (-e $f){unlink ($f);}\n    \
  }\n    foreach my $d (@TMPDIR_LIST)\n      {\n	c\
lean_dir ($d);\n      }\n    #No More Lock Release\
\n    #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\"\"); #\
release lock\n\n    if ( $debug_lock ){print STDER\
R \"\\n\\n\\n********** END PG: $PROGRAM ($$) ****\
*********\\n\";}\n    if ( $debug_lock ){print STD\
ERR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $\
$ *************\\n\";}\n  }\nEND\n  {\n\n    &clea\
nup();\n  }\n\nsub blast_com2new_blast_com\n    {\\
n      my $com=shift;\n      if ($com=~/formatdb/)\
\n	{\n	  $com=~s/formatdb/makeblastdb/;\n	  $com=~\
s/\\-i/\\-in/;\n	  if ($com =~/pF/){$com=~s/\\-pF/\
\\-dbtype nucl/;}\n	  if ($com =~/p F/){$com=~s/\\\
-p F/\\-dbtype nucl/;}\n	  $com=\"$com -logfile /d\
ev/null\";\n	  return $com;\n	}\n      else {retur\
n $com;}\n      \n    }\nsub safe_system\n{\n  my \
$com=shift;\n  my $ntry=shift;\n  my $ctry=shift;\\
n  my $pid;\n  my $status;\n  my $ppid=getppid();\\
n  if ($com eq \"\"){return 1;}\n \n  if ( ($com=~\
/^blast/) ||($com=~/^formatdb/)){$com=&blast_com2n\
ew_blast_com($com);}\n \n  if (($pid = fork ()) < \
0){return (-1);}\n  if ($pid == 0)\n    {\n      s\
et_lock($$, \" -SHELL- $com (tcg)\");\n      if( $\
debug_generic_method ) { printf \"~ exec: %s\\n\",\
 $com; }\n      exec ($com);\n      if( $debug_gen\
eric_method ) { printf \"~ exitcode: %s\\n\", $?; \
}\n    }\n  else\n    {\n      lock4tc ($$, \"LLOC\
K\", \"LSET\", \"$pid\\n\");#update parent\n      \
$PIDCHILD=$pid;\n    }\n  if ($debug_lock){printf \
STDERR \"\\n\\t .... safe_system (fasta_seq2hmm)  \
p: $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,\
WTERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\
\"LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAIL\
URE || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\")\
)\n    {\n      if ($ntry && $ctry <$ntry)\n	{\n\n\
	  add_warning ($$,$$,\"$com failed [retry: $ctry \
out of $ntry]\");\n	  lock4tc ($pid, \"LRELEASE\",\
 \"LERROR\", \"\");\n	  #if ($com=~/EBI/){$com=~s/\
EBI/NCBI/;}\n	  #elsif ($com=~/NCBI/){$com=~s/NCBI\
/EBI/;}\n\n	  return safe_system ($com, $ntry, ++$\
ctry);\n	}\n      elsif ($ntry == -1)\n	{\n	  if (\
!shift_lock ($pid, $$, \"LERROR\", \"LWARNING\", \\
"LSET\"))\n	    {\n	      add_warning ($$,$$,\"$co\
m failed\");\n	    }\n	  else\n	    {\n	      lock\
4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	    \
}\n	  return $?;}\n      else\n	{\n	  if (!shift_l\
ock ($pid,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	\
    {\n	      myexit(add_error ($EXIT_FAILURE,$$,$\
pid,getppid(), \"UNSPECIFIED system\", $com));\n	 \
   }\n	}\n    }\n  return $?;\n}\n\nsub check_conf\
iguration\n    {\n      my @l=@_;\n      my $v;\n \
     foreach my $p (@l)\n	{\n\n	  if   ( $p eq \"E\
MAIL\")\n	    {\n	      if ( !($EMAIL=~/@/))\n		{\\
n		add_warning($$,$$,\"Could Not Use EMAIL\");\n		\
myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"\
EMAIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p\
 eq \"INTERNET\")\n	    {\n	      if ( !&check_int\
ernet_connection())\n		{\n		  myexit(add_error ($E\
XIT_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"))\
;\n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\\
n	      if (!&pg_is_installed (\"wget\") && !&pg_i\
s_installed (\"curl\"))\n		{\n		  myexit(add_error\
 ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED\
:wget\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_i\
s_installed ($p)))\n	    {\n	      myexit(add_erro\
r ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLE\
D:$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n  \
  }\nsub nw\n      {\n	my($A,$B)=@_;\n	my ($i,$j, \
$s);\n	my $gep=-2;\n	my $match=+2;\n	my $mismatch=\
0;\n	my ($sub, $ins, $del);\n\n\n	if ($A eq $B){re\
turn ($A,$B);}\n	\n	$A=~s/[\\s\\d]//g;	\n	$B=~s/[\\
\s\\d]//g;	\n\n\n	my @rA=split ('',$A);\n	my @rB=s\
plit ('',$B);\n	\n	#evaluate substitutions\n	my $l\
enA=@rA;\n	my $lenB=@rB;\n	\n	for ($i=0; $i<=$lenA\
; $i++){$smat[$i][0]=$i*$gep;$tb[$i][0 ]= 1;}\n	fo\
r ($j=0; $j<=$lenB; $j++){$smat[0][$j]=$j*$gep;$tb\
[0 ][$j]=-1;}\n	\n	for ($i=1; $i<=$lenA; $i++)\n	 \
 {\n	    for ($j=1; $j<=$lenB; $j++)\n	      {\n		\
if ($rA[$i-1] eq $rB[$j-1]){$s=$match;}\n		else {$\
s=$mismatch;}\n		\n		$sub=$smat[$i-1][$j-1]+$s;\n	\
	$del=$smat[$i  ][$j-1]+$gep;\n		$ins=$smat[$i-1][\
$j  ]+$gep;\n		\n		if   ($sub>=$del && $sub>=$ins)\
{$smat[$i][$j]=$sub;$tb[$i][$j]=0;}\n		elsif($del>\
$ins){$smat[$i][$j]=$del;$tb[$i][$j]=-1;}\n		else \
{$smat[$i][$j]=$ins;$tb[$i][$j]=1;}\n		}\n	  }\n	#\
print \"\\n---- SCORE=$smat[$lenA][$lenB]\\n\";\n	\
\n	$i=$lenA;\n	$j=$lenB;\n	my $aln_len=0;\n\n	whil\
e (!($i==0 && $j==0))\n	  {\n	    if ($tb[$i][$j]=\
=0)\n	    {\n	      $aA[$aln_len]=$rA[--$i];\n	   \
   $aB[$aln_len]=$rB[--$j];\n	    }\n	  elsif ($tb\
[$i][$j]==-1)\n	    {\n	      $aA[$aln_len]='-';\n\
	      $aB[$aln_len]=$rB[--$j];\n	    }\n	  elsif \
($tb[$i][$j]==1)\n	    {\n	      $aA[$aln_len]=$rA\
[--$i];\n	      $aB[$aln_len]='-';\n	    }\n	  $al\
n_len++;\n	  }\n	\n	\n	@aA=reverse (@aA);\n	@aB=re\
verse (@aB);\n	my $sA=join('',@aA);\n	my $sB=join(\
'',@aB);\n	return ($sA,$sB);\n      }\n      \nsub\
 fasta2nseq\n	{\n	  \n	  my $f=@_[0];\n	  my $nseq\
;\n\n	  open (F, \"$f\") or return 0;\n	  while (<\
F>)\n	    {\n	      if ($_=~/\\>/){$nseq++;}\n	   \
 }\n	  close (F);\n	  return $nseq;\n	}\n	\nsub co\
mpress\n	  {\n	    my ($f, $mode)=@_;\n	    if    \
($mode eq \"gz\"){system (\"gzip $f\");}\n	    els\
if ($mode eq \"zip\" ){system (\"zip $f\");}\n	   \
 return;\n	  }\nsub uncompress \n	  {\n	    my $f=\
@_[0];\n	    if ( -e $f) {return \"\";}\n\n	    my\
 $gz=$f.\".gz\";\n	    if ( -e $gz)\n	      {\n		s\
ystem (\"gunzip $gz\");\n		return ($f, \"gz\");\n	\
      }\n	    my $gz=$f.\".zip\";\n	    \n	    if \
( -e $gz)\n	      {\n		system (\"unzip $gz\");\n		\
return ($f, \"zip\");\n	      }\n	    return \"\";\
\n	  }\nmy $program=\"T-COFFEE (Version_13.46.0.91\
9e8c6b)\";\n","use Env;\nuse strict;\nuse FileHand\
le;\nuse DirHandle;\nuse Cwd;\nuse File::Path;\nus\
e Sys::Hostname;\nuse File::Temp qw/ tempfile temp\
dir /;\n\nmy $QUIET=\">/dev/null 2>/dev/null\";\nm\
y $VERBOSE=$ENV{VERBOSE_4_DYNAMIC};\nmy $FATAL=\"[\
FATAL:mmseqs2prf]\";\nour $EXIT_FAILURE=1;\nour $E\
XIT_SUCCESS=0;\nour $LAST_COM=\"\";\n\nmy $tmpdir \
= File::Temp->newdir();\n\nmy $doquiet=0;\nmy ($ou\
tdir,$cachedb,$cacheq);\n\nmy $ff=new FileHandle;\\
nmy (%db, %q, %P, %H, %T);\nmy ($dbf, $qf, $out);\\
nmy $prot_min_cov=50;\nmy $prot_min_sim=0;\nmy $pr\
ot_max_sim=100;\nmy $psitrim_mode=\"sorttrim\";\n\\
nmy $psitrim_tree=\"codnd\";\nmy $psitrim=100;\nmy\
 $psiJ=1;\nmy $TF;\nmy $S;\n\nmy $updatedb=0;\nmy \
$updateq=0;\nmy $update=0;\nmy $qff;\nmy %CIRCULAR\
;\nmy $mmseqsR;\nmy %R;\nmy $split=1000000;\n\nfor\
 ($a=0; $a<=$#ARGV; $a++)\n  {\n    if    ($ARGV[$\
a] eq \"-protein_db\" || $ARGV[$a] eq \"-db\"){$db\
f=file2abs($ARGV[++$a]);}\n    elsif ($ARGV[$a] eq\
 \"-q\" || $ARGV[$a] eq \"-i\") {$qff =file2abs($A\
RGV[++$a]);}\n    elsif ($ARGV[$a] eq \"-update\")\
{$update=1;}\n    elsif ($ARGV[$a] eq \"-quiet\") \
{$doquiet=1;}\n    \n    elsif ($ARGV[$a] eq \"-od\
ir\") {$outdir=file2abs($ARGV[++$a]);}\n    elsif \
($ARGV[$a] eq \"-o\")    {$mmseqsR=file2abs($ARGV[\
++$a]);}\n    \n    elsif ($ARGV[$a] eq \"-templat\
e_file\" || $ARGV[$a] eq \"-tf\") {$TF=($ARGV[++$a\
]);}\n    elsif ($ARGV[$a] eq \"-cachedb\") {$cach\
edb=file2abs($ARGV[++$a]);}\n    elsif ($ARGV[$a] \
eq \"-updatedb\") {$updatedb=1;}\n    \n    \n    \
elsif ($ARGV[$a] eq \"-cacheq\")  {$cacheq=file2ab\
s($ARGV[++$a]);}\n    elsif ($ARGV[$a] eq \"-updat\
eq\") {$updateq=1;}\n    \n    elsif ($ARGV[$a] eq\
 \"-prot_min_sim\") {$prot_min_sim=($ARGV[++$a]);}\
\n    elsif ($ARGV[$a] eq \"-prot_max_sim\") {$pro\
t_max_sim=($ARGV[++$a]);}\n    elsif ($ARGV[$a] eq\
 \"-prot_min_cov\") {$prot_min_cov=($ARGV[++$a]);}\
\n    \n    elsif ($ARGV[$a] eq \"-psitrim_mode\")\
 {$psitrim_mode=$ARGV[++$a];}\n    elsif ($ARGV[$a\
] eq \"-psiJ\")         {$psiJ=$ARGV[++$a];}\n    \
elsif ($ARGV[$a] eq \"-psitrim_tree\") {$psitrim_t\
ree=$ARGV[++$a];}\n    elsif ($ARGV[$a] eq \"-psit\
rim\")      {$psitrim=$ARGV[++$a];}\n    elsif ($A\
RGV[$a] eq \"-s\")            {$S=\"-s \".$ARGV[++\
$a];}\n    elsif ($ARGV[$a] eq \"-split\")        \
{$split=$ARGV[++$a];}\n    else {die \"$ARGV[$a] i\
s an unknown argument $FATAL\";}\n  }\n\nif (!$dbf\
){print STDERR \"ERROR: mmseqs2prf required a data\
base via -protein_db $FATAL\\n\";exit ($EXIT_FAILU\
RE);}\nif (!$qff){print STDERR \"ERROR: mmseqs2prf\
 required a query    -q $FATAL\";exit ($EXIT_FAILU\
RE);}\nif (!is_installed(\"mmseqs\")){print STDERR\
 \"ERROR: mmseqs2prf required mmseqs to be install\
ed $FATAL\\n\";exit ($EXIT_FAILURE);}\n\n\nif   ($\
cachedb eq \"TMP\"){$cachedb=$tmpdir;}\nelsif(!$ca\
chedb){$cachedb=file2path($dbf);}\n\nif   (!$mmseq\
sR){$mmseqsR=\"./out.mmseqs\";}\nif   (!$cacheq ||\
 $cacheq eq \"TMP\"){$cacheq=\"$tmpdir/query\";}\n\
if   (!$outdir  ){$outdir=\"./R_dir\";}\nif   (!$d\
oquiet) {$QUIET=\"\";}\n\nmymkdir ($outdir,$cached\
b,$cacheq,$tmpdir);\nmy ($qf,%H)=filelist2h($outdi\
r,string2fasta_list($qff));\nmake_output_structure\
($outdir,%H);\n\n\nif ( ! -e $mmseqsR || $update){\
split_mmseqs($qf, $cacheq, $updateq, $dbf, $cached\
b, $updatedb,$S,$mmseqsR,$split);}\n\nmmseqs2prf (\
$mmseqsR,$outdir,$prot_min_sim,$prot_max_sim, $pro\
t_min_cov,%H);\nprf2trimprf($outdir,$psitrim_mode,\
 $psitrim_tree, $psitrim, %H);\n\nif ($TF){h2templ\
ate_file ($TF,%H);}\n\n\nsub h2template_file\n  {\\
n    my ($tf,%h)=@_;\n    my $fh =new FileHandle;\\
n    my %lu=my %lu=h2lu(%h);\n\n    open ($fh, \">\
$tf\");\n    foreach my $s (keys{%lu})\n      {\n	\
my $f=$lu{$s}{0};\n	print $fh \">$s _R_ $h{$f}{tpr\
f}{$s}{absolute}\\n\";\n      }\n    close ($fh);\\
n    return $tf;\n  }\n    \nsub prf2trimprf\n  {\\
n    my ($outdir,$psitrim_mode,$psitrim_tree,$psit\
rim, %h)=@_;\n    my $template =new FileHandle;\n \
   my $template_file;\n    my $qf=abs2file($qf);\n\
    \n    my %lu=h2lu(%h);\n    foreach my $s (key\
s(%lu))\n      {\n	my $f=$lu{$s}{0};\n	system (\"t\
_coffee -other_pg seq_reformat -in  $h{$f}{prf}{$s\
}{absolute} -treemode=$psitrim_tree -keep 1 -actio\
n +$psitrim_mode $psitrim -output fasta_aln -out $\
h{$f}{tprf}{$s}{absolute}\");\n	foreach my $i (key\
s%{$lu{$s}})\n	  {\n	   \n	    my $f2=$lu{$s}{$i};\
\n	    if ( !-e $h{$f2}{tprf}{$s}{absolute})\n	   \
   {\n		system (\"cp $h{$f}{tprf}{$s}{absolute} $h\
{$f2}{tprf}{$s}{absolute}\");\n	      }\n	  }\n   \
   }\n    return %P;\n  }\n\n\n\n\nsub h2lu\n  {\n\
    my %h=@_;\n    my %lu;\n    my %count;\n    fo\
reach my $f (keys (%h))\n      {\n	foreach my $s (\
keys(%{$h{$f}{name}}))\n	  {\n	    \n	    $lu{$s}{\
$count{$s}++}=$f;\n	  }\n      }\n    return %lu;\\
n  }\nsub mmseqs2prf\n  {\n    #\"query[0],target[\
1],qaln[2],taln[3],qstart[4],qend[5],pident[6],qco\
v[7],qlen[8]\\\" $QUIET\");\n    my ($out,$outdir,\
$min_id, $max_id,$min_cov, %h)=@_;\n    my $ff  =n\
ew FileHandle;\n    my $prf =new FileHandle;\n    \
my $nn;\n    my $tot;\n    my $psn;\n    my %lu=h2\
lu(%h);\n    my %luf;\n    \n    open ($ff,$out);\\
n    while (<$ff>)\n      {\n	my $l=$_;\n	my @ll=s\
plit (/\\s/, $l);\n	my $sn=$ll[0];\n	my $f =$lu{$s\
n}{0};\n	my $cf=$h {$f}{prf}{$sn}{absolute};\n	\n	\
if ($sn ne $psn)\n	  {\n	    close $prf;\n	    if \
($luf{$cf}){open ($prf, \">>$cf\");}\n	    else\n	\
      {		\n		open ($prf, \">$cf\");\n		print $prf \
\">$sn\\n$h{$f}{seq}{$sn}\\n\";\n	      }\n	  }\n	\
$nn=++$luf{$cf};\n	$psn=$sn;\n	\n	my $id=$ll[6]*10\
0;\n	my $cov=$ll[7]*100;\n	my $len=$ll[8];\n	print\
 \"$id $max_id $min_id\\n\";\n	if ($id<=$max_id &&\
 $id>=$min_id && $cov>$min_cov)\n	  {\n	    print \
\"Keep\";\n	    print $prf \">$sn\\_$nn\\n\";\n	  \
  for (my $a=1; $a<$ll[4]; $a++){print $prf \"-\"}\
\n	    \n	    my @ql=split (//,$ll[2]);\n	    my @\
tl=split (//,$ll[3]);\n	    my $qlen=length($ll[2]\
);\n	    for (my $a=0; $a<$qlen; $a++)\n	      {\n\
		if ($ql[$a] ne \"-\"){print $prf \"$tl[$a]\";}\n\
	      }\n	    for (my $a=$ll[5]; $a<$len; $a++){p\
rint $prf \"-\"}\n	    print $prf \"\\n\";\n	  }\n\
      }\n    close($prf);\n    close($ff);\n\n    \
# checkout the un-used ones\n    foreach my $sn (k\
eys(%lu))\n      {\n	my $f=$lu{$sn}{0};\n\n	if (!-\
e $h{$f}{prf}{$sn}{absolute})\n	  {\n	    open ($p\
rf,\">$h{$f}{prf}{$sn}{absolute}\");\n	    print $\
prf \">$sn\\n$h{$f}{seq}{$sn}\\n\";\n	    close (p\
rf);\n	  }\n      }\n    #duplicate prf files that\
 are shared by different input datasets\n    forea\
ch my $sn (keys (%lu))\n      {\n	my $f0=$lu{$sn}{\
0};\n	\n	foreach my $i (keys(%{$lu{$sn}}))\n	  {	\\
n	    my $f=$lu{$sn}{$i};\n	    if (! -e $h{$f}{pr\
f}{$sn}{absolute}){system (\"cp $h{$f0}{prf}{$sn}{\
absolute} $h{$f}{prf}{$sn}{absolute}\");}\n	  }\n \
     }\n  }\n\n\nsub file2db\n  {\n    my ($in,$di\
r, $update)=(@_);\n    my %f;\n    my $out;\n    \\
n    if ( !-e $in)\n      {\n	print \"$in does not\
 exists $FATAL \\n\";\n      }\n    if ($dir)\n   \
   {\n	if (!-d $dir){mkdir ($dir) or die \"Could n\
ot create $dir $FATAL\"; }\n	$out=$dir.\"/\".abs2f\
ile($in);}\n    else {$out=$in;}\n    \n    \n    \
$f{name}=$in;\n    $f{db}=\"$out\\.MMSEQSDB\";\n  \
  $f{index}=\"$out\\.MMSEQSINDEX\";\n    \n    \n \
   #Trigger automated update when source db younge\
r than mmseqs file\n    if (-e $f{db} && ((-M $f{d\
b})>(-M $f{name}))){$update=1;}\n    \n    if (!-e\
 $f{db} || $update)\n      {\n	system (\"mmseqs cr\
eatedb  $f{name}  $f{db} $QUIET\");\n      }\n    \
\n    if (!-d $f{index} || $update)\n      {\n	sys\
tem (\"mmseqs createindex $f{name} $f{index} $QUIE\
T\");\n      }\n\n    return %f;\n    \n  }\n\nsub\
 file2path\n    {\n      my ($f)=@_;\n      $f=fil\
e2abs($f);\n      $f=~/(.*\\/)[^\\/]*$/;\n      my\
 $cdir=$1;\n      return $cdir;\n    }\nsub file2a\
bs\n     {\n       my ($f, $mode)=@_;\n       my $\
cdir=getcwd();\n       if ($f=~/^\\//){return $f;}\
\n       return \"$cdir/$f\";\n   }\nsub abs2file\\
n    {\n      my $in=shift @_;\n      my $out;\n  \
    \n      if ( $in=~/\\//)\n	  {\n	    $in=~/.*\\
\/([^\\/]*)$/;\n	    $out=$1;\n	  }\n	else\n	  {\n\
	    $out=$in;\n	  }\n      \n    return $out;\n  \
  }\n\nsub read_fasta_seq\n  {\n    my $f=@_[0];\n\
    my %hseq;\n    my (@seq, @com, @name);\n    my\
 ($a, $s,$nseq);\n    my $fh=new FileHandle;\n    \
\n    open ($fh, $f);\n    while (<$fh>)\n      {\\
n	$s.=$_;\n      }\n    close ($fh);\n\n\n    @nam\
e=($s=~/>(\\S*).*\\n[^>]*/g);\n\n    @seq =($s=~/>\
.*.*\\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([\
^>]*)/g);\n\n\n    $nseq=$#name+1;\n\n    for ($a=\
0; $a<$nseq; $a++)\n      {\n	my $s;\n	my $n=$name\
[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{cname}=cle\
an_file_name($n);\n	\n	$seq[$a]=~s/[^A-Za-z]//g;\n\
	$hseq{$n}{order}=$a;\n	$hseq{$n}{seq}=$seq[$a];\n\
	$hseq{$n}{com}=$com[$a];\n\n      }\n    return %\
hseq;\n  }\nsub mymkdir\n    {\n      my @l=@_;\n \
     foreach my $a (@l)\n	{\n	  if ( $a && !-d $a)\
\n	    {\n	      system (\"mkdir -p $a\");\n	     \
 if ( !-d $a)\n		{\n		  die \"Could not Create $a \
$FATAL\\n\";\n		}\n	    }\n	}\n\n      return 1;\n\
    }\nsub clean_file_name\n  {\n    my $name=@_[0\
];\n\n    $name=~s/[^A-Za-z1-9.-]/_/g;\n    return\
 $name;\n  }\n\nsub string2fasta_list\n    {\n    \
  my $string=@_[0];\n      if (!-f $string && !-d \
$string && !($string=~/\\*/)){return ();}\n      i\
f ($CIRCULAR{$string}){print STDERR \"ERROR: CIRCU\
LAR REFERENCE $string   $FATAL\\n\";exit ($EXIT_FA\
ILURE);}\n      $CIRCULAR{$string}=1;\n      \n\n \
     my @l1=string2list($string);\n      my @l2;\n\
      \n      foreach my $f (@l1)\n	{\n\n	  if (is\
fasta($f)){push (@l2, $f);}\n	  else \n	    {\n	  \
    foreach my $string2 (file2list($f))\n		{\n		  \
my @l3=string2fasta_list($string2);\n		  foreach m\
y $string3 (@l3)\n		    {\n		      push (@l2, $str\
ing3);\n		    }\n		}\n	    }\n	}\n      return shr\
inklist(@l2);\n    }\nsub shrinklist\n      {\n	my\
 @l=@_;\n	my @l2;\n\n	foreach my $e (@l)\n	  {\n	 \
   if ($e)\n	      {\n		print \"PUSH [$e]\\n\";\n	\
	push (@l2,$e);\n	      }\n	  }\n	return @l2;\n   \
   }\nsub string2list\n    {\n      my $string=@_[\
0];\n      my @list;\n\n      if    (-d $string   \
   ){@list= dir2list($string);}\n      elsif (   $\
string=~/\\*/){@list= glob    ($string);}\n      e\
lsif (-f $string      ){@list=         ($string);}\
\n      return @list;\n    }\nsub dir2list\n  {\n \
   my $dir=shift;\n    my @list;\n    my $cdir=get\
cwd;\n    my $dh  =new DirHandle;\n    \n    opend\
ir ($dh, $dir);\n    my @dlist=readdir($dh);\n    \
closedir($dh);\n    \n    foreach my $f (@dlist)\n\
      {\n	if ($f eq \".\" || $f eq \"..\"){;}\n	el\
se {push (@list, string2list(\"$dir/$f\"));}\n    \
  }\n    return @list;\n  }\nsub file2list\n{\n  m\
y $file=shift;\n  my @list;\n  my $fh=new FileHand\
le;\n  \n  open ($fh, \"$file\");\n  while (<$fh>)\
\n    {\n      my $l=$_;\n      chomp ($l);\n     \
 if ($l){push(@list, $l);}\n    }\n  close ($fh);\\
n  return @list;\n}\nsub file2string\n  {\n    my \
$f=@_[0];\n    my ($string, $l);\n    my $fh= new \
FileHandle;\n    open ($fh,\"$f\");\n    while (<$\
fh>)\n      {\n\n	$l=$_;\n	$string.=$l;\n      }\n\
    close ($fh);\n    $string=~s/\\r\\n/\\n/g;\n  \
  return $string;\n  }\nsub isfasta\n  {\n    my $\
file=shift;\n    my $fh=new FileHandle;\n    \n   \
 open ($fh, \"$file\");\n    while (<$fh>)\n      \
{\n	my $l=$_;\n	close ($fh);\n	if ($l=~/^>/){retur\
n 1;}\n	return 0;\n      }\n  }\n\n  \n\n\n\nsub f\
ilelist2h\n  {\n    my ($outdir,@list)=@_;\n    my\
 $infile=\"$outdir/fullseq.fa\";\n    my %h;\n    \
my %lu;\n    my $fh=new FileHandle;\n    \n    ope\
n ($fh, \">$infile\");\n    foreach my $f (@list)\\
n      {\n	if (!$QUIET){print \"! Process $f\\n\";\
}\n	my %s=read_fasta_seq($f);\n	$f=abs2file ($f);\\
n\n	$h{$f}{template_dir }=\"$outdir/$f.template_di\
r\";\n	$h{$f}{template_file}{relative}=\"$f.R.temp\
late_file\";\n	$h{$f}{template_file}{absolute}=\"$\
h{$f}{template_dir}/$h{$f}{template_file}{relative\
}\";\n		\n	foreach my $s (keys(%s))\n	  {\n	    my\
 $cname=$s{$s}{cname};\n	    my $name =$s{$s}{ nam\
e};\n	    \n	    my $seq  =$s{$s}{seq};\n	    my $\
com  =$s{$s}{com};\n	    \n	    #clean \n	    $seq\
=~s/\\-//g;\n	    $seq=~s/\\.//g;\n	    \n	    if \
($lu{$s}{seq} && $lu{$s}{seq} ne $seq)\n	      {\n\
		#print \"$lu{$s}{seq} ---> \\n$lu{$s}{source}/$s\
 \\n\";\n		#print \"$seq ---> \\n$f/$s \\n\";\n		\\
n		#die;\n		print \"ERROR: two different sequences\
 where provided for $name: [$seq/$f] and [$lu{$s}{\
seq}/$lu{$s}{source}]$FATAL]\\n\";\n		close (F);\n\
		die;\n	      }\n	    else \n	      {\n		$lu{$s}{\
seq   }=$seq;\n		$lu{$s}{source}=$f;\n		\n		print \
$fh \">$name\\n$seq\\n\";\n	      }\n	    $h{$f}{s\
eq  }{$name}=$seq;\n	    $h{$f}{name }{$name}=$nam\
e;\n	    $h{$f}{cname}{$name}=$cname;\n	    $h{$f}\
{prf}{$name}{relative}=\"$cname.R.prf\";\n	    $h{\
$f}{prf}{$name}{absolute}=\"$tmpdir/$h{$f}{prf}{$n\
ame}{relative}\";\n	    \n	    $h{$f}{tprf}{$name}\
{relative}=\"$cname.R.prf\";\n	    $h{$f}{tprf}{$n\
ame}{absolute}=\"$h{$f}{template_dir}/$h{$f}{tprf}\
{$name}{relative}\";\n	    \n	    $h{$f}{templates\
}.=\">$name _R_ $h{$f}{tprf}{$name}{relative}\\n\"\
;\n	  }\n      }\n    close ($fh);\n    return $in\
file, %h;\n  }\nsub make_output_structure\n  {\n  \
  my ($outdir,%h)=@_;\n    my $fh=new FileHandle;\\
n    mymkdir ($outdir) || die \"Could not create $\
outdir\\n\";\n    \n    foreach my $f (keys (%h))\\
n	{\n	  mymkdir ($h{$f}{template_dir}) or die \"1 \
Could not create $h{$f}{template_dir}\\n\";\n	  op\
en  ($fh, \">$h{$f}{template_file}{absolute}\") or\
 die \"Could not open $h{$f}{template_file}{absolu\
te}\";\n	  print $fh \"$h{$f}{templates}\";\n	  cl\
ose $fh;\n	}\n    return 1;\n  }\nsub is_installed\
\n  {\n    my $p=@_[0];\n    my $r=0;\n    my $cwh\
ich=\"$tmpdir/which\";\n    if (-e $cwhich){unlink\
 (\"$cwhich\");}\n    \n    system (\"which $p >$c\
which 2>/dev/null\");\n    my $w=file2string ($cwh\
ich);\n    if (($w=~/mmseqs/)){$r=1;}\n    return \
$r;\n  }\n\n  \nsub split_mmseqs\n    {\n      my \
($qf,$cacheq,$updateq, $db,$cachedb, $updatedb, $s\
,$out,$split)=@_;\n      \n      my @dbl=splitfast\
a($split,(string2fasta_list($db)));\n      \n\n   \
   if ( -e $out){unlink ($out)}\n\n\n     \n      \
foreach my $d (@dbl)\n	{\n	  my $uid=getuid();\n	 \
 my $lcacheo=\"$tmpdir/$uid/search/\";\n	  my $lca\
chedb=(($d =~/$tmpdir/))?\"$tmpdir/$uid/db/\":$cac\
hedb;\n	  my $lcacheq=$cacheq;\n	  mymkdir ($lcach\
eo, $lcachedb);\n	  \n	  print \"! Process Databas\
e $d\\n\";\n	  \n	  my %db=file2db($d ,$lcachedb,$\
updatedb);\n	  my %q =file2db($qf,$lcacheq ,$updat\
eq );\n	  my $ld=abs2file ($d);\n\n	  if (! -d $lc\
acheo){die \"NO CACHE\";}\n\n	  $q{search }=\"$lca\
cheo/$ld\\.MMSEQSSEARCH\";\n	  $q{convert}=\"$lcac\
heo/$ld\\.MMSEQSCONVERT\";\n	 \n	  system (\"mmseq\
s search $q{db} $db{db} $q{search} $db{index} $s -\
a $QUIET\");\n	  system (\"mmseqs convertalis $q{d\
b} $db{db} $q{search} $q{convert} --format-output \
\\\"query,target,qaln,taln,qstart,qend,pident,qcov\
,qlen\\\" $QUIET\");\n	  system (\"cat $q{convert}\
 >> $out\");\n	  \n	}\n      return $out;\n    }\n\
    \nsub splitfasta \n      {\n	my ($split,@list)\
=@_;\n	my @fl;\n\n		\n	if (!$split){return @list;}\
\n	\n	foreach my $e (@list)\n	  {\n	    \n	    my \
$n=`grep -c \">\" $e`;\n	    \n	    if ($n>$split)\
\n	      {\n		my $uid=getuid();\n		my $odir=\"$tmp\
dir/$uid/\";\n		mymkdir ($odir);\n		system (\"t_co\
ffee -other_pg seq_reformat -action +odir $odir +s\
plit $e $split\");\n		push (@fl,string2list (\"$od\
ir/*.split\"));\n	      }\n	    else\n	      {\n		\
push (@fl, $e);\n	      }\n	  }\n	return @fl;\n   \
   }\nsub getuid\n	{\n	  my $n;\n	  my $l=3;\n	  m\
y $string=randomstring ($l);\n	  while ($R{$string\
})\n	    {\n	      $n++;\n	      \n	      if ($n==\
10){$l++;}\n	      $string=randomstring($l);\n	   \
 }\n	  $R{$string}=1;\n	  return $string;\n	}\n		\\
nsub randomstring\n	  {\n	    my $l=shift;\n	    m\
y @s;\n	    my @alp=split (//, 'abcdefghijklmnopqr\
stuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_');\n\
	    \n	    my $lalp=@alp;\n	  \n	    for ( my $i=\
0; $i<$l; $i++)\n	      {\n		$s[$i]=$alp[rand($lal\
p)];\n	      }\n	    return join ('',@s);\n	  }\n	\
  \n		  \nmy $program=\"T-COFFEE (Version_13.46.0.\
919e8c6b)\";\n","use Env;\nuse strict;\nuse FileHa\
ndle;\nuse Cwd;\nuse File::Path;\nuse Sys::Hostnam\
e;\nuse File::Temp qw/ tempfile tempdir /;\nmy $QU\
IET=\"2>/dev/null\";\nmy $VERBOSE=$ENV{VERBOSE_4_D\
YNAMIC};\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=\
0;\nour $LAST_COM=\"\";\nour $NOEXIT=0;\nmy %metho\
d;\nmy $method2use;\nmy $treeF;\nmy $tree=$ENV{\"c\
hild_tree_4_TCOFFEE\"};\nmy $thread=$ENV{\"child_t\
hread_4_TCOFFEE\"};\nmy $dynamic=$ENV{dynamic_conf\
ig_4_TCOFFEE};\nmy $clean;\nmy $treeFlag;\nmy $bla\
stFlag;\nmy $infile;\nmy $outfile;\nmy $flush;\nmy\
 $do_exit=0;\nmy ($h1, $h2);\nmy @tmpL;\nmy $tmpdi\
r = File::Temp->newdir();\nmy $stderrF=\"$tmpdir/s\
tderr\";\n$QUIET=\"2>$stderrF\";\nmy $CDIR=getcwd(\
);\nmy $threadFlag4tc;\nmy $threadFlag4famsa;\nmy \
$threadFlag;\nmy $tcarg;\n\nmy $QUIET_ENV=$ENV{QUI\
ET_ENV};\n\nif ($QUIET_ENV==1){$QUIET=\"\";}\n\nfo\
r ($a=0; $a<=$#ARGV; $a++)\n  {\n    if    ($ARGV[\
$a] eq \"-seq\"){$infile=file2abs($ARGV[++$a]);}\n\
    elsif ($ARGV[$a] eq \"-outfile\"){$outfile=fil\
e2abs($ARGV[++$a], \"new\");}\n    elsif ($ARGV[$a\
] eq \"-dynamic_config\"){\n    	$dynamic=file2abs\
($ARGV[++$a]);\n    	if ($VERBOSE){print \"\\n![dy\
namic.pl] --- -dynamic_config flag if/else--- $dyn\
amic\\n\";}\n	}\n    \n    elsif ($ARGV[$a] eq \"-\
tree\") {$tree=$ARGV[++$a];}\n    elsif ($ARGV[$a]\
 eq \"-method\") {$method2use=$ARGV[++$a];}\n    e\
lsif ($ARGV[$a] eq \"-verbose\"){$VERBOSE=1; $QUIE\
T=\"\";}\n    elsif ($ARGV[$a] eq \"-clean\"){$cle\
an=1;}\n  \n    elsif ($ARGV[$a] eq \"-thread\"){$\
thread=$ARGV[++$a]}\n    elsif ($ARGV[$a] eq \"-tc\
arg\") {$tcarg=file2string($ARGV[++$a]);}\n    els\
if ($ARGV[$a] eq \"-noexit\") {$NOEXIT=1;}\n    el\
se \n      {\n	add2tcenv($a++,@ARGV);\n      }\n  \
}\n\n\n\nif ($tree eq \"list\")\n  {\n    my $f=\"\
$tmpdir/f\";\n    open (F, \">$f\");\n    print F \
\">a\\nxxx\\n>b\\nyyyyy\\n\";\n    close (F);\n   \
 print STDOUT (\"**** Supported Guide tree modes:\\
\n\");\n    my_system (\"t_coffee -other_pg seq_re\
format -in $f -action +seq2dnd list \");\n    $do_\
exit=1;\n  }\n\nif ($method2use eq \"list\")\n  {\\
n    my %ml;\n    my $listfile=\"$tmpdir/list\";\n\
    \n    $ml{tcoffee}=1;\n    $ml{psicoffee}=1;\n\
    $ml{accurate}=1;\n    $ml{'3dcoffee'}=1;\n    \
$ml{expresso}=1;\n    $ml{clustalo}=1;\n    $ml{ma\
fft}=1;\n    $ml{famsa}=1;\n      print STDOUT (\"\
**** Supported MSA mode:\\n\");\n    my_system (\"\
t_coffee 2>/dev/null | grep _msa > $listfile\");\n\
    open (F, $listfile);\n    while (<F>)\n      {\
\n	my $l=$_;\n	$l=~/(.*_msa)\\s+(.*)/;\n	my $m=$1;\
\n	my $i=\"$2\\n\";\n	if ($m=~/mafftsparsescore/)\\
n	  {\n	   printf STDOUT \"%-20s DOES NOT Support \
[-tree] -- $i\", $m;\n	  }\n	elsif ($m=~/tcoffee/)\
{;}\n	elsif ($m=~/mafft/){;}\n	elsif (!$ml{$m})\n	\
  {\n	    printf STDOUT \"%-20s DOES     Support [\
-tree] -- $i\", $m;\n	  }\n      }\n    $do_exit=1\
;\n  }\nif ($do_exit){my_exit ($CDIR,$EXIT_SUCCESS\
);}\n\nmy $NSEQ=file2nseq($infile);\n\n\nif ($NSEQ\
==0)\n  {\n    print \"ERROR - No sequences provid\
ed [FATAL:dynamic.pl]\\n\";\n    exit ($EXIT_FAILU\
RE);\n  }\nif (!$outfile)\n  {\n    ($h1,$outfile)\
=tempfile();\n    push (@tmpL,$outfile);\n    $flu\
sh=1;\n  }\n\nif (!($method2use=~/dynamic/)){;}\ne\
lse \n  {\n    if (-e $dynamic)\n      {\n       i\
f ($VERBOSE){print \"\\n![dynamic.pl] --- -dynamic\
_config FILE: \\n\";}\n        my @dynamicFile;\n	\
open (F, $dynamic);\n	while (<F>)\n	  {\n	    my $\
f=$_;\n	    if ($VERBOSE){print \"\\n![dynamic.pl]\
 --- FILE content: $f\\n\";}\n	    ## $f=~/(\\W)+ \
(\\d)+/;\n	    @dynamicFile = split ' ', $f;\n	   \
 if ($VERBOSE){print \"\\n![dynamic.pl] --- -dynam\
ic_config --- $dynamicFile[0] :: $dynamicFile[1]\\\
n\";}\n	    $method{$dynamicFile[0]} = $dynamicFil\
e[1];\n	  }\n	close(F);\n      }\n    else\n      \
{\n	$method{\"psicoffee_msa\"}=50;\n	$method{\"fam\
sa_msa\"}=1000000000;\n      }\n    \n    foreach \
my $name (sort { $method{$a} <=> $method{$b} } key\
s %method) \n      {\n	if ($NSEQ<=$method{$name})\\
n	  {\n	    $method2use=$name;\n	    last;\n	  }\n\
      }\n  }\n\n\nif ($tree)\n  {\n    ($h2,$treeF\
)=tempfile();\n    my ($h2,$tmptree)=tempfile();\n\
    push (@tmpL,$treeF);\n    if ( $tree eq \"defa\
ult\"){$treeF=0;}\n    elsif ( -e $tree)\n      {\\
n	my_system (\"cp $tree $tmptree\");\n      }\n   \
 elsif ($tree eq \"master\" || $tree eq \"main\" |\
| $tree eq \"parent\")\n      {\n	if ($ENV{CHILD_T\
REEF_4_TCOFFEE} && -e $ENV{CHILD_TREEF_4_TCOFFEE})\
\n	  {\n	    my $ctree=$ENV{CHILD_TREEF_4_TCOFFEE}\
;\n	    my_system (\"mv $ctree $tmptree\");\n	  }\\
n	else\n	  {\n	    my $master_tree=$ENV{CHILD_TREE\
_4_TCOFFEE};\n	    my_system (\"t_coffee -other_pg\
 seq_reformat -in $master_tree -in2 $infile -actio\
n +prune_tree -output newick > $tmptree\");\n	  }\\
n      }\n    else \n      {\n	my_system (\"t_coff\
ee -other_pg seq_reformat -in $infile -action +seq\
2dnd $tree -output newick> $tmptree\");\n      }\n\
    \n    if ($method2use=~/mafft/)\n      {\n	#pr\
int \"cp $tmptree /Users/cnotredame/.Trash/$$.tmpt\
ree\\n\";\n	#system (\"cp $tmptree /Users/cnotreda\
me/.Trash/$$.tmptree\");\n	\n	my_system (\"t_coffe\
e -other_pg seq_reformat -in $tmptree -output maff\
tdndmatrix> $treeF\");\n      }\n    else\n      {\
\n	my_system (\"mv $tmptree $treeF\");\n      }\n \
 }\nchdir ($tmpdir);\n\nmy $CL4tc=get_cl4tc();#wil\
l collect from env every CLTCOFEE env variable\n\n\
if (!$treeF || $NSEQ<=2){$treeFlag=\"\";}\nelsif (\
 $method2use=~/coffee/ || $method2use=~/accurate/)\
{$treeFlag=\"-usetree $treeF \";}\nelsif ( $method\
2use=~/clustalo/){$treeFlag=\"--guidetree-in=$tree\
F \";}\nelsif ( $method2use=~/mafftsparsecore/){;}\
\nelsif ( $method2use=~/mafft/){$treeFlag=\"--tree\
in $treeF \";}\nelsif ( $method2use=~/famsa/){$tre\
eFlag=\"-gt import $treeF \";}\n$CL4tc.=\" $treeFl\
ag \";\n\n$threadFlag=($thread)?\"--thread $thread\
 \":\"--thread 1 \";\n$threadFlag4tc=($thread)?\"-\
thread $thread \":\"-thread 1 \";\n$threadFlag4fam\
sa=($thread)?\"-t $thread \":\"-t 1 \";\n$CL4tc.=\\
" $threadFlag4tc \";\n\nif ($VERBOSE){print \"\\n!\
[dynamic.pl] --- CL4tc == $CL4tc\\n\";}\n\nmy $cme\
thod=$method2use;\n$cmethod=~s/_pair/_msa/;\n$cmet\
hod=~s/_msa//;\n\nif ($VERBOSE){print \"\\n![dynam\
ic.pl] --- cmethod == $cmethod\\n\";}\n\nif (-e $m\
ethod2use)\n     {\n    \n      my $com=file2strin\
g($method2use);\n      \n      if (!($com=~s/\\$in\
put/$infile/))\n	  {\n	      print \"ERROR - provi\
ded command [$com] should specify <input> [FATAL:d\
ynamic.pl]\\n\";\n	      exit ($EXIT_FAILURE);\n	 \
 }\n      if (!($com=~s/\\$output/$outfile/))\n	  \
{\n	      print \"ERROR - provided command [$com] \
should specify <input> [FATAL:dynamic.pl]\\n\";\n	\
      exit ($EXIT_FAILURE);\n	  }\n            \n \
     my_system ($com);\n  }\n\nelsif ($cmethod eq \
\"tcoffee\"|| $cmethod eq \"t_coffee\" )\n  {\n   \
 my_system (\"t_coffee -seq $infile -outfile $outf\
ile -output fasta_aln $CL4tc>/dev/null  $QUIET\");\
    \n  }\nelsif ($cmethod=~/(.*coffee)/ || $cmeth\
od=~/(accurate)/ || $cmethod=~/(expresso)/)\n  {\n\
    my $mode=$1;\n    my_system (\"t_coffee  -mode\
 $mode -seq $infile -outfile $outfile -output fast\
a_aln $CL4tc >/dev/null  $QUIET\");    \n  }\nelsi\
f ($cmethod eq \"clustalo\")\n  {\n    my_system (\
\"clustalo -i $infile $treeFlag -o $outfile  --for\
ce $threadFlag $QUIET\");\n    }\nelsif ($cmethod \
=~/sparsecore/)\n{\n    my_system (\"mafft-sparsec\
ore.rb -i $infile -C \\\"--anysymbol\\\"> $outfile\
 $QUIET\");\n  }\nelsif (($cmethod =~/mafft/))\n  \
{\n    my $mm;\n    my $retree;\n    \n    if ( $c\
method eq \"mafft\" || $cmethod=~/\\-/ )\n      {\\
n	$mm=$cmethod;\n      }\n    elsif (($cmethod=~/m\
afft(.*)/))\n      {\n	$mm=\"mafft-\".$1;\n      }\
\n\n    if ($mm =~/1/)\n      {\n	$mm=~s/1/i/;\n	$\
retree=\"--retree 1 \"\n      }\n    \n    my_syst\
em (\"$mm --anysymbol $threadFlag $treeFlag $retre\
e $infile > $outfile $QUIET\");\n  }\n\nelsif ($me\
thod2use=~/famsa/)\n  {\n    print \"\\n![dynamic.\
pl] --- FAMSA DEFAULT\\n\";    \n    my_system (\"\
famsa $treeFlag $threadFlag4famsa $infile $outfile\
 >/dev/null $QUIET\");\n  }\nelsif ($method2use=~/\
famsaUpgma/)\n  {\n    print \"\\n![dynamic.pl] --\
- FAMSA Upgma\\n\";\n    print \"\\n![dynamic.pl] \
--- Command: famsa -gt upgma $treeFlag $threadFlag\
4famsa $infile $outfile >/dev/null $QUIET\\n\";\n \
   my_system (\"famsa -gt upgma $treeFlag $threadF\
lag4famsa $infile $outfile >/dev/null $QUIET\");\n\
  }\nelse\n  {\n    if ($treeF)\n      {\n	printf \
(STDERR \"WARNING: Method $method2use CANNOT use p\
re-sepecified guide tree [dynamic.pl]\\n\");\n    \
  }\n    my_system (\"t_coffee -in $infile -method\
 $method2use -outfile $outfile -output fasta_aln $\
tcarg -quiet $QUIET\");\n  }\n\nif ( ! -e $outfile\
 || -s $outfile <1)\n{\n    if ($NOEXIT==0)\n     \
 {  \n	print \"ERROR - No MSA computed - $LAST_COM\
 -- [FATAL:dynamic.pl]\\n\";\n	my_exit ($CDIR,$EXI\
T_FAILURE);\n      }\n    else\n      {\n	print \"\
WARNING - No MSA computed - $LAST_COM -- [WARNING:\
dynamic.pl]\\n\";\n      }\n}\nelsif ( $flush)\n {\
\n   open (F, $outfile);\n   while (<F>){print $_;\
}\n   close (F);\n }\nforeach my $f (@tmpL){unlink\
($f);}\n\n\nif ($VERBOSE!=-1)\n  {\n    open (F, \\
"$stderrF\");\n    while (<F>)\n      {\n	my $l=$_\
;\n	print (\"------  $l\");\n	if ( $VERBOSE || $l=\
~/WARNING/ || $l=~/ERROR/ || $l=~/INFORMATION/){pr\
int stderr \"$l\";}\n      }\n    close (F);\n    \
if (-e $stderrF){unlink ($stderrF);}\n  }\n\nmy_ex\
it ($CDIR,$EXIT_SUCCESS);\n\n\nsub file2nseq\n  {\\
n    my ($f)=@_;\n    my $n=`grep -c '>' $f`; \n  \
  \n    return $n;\n  }\nsub file2abs\n    {\n    \
  my ($f, $mode)=@_;\n      \n      if (!$f || $f=\
~/^\\//){return $f;}\n      elsif (!-e $f && $mode\
 eq \"new\"){return \"$CDIR/$f\";}\n      elsif (!\
-e $f){return $f;}\n    \n      return \"$CDIR/$f\\
";\n    }\nsub file2string \n    {\n      my ($f)=\
@_;\n      my $s;\n\n      open (F, $f) || return \
0;\n      while (<F>)\n	{\n	  $s.=$_;\n	}\n      c\
lose (F);\n      chomp($s);\n      return $s;\n   \
 }   \n\nsub get_psicl\n      {\n	my ($psitrim, $p\
sitrim_mode, $pisN);\n	my $cl;\n	\n	if ($ENV{psitr\
im_tree_4_TCOFFEE}){$cl.=\" -psitrim_tree=\".$ENV{\
psitrim_tree_4_TCOFFEE}.\" \";}\n	if ($ENV{psitrim\
_mode_4_TCOFFEE}){$cl.=\" -psitrim_mode=\".$ENV{ps\
itrim_mode_4_TCOFFEE}.\" \";}\n	if ($ENV{psitrim_4\
_TCOFFEE}){$cl.=\" -psitrim=\".$ENV{psitrim_4_TCOF\
FEE}.\" \";}\n	if ($ENV{psiJ_4_TCOFFEE}){$cl.=\" -\
psiJ=\".$ENV{psiJ_4_TCOFFEE}.\" \";}\n	\n\n	return\
 $cl;\n      }\n      \nsub get_cl4tc\n	{\n	  my $\
cl;\n	  \n	  foreach my $arg (keys(%ENV))\n	    {\\
n	      if ($arg=~/(.*)_4_CLTCOFFEE/)\n		{\n		  my\
 $name=$1;\n		  my $val=$ENV{$arg};\n		  if (-e $v\
al){$val=file2abs($val);}\n		  \n\n		  if ($val eq\
 \"FLAGSET\"){$val=\"\";}\n		  $cl.=\"-$name $val \
\";\n		}\n	    }\n	  if ($VERBOSE){print \"\\n![dy\
namic.pl] --- get_psicl --- $cl\\n\";}\n	  return \
$cl;\n	}\n\nsub add2tcenv\n	    {\n	      my ($p, \
@argv)=@_;\n\n	      my $flag=$argv[$p];\n	      $\
flag =~s/^-//;\n	      my $val =file2abs($argv[$p+\
1]);\n	      my $envv=\"$flag\\_4_CLTCOFFEE\";\n	 \
     $ENV{$envv}=$val;\n	    }\n	      \nsub my_ex\
it\n    {\n      my ($dir,$ec)=@_;\n      my $a;\n\
      if ($VERBOSE)\n	{\n	  print \"\\n![dynamic.p\
l] --- CDIR: $CDIR\\n\";\n	  print \"\\n![dynamic.\
pl] --- Processed $NSEQ\\n\";\n	  print \"\\n![dyn\
amic.pl] --- \";\n	  foreach my $arg (@ARGV)\n	   \
 {\n	      print \"$arg \";\n	    }\n	  \n	  print\
 \"\\n![dynamic.pl] --- EXIT: $ec ($EXIT_SUCCESS:s\
uccess, $EXIT_FAILURE:failure)-- Verbose mode -- u\
nset VERBOSE_4_DYNAMIC to turn verbose mode off\\n\
\";\n	}\n      chdir ($dir);\n      exit ($ec);\n \
   }\n	      \nsub my_system \n  {\n    my ($com)=\
@_;\n    $LAST_COM=$com;\n    \n    if ($VERBOSE){\
print \"\\n![dynamic.pl] --- SysCall --- $com\\n\"\
;}\n    \n    system ($com);\n  }\nmy $program=\"T\
-COFFEE (Version_13.46.0.919e8c6b)\";\n","use Env;\
\nuse FileHandle;\nuse Cwd;\nuse File::Path;\nuse \
Sys::Hostname;\nmy $f = new FileHandle;\n\nopen ($\
f, $ARGV[1]);\n$atom=$ARGV[0];\n\n$atom=~s/PRIME/\\
\'/;\nwhile (<$f>)\n  {\n    my $l=$_;\n\n    $l=~\
s/$atom/CA /;\n    \n    \n    $l=~s/  G /GLY /g;\\
n    $l=~s/  C /CYS /g;\n    $l=~s/  T /THR /g;\n \
   $l=~s/  A /ALA /g;\n    $l=~s/  U /THR /g;\n   \
 \n    $l=~s/ DG /GLY /g;\n    $l=~s/ DC /CYS /g;\\
n    $l=~s/ DT /THR /g;\n    $l=~s/ DA /ALA /g;\n \
   $l=~s/ DU /THR /g;\n    \n    print $l;\n  }\n\\
n\n\n","*TC_METHOD_FORMAT_01\n******************ge\
neric_method.tc_method*************\n*\n*       In\
corporating new methods in T-Coffee\n*       Cedri\
c Notredame 26/08/08\n*\n*************************\
******************************\n*This file is a me\
thod file\n*Copy it and adapt it to your need so t\
hat the method \n*you want to use can be incorpora\
ted within T-Coffee\n*****************************\
**************************\n*                  USA\
GE                              *\n***************\
****************************************\n*This fi\
le is passed to t_coffee via -in:\n*\n*	t_coffee -\
in Mgeneric_method.method\n*\n*	The method is pass\
ed to the shell using the following\n*call:\n*<EXE\
CUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FL\
AG><outname><PARAM>\n*\n*Conventions:\n*<FLAG_NAME\
> 	<TYPE>		<VALUE>\n*<VALUE>:	no_name 	<=> Replace\
d with a space\n*<VALUE>:	&nbsp	<=> Replaced with \
a space\n*\n**************************************\
*****************\n*                  ALN_MODE    \
                       *\n************************\
*******************************\n*pairwise   ->all\
 Vs all (no self )[(n2-n)/2aln]\n*m_pairwise ->all\
 Vs all (no self)[n^2-n]^2\n*s_pairwise ->all Vs a\
ll (self): [n^2-n]/2 + n\n*multiple   ->All the se\
quences in one go\n*\nALN_MODE		pairwise\n*\n*****\
**************************************************\
\n*                  OUT_MODE                     \
      *\n*****************************************\
**************\n* mode for the output:\n*External \
methods: \n* aln -> alignmnent File (Fasta or Clus\
talW Format)\n* lib-> Lib file (TC_LIB_FORMAT_01)\\
n*Internal Methods:\n* fL -> Internal Function ret\
urning a List (Librairie)\n* fA -> Internal Functi\
on returning an Alignmnent\n*\nOUT_MODE		aln\n****\
**************************************************\
*\n*                  SEQ_TYPE                    \
       *\n****************************************\
***************\n*G: Genomic, S: Sequence, P: PDB,\
 R: Profile\n*Examples:\n*SEQTYPE	S	sequences agai\
nst sequences (default)\n*SEQTYPE	S_P	sequence aga\
inst structure\n*SEQTYPE	P_P	structure against str\
ucture\n*SEQTYPE	PS	mix of sequences and structure\
	\n*\nSEQ_TYPE	S\n*\n\n***************************\
****************************\n*                COM\
MAND LINE                         *\n*EXECUTABLE P\
ARAM1 IN_FLAG OUT_FLAG PARAM             *\n******\
*************************************************\\
n*************************************************\
******\n*                  EXECUTABLE             \
            *\n***********************************\
********************\n*name of the executable\n*pa\
ssed to the shell: executable\n*	\nEXECUTABLE	tc_g\
eneric_method.pl\n*\n*****************************\
**************************\n*                  IN_\
FLAG                             *\n**************\
*****************************************\n*IN_FLA\
G\n*flag indicating the name of the in coming sequ\
ences\n*IN_FLAG S no_name ->no flag\n*IN_FLAG S &b\
nsp-in&bnsp -> \" -in \"\n*\nIN_FLAG		-infile=\n*\\
n*************************************************\
******\n*                  OUT_FLAG               \
            *\n***********************************\
********************\n*OUT_FLAG\n*flag indicating \
the name of the out-coming data\n*same conventions\
 as IN_FLAG\n*OUT_FLAG	S no_name ->no flag\n*if yo\
u want to redirect, pass the parameters via PARAM1\
\n*set OUT_FLAG to >\n*\nOUT_FLAG		-outfile=\n*\n*\
**************************************************\
****\n*                  PARAM_1                  \
            *\n***********************************\
********************\n*<EXECUTABLE><PARAM1><IN_FLA\
G><seq_file><PARAM2><OUT_FLAG><outname><PARAM>\n*P\
arameters sent to the EXECUTABLE and specified *be\
fore* IN_FLAG \n*If there is more than 1 PARAM lin\
e, the lines are\n*concatenated\n*Command_line: @E\
P@PARAM@-gapopen%e10%s-gapext%e20\n*	%s white spac\
e\n*	%e equal sign\n*\n*PARAM1	\n*\n*\n*\n********\
***********************************************\n*\
                  PARAM_2                         \
     *\n******************************************\
*************\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_\
file><PARAM2><OUT_FLAG><outname><PARAM>\n*Paramete\
rs sent to the EXECUTABLE and specified \n*after* \
IN_FLAG and *before* OUT_FLAG\n*If there is more t\
han 1 PARAM line, the lines are\n*concatenated\n*\\
n*PARAM1	\n*\n*\n*********************************\
**********************\n*                  PARAM  \
                            *\n*******************\
************************************\n*<EXECUTABLE\
><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><out\
name><PARAM>\n*Parameters sent to the EXECUTABLE a\
nd specified *after* OUT_FLAG\n*If there is more t\
han 1 PARAM line, the lines are\n*concatenated\n*\\
nPARAM	-mode=seq_msa -method=clustalw2\nPARAM   -O\
UTORDER=INPUT -NEWTREE=core -align -gapopen=-15\n*\
\n************************************************\
*******\n*                  END                   \
             *\n**********************************\
*********************\n","*TC_METHOD_FORMAT_01\n**\
*************clustalw_method.tc_method*********\nE\
XECUTABLE	clustalw\nALN_MODE		pairwise\nIN_FLAG		-\
INFILE=\nOUT_FLAG		-OUTFILE=\nOUT_MODE		aln\nPARAM\
		-gapopen=-10\nSEQ_TYPE		S\n*********************\
****************************\n","$VersionTag =    \
                                                  \
                                                  \
                           2.43;\nuse Env;\nuse Fi\
leHandle;\nuse Cwd;\nuse File::Path;\nuse Sys::Hos\
tname;\n\nour $PIDCHILD;\nour $ERROR_DONE;\nour @T\
MPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCC\
ESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCESS=0\
;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\"extract_\
from_pdb\";\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT_\
STARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\"};\no\
ur $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LO\
CKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=$ENV{\"E\
RRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\"ERROR\
FILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (isshellpi\
d(getppid())){lock4tc(getppid(), \"LLOCK\", \"LSET\
\", \"$$\\n\");}\n      \nour $SILENT=\" >/dev/nul\
l 2>/dev/null\";\nour $INTERNET=-1;\n\n\n\n\n\n\n\\
nour $BLAST_MAX_NRUNS=2;\nour $EXIT_SUCCESS=0;\nou\
r $EXIT_FAILURE=1;\nour $CONFIGURATION=-1;\nour $R\
EF_EMAIL=\"\";\nour $PROGRAM=\"extract_from_pdb\";\
\n\n\nmy %onelett_prot=&fill_onelett_prot();\nmy %\
threelett_prot=&fill_threelett_prot();\nmy %onelet\
t_RNA=&fill_onelett_RNA();\nmy %threelett_RNA=&fil\
l_threelett_RNA();\nmy %onelett_DNA=&fill_onelett_\
DNA();\nmy %threelett_DNA=&fill_threelett_DNA();\n\
\n\n\n\n\nmy %onelett = (\n'P' => \\%onelett_prot,\
\n'D' => \\%onelett_DNA,\n'R' => \\%onelett_RNA\n)\
;\n\n\nmy %threelett = (\n'P' => \\%threelett_prot\
,\n'D' => \\%threelett_DNA,\n'R' => \\%threelett_R\
NA\n);\n\n\n\n\n\n\n\nif($ARGV[0]=~/help/ ||$ARGV[\
0]=~/man/ || $ARGV[0]=~/HELP/ || $ARGV[0]=~/Man/ |\
| $ARGV[0] eq \"-h\"  || $ARGV[0] eq \"-H\"  )\n{d\
ie \"SYNTAX: extract_from_pdb Version $VersionTag	\
\n	Minimum:             [extract_from_pdb file] \n\
			   OR \n			     [... | extract_from_pdb]\n 	Fla\
gs (Default setting on the first line)\n	   -versi\
on...................[Returns the Version Number]\\
n           -force.....................[Forces the\
 file to be treated like a PDB file]\n            \
                          [Regenerates the header \
and SEQRES fields]\n           -force_name........\
........[Forces the file to be named after name]]\\
n           -infile.....file...........[Flag can b\
e omited]\n			              [File must be pdb or f\
ro pgm]\n                                      [Fi\
le can also be compressed Z or gz]\n              \
                        [In the case of a compress\
ed file, you can omit the gz|Z extension]\n       \
    -netfile...................[File will be fetch\
 from the net using wget]\n                       \
               [wget or curl must be installed]\n \
                                     [ftp://ftp.gn\
u.org/pub/gnu/wget/]\n                            \
          [http://curl.haxx.se/]\n                \
                      [Must also be used to retrie\
ve the file from a local pdb copy (cf netaddress)]\
\n           -netaddress................[Address u\
sed for the retrieving the netfile]\n             \
                         [http://www.rcsb.org/pdb/\
cgi/export.cgi/%%.pdb.gz?format=PDB&pdbId=%%&compr\
ession=gz]\n                                      \
[http://www.expasy.ch/cgi-bin/get-pdb-entry.pl?%%]\
\n                                      [local -> \
will get the file from pdb_dir (see pdb_dir)]\n   \
        -netcompression............[Extension if t\
he netfile comes compressed]\n                    \
                  [gz]\n           -pdb_dir.......\
............[address of the repertory where the pd\
b is installed]\n                                 \
     [Supports standard ftp style installation OR \
every stru in DIR]\n                              \
        [Give the ..../pdb/structure/ dir]\n      \
                                [If value omitted,\
 the pg gets it from the env variable PDB_DIR]\n  \
         -netcompression_pg.........[gunzip]\n    \
       -is_pdb_name..........name.[Returns 1 if th\
e name is a PDB ID, 0 otherwise]\n           -mode\
l_type...........name.[Returns the model type if v\
alid PDB name]\n           -is_released_pdb_name n\
ame.[Returns 1 if the name corresponds to a releas\
ed PDB file]\n           -get_pdb_chains.....name.\
..[Returns the list of chains corresponding to the\
 entry]\n           -get_pdb_id.........name...[Re\
turns the PDB id within the provided pdb file]\n  \
         -get_fugue_name.....name...[Turns a name \
into a name valid for fugue]\n                    \
                  [Uses the netaddress to do so]\n\
	   -chain......FIRST..........[Extract the first \
chain only]\n		       A B C..........[Extract Seve\
ral chains if needed]\n		       ALL............[Ex\
tract all the chains]	\n           -ligand.....ALL\
............[Extract the ligands in the chain (HET\
ATM)]\n                       <name1>,<name2>[Extr\
act All the named lignds]\n	   -ligand_only.......\
........[Extract only the ligands]\n           -li\
gand_list...............[Extract the list of ligan\
ds]\n	   -coor.......<start>..<end>.[Coordinates o\
f the fragment to extract]\n			              [Omit\
 end to include the Cter]\n           -num........\
absolute.......[absolute: relative to the seq] \n \
                      file...........[file: relati\
ve to file]\n           -num_out....new...........\
.[new: start 1->L]\n                       old....\
........[old: keep the file coordinates]\n        \
   -delete.....<start>..<end>.[Delete from residue\
 start to residue end]\n	   -atom.......CA........\
.....[Atoms to include, ALL for all of them]\n		  \
     CA O N.........[Indicate several atoms if nee\
ded]\n	   -code.......3..............[Use the 1 le\
tter code or the 3 letters code]\n	   -mode.......\
raw............[Output original pdb file]\n       \
                pdb............[Output something t\
hat looks like pdb]\n		       fasta..........[Outp\
ut the sequences in fasta format]\n		       simple\
.........[Output a format easy to parse in C ]\n  \
          -seq_field..ATOM...........[Field used t\
o extract the sequence]\n		       SEQRES.........[\
Use the complete sequence]\n	   -seq..............\
.........[Equivalent to  -mode fasta]\n	   -model.\
.....1..............[Chosen Model in an NMR file]\\
n           -nodiagnostic..............[Switches E\
rror Messages off]\n           -debug.............\
........[Sets the DEBUG ON]\n           -no_remote\
_pdb_dir.........[Do not look for a remote file]\n\
           -cache_pdb.................[Cache Value\
, default is $HOME/.t_coffee/cache, other values: \
NO<=> No cache]\n\n      Environement Variables\n \
          These variables can be set from the envi\
ronement\n           Command line values with the \
corresponding flag superseed evironement value\n  \
         NO_REMOTE_PDB_DIR..........[Prevents the \
program from searching remote file: faster]\n     \
      PDB_DIR....................[Indicates where \
PDB file must be fetched (localy)]\n\n	 PROBLEMS: \
please contact cedric.notredame\\@europe.com\\n\";\
\n	 exit ($EXIT_SUCCESS);\n}\n\n$np=0;\n$n_para=$#\
ARGV;\n$model=1;\n$pdb_dir=$ENV{'PDB_DIR'};if ($pd\
b_dir){$pdb_dir.=\"/\";}\n$debug=$ENV{'DEBUG_EXTRA\
CT_FROM_PDB'};\n\n$no_remote_pdb_dir=$ENV{NO_REMOT\
E_PDB_DIR};\n$HOME=$ENV{'HOME'};\nif ( $ENV{CACHE_\
4_TCOFFEE})\n{$cache=$ENV{CACHE_4_TCOFFEE};}\nelse\
\n{\n    $cache=\"$HOME/.t_coffee/cache/\";\n}\n\n\
   \n$netaddress=\"https://files.rcsb.org/download\
/%%.pdb.gz\";\n$netcompression_pg=\"gunzip\";\n$ne\
tcompression=\"gz\";\n\nforeach ($np=0; $np<=$n_pa\
ra; $np++)\n  {        \n    $value=$ARGV[$np];\n \
  \n    if  ($np==0 && !($value=~/^-.*/))\n      {\
 \n	$pdb_file= $ARGV[$np];\n      }\n    elsif ( !\
($value=~/^-.*/))\n      {\n	print \"@ARGV\";\n	di\
e;\n      } \n    \n    elsif ($value eq \"-nodiag\
nostic\"){$nodiagnostic=1;}\n    elsif ($value eq \
\"-force\")\n      {\n	$force_pdb=1;\n      }\n   \
 elsif ($value eq \"-force_name\")\n      {\n	$for\
ce_name=$ARGV[++$np];\n	$force_pdb=1;\n      }\n  \
  \n    elsif ($value eq \"-is_pdb_name\")\n      \
{\n	$pdb_file= $ARGV[++$np];	\n	$is_pdb_name=1;	\n\
      } \n    elsif ($value eq \"-is_released_pdb_\
name\")\n      {\n	$pdb_file= $ARGV[++$np];\n	\n	i\
f (!$pdb_file){print \"0\";exit (EXIT_SUCCESS);}\n\
	$is_released_pdb_name=1;\n      }\n    elsif ($va\
lue eq \"-model_type\")\n      {\n	$pdb_file= $ARG\
V[++$np];	\n	$model_type=1;\n      }\n    elsif ($\
value eq \"-debug\")\n{\n	$debug=1;\n}\n    elsif \
($value eq \"-get_pdb_chains\")\n{\n	$pdb_file= $A\
RGV[++$np];\n	$get_pdb_chains=1;\n}\n    elsif ($v\
alue eq \"-get_pdb_ligands\")\n{\n	$get_pdb_ligand\
s=1;\n}\n    \n    elsif ($value eq \"-get_pdb_id\\
")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb_id=1;\\
n	\n}\n    \n    elsif ( $value eq \"-get_fugue_na\
me\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_fugue_n\
ame=1;\n}\n    elsif ( $value eq \"-infile\")\n{\n\
       $pdb_file= $ARGV[++$np];\n} \n    elsif ($v\
alue eq \"-netfile\")\n{\n	$netfile=1;\n	if ( !($A\
RGV[$np+1]=~/^-.*/)){$pdb_file= $ARGV[++$np];}\n}\\
n    elsif (  $value eq \"-num\")\n{\n       $numb\
ering= $ARGV[++$np];\n}\n    elsif (  $value eq \"\
-num_out\")\n{\n       $numbering_out= $ARGV[++$np\
];\n}\n    elsif ( $value eq \"-netaddress\")\n{\n\
	$netadress=$ARGV[++$np];\n}\n     \n    elsif ( $\
value eq \"-netcompression\")\n{\n	 $netcompressio\
n=$ARGV[++$np];\n}\n    elsif ( $value eq \"-pdb_d\
ir\")\n{\n	 if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir\
= \"$ARGV[++$np]/\";}\n}\n     elsif ( $value eq \\
"-no_remote_pdb_dir\")\n{\n	$no_remote_pdb_dir=1;\\
n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir= \"$ARGV[\
++$np]/\";}\n}\n    elsif ( $value eq \"-cache\")\\
n{\n	$cache=$ARGV[++$np];\n}\n    \n    elsif ($va\
lue eq \"-netcompression_pg\")\n{\n	  $netcompress\
ion_pg=$ARGV[++$np];\n}\n     elsif ($value eq \"-\
mode\")\n{\n       $MODE=$ARGV[++$np];\n}\n\n    e\
lsif ( $value eq \"-model\")\n{\n       $model= $A\
RGV[++$np];\n}\n    elsif ($value eq \"-seq_field\\
" )\n{\n       $seq_field= $ARGV[++$np];\n}   \n  \
  elsif ($value eq \"-coor\" )\n{\n       $start= \
$ARGV[++$np];\n  \n       if (($ARGV[$np+1] eq \"\\
") ||($ARGV[$np+1]=~/^-.*/)){$end=\"*\";} \n      \
 else {$end=   $ARGV[++$np];}     \n       $coor_s\
et=1;\n}\n    elsif ($value eq \"-delete\" )\n{\n \
      $delete_start= $ARGV[++$np];\n       $delete\
_end= $ARGV[++$np];\n       $delete_set=1;\n}\n   \
 elsif  ($value eq \"-code\")\n{\n       $code= $A\
RGV[++$np];\n}\n    elsif  ($value eq \"-no_hetatm\
\")\n{\n       $no_hetatm=1;\n}\n    elsif ($value\
 eq \"-chain\")\n{\n       while (!($ARGV[$np+1] e\
q \"\") &&!($ARGV[$np+1]=~/^-.*/))\n{\n	      ++$n\
p;\n	      @c_chain=(@chain,  $ARGV[$np]);\n	     \
 $hc_chain{$ARGV[$np]}=$#c_chain+1;\n}           \\
n}\n    elsif ($value eq \"-atom\")\n{\n\n       w\
hile (!($ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/\
^-.*/))\n{\n	      ++$np;\n	      $atom[$n_atom++]\
=  $ARGV[$np];\n	      $atom_list{$ARGV[$np]}=1;	 \
     \n} \n       \n}\n    elsif ( $value eq \"-un\
fold\")\n{\n	$unfold=1;\n}\n    elsif ($value eq \\
"-seq\" ||$value eq \"-fasta\" )\n{\n       $MODE=\
\"fasta\";\n}\n    elsif ( $value eq \"-version\")\
\n{\n	print STDERR  \"\\nextract_from_pdb: Version\
 $VersionTag\\n\";\n	&myexit ($EXIT_SUCCESS);\n}\n\
    elsif ( $value eq \"-ligand\")\n{\n	while (!($\
ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/^-.*/))\n\
{\n	    ++$np;\n	    $ligand=1;\n	    $ligand_list\
{$ARGV[$np]}=1;	      \n} \n	$hc_chain{'LIGAND'}=1\
;\n}\n    elsif ( $value eq \"-ligand_only\")\n{\n\
	$ligand_only=1;\n}\n}\nif ( $debug)\n{\n    print\
 STDERR \"\\n[DEBUG:extract_from_pdb] NO_REMOTE_PD\
B_DIR: $no_remote_pdb_dir\\n\";\n    print STDERR \
\"\\n[DEBUG:extract_from_pdb] PDB_DIR: $pdb_dir\\n\
\";\n}\n\n\nif ( $is_pdb_name)\n  {\n    if (&remo\
te_is_pdb_name($pdb_file))\n      {\n	print \"1\";\
\n      }\n    else\n      {\n	print \"0\";\n     \
 }\n    exit ($EXIT_SUCCESS);\n  }\n\nif ( $is_rel\
eased_pdb_name)\n  {\n    \n    if (&is_released($\
pdb_file))\n      {\n	print \"1\";\n      }\n    e\
lse\n      {\n	print \"0\";\n      }\n    exit ($E\
XIT_SUCCESS);\n  }\nif ($model_type)\n  {\n   \n  \
  printf \"%s\", &pdb2model_type($pdb_file);\n    \
exit ($EXIT_SUCCESS);\n    \n  }\n    \n\nif (!$fo\
rce_name)\n{\n    $pdb_file=~/([^\\/]*)$/;\n    $f\
orce_name=$1;\n}\n\n$local_pdb_file=$pdb_file;\n\n\
if ( $debug){print STDERR \"\\n[DEBUG: extract_fro\
m_pdb] Scan For $local_pdb_file\\n\";}\n\n$mem=$no\
_remote_pdb_dir;\n$no_remote_pdb_dir=1;\n$tmp_pdb_\
file=get_pdb_file ($local_pdb_file);\n\nif ( !-e $\
tmp_pdb_file || $tmp_pdb_file eq \"\")\n  {\n    $\
local_pdb_file=$pdb_file;\n    ($local_pdb_file, $\
suffix_chain)=&pdb_name2name_and_chain($local_pdb_\
file);\n\n    if ($local_pdb_file)\n      {\n	if (\
 $debug){print STDERR \"\\nSplit $pdb_file into $l\
ocal_pdb_file and $suffix_chain \\n\";}\n	$tmp_pdb\
_file=get_pdb_file ($local_pdb_file);\n	if ( $tmp_\
pdb_file ne \"\")\n	  {\n	    @c_chain=();\n	    @\
c_chain=($suffix_chain);\n	    %hc_chain=();\n	   \
 $hc_chain{$suffix_chain}=1;\n	  }\n      }\n  }\n\
\n$no_remote_pdb_dir=$mem;\nif ($no_remote_pdb_dir\
==0)\n  {\n    \n    if ( !-e $tmp_pdb_file || $tm\
p_pdb_file eq \"\")\n      {\n	\n	$local_pdb_file=\
$pdb_file;\n	($local_pdb_file, $suffix_chain)=&pdb\
_name2name_and_chain($local_pdb_file);\n	if ($loca\
l_pdb_file)\n	  {\n	    \n	    if ( $debug){print \
STDERR \"\\nSplit $pdb_file into $local_pdb_file a\
nd $suffix_chain \\n\";}\n	    \n	    $tmp_pdb_fil\
e=get_pdb_file ($local_pdb_file);    \n	    \n	   \
 if ( $tmp_pdb_file ne \"\")\n	      {\n		@c_chain\
=();\n		@c_chain=($suffix_chain);\n		%hc_chain=();\
\n		$hc_chain{$suffix_chain}=1;\n	      }\n	  }\n \
     }\n  }\n\nif ( $debug){print STDERR \"\\n$pdb\
_file copied into ##$tmp_pdb_file##\\n\";}\n\nif (\
 !-e $tmp_pdb_file || $tmp_pdb_file eq \"\")\n{\n	\
if ($is_pdb_name)\n{\n	    print \"0\\n\"; exit ($\
EXIT_SUCCESS);\n}\n	else\n{\n  \n	    print \"\\nE\
XTRACT_FROM_PDB: NO RESULT for $pdb_file\\n\";\n	 \
   &myexit ($EXIT_SUCCESS);	\n}\n}\n\n\n\n\n%molec\
ule_type=&pdbfile2chaintype($tmp_pdb_file);\nif ( \
$debug){print STDERR \"\\n\\tSequence Type determi\
ned\\n\";}\n\n$pdb_id=&get_pdb_id ($tmp_pdb_file);\
\nif ( $debug){print STDERR \"\\n\\tID: $pdb_id (f\
or $tmp_pdb_file)\\n\";}\n\nif ( $pdb_id eq \"\"){\
$pdb_id=$force_name;}\n\n@f_chain=&get_chain_list \
($tmp_pdb_file);\nif ( $debug){print STDERR \"\\n\\
\tChain_list:@f_chain\\n\";}\n\nif ( $get_pdb_chai\
ns)\n{\n    print \"@f_chain\\n\";\n    &myexit ($\
EXIT_SUCCESS);\n}\nif ( $get_pdb_ligands)\n{\n    \
%complete_ligand_list=&get_ligand_list ($tmp_pdb_f\
ile);\n    print $complete_ligand_list{\"result\"}\
;\n    &myexit ($EXIT_SUCCESS);\n}\n\nelsif ( $get\
_pdb_id ||$get_fugue_name )\n{\n    if    (@c_chai\
n && $c_chain[0] eq \"FIRST\"){$pdb_id=$pdb_id.$f_\
chain[0];}\n    elsif (@c_chain && $c_chain[0] ne \
\" \"){$pdb_id=$pdb_id.$c_chain[0];}\n    \n    pr\
int \"$pdb_id\\n\";\n    &myexit ($EXIT_SUCCESS);\\
n    \n}\nelsif ( $is_pdb_name)\n{\n    printf \"1\
\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\n\n\n\n$st\
ructure_file=vtmpnam();\n\nif ( $debug){print STDE\
RR \"\\n\\tCheck_point #1: $tmp_pdb_file  $structu\
re_file\\n\";}\n\n$INFILE=vfopen (\"$tmp_pdb_file\\
", \"r\"); \nmy $TMP=vfopen (\"$structure_file\", \
\"w\");\n\n$print_model=1;\n$in_model=0;\n\nif ( $\
debug){print STDERR \"\\n\\tCheck_point #2\\n\";}\\
nwhile ( <$INFILE>)\n{\n  my $first_model=0;\n  $l\
ine=$_;\n\n  if ( !$first_model && ($line =~/^MODE\
L\\s*(\\d*)/))\n    {\n      $first_model=$1;\n   \
   if ($model==1){$model=$first_model;}\n    }\n  \
\n  if (($line =~/^MODEL\\s*(\\d*)/))\n    {\n    \
  if ($1==$model)\n	{\n	  $in_model=1;\n	  $print_\
model=1;\n	  $is_nmr=1;\n	}\n      elsif ( $in_mod\
el==0)\n	{\n	  $print_model=0;\n	}\n      elsif ( \
$in_model==1)\n	{\n	  last;\n	}\n    }\n  if ($pri\
nt_model){print $TMP $line;}  \n}\nclose ($TMP);\n\
close ($INFILE);\n\nif ( $debug){print STDERR \"\\\
n\\tCheck_point #3\\n\";}	\n\n  if ($numbering eq \
\"\"){$numbering=\"absolute\";}\n  if ($numbering_\
out eq \"\"){$numbering_out=\"new\";}\n\n  if ( $d\
elete_set && $coor_set) {die \"-delete and -coor a\
re mutually exclusive, sorry\\n\";}\n  if ( $n_ato\
m==0){$atom_list[$n_atom++]=\"ALL\";$atom_list{$at\
om_list[0]}=1;}\n  if ( $seq_field eq \"\"){$seq_f\
ield=\"ATOM\";}\n  \n  if ( $MODE eq \"\"){$MODE=\\
"pdb\";}\n  elsif ( $MODE eq \"simple\" && $code==\
0){$code=1;}\n\n  if ( $code==0){$code=3;}\n\n\nif\
 ($f_chain[0] eq \" \"){$hc_chain{' '}=1;$c_chain[\
0]=\" \";}\nelsif (!@c_chain){$hc_chain{FIRST}=1;$\
c_chain[0]=\"FIRST\";}#make sure the first chain i\
s taken by default\n\nif    ($hc_chain{ALL}) \n{\n\
      @c_chain=@f_chain;\n      foreach $e (@c_cha\
in){$hc_chain{$e}=1;}\n}\nelsif($hc_chain{FIRST})\\
n{\n	@c_chain=($f_chain[0]);\n	$hc_chain{$f_chain[\
0]}=1;\n}\n\n\n$MAIN_HOM_CODE=&get_main_hom_code (\
$structure_file);\n$INFILE=vfopen ($structure_file\
, \"r\");\n\n\nif ( $MODE eq \"raw_pdb\" || $MODE \
eq \"raw\")\n{\n    while (<$INFILE>)\n{	print \"$\
_\";}\n    close ( $INFILE);\n    &myexit($EXIT_SU\
CCESS);\n}    \nif ( $MODE eq \"raw4fugue\" )\n{\n\
    while (<$INFILE>)\n{	\n	$l=$_;\n	if ($l=~/^SEQ\
RES/)\n{\n	    \n	    $c= substr($l,11,1);\n	    i\
f ($hc_chain {$c}){print \"$l\";}\n}\n	elsif ( $l=\
~/^ATOM/)\n{\n	    $c=substr($l,21,1);\n	    if ($\
hc_chain {$c}){print \"$l\";}\n}\n}\n    close ( $\
INFILE);\n    &myexit($EXIT_SUCCESS);\n}    \n\nif\
 ( $MODE eq \"pdb\")\n{\n\n    $read_header=0;\n  \
  while (<$INFILE>) \n{\n	    $line=$_;\n	    if  \
  ($line =~ /^HEADER/){print \"$line\";$read_heade\
r=1;}\n}\n    close ($INFILE);\n\n    if (!$read_h\
eader)\n{\n	print \"HEADER    UNKNOWN             \
                    00-JAN-00   $force_name\\n\";\\
n}\n\n    $INFILE=vfopen ($structure_file, \"r\");\
\n    \n    print \"COMPND   1 CHAIN:\";\n    $las\
t=pop(@c_chain);\n    foreach $c ( @c_chain){ prin\
t \" $c,\";}\n    if ( $last eq \" \"){print \" NU\
LL;\\n\";}\n    else \n{\n      print \" $last;\\n\
\";\n}\n    @c_chain=(@c_chain, $last);\n    \n   \
 print \"REMARK Output of the program extract_from\
_pdb (Version $VersionTag)\\n\";\n    print \"REMA\
RK Legal PDB format not Guaranteed\\n\";\n    prin\
t \"REMARK This format is not meant to be used in \
place of the PDB format\\n\";\n    print \"REMARK \
The header refers to the original entry\\n\";\n   \
 print \"REMARK The sequence from the original fil\
e has been taken in the field: $seq_field\\n\";\n \
   print \"REMARK extract_from_pdb, 2001, 2002, 20\
03, 2004, 2005 2006 (c) CNRS and Cedric Notredame\\
\n\";   \n    if ( $coor_set)\n{\n       print \"R\
EMARK Partial chain: Start $start End $end\\n\";\n\
}\n    if ( $is_nmr)\n{\n       print \"REMARK NMR\
 structure: MODEL $model\\n\";\n}\n   \n    if ( $\
n_atom!=0)\n{\n       print \"REMARK Contains Coor\
dinates of: \";\n       foreach $a (@atom){print \\
"$a \";}\n       print \"\\n\";\n}  \n}\n\n\n\n\nm\
y $residue_index = -999;\nmy $old_c = \"TemporaryC\
hain\";\n\nwhile (<$INFILE>) \n{\n	$line=$_;\n\n\n\
	if ($line =~ /^SEQRES/)\n{\n\n		@field=/(\\S*)\\s\
*/g;\n\n		$c= substr($_,11,1);\n\n		\n		$l=$#field\
;\n		for ($a=4; $a<$#field ;)\n{\n			if (!$onelett\
{$molecule_type{$c}}->{$field[$a]})\n{\n				splice\
 @field, $a, 1;\n}\n			else \n{\n				$a++;\n}\n}\n\
	\n		if ( $c ne $in_chain)\n{\n			$pdb_chain_list[\
$n_pdb_chains]=$c;\n			$pdb_chain_len [$n_pdb_chai\
ns]=$len;\n			$in_chain=$c;\n			$n_pdb_chains++;\n\
}\n	\n		for ( $a=4; $a<$#field;$a++)\n{\n			$compl\
ete_seq{$c}[$complete_seq_len{$c}++]=$field[$a];\n\
}\n}\n    elsif ( $line=~/^ATOM/ || ($line=~/^HETA\
TM/ && &is_aa(substr($line,17,3),substr($line,21,1\
)) && !$no_hetatm))\n{\n\n	 \n    $RAW_AT_ID=$AT_I\
D=substr($line,12,4);\n	$RES_ID=&is_aa(substr($lin\
e,17,3),substr($line,21,1));\n	$CHAIN=substr($line\
,21,1);\n\n    $RES_NO=substr($line,22,4);\n	$HOM_\
CODE=substr ($line, 26, 1);\n	$TEMP=substr($line,6\
0,6);\n	\n	$TEMP=~s/\\s//g;\n        $AT_ID=~s/\\s\
//g;\n	$RES_ID=~s/\\s//g;\n        $RES_NO=~s/\\s/\
/g;\n		\n	if ( $HOM_CODE ne $MAIN_HOM_CODE){next;}\
\n	elsif ( $already_read2{$CHAIN}{$RES_ID}{$AT_ID}\
{$RES_NO}){next;}\n	else{$already_read2{$CHAIN}{$R\
ES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	\n	if ($coor_set \
&& $numbering eq \"file\" && $residue_index ne $RE\
S_NO)\n{\n	    \n	    if ( $RES_NO<=$start){$real_\
start{$CHAIN}++;}\n	    if ( $RES_NO<=$end){$real_\
end{$CHAIN}++;}\n}\n	elsif ($numbering eq \"absolu\
te\")\n{\n	    $real_start{$CHAIN}=$start;\n	    $\
real_end{$CHAIN}=$end;\n}\n\n        $KEY=\"ALL\";\
\n        if ( $CHAIN ne $in_atom_chain)\n{\n	    \
\n	  $pdb_atom_chain_list[$n_pdb_atom_chains]=$c;\\
n	  $pdb_atom_chain_len [$n_pdb_atom_chains]=$len;\
\n	  $in_atom_chain=$c;\n	  $n_pdb_atom_chains++;\\
n}\n	\n	if ( $residue_index ne $RES_NO)\n{\n	     \
$residue_index = $RES_NO;\n	     $atom_seq{$CHAIN}\
[$atom_seq_len{$CHAIN}++]=$RES_ID;;\n}\n}\n\n}\ncl\
ose ($INFILE);\n\n\n\n\n\n\n$INFILE=vfopen ($struc\
ture_file, \"r\");\nforeach $c (@c_chain)\n{\n\n	i\
f    ( $seq_field eq \"SEQRES\"){@pdb_seq=@{$compl\
ete_seq{$c}};}\n	elsif ( $seq_field eq \"ATOM\")  \
{@pdb_seq=@{$atom_seq{$c}};}\n	\n\n	$full_length=$\
l=$#pdb_seq+1;\n		\n	if ( $real_end{$c}==\"*\"){$r\
eal_end{$c}=$full_length;}\n	if ( $coor_set)\n{	  \
 \n\n	   if ( $real_end{$c} < $l){splice @pdb_seq,\
 $real_end{$c}, $l;}\n	   if ( $real_start{$c} < $\
l){splice @pdb_seq, 0, $real_start{$c}-1;}	  	   \\
n	   $l=$#pdb_seq;\n}\n\n	elsif ( $delete_set)\n{\\
n	   splice @pdb_seq, $delete_start, $delete_end-$\
delete_start+1;\n	   $l=$#pdb_seq;\n}\n	\n	$new_fa\
sta_name=\"$pdb_id$c\";\n	if ( $coor_set)\n{\n	   \
if ( $n_pdb_chains==0){$new_fasta_name=\"$new_fast\
a_name$c\";}\n	   $new_fasta_name= $new_fasta_name\
.\"\\_$start\\_$end\";\n}\n	   \n	if ( $MODE eq \"\
pdb\")\n{\n	   $nl=1;\n	   $n=0;\n	   \n	   foreac\
h $res ( @pdb_seq)\n		{\n		if ( !$n)\n		{\n		\n		 \
printf \"SEQRES %3d %1s %4d  \", $nl,$c, $l;\n		 $\
nl++;\n	}\n	     $res=~s/\\s//g;\n	     \n	     if\
 ($code==1){ printf \"%3s \",$onelett{$molecule_ty\
pe{$c}}->{$res};}\n	     elsif  ($code==3){ printf\
 \"%3s \",$res};\n	     \n	     $n++;		  \n	     i\
f ( $n==13){$n=0;print \"\\n\";}\n}\n	  if ( $n!=0\
){print \"\\n\"; $n=0;}\n}\n	elsif ( $MODE eq \"si\
mple\")\n{\n	  print \"# SIMPLE_PDB_FORMAT\\n\";\n\
	  if ( $new_fasta_name eq \" \"){$new_fasta_name=\
\"dummy_name\";}\n	  print \">$new_fasta_name\\n\"\
;\n\n	  foreach $res ( @pdb_seq)\n{\n	      print \
\"$onelett{$molecule_type{$c}}->{$res}\";\n}\n	  p\
rint \"\\n\";\n}\n	elsif ( $MODE eq \"fasta\")\n{\\
n	  $n=0;\n	  print \">$new_fasta_name\\n\";\n	  \\
n	  foreach $res ( @pdb_seq)\n{\n\n	      print \"\
$onelett{$molecule_type{$c}}->{$res}\";\n         \
     $n++;\n	      if ( $n==60){print \"\\n\"; $n=\
0;}\n}\n	  print \"\\n\"; \n}\n}\n\nif ( $MODE eq \
\"fasta\")\n{\n     &myexit($EXIT_SUCCESS);\n  \n}\
\n\n  \n  $charcount=0;\n  $inchain=\"BEGIN\";\n  \
$n=0;\n  while (<$INFILE>) \n{\n    $line=$_;\n   \
  \n    if ($line =~/^ATOM/  ||  ($line=~/^HETATM/\
))\n{\n	$line_header=\"UNKNWN\";\n	$RES_ID=substr(\
$line,17,3);\n	$chain = substr($line,21,1);\n\n	if\
 ($line =~/^ATOM/)\n{\n	    $line_header=\"ATOM\";\
\n	    $RES_ID=(&is_aa($RES_ID,$chain))?&is_aa($RE\
S_ID,$chain):$RES_ID;\n}\n	elsif ($line=~/^HETATM/\
 && ($ligand_list {$RES_ID} ||$ligand_list {'ALL'}\
 || !&is_aa($RES_ID,$chain)))\n{\n	    $line_heade\
r=\"HETATM\";\n}\n	elsif ($line=~/^HETATM/ && (&is\
_aa($RES_ID,$chain) && !$no_hetatm))\n{\n	    $lin\
e_header=\"ATOM\";\n	    $RES_ID=&is_aa($RES_ID,$c\
hain);\n}\n	else\n{\n	    next;\n}\n\n	\n\n	$X=sub\
str($line,30,8);     \n	$Y=substr($line,38,8);\n	$\
Z=substr($line,46,8);\n	$TEMP=substr($line,60,6);\\
n	\n	$RAW_AT_ID=$AT_ID=substr($line,12,4);\n	$CHAI\
N=substr($line,21,1);\n	$RES_NO=substr($line,22,4)\
;\n	$HOM_CODE=substr ($line, 26, 1);\n	\n	$X=~s/\\\
s//g;\n	$Y=~s/\\s//g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\\
\s//g;\n	\n	$AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g;\
\n	$RES_NO=~s/\\s//g;\n\n	\n	if ( $HOM_CODE ne $MA\
IN_HOM_CODE){next;}\n	elsif ( $already_read{$CHAIN\
}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n	else{$alread\
y_read{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	\
$KEY=\"ALL\";\n\n      	if ( $RES_NO ==0){$start_a\
t_zero=1;}\n\n	$RES_NO+=$start_at_zero;    \n	\n	i\
f ( $current_chain ne $CHAIN)\n{\n	    $current_ch\
ain=$CHAIN;\n	    $pos=$current_residue=0;\n	    $\
offset=($coor_set)?($real_start{$CHAIN}-1):0;\n	  \
  if    ( $seq_field eq \"SEQRES\"){@ref_seq=@{$co\
mplete_seq{$CHAIN}};}\n	    elsif ( $seq_field eq \
\"ATOM\")  {@ref_seq=@{$atom_seq{$CHAIN}};}\n}\n	\\
n	if ($current_residue != $RES_NO)\n{\n	    $curre\
nt_residue=$RES_NO;\n	    if    ( $seq_field eq \"\
SEQRES\"){$pos=$current_residue;}\n	    elsif ( $s\
eq_field eq \"ATOM\"){$pos++;}\n}\n	\n	\n	if ($n_a\
tom==0 || $atom_list{$AT_ID}==1 || $atom_list{$KEY\
}==1)\n{ 	\n	    \n	    $do_it=(!@c_chain || $hc_c\
hain{$CHAIN} ||$hc_chain{'LIGAND'} );\n	    \n	   \
 $do_it= ($do_it==1) && ($coor_set==0 ||($pos>=$re\
al_start{$CHAIN} && $pos<=$real_end{$CHAIN}));\n	 \
   $do_it= ($do_it==1) && ($delete_set==0 || $pos<\
$delete_start ||$pos>$delete_end );\n	    if ($lig\
and==0 && $line_header eq \"HETATM\" ){$do_it=0;}\\
n	    if ($ligand_only==1 && $line_header eq \"ATO\
M\" ){$do_it=0;}\n	    if ($ligand==1 && $line_hea\
der eq \"HETATM\" && $ligand_list{$RES_ID}==0 && $\
ligand_list{\"ALL\"}==0){$do_it=0;} \n	    \n	    \
\n	    if ( $do_it)\n{\n		$n++;\n		$out_pos=$pos;\\
n		\n	       if ( $delete_set)\n{\n		  if ( $out_p\
os< $delete_start){;}\n		  else {$offset=$delete_e\
nd-$delete_start;}\n}       \n	       \n	       if\
 ( $numbering_out eq \"new\"){$out_pos-=$offset;}\\
n	       elsif ( $numbering_out eq \"old\"){$out_p\
os=$RES_NO;}\n	       \n       \n	       \n	      \
 if ( $code==1){$RES_ID=$onelett{$molecule_type{$c\
}}->{$RES_ID};}\n	    \n	       if ($unfold)\n{\n	\
	   $unfolded_x+=5;\n		   $X=$unfolded_x;\n		   $Y\
=0;\n		   $Z=0;\n		   $float=1;\n}\n	       else\n\
{\n		   $float=3;\n}\n\n	       if ( $MODE eq \"pd\
b\")\n{\n		   printf \"%-6s%5d %-4s %3s %s%4d    %\
8.3f%8.3f%8.3f  1.00 %5.2f\\n\",$line_header, $n, \
$RAW_AT_ID,$RES_ID,$CHAIN,$out_pos, $X, $Y, $Z,$TE\
MP;		  \n}\n	       elsif ( $MODE eq \"simple\")\n\
{\n		    if ( $RES_ID eq \"\"){$RES_ID=\"X\";}\n		\
  printf \"%-6s %5s %s %2s %4d    %8.3f %8.3f %8.3\
f\\n\",$line_header, $AT_ID, $RES_ID,($CHAIN eq\"\\
" || $CHAIN eq \" \")?\"A\":$CHAIN,$out_pos, $X, $\
Y, $Z,$TEMP;\n}\n\n}\n}\n}\n}\nprint \"\\n\";\nclo\
se($INFILE);\n\n\nif ( $error ne \"\") \n{$error=$\
error.\"\\nDiagnostic:    SEQRES and the residues \
in ATOM are probably Incompatible\\n\";\n    $erro\
r=$error.  \"Recomendation: Rerun with '-fix 1' in\
 order to ignore the SEQRES sequences\\n\";\n}\nif\
 (!$nodiagnostic){print STDERR $error;}\n&myexit (\
 $EXIT_SUCCESS);\n\nsub get_pdb_entry_type_file\n \
 {\n    my $cache_file=\"$cache/pdb_entry_type.txt\
\";\n    my $env_file  = $ENV{\"PDB_ENTRY_TYPE_FIL\
E\"};\n    my $pdb_file  =\"$ENV{'PDB_DIR'}/derive\
d_data/pdb_entry_type.txt\";\n    \n    \n    if (\
-z $cache_file){unlink ($cache_file);}#will get up\
dated\n    if (-z $env_file){$env_file=\"\";}    #\
cannot update\n    if (-z $pdb_file){$pdb_file=\"\\
";}    #cannot update\n    \n    if    (-e $env_fi\
le){return $env_file;} #env wins: user decides\n  \
  elsif (-e $pdb_file){return $pdb_file;} #local d\
atabase wins: network file may be out of sync\n   \
 elsif ($no_remote_pdb_dir==1)\n      {\n	if (-e $\
cache_file){return $cache_file;}\n	else\n	  {add_w\
arning($$,$$,\"PDB_ENTRY_TYPE_FILE must be set to \
the location of <pdb>/derived_data/pdb_entry_type.\
txt when using NO_REMOTE_PDB_DIR=1\");\n	   return\
 \"\";\n	 }\n      }\n    else #update can only ta\
ke place if the file lives in cache\n      {\n	my \
$new_file;\n	if (!-e $cache_file || (-M $cache_fil\
e)>1)\n	  {\n	    $new_file=vtmpnam();\n	    &url2\
file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_data/pd\
b_entry_type.txt\", $new_file);\n	    if ( !-z $ne\
w_file){system (\"mv $new_file $cache_file\"); unl\
ink ($new_file); $new_file=$cache_file;}\n	    els\
e {unlink($new_file);}\n	  }\n	else\n	  {\n	    $n\
ew_file=$cache_file;\n	  }\n	\n	if (!-e $cache_fil\
e && !-e $new_file)\n	  {\n	    add_warning($$,$$,\
\"Could not download ftp://ftp.wwpdb.org/pub/pdb/d\
erived_data/pdb_entry_type.txt\");\n	    return \"\
\";\n	  }\n	elsif (-e $cache_file && !-e $new_file\
)\n	  {\n	    my $m=(-M $cache_file);\n	    add_wa\
rning($$,$$,\"Could not update file ftp://ftp.wwpd\
b.org/pub/pdb/derived_data/pdb_entry_type.txt. Old\
er Version [$cache_file]($m Month(s) old) will be \
used instead\");\n	    return $cache_file;\n	  }\n\
	else\n	  {\n	    return $new_file;\n	  }\n      }\
\n  }\n\n\n\nsub get_unrealeased_file\n  {\n    my\
 $cache_file=\"$cache/unreleased_entries.json.gz\"\
;\n    my $env_file  = $ENV{\"PDB_UNREALEASED_FILE\
\"};\n    my $pdb_file  =\"$ENV{'PDB_DIR'}/derived\
_data/unreleased_entries.json.gz\";\n    \n    \n \
   if ($env_file eq \"NO\" || $env_file eq \"No\" \
|| $env_file eq \"no\" || $env_file eq \"0\"){retu\
rn \"NO\";}\n\n    if (-z $cache_file){unlink ($ca\
che_file);}#will get updated\n    if (-z $env_file\
){unlink($env_file);}     #will update\n    if (-z\
 $pdb_file){$pdb_file=\"\";}          #cannot upda\
te\n    \n    if    (-e $env_file){return $env_fil\
e;} #env wins: user decides\n    elsif (-e $pdb_fi\
le){return $pdb_file;} #local database wins: netwo\
rk file may be out of sync\n    elsif ($no_remote_\
pdb_dir==1)        \n      {\n	if (-e $cache_file)\
{return $cache_file;}\n	elsif ( $env_file && ! -e \
$env_file)\n	  {\n	    &url2file(\"https://ftp.rcs\
b.org/pub/pdb/holdings/unreleased_entries.json.gz\\
",$env_file);\n	    if ( -e $env_file && !-z $env_\
file){return $env_file;}\n	  }\n	else\n	  {\n	    \
add_warning($$,$$,\"UNREALEASED_FILE must be set t\
o the location of your unreleased_entries.json.gz \
file as downloaded from https://ftp.rcsb.org/pub/p\
db/holdings/unreleased_entries.json.gz when using \
NO_REMOTE_PDB_DIR=1\");\n	    return \"\";\n	  }\n\
      }\n    else #update can only take place if t\
he file lives in cache\n      {\n	my $new_file=vtm\
pnam ();\n	if (!-e $cache_file || (-M $cache_file)\
>1)\n	  {\n	    &url2file(\"https://ftp.rcsb.org/p\
ub/pdb/holdings/unreleased_entries.json.gz\",$new_\
file);\n	    if ( !-z $new_file){system (\"mv $new\
_file $cache_file\"); unlink ($new_file); $new_fil\
e=$cache_file;}\n	    else {unlink($new_file);}\n	\
  }\n	else\n	  {\n	    $new_file=$cache_file;\n	  \
}\n	\n	if (!-e $cache_file && !-e $new_file)\n	  {\
\n	    add_warning($$,$$,\"Could not download http\
s://ftp.rcsb.org/pub/pdb/holdings/unreleased_entri\
es.json.gz\");\n	    return \"\";\n	  }\n	elsif (-\
e $cache_file && !-e $new_file)\n	  {\n	    my $m=\
(-M $cache_file);\n	    add_warning($$,$$,\"Could \
not update file https://ftp.rcsb.org/pub/pdb/holdi\
ngs/unreleased_entries.json.gz. Older Version [$ca\
che_file]($m Month(s) ) will be used\");\n	    ret\
urn $cache_file;\n	  }\n	else\n	  {\n	    return $\
new_file;\n	  }\n      }\n  }\n\nsub is_released \\
n  {\n    my ($r);\n    my $in=@_[0];\n    my $nam\
e=&remote_is_pdb_name ($in);\n    my $hold=&remote\
_is_on_hold($in);\n    \n    $r=($name && !$hold)?\
1:0;\n    return $r;\n  }\n\nsub remote_is_pdb_nam\
e \n  {\n    my $in=@_[0];\n    my ($pdb);\n    my\
 ($value,$value1,$value2);\n    my $max=2;\n    \n\
    \n    \n    my $ref_file=&get_pdb_entry_type_f\
ile();\n    \n    if ( $in=~/[^\\w\\d\\:\\_]/){ret\
urn 0;}\n    elsif (!-e $ref_file)\n      {\n	add_\
warning ($$,$$,\"Cannot find pdb_entry_type.txt;  \
$in is assumed to be valid; add ftp://ftp.wwpdb.or\
g/pub/pdb/derived_data/pdb_entry_type.txt in $cach\
e to automatically check name status\");\n	return \
1;\n      }\n    else\n      {\n	$pdb=substr ($in,\
0, 4);\n	chomp(($value1=`grep -c $pdb $ref_file`))\
;\n	$pdb=lc($pdb);\n	chomp(($value2=`grep -c $pdb \
$ref_file`));\n	$value=($value1 || $value2)?1:0;\n\
	$value=($value>0)?1:0;\n	\n	return $value;\n     \
 }\n  }\n\n\n\nsub pdb2model_type\n{\n    my $in=@\
_[0];\n    my ($ref_file, $pdb);\n    my ($value, \
$ret);\n\n    if ( $in=~/[^\\w\\d\\:\\_]/){return \
0;}\n    $ref_file=&get_pdb_entry_type_file();\n  \
  if (!-e $ref_file)\n      {\n	add_warning ($$,$$\
,\"Cannot find pdb_entry_type.txt;  $in is assumed\
 to be diffraction; add ftp://ftp.wwpdb.org/pub/pd\
b/derived_data/pdb_entry_type.txt in $cache to che\
ck name status\");\n	return \"diffraction\";\n    \
  }\n    else\n      {\n	$pdb=substr ($in,0, 4);\n\
	$pdb=lc($pdb);\n	\n	chomp(($value=`grep $pdb $ref\
_file`));\n	\n	$value=~/^\\S+\\s+\\S+\\s+(\\S+)/;\\
n	$ret=$1;\n	if ( $ret eq\"\"){return \"UNKNOWN\";\
}\n	\n	return $ret;\n      }\n  }\nsub remote_is_o\
n_hold\n  {\n    my $in=@_[0];\n    my ($ref_file,\
 $pdb);\n    my ($value1, $value2,$value);\n    \n\
\n\n    \n    $ref_file=&get_unrealeased_file();\n\
    if ($ref_file eq \"NO\"){return 0;}\n\n\n    i\
f ($no_remote_pdb==1){return 0;}\n    if ( $in=~/[\
^\\w\\d\\:\\_]/){return 0;}\n    \n    $ref_file=&\
get_unrealeased_file();\n    if (!-e $ref_file)\n \
     {\n	add_warning ($$,$$,\"Cannot find https://\
ftp.rcsb.org/pub/pdb/holdings/unreleased_entries.j\
son.gz;  $in is assumed to be released;\");\n	retu\
rn 1;\n      }\n    \n    $pdb=substr ($in,0, 4);\\
n    chomp(($value1=`grep -c $pdb $ref_file`));\n \
   $pdb=lc($pdb);\n    chomp(($value2=`grep -c $pd\
b $ref_file`));\n    $value=($value1 || $value2)?1\
:0;\n    $value=($value>0)?1:0;\n    return $value\
;\n  }\n\nsub is_pdb_file\n  {\n    my @arg=@_;\n \
   \n    if ( !-e $arg[0]){return 0;}\n    \n    $\
F=vfopen ($arg[0], \"r\");\n    while ( <$F>)\n   \
   {\n	if (/^HEADER/)\n	  {\n	    close $F;\n	    \
return 1;\n	  }\n	elsif ( /^SEQRES/)\n	  {\n	    c\
lose $F;\n	    return 1;\n	  }\n	elsif ( /^ATOM/)\\
n	  {\n	    close $F;\n	    return 1;\n	  }\n     \
 }\n    return 0;\n  }\nsub get_pdb_id\n{\n    my \
$header_file=@_[0];\n    my $id;\n    my $F= new F\
ileHandle;\n    \n    \n    $F=vfopen (\"$header_f\
ile\", \"r\");\n\n    while ( <$F>)\n      {\n	if \
( /HEADER/)\n	  {\n	    if ($debug){print \"$_\";}\
\n	    $id=substr($_,62,4 );\n	    return $id;\n	 \
 }\n      }\n    close ($F);\n    \n    return \"\\
";\n}\n\nsub get_ligand_list\n{\n    my $pdb_file=\
@_[0];\n    my $chain;\n    my $ligand;\n    my %c\
omplete_ligand_list;\n    \n\n    $F=vfopen ($pdb_\
file, \"r\");\n    while ( <$F>)\n{\n	if ( /^HETAT\
M/)\n{\n	    $line=$_;\n	    $chain=substr($line,2\
1,1);\n	    $ligand=substr($line,17,3);\n	    \n	 \
   if (!$complete_ligand_list{$chain}{$ligand})\n{\
\n		\n		$complete_ligand_list{\"result\"}.=\"CHAIN\
 $chain LIGAND $ligand\\n\";\n		$complete_ligand_l\
ist{$chain}{$ligand}=1;\n}\n}\n}\n    close ($F);\\
n    return %complete_ligand_list;\n}\n\nsub get_c\
hain_list \n{\n    my $header_file;\n    my @chain\
_list;\n    my @list;\n    my $n_chains;\n    my %\
chain_hasch;\n    my $pdb_file=@_[0];\n    my $c;\\
n    my %hasch;\n    my $chain;\n  \n    \n    $F=\
vfopen ($pdb_file, \"r\");\n    while ( <$F>)\n{\n\
\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+)/)\n	  {\n	    $\
chain = substr($_,11,1);$chain=~s/\\s//g;if ( $cha\
in eq \"\"){$chain=\" \";}\n	    if (!$hasch{$chai\
n}){$hasch{$chain}=1;push @chain_list, $chain;}\n	\
  }\n	if (/^ATOM/ || /^HETATM/)\n	  {\n	    $chain\
 = substr($_,21,1); $chain=~s/\\s//g;if ( $chain e\
q \"\"){$chain=\" \";}\n	    if (!$hasch{$chain}){\
$hasch{$chain}=1;push @chain_list, $chain;}\n	  }\\
n      }\n\n\nclose ($F);\nif (!@chain_list)\n  {\\
n    @chain_list=(\"A\");\n  }\n\n\nreturn @chain_\
list;\n}\n\nsub token_is_in_list\n{\n\n    my @lis\
t=@_;\n    my $a;\n    \n    for ($a=1; $a<=$#list\
; $a++)\n{\n	if ( $list[$a] eq $list[0]){return $a\
;}\n}\n}\n\nsub pdb_name2name_and_chain \n{\n    m\
y $pdb_file=@_[0];\n    my $pdb_file_in;\n    my @\
array;\n    my $chain;\n    my $c;\n\n    $pdb_fil\
e_in=$pdb_file;\n\n    $pdb_file=~/^(.{4})/;$pdb_i\
d=$1;\n    @array=($pdb_file=~/([\\w])/g);\n  \n  \
\n    $chain=uc ($array[4]);\n    $chain=($chain e\
q \"\")?\"FIRST\":$chain;\n    \n    return ( $pdb\
_id, $chain);\n\n    if ( $#array==3){return ($pdb\
_id, \"FIRST\");}\n    elsif ( $#array<4){ return \
($pdb_id, \"\");}\n    else {return ( $pdb_id, $ch\
ain);}\n      \n    \n    \n}\nsub get_main_hom_co\
de \n{\n    my $pdb_file=@_[0];\n    my %hom, $n, \
$best, $best_h;\n    open (F, $pdb_file);\n    whi\
le (<F>)\n{\n	if ( $_=~/^ATOM/)\n{\n	    $h=substr\
 ($_,26, 1);\n	    $n=++$hom{$h};\n	    if ($n>$be\
st)\n{\n		$best=$n;\n		$best_h=$h;\n}\n}\n}\n    c\
lose (F);\n    return $best_h;\n}\n\n\nsub get_pdb\
_file \n{\n    my ($pdb_file_in)=(@_);\n    my $re\
sult;\n    my @letter;\n    my @chain;\n    my $v;\
\n    my $pdb_file=$pdb_file_in;\n\n    $pdb_file=\
($pdb_file_in=~/\\S+_S_(\\S+)/)?$1:$pdb_file_in;\n\
    \n    if ($no_remote_pdb_dir==0)\n      {\n	$n\
o_remote_pdb_dir=1;\n	$result=get_pdb_file3 ($pdb_\
file);\n	$no_remote_pdb_dir=0;\n	if ( $result){ret\
urn $result;}\n	else\n	  {\n	    \n	    lc ($pdb_f\
ile);\n	    $result=get_pdb_file3($pdb_file);\n	  \
  return  $result;\n	  }\n      }\n    else\n     \
 {\n	return get_pdb_file3 ($pdb_file);\n      }\n \
   \n  }\n\nsub get_pdb_file3 \n{\n    my $pdb_fil\
e_in=@_[0];\n    my $result;\n    my @letter;\n   \
 my @chain;\n    my $lcfile;\n    my $ucfile;\n   \
 my $pdb_file=$pdb_file_in;\n    \n    $lcfile=lc \
$pdb_file;\n    $ucfile=uc $pdb_file;\n\n    if ( \
($result=get_pdb_file2 ($pdb_file))){return $resul\
t;}\n    \n\n    if ($lcfile ne $pdb_file && ($res\
ult=get_pdb_file2 ($lcfile))){return $result;}\n  \
  if ($ucfile ne $pdb_file && ($result=get_pdb_fil\
e2 ($ucfile))){return $result;}\n    \n   \n    \n\
    return \"\";\n}\nsub get_pdb_file2\n{\n    my \
$pdb_file=@_[0];\n    my $return_value;\n    \n   \
 $return_value=\"\";\n    \n    if ( ($result=get_\
pdb_file1 ($pdb_file))){$return_value=$result;}\n \
   elsif ( !($pdb_file=~/\\.pdb/) && !($pdb_file=~\
/\\.PDB/))\n{\n	if ( ($result=get_pdb_file1 (\"$pd\
b_file.pdb\"))){$return_value=$result;}\n	elsif ( \
($result=get_pdb_file1 (\"$pdb_file.PDB\"))){$retu\
rn_value=$result;}\n\n	elsif ( ($result=get_pdb_fi\
le1 (\"pdb$pdb_file.pdb\"))){$return_value=$result\
;}	\n	elsif ( ($result=get_pdb_file1 (\"pdb$pdb_fi\
le.PDB\"))){$return_value=$result;}\n	elsif ( ($re\
sult=get_pdb_file1 (\"PDB$pdb_file.PDB\"))){$retur\
n_value=$result;}\n	elsif ( ($result=get_pdb_file1\
 (\"PDB$pdb_file.pdb\"))){$return_value=$result;}\\
n	\n	\n	elsif ( ($result=get_pdb_file1 (\"$pdb_fil\
e.ent\"))){$return_value=$result;}\n	elsif ( ($res\
ult=get_pdb_file1 (\"pdb$pdb_file.ent\"))){$return\
_value=$result;}\n	elsif ( ($result=get_pdb_file1 \
(\"PDB$pdb_file.ent\"))){$return_value=$result;}\n\
\n	elsif ( ($result=get_pdb_file1 (\"$pdb_file.ENT\
\"))){$return_value=$result;}\n	elsif ( ($result=g\
et_pdb_file1 (\"pdb$pdb_file.ENT\"))){$return_valu\
e=$result;}\n	elsif ( ($result=get_pdb_file1 (\"PD\
B$pdb_file.ENT\"))){$return_value=$result;}\n	\n	\\
n	\n}\n    return $return_value;\n}\n    \nsub get\
_pdb_file1\n{\n    my ($pdb_file)=(@_);\n    my $r\
eturn_value;\n    \n\n    $return_value=\"\";\n   \
 if ( ($result=get_pdb_file0 ($pdb_file))){$return\
_value=$result;}\n    elsif ( ($result=get_pdb_fil\
e0 (\"$pdb_file.Z\"))){$return_value=$result;}\n  \
  elsif ( ($result=get_pdb_file0 (\"$pdb_file.gz\"\
))){$return_value=$result;}\n    elsif ( ($result=\
get_pdb_file0 (\"$pdb_file.GZ\"))){$return_value=$\
result;}\n    return $return_value;\n}\nsub get_pd\
b_file0 \n{ \n    my ($pdb_file, $attempt)=(@_);\n\
    my $pdb_file=@_[0];\n    my $tmp_pdb_file;    \
\n    my $return_value;\n\n    if ( !$attempt){$at\
tempt=1;}\n    \n    $local_pdb_file=\"$pdb_file\"\
;\n    if ( $local_pdb_file eq \"\")\n{\n	$tmp_pdb\
_file=vtmpnam();\n	open F, \">$tmp_pdb_file\";\n	\\
n	while (<STDIN>){print F \"$_\";}\n	close (F);\n	\
\n	if (-e $tmp_pdb_file && &is_pdb_file ( $local_p\
db_file))\n{return $tmp_pdb_file;}\n}\n\n    $loca\
l_pdb_file=\"$pdb_file\";\n    &debug_print (\"\\n\
Try access local file: $local_pdb_file\");\n    \n\
    $local_pdb_file=&check_pdb_file4compression ($\
local_pdb_file);\n    if ( -e $local_pdb_file && (\
&is_pdb_file ($local_pdb_file) || $force_pdb))\n{\\
n	&debug_print ( \"\\n\\tIs in Current Dir\");\n	$\
tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_file $tmp\
_pdb_file`;\n	return $tmp_pdb_file;\n}\n    else\n\
{\n	&debug_print (\"\\n\\tFile Not in Current Dir\\
");\n}\n\n    if ($pdb_file=~/^pdb/||$pdb_file=~/^\
PDB/){$pdb_div=substr ($pdb_file, 4, 2);}\n    els\
e\n{\n	  $pdb_div=substr ($pdb_file, 1, 2);\n}\n  \
  $local_pdb_file=\"$pdb_dir/$pdb_div/$pdb_file\";\
\n    $local_pdb_file=&check_pdb_file4compression \
( $local_pdb_file);\n    &debug_print (\"\\nTry ac\
cess file From PDB_DIR: $local_pdb_file\");\n    i\
f ($pdb_dir && -e $local_pdb_file && &is_pdb_file \
($local_pdb_file))\n{\n	&debug_print ( \"\\n\\tIs \
in Local PDB DIR\");\n	$tmp_pdb_file=vtmpnam();\n	\
`cp $local_pdb_file $tmp_pdb_file`;\n	return $tmp_\
pdb_file;\n}\n\n    $local_pdb_file=\"$pdb_dir/$pd\
b_file\";\n    $local_pdb_file=&check_pdb_file4com\
pression ( $local_pdb_file);\n    &debug_print (\"\
\\nTry access file From PDB_DIR: local_pdb_file\")\
;\n    if ($pdb_dir && -e $local_pdb_file && &is_p\
db_file ($local_pdb_file))\n{\n	&debug_print ( \"\\
\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_file=vtmpn\
am();\n	`cp $local_pdb_file $tmp_pdb_file`;\n	retu\
rn $tmp_pdb_file;\n}\n\n    $local_pdb_file=\"$pdb\
_dir$pdb_file\";\n    $local_pdb_file=&check_pdb_f\
ile4compression ( $local_pdb_file);\n    &debug_pr\
int (\"\\nTry access file From PDB_DIR: $local_pdb\
_file\");\n    if ($pdb_dir && -e $local_pdb_file \
&& &is_pdb_file ($local_pdb_file))\n{\n	&debug_pri\
nt ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_fi\
le=vtmpnam();\n	`cp $local_pdb_file $tmp_pdb_file`\
;\n	return $tmp_pdb_file;\n}\n    else\n{&debug_pr\
int ( \"\\n\\tNot In Local Pdb Dir\");}\n\n    if \
($cache ne \"NO\" && $cache ne \"no\")\n{\n\n	$loc\
al_pdb_file=\"$cache/$pdb_file\";\n	$local_pdb_fil\
e=&check_pdb_file4compression ( $local_pdb_file);\\
n	&debug_print(\"\\nTry access file From Cache: $l\
ocal_pdb_file\");\n	if (-e $local_pdb_file && &is_\
pdb_file ($local_pdb_file))\n{\n	    &debug_print \
( \"\\n\\tIs in T-Coffee Cache\");\n	    $tmp_pdb_\
file=vtmpnam();\n	    `cp $local_pdb_file $tmp_pdb\
_file`;\n	    return $tmp_pdb_file;\n}\n	else{&deb\
ug_print ( \"\\n\\tNot in Cache Dir\");}\n}\n\nif \
(!$no_remote_pdb_dir) \n  {\n    my $value=&is_rel\
eased ($pdb_file);\n    my $return_value=\"\";\n  \
  if ($value==1)\n      {\n	\n	&debug_print (\"\\n\
**************************************************\
***\\nTry Remote Access for $pdb_file\");\n	$tmp_p\
db_file=vtmpnam();\n	$netcommand=$netaddress;\n	$n\
etcommand=~s/%%/$pdb_file/g;\n	&url2file(\"$netcom\
mand\", \"$tmp_pdb_file.$netcompression\");\n	&deb\
ug_print(\"\\nREMOTE: $netcommand\\n\");\n	\n	$com\
pressed_tmp_file_name=\"$tmp_pdb_file.$netcompress\
ion\";\n	\n	if ($netcompression && -B $compressed_\
tmp_file_name && $attempt<5)\n	  {\n	    my $r;\n	\
    &debug_print (\"\\n\\tFile Found Remotely\");\\
n	    if (($r=safe_system ( \"$netcompression_pg $\
compressed_tmp_file_name\")!=$EXIT_SUCCESS) && $at\
tempts<5)\n	      {\n		&debug_print (\"\\n\\tPrope\
r Download Failed Try again\");\n		unlink $compres\
sed_tmp_file_name;\n		print \"\\nFailed to Downloa\
d $compressed_tmp_file_name. New Attempt $attempt/\
5\\n\";\n		return &get_pdb_file0($pdb_file, $attem\
pt+1);\n	      }\n	    elsif ($r== $EXIT_SUCCESS)\\
n	      {\n		&debug_print (\"\\n\\tProper Download\
 Succeeded \");\n		$return_value=$tmp_pdb_file;\n	\
      }\n	    else\n	      {\n		&debug_print (\"\\\
n\\tProper Download Failed \");\n		&debug_print (\\
"\\nFile Not Found Remotely\");\n		unlink $compres\
sed_tmp_file_name;\n	      }\n	  }\n	else\n	  {\n\\
n	    &debug_print (\"\\nFile Not Found Remotely\"\
);\n	    unlink $compressed_tmp_file_name;\n	  }\n\
	#Update cache if required\n	if ($cache ne \"no\" \
&& $cache ne \"update\" && -e $return_value)\n	  {\
\n	    `cp $return_value $cache/$pdb_file.pdb`;\n	\
    #`t_coffee -other_pg clean_cache.pl -file $pdb\
_file.pdb -dir $cache`;\n	  }\n      }\n    &debug\
_print (\"\\nRemote Download Finished\");\n    ret\
urn $return_value;\n  }\nreturn \"\";\n}\n\nsub ch\
eck_pdb_file4compression \n{\n    my $file=@_[0];\\
n    my $tmp;\n    my $r;\n    \n    $tmp=&vtmpnam\
();\n    if (-e $tmp){unlink $tmp;}\n    \n    $fi\
le=~s/\\/\\//\\//g;\n    if    (-B $file && ($file\
=~/\\.Z/)) {`cp $file $tmp.Z`;`rm $tmp`;`gunzip $t\
mp.Z $SILENT`;$r=$tmp;}\n    elsif (-B $file && ($\
file=~/\\.gz/)){`cp $file $tmp.gz`;`gunzip $tmp.gz\
 $SILENT`;return $r=$tmp;}\n    elsif (-B $file ){\
`cp $file $tmp.gz`;`gunzip $tmp.gz $SILENT`;$r=$tm\
p;}\n    elsif ( -e $file ) {$r= $file;}\n    elsi\
f ( -e \"$file.gz\" ){ `cp $file.gz $tmp.gz`;`gunz\
ip     $tmp.gz $SILENT`;$r=$tmp;}    \n    elsif (\
 -e \"$file.Z\") {`cp $file.Z  $tmp.Z`; `gunzip $t\
mp.Z $SILENT`;$r=$tmp;}\n    else  {$r= $file;}\n\\
n    if ( -e \"$tmp.Z\"){unlink \"$tmp.Z\";}\n    \
if ( -e \"$tmp.gz\"){unlink \"$tmp.gz\";}\n    \n \
   return $r;\n    \n}\n\n\n\n\n\n    \n\n\n\n\n\n\
\n\nsub vfopen \n{\n    my $file=@_[0];\n    my $m\
ode=@_[1];\n    my $tmp;\n    my $F = new FileHand\
le;\n    \n    \n    $tmp=$file;\n	\n    \n    if \
( $mode eq \"r\" && !-e $file){ myexit(flush_error\
 (\"Cannot open file $file\"));}\n    elsif ($mode\
 eq \"w\"){$tmp=\">$file\";}\n    elsif ($mode eq \
\"a\"){$tmp=\">>$file\";}\n    \n    \n    open ($\
F,$tmp);\n    return $F;\n}\nsub debug_print\n{\n \
   my $message =@_[0];\n    if ($debug){print STDE\
RR \"NO_REMOTE_PDB_DIR: $no_remote_pdb_dir - $mess\
age [DEBUG:extract_from_pdb]\";}\n    return;\n}\n\
sub is_aa \n{\n    my ($aa, $chain) =@_;\n\n    my\
 $one;\n    my $trhee;\n    \n    if ( $onelett{$m\
olecule_type{$chain}}->{$aa} eq 'X' || !$onelett{$\
molecule_type{$chain}}->{$aa} ){return '';}\n    e\
lse\n      {\n	$one=$onelett{$molecule_type{$chain\
}}->{$aa};\n\n	$three=$threelett{$molecule_type{$c\
hain}}->{$one};\n	\n\n	return $three;\n      }\n  \
}\n\n\n\n\n\nsub url2file\n{\n    my ($address, $o\
ut, $wget_arg, $curl_arg)=(@_);\n    my ($pg, $fla\
g, $r, $arg, $count);\n    \n    if (!$CONFIGURATI\
ON){&check_configuration (\"wget\", \"INTERNET\", \
\"gzip\");$CONFIGURATION=1;}\n    \n\n    if (&pg_\
is_installed (\"curl\")){$pg=\"curl\"; $flag=\"-o\\
";$arg=$curl_arg;}\n    elsif (&pg_is_installed (\\
"wget\"))   {$pg=\"wget\"; $flag=\"-O\";$arg=$wget\
_arg;}\n    \n    return safe_system (\"$pg $flag$\
out $address >/dev/null 2>/dev/null\");\n\n}\n\n\n\
\n\nsub pdbfile2chaintype\n  {\n    my $file=@_[0]\
;\n    my %ct;\n    my $F;\n    \n    $F=vfopen ($\
file, \"r\");\n    while (<$F>)\n      {\n	my $lin\
e=$_;\n	if ($line =~/^ATOM/)\n	  {\n	    my $C=sub\
str($line,21,1);\n	    if (!$ct{$C})\n	      {	\n	\
	my $r=substr($line,17,3);\n		$r=~s/\\s+//;\n		if \
(length ($r)==1){$ct{$C}=\"R\";}\n		elsif (length \
($r)==2){$ct{$C}=\"D\";}\n		elsif (length ($r)==3)\
{$ct{$C}=\"P\";}\n		else \n		  {\n		    myexit(flu\
sh_error(\"ERROR: Could not read RES_ID field in f\
ile $file\"));\n		  }\n	      }\n	  }\n      }\n  \
  close ($F);\n    return %ct;\n  }\n   \n    \n\n\
\n\nsub fill_threelett_RNA\n{\n\n	my %threelett=(\\
n	'A', '  A',\n	'T', '  T',\n	'U', '  U',\n	'C', '\
  C',\n	'G', '  G',\n	'I', '  I', #Inosine\n	);\n	\
\n	return %threelett;\n\n}\n\n\nsub fill_onelett_R\
NA\n{\n	my   %onelett=(\n	'  A' => 'A',\n	'  T' =>\
 'T',\n	'  U' => 'U',\n	'  C' => 'C',\n	'  G' => '\
G',\n	'CSL' => 'X',\n	'UMS' => 'X',\n	'  I' => 'I'\
,\n	'A' => 'A',\n	'T' => 'T',\n	'U' => 'U',\n	'C' \
=> 'C',\n	'G' => 'G',\n	'I' => 'I',\n	);\n\n	retur\
n %onelett;\n\n}\n\n\nsub fill_onelett_DNA\n{\n	my\
   %onelett=(\n	' DA', 'A',\n	' DT', 'T',\n	' DC',\
 'C',\n	' DG', 'G',\n	'DA', 'A',\n	'DT', 'T',\n	'D\
C', 'C',\n	'DG', 'G',\n	);\n\n	return %onelett;\n\\
n}\n\nsub fill_threelett_DNA\n{\n\n	my %threelett=\
(\n	'A', ' DA',\n	'T', ' DT',\n	'C', ' DC',\n	'G',\
 ' DG',\n	);\n\n	return %threelett;\n\n}\n\n\n\n\n\
sub fill_threelett_prot\n{  \n  my %threelett;\n\n\
  %threelett=(\n'A', 'ALA',\n'C', 'CYS',\n'D', 'AS\
P',\n'E', 'GLU',\n'F', 'PHE',\n'G', 'GLY',\n'H', '\
HIS',\n'I', 'ILE',\n'K', 'LYS',\n'L', 'LEU',\n'N',\
 'ASN',\n'M', 'MET',\n'P', 'PRO',\n'Q', 'GLN',\n'R\
', 'ARG',\n'S', 'SER',\n'T', 'THR',\n'V', 'VAL',\n\
'W', 'TRP',\n'Y', 'TYR',\n);\n\nreturn %threelett;\
\n\n\n}\n\nsub fill_onelett_prot\n{\n    my %onele\
tt;\n    \n    %onelett=(\n\n'10A', 'X',\n'11O', '\
X',\n'12A', 'X',\n'13P', 'X',\n'13R', 'X',\n'13S',\
 'X',\n'14W', 'X',\n'15P', 'X',\n'16A', 'X',\n'16G\
', 'X',\n'1AN', 'X',\n'1AP', 'X',\n'1AR', 'X',\n'1\
BH', 'X',\n'1BO', 'X',\n'1C5', 'X',\n'1CU', 'X',\n\
'1DA', 'X',\n'1GL', 'X',\n'1GN', 'X',\n'1IN', 'X',\
\n'1LU', 'L',\n'1MA', 'X',\n'1MC', 'X',\n'1MG', 'X\
',\n'1MZ', 'X',\n'1NA', 'X',\n'1NB', 'X',\n'1NI', \
'X',\n'1PA', 'A',\n'1PC', 'X',\n'1PE', 'X',\n'1PG'\
, 'X',\n'1PI', 'A',\n'1PM', 'X',\n'1PN', 'X',\n'1P\
U', 'X',\n'1PY', 'X',\n'1UN', 'X',\n'24T', 'X',\n'\
25T', 'X',\n'26P', 'X',\n'2AB', 'X',\n'2AM', 'X',\\
n'2AN', 'X',\n'2AP', 'X',\n'2AR', 'X',\n'2AS', 'D'\
,\n'2BL', 'X',\n'2BM', 'X',\n'2CP', 'X',\n'2DA', '\
X',\n'2DG', 'X',\n'2DP', 'X',\n'2DT', 'X',\n'2EP',\
 'X',\n'2EZ', 'X',\n'2FG', 'X',\n'2FL', 'X',\n'2FP\
', 'X',\n'2FU', 'X',\n'2GL', 'X',\n'2GP', 'X',\n'2\
HP', 'X',\n'2IB', 'X',\n'2IP', 'X',\n'2LU', 'L',\n\
'2MA', 'X',\n'2MD', 'X',\n'2ME', 'X',\n'2MG', 'X',\
\n'2ML', 'L',\n'2MO', 'X',\n'2MR', 'R',\n'2MU', 'X\
',\n'2MZ', 'X',\n'2NO', 'X',\n'2NP', 'X',\n'2OG', \
'X',\n'2PA', 'X',\n'2PC', 'X',\n'2PE', 'X',\n'2PG'\
, 'X',\n'2PH', 'X',\n'2PI', 'X',\n'2PL', 'X',\n'2P\
P', 'X',\n'2PU', 'X',\n'2SI', 'X',\n'2TB', 'X',\n'\
34C', 'X',\n'35G', 'X',\n'3AA', 'X',\n'3AD', 'X',\\
n'3AH', 'H',\n'3AN', 'X',\n'3AP', 'X',\n'3AT', 'X'\
,\n'3BT', 'X',\n'3CH', 'X',\n'3CN', 'X',\n'3CO', '\
X',\n'3CP', 'X',\n'3DR', 'X',\n'3EP', 'X',\n'3FM',\
 'X',\n'3GA', 'X',\n'3GP', 'X',\n'3HB', 'X',\n'3HC\
', 'X',\n'3HP', 'X',\n'3IB', 'X',\n'3ID', 'X',\n'3\
IN', 'X',\n'3MA', 'X',\n'3MB', 'X',\n'3MC', 'X',\n\
'3MD', 'D',\n'3MF', 'X',\n'3MP', 'X',\n'3MT', 'X',\
\n'3OL', 'X',\n'3PA', 'X',\n'3PG', 'X',\n'3PO', 'X\
',\n'3PP', 'X',\n'3PY', 'X',\n'49A', 'X',\n'4AB', \
'X',\n'4AM', 'X',\n'4AN', 'X',\n'4AP', 'X',\n'4BA'\
, 'X',\n'4BT', 'X',\n'4CA', 'X',\n'4CO', 'X',\n'4H\
P', 'X',\n'4IP', 'X',\n'4MO', 'X',\n'4MV', 'X',\n'\
4MZ', 'X',\n'4NC', 'X',\n'4NP', 'X',\n'4OX', 'X',\\
n'4PB', 'X',\n'4PN', 'X',\n'4PP', 'X',\n'4SC', 'X'\
,\n'4SU', 'X',\n'4TB', 'X',\n'55C', 'X',\n'5AD', '\
X',\n'5AN', 'X',\n'5AT', 'X',\n'5CM', 'X',\n'5GP',\
 'X',\n'5HP', 'E',\n'5HT', 'X',\n'5IT', 'X',\n'5IU\
', 'X',\n'5MB', 'X',\n'5MC', 'X',\n'5MD', 'X',\n'5\
MP', 'X',\n'5MU', 'X',\n'5NC', 'X',\n'5OB', 'X',\n\
'5PA', 'X',\n'5PV', 'X',\n'6AB', 'X',\n'6CT', 'X',\
\n'6HA', 'X',\n'6HC', 'X',\n'6HG', 'X',\n'6HT', 'X\
',\n'6IN', 'X',\n'6MO', 'X',\n'6MP', 'X',\n'6PG', \
'X',\n'6WO', 'X',\n'70U', 'X',\n'7DG', 'X',\n'7HP'\
, 'X',\n'7I2', 'X',\n'7MG', 'X',\n'7MQ', 'X',\n'7N\
I', 'X',\n'87Y', 'X',\n'8AD', 'X',\n'8BR', 'X',\n'\
8IG', 'X',\n'8IN', 'X',\n'8OG', 'X',\n'95A', 'X',\\
n'9AD', 'X',\n'9AM', 'X',\n'9AP', 'X',\n'9DG', 'X'\
,\n'9DI', 'X',\n'9HX', 'X',\n'9OH', 'X',\n'9TA', '\
X',\n'A12', 'X',\n'A15', 'X',\n'A23', 'X',\n'A24',\
 'X',\n'A26', 'X',\n'A2G', 'X',\n'A2P', 'X',\n'A32\
', 'X',\n'A3P', 'X',\n'A4P', 'X',\n'A5P', 'X',\n'A\
70', 'X',\n'A76', 'X',\n'A77', 'X',\n'A78', 'X',\n\
'A79', 'X',\n'A80', 'X',\n'A85', 'X',\n'A88', 'X',\
\n'A9A', 'X',\n'AA3', 'X',\n'AA4', 'X',\n'AA6', 'X\
',\n'AAA', 'X',\n'AAB', 'X',\n'AAC', 'X',\n'AAE', \
'X',\n'AAG', 'R',\n'AAH', 'X',\n'AAM', 'X',\n'AAN'\
, 'X',\n'AAP', 'X',\n'AAR', 'R',\n'AAS', 'X',\n'AA\
T', 'X',\n'ABA', 'X',\n'ABC', 'X',\n'ABD', 'X',\n'\
ABE', 'X',\n'ABH', 'X',\n'ABI', 'X',\n'ABK', 'X',\\
n'ABM', 'X',\n'ABN', 'X',\n'ABP', 'X',\n'ABR', 'X'\
,\n'ABS', 'X',\n'ABU', 'X',\n'AC1', 'X',\n'AC2', '\
X',\n'ACA', 'X',\n'ACB', 'D',\n'ACC', 'C',\n'ACD',\
 'X',\n'ACE', 'X',\n'ACH', 'X',\n'ACI', 'X',\n'ACL\
', 'R',\n'ACM', 'X',\n'ACN', 'X',\n'ACO', 'X',\n'A\
CP', 'X',\n'ACQ', 'X',\n'ACR', 'X',\n'ACS', 'X',\n\
'ACT', 'X',\n'ACV', 'V',\n'ACX', 'X',\n'ACY', 'X',\
\n'AD2', 'X',\n'AD3', 'X',\n'ADC', 'X',\n'ADD', 'X\
',\n'ADE', 'X',\n'ADH', 'X',\n'ADI', 'X',\n'ADM', \
'X',\n'ADN', 'X',\n'ADP', 'X',\n'ADQ', 'X',\n'ADR'\
, 'X',\n'ADS', 'X',\n'ADT', 'X',\n'ADU', 'X',\n'AD\
W', 'X',\n'ADX', 'X',\n'AE2', 'X',\n'AEA', 'X',\n'\
AEB', 'X',\n'AEI', 'D',\n'AEN', 'X',\n'AET', 'T',\\
n'AF1', 'X',\n'AF3', 'X',\n'AFA', 'D',\n'AFP', 'X'\
,\n'AG7', 'X',\n'AGB', 'X',\n'AGF', 'X',\n'AGL', '\
X',\n'AGM', 'R',\n'AGN', 'X',\n'AGP', 'X',\n'AGS',\
 'X',\n'AGU', 'X',\n'AH0', 'X',\n'AH1', 'X',\n'AHA\
', 'X',\n'AHB', 'D',\n'AHC', 'X',\n'AHF', 'X',\n'A\
HG', 'X',\n'AHH', 'X',\n'AHM', 'X',\n'AHO', 'X',\n\
'AHP', 'X',\n'AHS', 'X',\n'AHT', 'Y',\n'AHU', 'X',\
\n'AHX', 'X',\n'AI1', 'X',\n'AI2', 'X',\n'AIB', 'X\
',\n'AIC', 'X',\n'AIM', 'X',\n'AIP', 'X',\n'AIQ', \
'X',\n'AIR', 'X',\n'AJ3', 'X',\n'AKB', 'X',\n'AKG'\
, 'X',\n'AKR', 'X',\n'AL1', 'X',\n'AL2', 'X',\n'AL\
3', 'X',\n'AL4', 'X',\n'AL5', 'X',\n'AL6', 'X',\n'\
AL7', 'X',\n'AL8', 'X',\n'AL9', 'X',\n'ALA', 'A',\\
n'ALB', 'X',\n'ALC', 'X',\n'ALD', 'L',\n'ALE', 'X'\
,\n'ALF', 'X',\n'ALG', 'X',\n'ALL', 'X',\n'ALM', '\
A',\n'ALN', 'A',\n'ALO', 'T',\n'ALP', 'X',\n'ALQ',\
 'X',\n'ALR', 'X',\n'ALS', 'X',\n'ALT', 'A',\n'ALY\
', 'K',\n'ALZ', 'X',\n'AMA', 'X',\n'AMB', 'X',\n'A\
MC', 'X',\n'AMD', 'X',\n'AMG', 'X',\n'AMH', 'X',\n\
'AMI', 'X',\n'AML', 'X',\n'AMN', 'X',\n'AMO', 'X',\
\n'AMP', 'X',\n'AMQ', 'X',\n'AMR', 'X',\n'AMS', 'X\
',\n'AMT', 'X',\n'AMU', 'X',\n'AMW', 'X',\n'AMX', \
'X',\n'AMY', 'X',\n'ANA', 'X',\n'ANB', 'X',\n'ANC'\
, 'X',\n'AND', 'X',\n'ANE', 'X',\n'ANI', 'X',\n'AN\
L', 'X',\n'ANO', 'X',\n'ANP', 'X',\n'ANS', 'X',\n'\
ANT', 'X',\n'AOE', 'X',\n'AOP', 'X',\n'AP1', 'X',\\
n'AP2', 'X',\n'AP3', 'X',\n'AP4', 'X',\n'AP5', 'X'\
,\n'AP6', 'X',\n'APA', 'X',\n'APB', 'X',\n'APC', '\
X',\n'APE', 'F',\n'APF', 'X',\n'APG', 'X',\n'APH',\
 'A',\n'API', 'X',\n'APL', 'X',\n'APM', 'X',\n'APN\
', 'G',\n'APP', 'X',\n'APQ', 'X',\n'APR', 'X',\n'A\
PS', 'X',\n'APT', 'X',\n'APU', 'X',\n'APX', 'X',\n\
'APY', 'X',\n'APZ', 'X',\n'AQS', 'X',\n'AR1', 'X',\
\n'AR2', 'X',\n'ARA', 'X',\n'ARB', 'X',\n'ARC', 'X\
',\n'ARD', 'X',\n'ARG', 'R',\n'ARH', 'X',\n'ARI', \
'X',\n'ARM', 'R',\n'ARN', 'X',\n'ARO', 'R',\n'ARP'\
, 'X',\n'ARQ', 'X',\n'ARS', 'X',\n'AS1', 'R',\n'AS\
2', 'X',\n'ASA', 'D',\n'ASB', 'D',\n'ASC', 'X',\n'\
ASD', 'X',\n'ASE', 'X',\n'ASF', 'X',\n'ASI', 'X',\\
n'ASK', 'D',\n'ASL', 'X',\n'ASM', 'N',\n'ASO', 'X'\
,\n'ASP', 'D',\n'ASQ', 'X',\n'ASU', 'X',\n'ATA', '\
X',\n'ATC', 'X',\n'ATD', 'X',\n'ATF', 'X',\n'ATG',\
 'X',\n'ATH', 'X',\n'ATM', 'X',\n'ATO', 'X',\n'ATP\
', 'X',\n'ATQ', 'X',\n'ATR', 'X',\n'ATT', 'X',\n'A\
TY', 'X',\n'ATZ', 'X',\n'AUC', 'X',\n'AUR', 'X',\n\
'AVG', 'X',\n'AXP', 'X',\n'AYA', 'A',\n'AZ2', 'X',\
\n'AZA', 'X',\n'AZC', 'X',\n'AZD', 'X',\n'AZE', 'X\
',\n'AZI', 'X',\n'AZL', 'X',\n'AZM', 'X',\n'AZR', \
'X',\n'AZT', 'X',\n'B12', 'X',\n'B1F', 'F',\n'B2A'\
, 'A',\n'B2F', 'F',\n'B2I', 'I',\n'B2V', 'V',\n'B3\
I', 'X',\n'B3P', 'X',\n'B7G', 'X',\n'B96', 'X',\n'\
B9A', 'X',\n'BA1', 'X',\n'BAA', 'X',\n'BAB', 'X',\\
n'BAC', 'X',\n'BAF', 'X',\n'BAH', 'X',\n'BAI', 'X'\
,\n'BAK', 'X',\n'BAL', 'A',\n'BAM', 'X',\n'BAO', '\
X',\n'BAP', 'X',\n'BAR', 'X',\n'BAS', 'X',\n'BAT',\
 'F',\n'BAY', 'X',\n'BAZ', 'X',\n'BB1', 'X',\n'BB2\
', 'X',\n'BBA', 'X',\n'BBH', 'X',\n'BBS', 'X',\n'B\
BT', 'X',\n'BBZ', 'X',\n'BCA', 'X',\n'BCB', 'X',\n\
'BCC', 'X',\n'BCD', 'X',\n'BCL', 'X',\n'BCN', 'X',\
\n'BCR', 'X',\n'BCS', 'C',\n'BCT', 'X',\n'BCY', 'X\
',\n'BCZ', 'X',\n'BDA', 'X',\n'BDG', 'X',\n'BDK', \
'X',\n'BDM', 'X',\n'BDN', 'X',\n'BDS', 'X',\n'BE1'\
, 'X',\n'BE2', 'X',\n'BEA', 'X',\n'BEF', 'X',\n'BE\
N', 'X',\n'BEO', 'X',\n'BEP', 'X',\n'BER', 'X',\n'\
BES', 'X',\n'BET', 'X',\n'BEZ', 'X',\n'BF2', 'X',\\
n'BFA', 'X',\n'BFD', 'X',\n'BFP', 'X',\n'BFS', 'X'\
,\n'BFU', 'X',\n'BG6', 'X',\n'BGF', 'X',\n'BGG', '\
X',\n'BGL', 'X',\n'BGN', 'X',\n'BGP', 'X',\n'BGX',\
 'X',\n'BH4', 'X',\n'BHA', 'X',\n'BHC', 'X',\n'BHD\
', 'D',\n'BHO', 'X',\n'BHS', 'X',\n'BIC', 'X',\n'B\
IN', 'X',\n'BIO', 'X',\n'BIP', 'X',\n'BIS', 'X',\n\
'BIZ', 'X',\n'BJH', 'X',\n'BJI', 'X',\n'BJP', 'X',\
\n'BLA', 'X',\n'BLB', 'X',\n'BLE', 'L',\n'BLG', 'P\
',\n'BLI', 'X',\n'BLM', 'X',\n'BLV', 'X',\n'BLY', \
'K',\n'BM1', 'X',\n'BM2', 'X',\n'BM5', 'X',\n'BM9'\
, 'X',\n'BMA', 'X',\n'BMD', 'X',\n'BME', 'X',\n'BM\
P', 'X',\n'BMQ', 'X',\n'BMS', 'X',\n'BMT', 'T',\n'\
BMU', 'X',\n'BMY', 'X',\n'BMZ', 'X',\n'BNA', 'X',\\
n'BNG', 'X',\n'BNI', 'X',\n'BNN', 'F',\n'BNO', 'L'\
,\n'BNS', 'X',\n'BNZ', 'X',\n'BO3', 'X',\n'BO4', '\
X',\n'BOC', 'X',\n'BOG', 'X',\n'BOM', 'X',\n'BOT',\
 'X',\n'BOX', 'X',\n'BOZ', 'X',\n'BPA', 'X',\n'BPB\
', 'X',\n'BPD', 'X',\n'BPG', 'X',\n'BPH', 'X',\n'B\
PI', 'X',\n'BPJ', 'X',\n'BPM', 'X',\n'BPN', 'X',\n\
'BPO', 'X',\n'BPP', 'X',\n'BPT', 'X',\n'BPY', 'X',\
\n'BRB', 'X',\n'BRC', 'X',\n'BRE', 'X',\n'BRI', 'X\
',\n'BRL', 'X',\n'BRM', 'X',\n'BRN', 'X',\n'BRO', \
'X',\n'BRS', 'X',\n'BRU', 'X',\n'BRZ', 'X',\n'BSB'\
, 'X',\n'BSI', 'X',\n'BSP', 'X',\n'BT1', 'X',\n'BT\
2', 'X',\n'BT3', 'X',\n'BTA', 'L',\n'BTB', 'X',\n'\
BTC', 'C',\n'BTD', 'X',\n'BTN', 'X',\n'BTP', 'X',\\
n'BTR', 'W',\n'BU1', 'X',\n'BUA', 'X',\n'BUB', 'X'\
,\n'BUC', 'X',\n'BUG', 'X',\n'BUL', 'X',\n'BUM', '\
X',\n'BUQ', 'X',\n'BUT', 'X',\n'BVD', 'X',\n'BX3',\
 'X',\n'BYS', 'X',\n'BZ1', 'X',\n'BZA', 'X',\n'BZB\
', 'X',\n'BZC', 'X',\n'BZD', 'X',\n'BZF', 'X',\n'B\
ZI', 'X',\n'BZM', 'X',\n'BZO', 'X',\n'BZP', 'X',\n\
'BZQ', 'X',\n'BZS', 'X',\n'BZT', 'X',\n'C02', 'X',\
\n'C11', 'X',\n'C1O', 'X',\n'C20', 'X',\n'C24', 'X\
',\n'C2F', 'X',\n'C2O', 'X',\n'C2P', 'X',\n'C3M', \
'X',\n'C3P', 'X',\n'C3X', 'X',\n'C48', 'X',\n'C4M'\
, 'X',\n'C4X', 'X',\n'C5C', 'X',\n'C5M', 'X',\n'C5\
P', 'X',\n'C5X', 'X',\n'C60', 'X',\n'C6C', 'X',\n'\
C6M', 'X',\n'C78', 'X',\n'C8E', 'X',\n'CA3', 'X',\\
n'CA5', 'X',\n'CAA', 'X',\n'CAB', 'X',\n'CAC', 'X'\
,\n'CAD', 'X',\n'CAF', 'C',\n'CAG', 'X',\n'CAH', '\
X',\n'CAL', 'X',\n'CAM', 'X',\n'CAN', 'X',\n'CAO',\
 'X',\n'CAP', 'X',\n'CAQ', 'X',\n'CAR', 'X',\n'CAS\
', 'C',\n'CAT', 'X',\n'CAV', 'X',\n'CAY', 'C',\n'C\
AZ', 'X',\n'CB3', 'X',\n'CB4', 'X',\n'CBA', 'X',\n\
'CBD', 'X',\n'CBG', 'X',\n'CBI', 'X',\n'CBL', 'X',\
\n'CBM', 'X',\n'CBN', 'X',\n'CBO', 'X',\n'CBP', 'X\
',\n'CBS', 'X',\n'CBX', 'X',\n'CBZ', 'X',\n'CC0', \
'X',\n'CC1', 'X',\n'CCC', 'X',\n'CCH', 'X',\n'CCI'\
, 'X',\n'CCM', 'X',\n'CCN', 'X',\n'CCO', 'X',\n'CC\
P', 'X',\n'CCR', 'X',\n'CCS', 'C',\n'CCV', 'X',\n'\
CCY', 'X',\n'CD1', 'X',\n'CDC', 'X',\n'CDE', 'X',\\
n'CDF', 'X',\n'CDI', 'X',\n'CDL', 'X',\n'CDM', 'X'\
,\n'CDP', 'X',\n'CDR', 'X',\n'CDU', 'X',\n'CE1', '\
X',\n'CEA', 'C',\n'CEB', 'X',\n'CEC', 'X',\n'CED',\
 'X',\n'CEF', 'X',\n'CEH', 'X',\n'CEM', 'X',\n'CEO\
', 'X',\n'CEP', 'X',\n'CEQ', 'X',\n'CER', 'X',\n'C\
ES', 'G',\n'CET', 'X',\n'CFC', 'X',\n'CFF', 'X',\n\
'CFM', 'X',\n'CFO', 'X',\n'CFP', 'X',\n'CFS', 'X',\
\n'CFX', 'X',\n'CGN', 'X',\n'CGP', 'X',\n'CGS', 'X\
',\n'CGU', 'E',\n'CH2', 'X',\n'CH3', 'X',\n'CHA', \
'X',\n'CHB', 'X',\n'CHD', 'X',\n'CHF', 'X',\n'CHG'\
, 'G',\n'CHI', 'X',\n'CHN', 'X',\n'CHO', 'X',\n'CH\
P', 'G',\n'CHR', 'X',\n'CHS', 'F',\n'CHT', 'X',\n'\
CHX', 'X',\n'CIC', 'X',\n'CIN', 'X',\n'CIP', 'X',\\
n'CIR', 'X',\n'CIT', 'X',\n'CIU', 'X',\n'CKI', 'X'\
,\n'CL1', 'X',\n'CL2', 'X',\n'CLA', 'X',\n'CLB', '\
A',\n'CLC', 'S',\n'CLD', 'A',\n'CLE', 'L',\n'CLF',\
 'X',\n'CLK', 'S',\n'CLL', 'X',\n'CLM', 'X',\n'CLN\
', 'X',\n'CLO', 'X',\n'CLP', 'X',\n'CLQ', 'X',\n'C\
LR', 'X',\n'CLS', 'X',\n'CLT', 'X',\n'CLX', 'X',\n\
'CLY', 'X',\n'CMA', 'R',\n'CMC', 'X',\n'CMD', 'X',\
\n'CME', 'C',\n'CMG', 'X',\n'CMK', 'X',\n'CMN', 'X\
',\n'CMO', 'X',\n'CMP', 'X',\n'CMR', 'X',\n'CMS', \
'X',\n'CMT', 'C',\n'CMX', 'X',\n'CNA', 'X',\n'CNC'\
, 'X',\n'CND', 'X',\n'CNH', 'X',\n'CNM', 'X',\n'CN\
N', 'X',\n'CNO', 'X',\n'CNP', 'X',\n'CO2', 'X',\n'\
CO3', 'X',\n'CO5', 'X',\n'CO8', 'X',\n'COA', 'X',\\
n'COB', 'X',\n'COC', 'X',\n'COD', 'X',\n'COE', 'X'\
,\n'COF', 'X',\n'COH', 'X',\n'COI', 'X',\n'COJ', '\
X',\n'COL', 'X',\n'COM', 'X',\n'CON', 'X',\n'COP',\
 'X',\n'COR', 'X',\n'COS', 'X',\n'COT', 'X',\n'COY\
', 'X',\n'CP1', 'G',\n'CP2', 'X',\n'CP4', 'X',\n'C\
PA', 'X',\n'CPB', 'X',\n'CPC', 'X',\n'CPD', 'X',\n\
'CPG', 'X',\n'CPH', 'X',\n'CPI', 'X',\n'CPM', 'X',\
\n'CPN', 'G',\n'CPO', 'X',\n'CPP', 'X',\n'CPQ', 'X\
',\n'CPR', 'X',\n'CPS', 'X',\n'CPT', 'X',\n'CPU', \
'X',\n'CPV', 'X',\n'CPY', 'X',\n'CR1', 'X',\n'CR6'\
, 'X',\n'CRA', 'X',\n'CRB', 'X',\n'CRC', 'X',\n'CR\
G', 'X',\n'CRH', 'X',\n'CRO', 'T',\n'CRP', 'X',\n'\
CRQ', 'X',\n'CRS', 'X',\n'CRT', 'X',\n'CRY', 'X',\\
n'CSA', 'C',\n'CSB', 'X',\n'CSD', 'C',\n'CSE', 'C'\
,\n'CSH', 'X',\n'CSI', 'X',\n'CSN', 'X',\n'CSO', '\
C',\n'CSP', 'C',\n'CSR', 'C',\n'CSS', 'C',\n'CST',\
 'X',\n'CSW', 'C',\n'CSX', 'C',\n'CSY', 'X',\n'CSZ\
', 'C',\n'CT3', 'X',\n'CTA', 'X',\n'CTB', 'X',\n'C\
TC', 'X',\n'CTD', 'X',\n'CTH', 'T',\n'CTO', 'X',\n\
'CTP', 'X',\n'CTR', 'X',\n'CTS', 'X',\n'CTT', 'X',\
\n'CTY', 'X',\n'CTZ', 'X',\n'CU1', 'X',\n'CUA', 'X\
',\n'CUC', 'X',\n'CUL', 'X',\n'CUO', 'X',\n'CUZ', \
'X',\n'CVI', 'X',\n'CXF', 'X',\n'CXL', 'X',\n'CXM'\
, 'M',\n'CXN', 'X',\n'CXP', 'X',\n'CXS', 'X',\n'CY\
1', 'C',\n'CY3', 'X',\n'CYB', 'X',\n'CYC', 'X',\n'\
CYF', 'C',\n'CYG', 'C',\n'CYH', 'X',\n'CYL', 'X',\\
n'CYM', 'C',\n'CYN', 'X',\n'CYO', 'X',\n'CYP', 'X'\
,\n'CYQ', 'C',\n'CYS', 'C',\n'CYU', 'X',\n'CYY', '\
X',\n'CYZ', 'X',\n'CZH', 'X',\n'CZZ', 'C',\n'D12',\
 'X',\n'D13', 'X',\n'D16', 'X',\n'D18', 'X',\n'D19\
', 'X',\n'D1P', 'X',\n'D24', 'X',\n'D34', 'X',\n'D\
35', 'X',\n'D4D', 'X',\n'D4T', 'X',\n'D6G', 'X',\n\
'DA2', 'R',\n'DA3', 'X',\n'DA6', 'X',\n'DA7', 'X',\
\n'DAA', 'X',\n'DAB', 'X',\n'DAC', 'X',\n'DAD', 'X\
',\n'DAE', 'X',\n'DAF', 'X',\n'DAG', 'X',\n'DAH', \
'A',\n'DAJ', 'X',\n'DAK', 'X',\n'DAL', 'A',\n'DAM'\
, 'A',\n'DAN', 'X',\n'DAO', 'X',\n'DAP', 'X',\n'DA\
Q', 'X',\n'DAR', 'R',\n'DAS', 'D',\n'DAT', 'X',\n'\
DAU', 'X',\n'DAV', 'X',\n'DBA', 'X',\n'DBD', 'X',\\
n'DBF', 'X',\n'DBG', 'X',\n'DBI', 'X',\n'DBV', 'X'\
,\n'DBY', 'Y',\n'DCA', 'X',\n'DCB', 'X',\n'DCE', '\
X',\n'DCF', 'X',\n'DCG', 'X',\n'DCH', 'X',\n'DCI',\
 'I',\n'DCL', 'X',\n'DCM', 'X',\n'DCP', 'X',\n'DCS\
', 'X',\n'DCT', 'X',\n'DCY', 'C',\n'DCZ', 'X',\n'D\
DA', 'X',\n'DDB', 'X',\n'DDC', 'X',\n'DDF', 'X',\n\
'DDG', 'X',\n'DDH', 'X',\n'DDL', 'X',\n'DDM', 'X',\
\n'DDO', 'L',\n'DDP', 'X',\n'DDQ', 'X',\n'DDT', 'Y\
',\n'DDU', 'X',\n'DEA', 'X',\n'DEB', 'X',\n'DEC', \
'X',\n'DEF', 'X',\n'DEL', 'X',\n'DEM', 'X',\n'DEN'\
, 'X',\n'DEP', 'X',\n'DEQ', 'X',\n'DES', 'X',\n'DE\
T', 'X',\n'DFC', 'X',\n'DFG', 'X',\n'DFI', 'X',\n'\
DFL', 'X',\n'DFO', 'X',\n'DFP', 'X',\n'DFR', 'X',\\
n'DFT', 'X',\n'DFV', 'X',\n'DFX', 'X',\n'DG2', 'X'\
,\n'DG3', 'X',\n'DG6', 'X',\n'DGA', 'X',\n'DGD', '\
X',\n'DGG', 'X',\n'DGL', 'E',\n'DGN', 'Q',\n'DGP',\
 'X',\n'DGT', 'X',\n'DGX', 'X',\n'DH2', 'X',\n'DHA\
', 'A',\n'DHB', 'X',\n'DHC', 'X',\n'DHD', 'X',\n'D\
HE', 'X',\n'DHF', 'X',\n'DHG', 'X',\n'DHI', 'H',\n\
'DHL', 'X',\n'DHM', 'X',\n'DHN', 'V',\n'DHP', 'X',\
\n'DHQ', 'X',\n'DHR', 'X',\n'DHS', 'X',\n'DHT', 'X\
',\n'DHU', 'X',\n'DHY', 'X',\n'DHZ', 'X',\n'DI2', \
'X',\n'DI3', 'G',\n'DI4', 'X',\n'DI5', 'X',\n'DIA'\
, 'X',\n'DIC', 'X',\n'DIF', 'X',\n'DIG', 'X',\n'DI\
I', 'X',\n'DIL', 'I',\n'DIM', 'X',\n'DIO', 'X',\n'\
DIP', 'X',\n'DIQ', 'X',\n'DIS', 'X',\n'DIT', 'X',\\
n'DIV', 'V',\n'DIX', 'X',\n'DIY', 'X',\n'DKA', 'X'\
,\n'DLA', 'X',\n'DLE', 'L',\n'DLF', 'X',\n'DLS', '\
K',\n'DLY', 'K',\n'DM1', 'X',\n'DM2', 'X',\n'DM3',\
 'X',\n'DM4', 'X',\n'DM5', 'X',\n'DM6', 'X',\n'DM7\
', 'X',\n'DM8', 'X',\n'DM9', 'X',\n'DMA', 'X',\n'D\
MB', 'X',\n'DMC', 'X',\n'DMD', 'X',\n'DME', 'X',\n\
'DMF', 'X',\n'DMG', 'G',\n'DMH', 'N',\n'DMI', 'X',\
\n'DMJ', 'X',\n'DML', 'X',\n'DMM', 'X',\n'DMN', 'X\
',\n'DMO', 'X',\n'DMP', 'X',\n'DMQ', 'X',\n'DMR', \
'X',\n'DMS', 'X',\n'DMT', 'X',\n'DMV', 'X',\n'DMY'\
, 'X',\n'DNC', 'X',\n'DND', 'X',\n'DNH', 'X',\n'DN\
J', 'X',\n'DNN', 'X',\n'DNP', 'X',\n'DNQ', 'X',\n'\
DNR', 'X',\n'DO2', 'X',\n'DO3', 'X',\n'DOA', 'X',\\
n'DOB', 'X',\n'DOC', 'X',\n'DOH', 'D',\n'DOM', 'X'\
,\n'DOS', 'X',\n'DOX', 'X',\n'DP5', 'X',\n'DP7', '\
X',\n'DPA', 'X',\n'DPC', 'X',\n'DPD', 'X',\n'DPE',\
 'X',\n'DPG', 'X',\n'DPH', 'F',\n'DPM', 'X',\n'DPN\
', 'F',\n'DPO', 'X',\n'DPP', 'X',\n'DPR', 'P',\n'D\
PS', 'X',\n'DPT', 'X',\n'DPX', 'X',\n'DPY', 'X',\n\
'DPZ', 'X',\n'DQH', 'X',\n'DQN', 'X',\n'DR1', 'X',\
\n'DRB', 'X',\n'DRC', 'X',\n'DRI', 'X',\n'DRP', 'X\
',\n'DRT', 'X',\n'DRU', 'X',\n'DSA', 'X',\n'DSB', \
'X',\n'DSC', 'X',\n'DSD', 'X',\n'DSE', 'S',\n'DSI'\
, 'X',\n'DSN', 'S',\n'DSP', 'D',\n'DSR', 'X',\n'DS\
S', 'X',\n'DSX', 'X',\n'DSY', 'X',\n'DTB', 'X',\n'\
DTD', 'X',\n'DTH', 'T',\n'DTN', 'X',\n'DTO', 'X',\\
n'DTP', 'X',\n'DTQ', 'X',\n'DTR', 'W',\n'DTT', 'X'\
,\n'DTY', 'Y',\n'DUD', 'X',\n'DUO', 'X',\n'DUR', '\
X',\n'DUT', 'X',\n'DVA', 'V',\n'DVR', 'X',\n'DX9',\
 'X',\n'DXA', 'X',\n'DXB', 'X',\n'DXC', 'X',\n'DXG\
', 'X',\n'DXX', 'X',\n'DZF', 'X',\n'E09', 'X',\n'E\
20', 'X',\n'E2P', 'X',\n'E3G', 'X',\n'E4N', 'X',\n\
'E4P', 'X',\n'E64', 'X',\n'E6C', 'X',\n'E96', 'X',\
\n'E97', 'X',\n'EA2', 'X',\n'EAA', 'X',\n'EAP', 'X\
',\n'EBP', 'X',\n'EBW', 'X',\n'ECO', 'X',\n'EDA', \
'X',\n'EDC', 'X',\n'EDE', 'X',\n'EDO', 'X',\n'EDR'\
, 'X',\n'EEB', 'X',\n'EEE', 'X',\n'EFC', 'X',\n'EF\
Z', 'X',\n'EG1', 'X',\n'EG2', 'X',\n'EG3', 'X',\n'\
EGC', 'X',\n'EGL', 'X',\n'EHP', 'A',\n'EIC', 'X',\\
n'EJT', 'X',\n'ELA', 'X',\n'EMB', 'X',\n'EMC', 'X'\
,\n'EMD', 'X',\n'EMM', 'X',\n'EMO', 'X',\n'EMP', '\
X',\n'EMR', 'X',\n'ENA', 'X',\n'ENC', 'X',\n'ENH',\
 'X',\n'ENO', 'X',\n'ENP', 'X',\n'EOA', 'X',\n'EOH\
', 'X',\n'EOT', 'X',\n'EOX', 'X',\n'EPA', 'X',\n'E\
PE', 'X',\n'EPH', 'X',\n'EPI', 'X',\n'EPN', 'X',\n\
'EPO', 'X',\n'EPT', 'X',\n'EPU', 'X',\n'EPX', 'X',\
\n'EPY', 'X',\n'EQI', 'X',\n'EQP', 'X',\n'EQU', 'X\
',\n'ERG', 'X',\n'ERI', 'X',\n'ERY', 'X',\n'ESC', \
'X',\n'ESD', 'X',\n'ESI', 'X',\n'ESO', 'X',\n'ESP'\
, 'X',\n'EST', 'X',\n'ESX', 'X',\n'ETA', 'X',\n'ET\
C', 'X',\n'ETD', 'X',\n'ETF', 'X',\n'ETH', 'X',\n'\
ETI', 'X',\n'ETN', 'X',\n'ETO', 'X',\n'ETP', 'X',\\
n'ETR', 'X',\n'ETS', 'X',\n'ETY', 'X',\n'EU3', 'X'\
,\n'EUG', 'X',\n'EYS', 'C',\n'F09', 'X',\n'F2B', '\
X',\n'F3S', 'X',\n'F42', 'X',\n'F43', 'X',\n'F4S',\
 'X',\n'F6B', 'X',\n'F6P', 'X',\n'F89', 'X',\n'FA1\
', 'X',\n'FA5', 'F',\n'FAA', 'X',\n'FAB', 'X',\n'F\
AC', 'X',\n'FAD', 'X',\n'FAF', 'X',\n'FAG', 'X',\n\
'FAM', 'X',\n'FAR', 'X',\n'FAS', 'X',\n'FAT', 'X',\
\n'FBA', 'X',\n'FBE', 'X',\n'FBI', 'X',\n'FBP', 'X\
',\n'FBQ', 'X',\n'FBS', 'X',\n'FBT', 'X',\n'FBU', \
'X',\n'FCA', 'X',\n'FCB', 'X',\n'FCI', 'X',\n'FCN'\
, 'X',\n'FCO', 'X',\n'FCR', 'X',\n'FCT', 'X',\n'FC\
X', 'X',\n'FCY', 'C',\n'FD1', 'F',\n'FD2', 'F',\n'\
FD3', 'F',\n'FD4', 'F',\n'FDA', 'X',\n'FDC', 'X',\\
n'FDI', 'X',\n'FDP', 'X',\n'FDS', 'X',\n'FE2', 'X'\
,\n'FEA', 'X',\n'FEL', 'X',\n'FEM', 'X',\n'FEN', '\
X',\n'FEO', 'X',\n'FEP', 'X',\n'FER', 'X',\n'FES',\
 'X',\n'FFB', 'X',\n'FFC', 'X',\n'FFF', 'X',\n'FFO\
', 'X',\n'FGL', 'G',\n'FHB', 'X',\n'FHC', 'X',\n'F\
HP', 'X',\n'FHU', 'X',\n'FID', 'X',\n'FII', 'X',\n\
'FIP', 'X',\n'FK5', 'X',\n'FKA', 'X',\n'FKI', 'X',\
\n'FKP', 'X',\n'FL2', 'X',\n'FL9', 'X',\n'FLA', 'A\
',\n'FLC', 'X',\n'FLD', 'X',\n'FLE', 'L',\n'FLF', \
'X',\n'FLO', 'X',\n'FLP', 'X',\n'FLT', 'Y',\n'FLU'\
, 'X',\n'FLX', 'X',\n'FM1', 'X',\n'FM2', 'X',\n'FM\
A', 'X',\n'FMB', 'X',\n'FMC', 'X',\n'FME', 'M',\n'\
FMN', 'X',\n'FMP', 'X',\n'FMR', 'X',\n'FMS', 'X',\\
n'FMT', 'X',\n'FNE', 'X',\n'FNP', 'X',\n'FNS', 'X'\
,\n'FOC', 'X',\n'FOE', 'X',\n'FOG', 'F',\n'FOH', '\
X',\n'FOK', 'X',\n'FOL', 'X',\n'FON', 'X',\n'FOP',\
 'X',\n'FOR', 'X',\n'FOS', 'X',\n'FPA', 'X',\n'FPC\
', 'X',\n'FPI', 'X',\n'FPO', 'X',\n'FPP', 'X',\n'F\
PT', 'X',\n'FQP', 'X',\n'FRA', 'X',\n'FRD', 'F',\n\
'FRU', 'X',\n'FS3', 'X',\n'FS4', 'X',\n'FSB', 'X',\
\n'FSO', 'X',\n'FSX', 'X',\n'FTC', 'X',\n'FTP', 'X\
',\n'FTR', 'W',\n'FTT', 'X',\n'FTY', 'Y',\n'FUA', \
'X',\n'FUC', 'X',\n'FUM', 'X',\n'FUP', 'X',\n'FVF'\
, 'X',\n'FXP', 'X',\n'FXV', 'X',\n'FYA', 'F',\n'G1\
6', 'X',\n'G1P', 'X',\n'G20', 'X',\n'G21', 'X',\n'\
G23', 'X',\n'G26', 'X',\n'G28', 'X',\n'G2F', 'X',\\
n'G37', 'X',\n'G39', 'X',\n'G3H', 'X',\n'G3P', 'X'\
,\n'G4D', 'X',\n'G6D', 'X',\n'G6P', 'X',\n'G6Q', '\
X',\n'G7M', 'X',\n'GA2', 'X',\n'GAA', 'X',\n'GAB',\
 'X',\n'GAC', 'X',\n'GAI', 'X',\n'GAL', 'X',\n'GAM\
', 'X',\n'GAN', 'X',\n'GAO', 'X',\n'GAP', 'X',\n'G\
AR', 'G',\n'GAS', 'X',\n'GAT', 'X',\n'GBC', 'X',\n\
'GBI', 'X',\n'GBP', 'X',\n'GBS', 'X',\n'GBX', 'X',\
\n'GC4', 'X',\n'GCA', 'X',\n'GCD', 'X',\n'GCG', 'G\
',\n'GCH', 'G',\n'GCK', 'X',\n'GCL', 'X',\n'GCM', \
'X',\n'GCN', 'X',\n'GCO', 'X',\n'GCP', 'X',\n'GCR'\
, 'X',\n'GCS', 'X',\n'GCU', 'X',\n'GD3', 'X',\n'GD\
B', 'X',\n'GDM', 'X',\n'GDN', 'X',\n'GDP', 'X',\n'\
GDS', 'X',\n'GDU', 'X',\n'GE1', 'X',\n'GE2', 'X',\\
n'GE3', 'X',\n'GEA', 'X',\n'GEL', 'X',\n'GEM', 'X'\
,\n'GEN', 'X',\n'GEP', 'X',\n'GER', 'X',\n'GFP', '\
X',\n'GGB', 'X',\n'GGL', 'E',\n'GGP', 'X',\n'GHP',\
 'G',\n'GIP', 'X',\n'GIS', 'X',\n'GKR', 'X',\n'GL2\
', 'X',\n'GL3', 'G',\n'GL4', 'X',\n'GL5', 'X',\n'G\
L7', 'X',\n'GL9', 'X',\n'GLA', 'X',\n'GLB', 'X',\n\
'GLC', 'X',\n'GLD', 'X',\n'GLE', 'X',\n'GLF', 'X',\
\n'GLG', 'X',\n'GLH', 'Q',\n'GLI', 'X',\n'GLL', 'X\
',\n'GLM', 'G',\n'GLN', 'Q',\n'GLO', 'X',\n'GLP', \
'X',\n'GLR', 'X',\n'GLS', 'X',\n'GLT', 'X',\n'GLU'\
, 'E',\n'GLV', 'X',\n'GLW', 'X',\n'GLY', 'G',\n'GL\
Z', 'X',\n'GM1', 'X',\n'GMA', 'X',\n'GMC', 'X',\n'\
GMH', 'X',\n'GMP', 'X',\n'GMY', 'X',\n'GN7', 'X',\\
n'GNA', 'X',\n'GNB', 'X',\n'GNH', 'X',\n'GNP', 'X'\
,\n'GNT', 'X',\n'GOA', 'X',\n'GOL', 'X',\n'GOX', '\
X',\n'GP1', 'X',\n'GP3', 'X',\n'GP4', 'X',\n'GP6',\
 'X',\n'GP8', 'X',\n'GPB', 'E',\n'GPC', 'X',\n'GPE\
', 'X',\n'GPG', 'X',\n'GPI', 'X',\n'GPJ', 'X',\n'G\
PL', 'K',\n'GPM', 'X',\n'GPN', 'G',\n'GPP', 'X',\n\
'GPR', 'X',\n'GPS', 'X',\n'GPX', 'X',\n'GR1', 'X',\
\n'GR3', 'X',\n'GR4', 'X',\n'GSA', 'X',\n'GSB', 'X\
',\n'GSC', 'G',\n'GSE', 'S',\n'GSH', 'X',\n'GSP', \
'X',\n'GSR', 'X',\n'GSS', 'X',\n'GT9', 'C',\n'GTA'\
, 'X',\n'GTB', 'X',\n'GTD', 'X',\n'GTE', 'X',\n'GT\
H', 'T',\n'GTN', 'X',\n'GTO', 'X',\n'GTP', 'X',\n'\
GTR', 'X',\n'GTS', 'X',\n'GTT', 'X',\n'GTX', 'X',\\
n'GTZ', 'X',\n'GU7', 'X',\n'GUA', 'X',\n'GUD', 'X'\
,\n'GUM', 'X',\n'GUN', 'X',\n'GUP', 'X',\n'GUR', '\
X',\n'GW3', 'X',\n'GZZ', 'X',\n'H2B', 'X',\n'H2P',\
 'H',\n'H2S', 'X',\n'H2U', 'X',\n'H4B', 'X',\n'H5M\
', 'P',\n'H5P', 'X',\n'HAA', 'X',\n'HAB', 'X',\n'H\
AC', 'A',\n'HAD', 'X',\n'HAE', 'X',\n'HAG', 'X',\n\
'HAI', 'X',\n'HAM', 'X',\n'HAP', 'X',\n'HAQ', 'X',\
\n'HAR', 'R',\n'HAS', 'X',\n'HAV', 'V',\n'HAX', 'X\
',\n'HAZ', 'X',\n'HBA', 'X',\n'HBC', 'X',\n'HBD', \
'X',\n'HBI', 'X',\n'HBO', 'X',\n'HBU', 'X',\n'HBY'\
, 'X',\n'HC0', 'X',\n'HC1', 'X',\n'HC4', 'X',\n'HC\
A', 'X',\n'HCC', 'X',\n'HCI', 'X',\n'HCS', 'X',\n'\
HDA', 'X',\n'HDD', 'X',\n'HDF', 'X',\n'HDN', 'X',\\
n'HDS', 'X',\n'HDZ', 'X',\n'HE1', 'X',\n'HE6', 'X'\
,\n'HEA', 'X',\n'HEB', 'X',\n'HEC', 'X',\n'HED', '\
X',\n'HEE', 'X',\n'HEF', 'X',\n'HEG', 'X',\n'HEM',\
 'X',\n'HEN', 'X',\n'HEO', 'X',\n'HEP', 'X',\n'HEU\
', 'X',\n'HEV', 'X',\n'HEX', 'X',\n'HEZ', 'X',\n'H\
F1', 'X',\n'HFA', 'X',\n'HFP', 'X',\n'HGA', 'Q',\n\
'HGB', 'X',\n'HGC', 'X',\n'HGI', 'X',\n'HGU', 'X',\
\n'HHO', 'X',\n'HHP', 'X',\n'HIB', 'X',\n'HIC', 'H\
',\n'HII', 'X',\n'HIN', 'X',\n'HIO', 'X',\n'HIP', \
'H',\n'HIS', 'H',\n'HLE', 'X',\n'HLT', 'X',\n'HMA'\
, 'A',\n'HMB', 'X',\n'HMC', 'X',\n'HMD', 'X',\n'HM\
F', 'A',\n'HMG', 'X',\n'HMH', 'X',\n'HMI', 'L',\n'\
HMM', 'X',\n'HMN', 'X',\n'HMO', 'X',\n'HMP', 'X',\\
n'HMR', 'R',\n'HNI', 'X',\n'HNP', 'X',\n'HOA', 'X'\
,\n'HOE', 'X',\n'HOH', 'X',\n'HOM', 'X',\n'HOP', '\
X',\n'HOQ', 'X',\n'HP1', 'A',\n'HP2', 'A',\n'HP3',\
 'X',\n'HPA', 'X',\n'HPB', 'X',\n'HPC', 'X',\n'HPD\
', 'X',\n'HPE', 'A',\n'HPG', 'X',\n'HPH', 'F',\n'H\
PP', 'X',\n'HPQ', 'F',\n'HPR', 'X',\n'HPT', 'X',\n\
'HPY', 'X',\n'HQO', 'X',\n'HQQ', 'X',\n'HQU', 'X',\
\n'HRG', 'R',\n'HRI', 'X',\n'HSA', 'X',\n'HSE', 'S\
',\n'HSF', 'X',\n'HSM', 'X',\n'HSO', 'H',\n'HSP', \
'X',\n'HT1', 'X',\n'HT2', 'X',\n'HTA', 'X',\n'HTL'\
, 'X',\n'HTO', 'X',\n'HTP', 'X',\n'HTR', 'W',\n'HU\
P', 'X',\n'HUX', 'X',\n'HV5', 'A',\n'HV7', 'X',\n'\
HV8', 'X',\n'HXA', 'X',\n'HXC', 'X',\n'HXP', 'X',\\
n'HY1', 'X',\n'HYA', 'X',\n'HYB', 'X',\n'HYD', 'X'\
,\n'HYG', 'X',\n'HYP', 'P',\n'I06', 'X',\n'I10', '\
X',\n'I11', 'X',\n'I17', 'X',\n'I2P', 'X',\n'I3N',\
 'X',\n'I3P', 'X',\n'I40', 'X',\n'I48', 'X',\n'I4B\
', 'X',\n'I52', 'X',\n'I5P', 'X',\n'I84', 'G',\n'I\
AG', 'G',\n'IAS', 'X',\n'IB2', 'X',\n'IBB', 'X',\n\
'IBP', 'X',\n'IBR', 'X',\n'IBS', 'X',\n'IBZ', 'X',\
\n'IC1', 'X',\n'ICA', 'X',\n'ICI', 'X',\n'ICL', 'X\
',\n'ICP', 'X',\n'ICT', 'X',\n'ICU', 'X',\n'ID2', \
'X',\n'IDC', 'X',\n'IDG', 'X',\n'IDH', 'X',\n'IDM'\
, 'X',\n'IDO', 'X',\n'IDP', 'X',\n'IDR', 'X',\n'ID\
S', 'X',\n'IDT', 'X',\n'IDU', 'X',\n'IFG', 'X',\n'\
IFP', 'X',\n'IGL', 'X',\n'IGN', 'X',\n'IGP', 'X',\\
n'IGU', 'X',\n'IH1', 'X',\n'IH2', 'X',\n'IH3', 'X'\
,\n'IHB', 'X',\n'IHN', 'X',\n'IHP', 'X',\n'IIC', '\
X',\n'IIL', 'I',\n'IIP', 'X',\n'IK2', 'X',\n'IKT',\
 'X',\n'ILA', 'I',\n'ILE', 'I',\n'ILG', 'X',\n'ILO\
', 'X',\n'ILX', 'I',\n'IM1', 'X',\n'IM2', 'X',\n'I\
MC', 'X',\n'IMD', 'X',\n'IME', 'X',\n'IMF', 'X',\n\
'IMG', 'X',\n'IMH', 'X',\n'IMI', 'X',\n'IML', 'I',\
\n'IMM', 'X',\n'IMN', 'X',\n'IMO', 'X',\n'IMP', 'X\
',\n'IMR', 'X',\n'IMU', 'X',\n'IN0', 'D',\n'IN1', \
'R',\n'IN2', 'K',\n'IN3', 'L',\n'IN4', 'X',\n'IN5'\
, 'A',\n'IN6', 'L',\n'IN7', 'X',\n'IN8', 'X',\n'IN\
9', 'X',\n'INA', 'L',\n'INB', 'X',\n'INC', 'X',\n'\
IND', 'X',\n'INE', 'X',\n'INF', 'F',\n'ING', 'F',\\
n'INH', 'R',\n'INI', 'X',\n'INJ', 'X',\n'INK', 'X'\
,\n'INL', 'X',\n'INM', 'X',\n'INN', 'A',\n'INO', '\
X',\n'INP', 'X',\n'INQ', 'X',\n'INR', 'X',\n'INS',\
 'X',\n'INT', 'V',\n'INU', 'X',\n'INV', 'X',\n'INW\
', 'X',\n'INX', 'X',\n'INY', 'X',\n'INZ', 'X',\n'I\
OA', 'X',\n'IOB', 'X',\n'IOC', 'X',\n'IOD', 'X',\n\
'IOE', 'X',\n'IOF', 'X',\n'IOH', 'X',\n'IOL', 'X',\
\n'IOP', 'X',\n'IP1', 'X',\n'IP2', 'X',\n'IP3', 'X\
',\n'IP4', 'X',\n'IPA', 'X',\n'IPB', 'X',\n'IPD', \
'X',\n'IPG', 'G',\n'IPH', 'X',\n'IPL', 'X',\n'IPM'\
, 'X',\n'IPN', 'X',\n'IPO', 'F',\n'IPP', 'X',\n'IP\
S', 'X',\n'IPT', 'X',\n'IPU', 'X',\n'IPY', 'A',\n'\
IQB', 'X',\n'IQP', 'X',\n'IQS', 'X',\n'IR3', 'X',\\
n'IRI', 'X',\n'IRP', 'X',\n'ISA', 'X',\n'ISF', 'X'\
,\n'ISO', 'X',\n'ISP', 'X',\n'ISQ', 'X',\n'ISU', '\
X',\n'ITM', 'X',\n'ITP', 'X',\n'ITR', 'W',\n'ITS',\
 'X',\n'ITU', 'X',\n'IU5', 'X',\n'IUM', 'X',\n'IUR\
', 'X',\n'IVA', 'X',\n'IYG', 'G',\n'IYR', 'Y',\n'J\
77', 'X',\n'J78', 'X',\n'J80', 'X',\n'JE2', 'X',\n\
'JEN', 'X',\n'JST', 'X',\n'K21', 'X',\n'KAH', 'X',\
\n'KAI', 'X',\n'KAM', 'X',\n'KAN', 'X',\n'KAP', 'X\
',\n'KCP', 'X',\n'KCX', 'K',\n'KDO', 'X',\n'KEF', \
'X',\n'KET', 'X',\n'KGR', 'X',\n'KH1', 'X',\n'KIF'\
, 'X',\n'KIV', 'V',\n'KNI', 'X',\n'KPH', 'K',\n'KT\
H', 'X',\n'KTN', 'X',\n'KTP', 'X',\n'KWT', 'X',\n'\
L04', 'X',\n'L1P', 'X',\n'L24', 'E',\n'L2P', 'X',\\
n'L34', 'E',\n'L37', 'E',\n'L3P', 'X',\n'L4P', 'X'\
,\n'L75', 'X',\n'LAC', 'X',\n'LAD', 'X',\n'LAK', '\
X',\n'LAM', 'X',\n'LAR', 'X',\n'LAT', 'X',\n'LAX',\
 'X',\n'LCO', 'X',\n'LCP', 'X',\n'LCS', 'X',\n'LDA\
', 'X',\n'LDO', 'L',\n'LDP', 'X',\n'LEA', 'X',\n'L\
EO', 'X',\n'LEU', 'L',\n'LG2', 'X',\n'LG6', 'X',\n\
'LGC', 'X',\n'LGP', 'X',\n'LHG', 'X',\n'LHY', 'F',\
\n'LI1', 'X',\n'LIG', 'X',\n'LIL', 'X',\n'LIM', 'X\
',\n'LIN', 'X',\n'LIO', 'X',\n'LIP', 'X',\n'LLA', \
'X',\n'LLP', 'K',\n'LLY', 'K',\n'LMG', 'X',\n'LML'\
, 'X',\n'LMT', 'X',\n'LMU', 'X',\n'LMZ', 'X',\n'LN\
K', 'X',\n'LNL', 'X',\n'LNO', 'X',\n'LOF', 'X',\n'\
LOL', 'L',\n'LOM', 'X',\n'LOR', 'X',\n'LOS', 'X',\\
n'LOV', 'L',\n'LOX', 'X',\n'LP1', 'X',\n'LP2', 'R'\
,\n'LPA', 'X',\n'LPC', 'X',\n'LPF', 'X',\n'LPL', '\
X',\n'LPM', 'X',\n'LPP', 'X',\n'LRB', 'X',\n'LRU',\
 'X',\n'LS1', 'X',\n'LS2', 'X',\n'LS3', 'X',\n'LS4\
', 'X',\n'LS5', 'X',\n'LTA', 'X',\n'LTL', 'X',\n'L\
TR', 'W',\n'LUM', 'X',\n'LVS', 'L',\n'LXC', 'X',\n\
'LY2', 'X',\n'LY3', 'X',\n'LYA', 'X',\n'LYB', 'X',\
\n'LYC', 'X',\n'LYD', 'X',\n'LYM', 'K',\n'LYN', 'X\
',\n'LYS', 'K',\n'LYT', 'X',\n'LYW', 'X',\n'LYZ', \
'K',\n'M1A', 'X',\n'M1G', 'X',\n'M2G', 'X',\n'M3L'\
, 'K',\n'M6P', 'X',\n'M6T', 'X',\n'M7G', 'X',\n'MA\
1', 'X',\n'MA2', 'X',\n'MA3', 'X',\n'MA4', 'X',\n'\
MA6', 'X',\n'MAA', 'A',\n'MAB', 'X',\n'MAC', 'X',\\
n'MAE', 'X',\n'MAG', 'X',\n'MAH', 'X',\n'MAI', 'R'\
,\n'MAK', 'X',\n'MAL', 'X',\n'MAM', 'X',\n'MAN', '\
X',\n'MAO', 'X',\n'MAP', 'X',\n'MAR', 'X',\n'MAS',\
 'X',\n'MAT', 'X',\n'MAU', 'X',\n'MAZ', 'X',\n'MBA\
', 'X',\n'MBD', 'X',\n'MBG', 'X',\n'MBH', 'X',\n'M\
BN', 'X',\n'MBO', 'X',\n'MBR', 'X',\n'MBS', 'X',\n\
'MBV', 'X',\n'MBZ', 'X',\n'MCA', 'X',\n'MCD', 'X',\
\n'MCE', 'X',\n'MCG', 'G',\n'MCI', 'X',\n'MCN', 'X\
',\n'MCP', 'X',\n'MCT', 'X',\n'MCY', 'X',\n'MD2', \
'X',\n'MDA', 'X',\n'MDC', 'X',\n'MDG', 'X',\n'MDH'\
, 'X',\n'MDL', 'X',\n'MDM', 'X',\n'MDN', 'X',\n'MD\
P', 'X',\n'ME6', 'X',\n'MEB', 'X',\n'MEC', 'X',\n'\
MEL', 'X',\n'MEN', 'N',\n'MEP', 'X',\n'MER', 'X',\\
n'MES', 'X',\n'MET', 'M',\n'MEV', 'X',\n'MF2', 'X'\
,\n'MF3', 'M',\n'MFB', 'X',\n'MFD', 'X',\n'MFU', '\
X',\n'MG7', 'X',\n'MGA', 'X',\n'MGB', 'X',\n'MGD',\
 'X',\n'MGG', 'R',\n'MGL', 'X',\n'MGN', 'Q',\n'MGO\
', 'X',\n'MGP', 'X',\n'MGR', 'X',\n'MGS', 'X',\n'M\
GT', 'X',\n'MGU', 'X',\n'MGY', 'G',\n'MHB', 'X',\n\
'MHF', 'X',\n'MHL', 'L',\n'MHM', 'X',\n'MHO', 'M',\
\n'MHS', 'H',\n'MHZ', 'X',\n'MIA', 'X',\n'MIC', 'X\
',\n'MID', 'X',\n'MIL', 'X',\n'MIM', 'X',\n'MIN', \
'G',\n'MIP', 'X',\n'MIS', 'S',\n'MIT', 'X',\n'MJI'\
, 'X',\n'MK1', 'X',\n'MKC', 'X',\n'MLA', 'X',\n'ML\
C', 'X',\n'MLE', 'L',\n'MLN', 'X',\n'MLT', 'X',\n'\
MLY', 'K',\n'MLZ', 'K',\n'MM3', 'X',\n'MM4', 'X',\\
n'MMA', 'X',\n'MMC', 'X',\n'MME', 'M',\n'MMO', 'R'\
,\n'MMP', 'X',\n'MMQ', 'X',\n'MMT', 'X',\n'MN1', '\
X',\n'MN2', 'X',\n'MN3', 'X',\n'MN5', 'X',\n'MN7',\
 'X',\n'MN8', 'X',\n'MNA', 'X',\n'MNB', 'X',\n'MNC\
', 'X',\n'MNG', 'X',\n'MNL', 'L',\n'MNO', 'X',\n'M\
NP', 'X',\n'MNQ', 'X',\n'MNS', 'X',\n'MNT', 'X',\n\
'MNV', 'V',\n'MO1', 'X',\n'MO2', 'X',\n'MO3', 'X',\
\n'MO4', 'X',\n'MO5', 'X',\n'MO6', 'X',\n'MOA', 'X\
',\n'MOB', 'X',\n'MOC', 'X',\n'MOE', 'X',\n'MOG', \
'X',\n'MOH', 'X',\n'MOL', 'X',\n'MOO', 'X',\n'MOP'\
, 'X',\n'MOR', 'X',\n'MOS', 'X',\n'MOT', 'X',\n'MO\
X', 'X',\n'MP1', 'X',\n'MP3', 'X',\n'MPA', 'X',\n'\
MPB', 'X',\n'MPC', 'X',\n'MPD', 'X',\n'MPG', 'X',\\
n'MPH', 'M',\n'MPI', 'X',\n'MPJ', 'M',\n'MPL', 'X'\
,\n'MPN', 'X',\n'MPO', 'X',\n'MPP', 'X',\n'MPQ', '\
G',\n'MPR', 'X',\n'MPS', 'X',\n'MQ0', 'X',\n'MQ7',\
 'X',\n'MQ8', 'X',\n'MQ9', 'X',\n'MQI', 'X',\n'MR2\
', 'X',\n'MRC', 'X',\n'MRM', 'X',\n'MRP', 'X',\n'M\
S2', 'X',\n'MSA', 'X',\n'MSB', 'X',\n'MSD', 'X',\n\
'MSE', 'M',\n'MSF', 'X',\n'MSI', 'X',\n'MSO', 'M',\
\n'MSQ', 'X',\n'MST', 'X',\n'MSU', 'X',\n'MTA', 'X\
',\n'MTB', 'X',\n'MTC', 'X',\n'MTD', 'X',\n'MTE', \
'X',\n'MTF', 'X',\n'MTG', 'X',\n'MTO', 'X',\n'MTS'\
, 'X',\n'MTT', 'X',\n'MTX', 'X',\n'MTY', 'Y',\n'MU\
G', 'X',\n'MUP', 'X',\n'MUR', 'X',\n'MVA', 'V',\n'\
MW1', 'X',\n'MW2', 'X',\n'MXA', 'X',\n'MXY', 'X',\\
n'MYA', 'X',\n'MYC', 'X',\n'MYG', 'X',\n'MYR', 'X'\
,\n'MYS', 'X',\n'MYT', 'X',\n'MZM', 'X',\n'N1T', '\
X',\n'N25', 'X',\n'N2B', 'X',\n'N3T', 'X',\n'N4B',\
 'X',\n'NA2', 'X',\n'NA5', 'X',\n'NA6', 'X',\n'NAA\
', 'X',\n'NAB', 'X',\n'NAC', 'X',\n'NAD', 'X',\n'N\
AE', 'X',\n'NAF', 'X',\n'NAG', 'X',\n'NAH', 'X',\n\
'NAI', 'X',\n'NAL', 'A',\n'NAM', 'A',\n'NAN', 'X',\
\n'NAO', 'X',\n'NAP', 'X',\n'NAQ', 'X',\n'NAR', 'X\
',\n'NAS', 'X',\n'NAU', 'X',\n'NAV', 'X',\n'NAW', \
'X',\n'NAX', 'X',\n'NAY', 'X',\n'NBA', 'X',\n'NBD'\
, 'X',\n'NBE', 'X',\n'NBG', 'X',\n'NBN', 'X',\n'NB\
P', 'X',\n'NBS', 'X',\n'NBU', 'X',\n'NCA', 'X',\n'\
NCB', 'A',\n'NCD', 'X',\n'NCH', 'X',\n'NCM', 'X',\\
n'NCN', 'X',\n'NCO', 'X',\n'NCR', 'X',\n'NCS', 'X'\
,\n'ND4', 'X',\n'NDA', 'X',\n'NDC', 'X',\n'NDD', '\
X',\n'NDO', 'X',\n'NDP', 'X',\n'NDT', 'X',\n'NEA',\
 'X',\n'NEB', 'X',\n'NED', 'X',\n'NEM', 'H',\n'NEN\
', 'X',\n'NEO', 'X',\n'NEP', 'H',\n'NEQ', 'X',\n'N\
ES', 'X',\n'NET', 'X',\n'NEV', 'X',\n'NFA', 'F',\n\
'NFE', 'X',\n'NFG', 'X',\n'NFP', 'X',\n'NFS', 'X',\
\n'NG6', 'X',\n'NGA', 'X',\n'NGL', 'X',\n'NGM', 'X\
',\n'NGO', 'X',\n'NGP', 'X',\n'NGT', 'X',\n'NGU', \
'X',\n'NH2', 'X',\n'NH3', 'X',\n'NH4', 'X',\n'NHD'\
, 'X',\n'NHE', 'X',\n'NHM', 'X',\n'NHP', 'X',\n'NH\
R', 'X',\n'NHS', 'X',\n'NI1', 'X',\n'NI2', 'X',\n'\
NIC', 'X',\n'NID', 'X',\n'NIK', 'X',\n'NIO', 'X',\\
n'NIP', 'X',\n'NIT', 'X',\n'NIU', 'X',\n'NIY', 'Y'\
,\n'NLA', 'X',\n'NLE', 'L',\n'NLG', 'X',\n'NLN', '\
L',\n'NLP', 'L',\n'NM1', 'X',\n'NMA', 'A',\n'NMB',\
 'X',\n'NMC', 'G',\n'NMD', 'X',\n'NME', 'X',\n'NMN\
', 'X',\n'NMO', 'X',\n'NMQ', 'X',\n'NMX', 'X',\n'N\
MY', 'X',\n'NNH', 'R',\n'NNO', 'X',\n'NO2', 'X',\n\
'NO3', 'X',\n'NOA', 'X',\n'NOD', 'X',\n'NOJ', 'X',\
\n'NON', 'X',\n'NOP', 'X',\n'NOR', 'X',\n'NOS', 'X\
',\n'NOV', 'X',\n'NOX', 'X',\n'NP3', 'X',\n'NPA', \
'X',\n'NPC', 'X',\n'NPD', 'X',\n'NPE', 'X',\n'NPF'\
, 'X',\n'NPH', 'C',\n'NPI', 'X',\n'NPL', 'X',\n'NP\
N', 'X',\n'NPO', 'X',\n'NPP', 'X',\n'NPT', 'X',\n'\
NPY', 'X',\n'NRG', 'R',\n'NRI', 'X',\n'NS1', 'X',\\
n'NS5', 'X',\n'NSP', 'X',\n'NTA', 'X',\n'NTB', 'X'\
,\n'NTC', 'X',\n'NTH', 'X',\n'NTM', 'X',\n'NTP', '\
X',\n'NTS', 'X',\n'NTU', 'X',\n'NTZ', 'X',\n'NU1',\
 'X',\n'NVA', 'V',\n'NVI', 'X',\n'NVP', 'X',\n'NW1\
', 'X',\n'NYP', 'X',\n'O4M', 'X',\n'OAA', 'X',\n'O\
AI', 'X',\n'OAP', 'X',\n'OAR', 'X',\n'OAS', 'S',\n\
'OBA', 'X',\n'OBN', 'X',\n'OC1', 'X',\n'OC2', 'X',\
\n'OC3', 'X',\n'OC4', 'X',\n'OC5', 'X',\n'OC6', 'X\
',\n'OC7', 'X',\n'OCL', 'X',\n'OCM', 'X',\n'OCN', \
'X',\n'OCO', 'X',\n'OCP', 'X',\n'OCS', 'C',\n'OCT'\
, 'X',\n'OCV', 'K',\n'OCY', 'C',\n'ODA', 'X',\n'OD\
S', 'X',\n'OES', 'X',\n'OET', 'X',\n'OF1', 'X',\n'\
OF2', 'X',\n'OF3', 'X',\n'OFL', 'X',\n'OFO', 'X',\\
n'OHE', 'X',\n'OHO', 'X',\n'OHT', 'X',\n'OIC', 'X'\
,\n'OIP', 'X',\n'OKA', 'X',\n'OLA', 'X',\n'OLE', '\
X',\n'OLI', 'X',\n'OLO', 'X',\n'OMB', 'X',\n'OMC',\
 'X',\n'OMD', 'X',\n'OME', 'X',\n'OMG', 'X',\n'OMP\
', 'X',\n'OMT', 'M',\n'OMU', 'X',\n'ONE', 'X',\n'O\
NL', 'L',\n'ONP', 'X',\n'OPA', 'X',\n'OPD', 'X',\n\
'OPE', 'X',\n'OPG', 'X',\n'OPH', 'X',\n'OPN', 'X',\
\n'OPP', 'X',\n'OPR', 'R',\n'ORN', 'X',\n'ORO', 'X\
',\n'ORP', 'X',\n'OSB', 'X',\n'OSS', 'X',\n'OTA', \
'X',\n'OTB', 'X',\n'OTE', 'X',\n'OTG', 'X',\n'OUT'\
, 'X',\n'OVA', 'X',\n'OWQ', 'X',\n'OXA', 'X',\n'OX\
E', 'X',\n'OXI', 'X',\n'OXL', 'X',\n'OXM', 'X',\n'\
OXN', 'X',\n'OXO', 'X',\n'OXP', 'X',\n'OXS', 'X',\\
n'OXY', 'X',\n'P11', 'A',\n'P24', 'X',\n'P28', 'X'\
,\n'P2P', 'X',\n'P2U', 'X',\n'P3M', 'X',\n'P4C', '\
X',\n'P4P', 'X',\n'P5P', 'X',\n'P6G', 'X',\n'PA1',\
 'X',\n'PA2', 'X',\n'PA3', 'X',\n'PA4', 'X',\n'PA5\
', 'X',\n'PAA', 'X',\n'PAB', 'X',\n'PAC', 'X',\n'P\
AD', 'X',\n'PAE', 'X',\n'PAG', 'X',\n'PAH', 'X',\n\
'PAI', 'X',\n'PAL', 'D',\n'PAM', 'X',\n'PAN', 'X',\
\n'PAO', 'X',\n'PAP', 'A',\n'PAQ', 'F',\n'PAR', 'X\
',\n'PAS', 'X',\n'PAT', 'W',\n'PBA', 'X',\n'PBB', \
'X',\n'PBC', 'X',\n'PBF', 'F',\n'PBG', 'X',\n'PBI'\
, 'X',\n'PBM', 'X',\n'PBN', 'X',\n'PBP', 'X',\n'PB\
R', 'X',\n'PBZ', 'X',\n'PC2', 'X',\n'PCA', 'E',\n'\
PCB', 'X',\n'PCD', 'X',\n'PCE', 'X',\n'PCG', 'X',\\
n'PCH', 'X',\n'PCL', 'X',\n'PCM', 'X',\n'PCP', 'X'\
,\n'PCR', 'X',\n'PCS', 'X',\n'PCU', 'X',\n'PCV', '\
X',\n'PCY', 'X',\n'PD1', 'X',\n'PDA', 'X',\n'PDC',\
 'X',\n'PDD', 'A',\n'PDE', 'A',\n'PDI', 'X',\n'PDL\
', 'A',\n'PDN', 'X',\n'PDO', 'X',\n'PDP', 'X',\n'P\
DT', 'X',\n'PDU', 'X',\n'PE2', 'X',\n'PE6', 'X',\n\
'PEA', 'X',\n'PEB', 'X',\n'PEC', 'X',\n'PED', 'X',\
\n'PEE', 'X',\n'PEF', 'X',\n'PEG', 'X',\n'PEL', 'X\
',\n'PEO', 'X',\n'PEP', 'X',\n'PEQ', 'X',\n'PER', \
'X',\n'PET', 'X',\n'PFB', 'X',\n'PFC', 'X',\n'PFG'\
, 'X',\n'PFL', 'X',\n'PFM', 'X',\n'PFZ', 'X',\n'PG\
4', 'X',\n'PG5', 'X',\n'PG6', 'X',\n'PGA', 'X',\n'\
PGC', 'X',\n'PGD', 'X',\n'PGE', 'X',\n'PGG', 'G',\\
n'PGH', 'X',\n'PGL', 'X',\n'PGO', 'X',\n'PGP', 'X'\
,\n'PGQ', 'X',\n'PGR', 'X',\n'PGS', 'X',\n'PGU', '\
X',\n'PGX', 'X',\n'PGY', 'G',\n'PH1', 'X',\n'PH2',\
 'X',\n'PH3', 'X',\n'PHA', 'F',\n'PHB', 'X',\n'PHC\
', 'X',\n'PHD', 'X',\n'PHE', 'F',\n'PHG', 'X',\n'P\
HH', 'X',\n'PHI', 'F',\n'PHL', 'F',\n'PHM', 'X',\n\
'PHN', 'X',\n'PHO', 'X',\n'PHP', 'X',\n'PHQ', 'X',\
\n'PHS', 'H',\n'PHT', 'X',\n'PHW', 'P',\n'PHY', 'X\
',\n'PI1', 'X',\n'PI2', 'X',\n'PI3', 'X',\n'PI4', \
'X',\n'PI5', 'X',\n'PI6', 'X',\n'PI7', 'X',\n'PI8'\
, 'X',\n'PI9', 'X',\n'PIA', 'X',\n'PIB', 'X',\n'PI\
C', 'X',\n'PID', 'X',\n'PIG', 'X',\n'PIH', 'X',\n'\
PIM', 'X',\n'PIN', 'X',\n'PIO', 'X',\n'PIP', 'X',\\
n'PIQ', 'X',\n'PIR', 'X',\n'PIV', 'X',\n'PKF', 'X'\
,\n'PL1', 'X',\n'PL9', 'X',\n'PLA', 'D',\n'PLC', '\
X',\n'PLE', 'L',\n'PLG', 'G',\n'PLH', 'X',\n'PLM',\
 'X',\n'PLP', 'X',\n'PLS', 'S',\n'PLT', 'W',\n'PLU\
', 'L',\n'PLY', 'X',\n'PMA', 'X',\n'PMB', 'X',\n'P\
MC', 'X',\n'PME', 'F',\n'PML', 'X',\n'PMM', 'X',\n\
'PMO', 'X',\n'PMP', 'X',\n'PMS', 'X',\n'PMY', 'X',\
\n'PN2', 'X',\n'PNA', 'X',\n'PNB', 'X',\n'PNC', 'G\
',\n'PND', 'X',\n'PNE', 'A',\n'PNF', 'X',\n'PNG', \
'X',\n'PNI', 'X',\n'PNL', 'X',\n'PNM', 'X',\n'PNN'\
, 'X',\n'PNO', 'X',\n'PNP', 'X',\n'PNQ', 'X',\n'PN\
S', 'X',\n'PNT', 'X',\n'PNU', 'X',\n'PO2', 'X',\n'\
PO4', 'X',\n'POB', 'X',\n'POC', 'X',\n'POL', 'X',\\
n'POM', 'P',\n'PON', 'X',\n'POP', 'X',\n'POR', 'X'\
,\n'POS', 'X',\n'PP1', 'X',\n'PP2', 'X',\n'PP3', '\
A',\n'PP4', 'X',\n'PP5', 'X',\n'PP6', 'X',\n'PP7',\
 'X',\n'PP8', 'N',\n'PP9', 'X',\n'PPB', 'X',\n'PPC\
', 'X',\n'PPD', 'X',\n'PPE', 'E',\n'PPG', 'X',\n'P\
PH', 'F',\n'PPI', 'X',\n'PPJ', 'V',\n'PPL', 'X',\n\
'PPM', 'X',\n'PPN', 'A',\n'PPO', 'X',\n'PPP', 'X',\
\n'PPQ', 'X',\n'PPR', 'X',\n'PPS', 'X',\n'PPT', 'X\
',\n'PPU', 'X',\n'PPX', 'F',\n'PPY', 'X',\n'PPZ', \
'X',\n'PQ0', 'X',\n'PQN', 'X',\n'PQQ', 'X',\n'PR1'\
, 'X',\n'PR2', 'X',\n'PR3', 'X',\n'PRA', 'X',\n'PR\
B', 'X',\n'PRC', 'X',\n'PRD', 'X',\n'PRE', 'X',\n'\
PRF', 'X',\n'PRH', 'X',\n'PRI', 'P',\n'PRL', 'X',\\
n'PRN', 'X',\n'PRO', 'P',\n'PRP', 'X',\n'PRR', 'A'\
,\n'PRS', 'P',\n'PRZ', 'X',\n'PS0', 'X',\n'PSA', '\
X',\n'PSD', 'X',\n'PSE', 'X',\n'PSF', 'S',\n'PSG',\
 'X',\n'PSI', 'X',\n'PSO', 'X',\n'PSQ', 'X',\n'PSS\
', 'X',\n'PST', 'X',\n'PSU', 'X',\n'PT1', 'X',\n'P\
T3', 'X',\n'PTA', 'X',\n'PTC', 'X',\n'PTD', 'X',\n\
'PTE', 'X',\n'PTH', 'Y',\n'PTL', 'X',\n'PTM', 'Y',\
\n'PTN', 'X',\n'PTO', 'X',\n'PTP', 'X',\n'PTR', 'Y\
',\n'PTS', 'X',\n'PTT', 'X',\n'PTU', 'X',\n'PTY', \
'X',\n'PUA', 'X',\n'PUB', 'X',\n'PUR', 'X',\n'PUT'\
, 'X',\n'PVA', 'X',\n'PVB', 'X',\n'PVH', 'H',\n'PV\
L', 'X',\n'PXA', 'X',\n'PXF', 'X',\n'PXG', 'X',\n'\
PXP', 'X',\n'PXY', 'X',\n'PXZ', 'X',\n'PY2', 'X',\\
n'PY4', 'X',\n'PY5', 'X',\n'PY6', 'X',\n'PYA', 'A'\
,\n'PYC', 'X',\n'PYD', 'X',\n'PYE', 'X',\n'PYL', '\
X',\n'PYM', 'X',\n'PYO', 'X',\n'PYP', 'X',\n'PYQ',\
 'X',\n'PYR', 'X',\n'PYS', 'X',\n'PYT', 'X',\n'PYX\
', 'X',\n'PYY', 'X',\n'PYZ', 'X',\n'PZQ', 'X',\n'Q\
82', 'X',\n'QNC', 'X',\n'QND', 'X',\n'QSI', 'Q',\n\
'QTR', 'X',\n'QUA', 'X',\n'QUE', 'X',\n'QUI', 'X',\
\n'QUO', 'X',\n'R11', 'X',\n'R12', 'X',\n'R13', 'X\
',\n'R18', 'X',\n'R1P', 'X',\n'R56', 'X',\n'R5P', \
'X',\n'RA2', 'X',\n'RAD', 'X',\n'RAI', 'X',\n'RAL'\
, 'X',\n'RAM', 'X',\n'RAN', 'X',\n'RAP', 'X',\n'RB\
F', 'X',\n'RBU', 'X',\n'RCA', 'X',\n'RCL', 'X',\n'\
RCO', 'X',\n'RDC', 'X',\n'RDF', 'W',\n'RE9', 'X',\\
n'REA', 'X',\n'RED', 'K',\n'REO', 'X',\n'REP', 'X'\
,\n'RET', 'X',\n'RFA', 'X',\n'RFB', 'X',\n'RFL', '\
X',\n'RFP', 'X',\n'RG1', 'X',\n'RGS', 'X',\n'RH1',\
 'X',\n'RHA', 'X',\n'RHC', 'X',\n'RHD', 'X',\n'RHM\
', 'X',\n'RHO', 'X',\n'RHQ', 'X',\n'RHS', 'X',\n'R\
IA', 'X',\n'RIB', 'X',\n'RIC', 'X',\n'RIF', 'X',\n\
'RIN', 'X',\n'RIP', 'X',\n'RIT', 'X',\n'RMB', 'X',\
\n'RMN', 'X',\n'RMP', 'X',\n'RNG', 'X',\n'RNS', 'X\
',\n'RNT', 'X',\n'RO2', 'X',\n'RO4', 'X',\n'ROC', \
'N',\n'ROI', 'X',\n'ROM', 'X',\n'RON', 'V',\n'ROP'\
, 'X',\n'ROS', 'X',\n'ROX', 'X',\n'RPA', 'X',\n'RP\
D', 'X',\n'RPH', 'X',\n'RPL', 'X',\n'RPP', 'X',\n'\
RPR', 'X',\n'RPX', 'X',\n'RQ3', 'X',\n'RR1', 'X',\\
n'RR6', 'X',\n'RRS', 'X',\n'RS1', 'X',\n'RS2', 'X'\
,\n'RS7', 'X',\n'RSS', 'X',\n'RTA', 'X',\n'RTB', '\
X',\n'RTC', 'X',\n'RTL', 'X',\n'RUB', 'X',\n'RUN',\
 'X',\n'RWJ', 'X',\n'RXP', 'X',\n'S02', 'X',\n'S11\
', 'X',\n'S1H', 'S',\n'S27', 'X',\n'S2C', 'C',\n'S\
3P', 'X',\n'S4U', 'X',\n'S57', 'X',\n'S58', 'X',\n\
'S5H', 'X',\n'S6G', 'X',\n'S80', 'X',\n'SAA', 'X',\
\n'SAB', 'X',\n'SAC', 'S',\n'SAD', 'X',\n'SAE', 'X\
',\n'SAF', 'X',\n'SAH', 'C',\n'SAI', 'C',\n'SAL', \
'X',\n'SAM', 'M',\n'SAN', 'X',\n'SAP', 'X',\n'SAR'\
, 'X',\n'SAS', 'X',\n'SB1', 'X',\n'SB2', 'X',\n'SB\
3', 'X',\n'SB4', 'X',\n'SB5', 'X',\n'SB6', 'X',\n'\
SBA', 'L',\n'SBB', 'X',\n'SBD', 'A',\n'SBI', 'X',\\
n'SBL', 'A',\n'SBN', 'X',\n'SBO', 'X',\n'SBR', 'X'\
,\n'SBS', 'X',\n'SBT', 'X',\n'SBU', 'X',\n'SBX', '\
X',\n'SC4', 'X',\n'SCA', 'X',\n'SCC', 'X',\n'SCD',\
 'X',\n'SCH', 'C',\n'SCI', 'X',\n'SCL', 'X',\n'SCM\
', 'X',\n'SCN', 'X',\n'SCO', 'X',\n'SCP', 'S',\n'S\
CR', 'X',\n'SCS', 'X',\n'SCV', 'C',\n'SCY', 'C',\n\
'SD8', 'X',\n'SDK', 'X',\n'SDZ', 'X',\n'SE4', 'X',\
\n'SEA', 'X',\n'SEB', 'S',\n'SEC', 'X',\n'SEG', 'A\
',\n'SEI', 'X',\n'SEL', 'S',\n'SEM', 'X',\n'SEO', \
'X',\n'SEP', 'S',\n'SER', 'S',\n'SES', 'X',\n'SET'\
, 'S',\n'SEU', 'X',\n'SF4', 'X',\n'SFG', 'X',\n'SF\
N', 'X',\n'SFO', 'X',\n'SGA', 'X',\n'SGC', 'X',\n'\
SGL', 'X',\n'SGM', 'X',\n'SGN', 'X',\n'SGP', 'X',\\
n'SHA', 'X',\n'SHC', 'X',\n'SHF', 'X',\n'SHH', 'X'\
,\n'SHP', 'G',\n'SHR', 'E',\n'SHT', 'T',\n'SHU', '\
X',\n'SI2', 'X',\n'SIA', 'X',\n'SIF', 'X',\n'SIG',\
 'X',\n'SIH', 'X',\n'SIM', 'X',\n'SIN', 'X',\n'SKD\
', 'X',\n'SKF', 'X',\n'SLB', 'X',\n'SLE', 'X',\n'S\
LZ', 'K',\n'SMA', 'X',\n'SMC', 'C',\n'SME', 'M',\n\
'SML', 'X',\n'SMM', 'M',\n'SMN', 'X',\n'SMP', 'X',\
\n'SMS', 'X',\n'SN1', 'X',\n'SN6', 'X',\n'SN7', 'X\
',\n'SNC', 'C',\n'SNN', 'X',\n'SNP', 'X',\n'SO1', \
'X',\n'SO2', 'X',\n'SO3', 'X',\n'SO4', 'X',\n'SOA'\
, 'X',\n'SOC', 'C',\n'SOM', 'X',\n'SOR', 'X',\n'SO\
T', 'X',\n'SOX', 'X',\n'SPA', 'X',\n'SPB', 'X',\n'\
SPC', 'X',\n'SPD', 'X',\n'SPE', 'X',\n'SPG', 'X',\\
n'SPH', 'X',\n'SPI', 'X',\n'SPK', 'X',\n'SPM', 'X'\
,\n'SPN', 'X',\n'SPO', 'X',\n'SPP', 'X',\n'SPS', '\
X',\n'SPY', 'X',\n'SQU', 'X',\n'SRA', 'X',\n'SRB',\
 'X',\n'SRD', 'X',\n'SRL', 'X',\n'SRM', 'X',\n'SRS\
', 'X',\n'SRY', 'X',\n'SSA', 'X',\n'SSB', 'X',\n'S\
SG', 'X',\n'SSP', 'X',\n'ST1', 'X',\n'ST2', 'X',\n\
'ST3', 'X',\n'ST4', 'X',\n'ST5', 'X',\n'ST6', 'X',\
\n'STA', 'X',\n'STB', 'X',\n'STE', 'X',\n'STG', 'X\
',\n'STI', 'X',\n'STL', 'X',\n'STN', 'X',\n'STO', \
'X',\n'STP', 'X',\n'STR', 'X',\n'STU', 'X',\n'STY'\
, 'Y',\n'SU1', 'X',\n'SU2', 'X',\n'SUC', 'X',\n'SU\
I', 'X',\n'SUL', 'X',\n'SUR', 'X',\n'SVA', 'S',\n'\
SWA', 'X',\n'T16', 'X',\n'T19', 'X',\n'T23', 'X',\\
n'T29', 'X',\n'T33', 'X',\n'T3P', 'X',\n'T42', 'A'\
,\n'T44', 'X',\n'T5A', 'X',\n'T6A', 'T',\n'T6P', '\
X',\n'T80', 'X',\n'T87', 'X',\n'TA1', 'X',\n'TAA',\
 'X',\n'TAB', 'X',\n'TAC', 'X',\n'TAD', 'X',\n'TAF\
', 'X',\n'TAM', 'X',\n'TAP', 'X',\n'TAR', 'X',\n'T\
AS', 'X',\n'TAU', 'X',\n'TAX', 'X',\n'TAZ', 'X',\n\
'TB9', 'X',\n'TBA', 'X',\n'TBD', 'X',\n'TBG', 'G',\
\n'TBH', 'X',\n'TBM', 'T',\n'TBO', 'X',\n'TBP', 'X\
',\n'TBR', 'X',\n'TBS', 'X',\n'TBT', 'X',\n'TBU', \
'X',\n'TBZ', 'X',\n'TC4', 'X',\n'TCA', 'X',\n'TCB'\
, 'X',\n'TCH', 'X',\n'TCK', 'X',\n'TCL', 'X',\n'TC\
M', 'X',\n'TCN', 'X',\n'TCP', 'X',\n'TCR', 'W',\n'\
TCS', 'X',\n'TCZ', 'X',\n'TDA', 'X',\n'TDB', 'X',\\
n'TDG', 'X',\n'TDP', 'X',\n'TDR', 'X',\n'TDX', 'X'\
,\n'TEA', 'X',\n'TEM', 'X',\n'TEN', 'X',\n'TEO', '\
X',\n'TEP', 'X',\n'TER', 'X',\n'TES', 'X',\n'TET',\
 'X',\n'TFA', 'X',\n'TFB', 'X',\n'TFH', 'X',\n'TFI\
', 'X',\n'TFK', 'X',\n'TFP', 'X',\n'THA', 'X',\n'T\
HB', 'X',\n'THC', 'T',\n'THD', 'X',\n'THE', 'X',\n\
'THF', 'X',\n'THJ', 'X',\n'THK', 'X',\n'THM', 'X',\
\n'THN', 'X',\n'THO', 'T',\n'THP', 'X',\n'THQ', 'X\
',\n'THR', 'T',\n'THS', 'X',\n'THT', 'X',\n'THU', \
'X',\n'THX', 'X',\n'THZ', 'X',\n'TI1', 'X',\n'TI2'\
, 'X',\n'TI3', 'P',\n'TIA', 'X',\n'TIH', 'A',\n'TK\
4', 'X',\n'TLA', 'X',\n'TLC', 'X',\n'TLM', 'X',\n'\
TLN', 'X',\n'TLX', 'X',\n'TM5', 'X',\n'TM6', 'X',\\
n'TMA', 'X',\n'TMB', 'T',\n'TMC', 'X',\n'TMD', 'T'\
,\n'TME', 'X',\n'TMF', 'X',\n'TML', 'K',\n'TMM', '\
X',\n'TMN', 'X',\n'TMP', 'X',\n'TMQ', 'X',\n'TMR',\
 'X',\n'TMT', 'X',\n'TMZ', 'X',\n'TNB', 'C',\n'TND\
', 'X',\n'TNK', 'X',\n'TNP', 'X',\n'TNT', 'X',\n'T\
OA', 'X',\n'TOB', 'X',\n'TOC', 'X',\n'TOL', 'X',\n\
'TOP', 'X',\n'TOS', 'X',\n'TOT', 'X',\n'TP1', 'G',\
\n'TP2', 'P',\n'TP3', 'E',\n'TP4', 'E',\n'TP7', 'T\
',\n'TPA', 'X',\n'TPE', 'X',\n'TPF', 'X',\n'TPI', \
'X',\n'TPL', 'W',\n'TPM', 'X',\n'TPN', 'G',\n'TPO'\
, 'T',\n'TPP', 'X',\n'TPQ', 'A',\n'TPR', 'P',\n'TP\
S', 'X',\n'TPT', 'X',\n'TPV', 'X',\n'TPX', 'X',\n'\
TPY', 'X',\n'TQ3', 'X',\n'TQ4', 'X',\n'TQ5', 'X',\\
n'TQ6', 'X',\n'TR1', 'X',\n'TRA', 'X',\n'TRB', 'X'\
,\n'TRC', 'X',\n'TRD', 'X',\n'TRE', 'X',\n'TRF', '\
W',\n'TRG', 'K',\n'TRH', 'X',\n'TRI', 'X',\n'TRJ',\
 'X',\n'TRM', 'X',\n'TRN', 'W',\n'TRO', 'W',\n'TRP\
', 'W',\n'TRQ', 'X',\n'TRS', 'X',\n'TRX', 'W',\n'T\
RZ', 'X',\n'TS2', 'X',\n'TS3', 'X',\n'TS4', 'X',\n\
'TS5', 'X',\n'TSA', 'X',\n'TSB', 'X',\n'TSI', 'X',\
\n'TSM', 'X',\n'TSN', 'X',\n'TSP', 'X',\n'TSU', 'X\
',\n'TTA', 'X',\n'TTE', 'X',\n'TTN', 'X',\n'TTO', \
'X',\n'TTP', 'X',\n'TTX', 'X',\n'TXL', 'X',\n'TYA'\
, 'Y',\n'TYB', 'Y',\n'TYD', 'X',\n'TYI', 'Y',\n'TY\
L', 'X',\n'TYM', 'W',\n'TYN', 'Y',\n'TYQ', 'Y',\n'\
TYR', 'Y',\n'TYS', 'Y',\n'TYV', 'X',\n'TYY', 'A',\\
n'TZB', 'X',\n'TZC', 'X',\n'TZE', 'X',\n'TZL', 'X'\
,\n'TZO', 'X',\n'TZP', 'X',\n'U01', 'X',\n'U02', '\
X',\n'U03', 'X',\n'U04', 'X',\n'U05', 'X',\n'U0E',\
 'X',\n'U10', 'X',\n'U18', 'X',\n'U2G', 'X',\n'U3P\
', 'X',\n'U49', 'X',\n'U55', 'X',\n'U5P', 'X',\n'U\
66', 'X',\n'U89', 'X',\n'U8U', 'X',\n'UAA', 'X',\n\
'UAG', 'A',\n'UAP', 'X',\n'UAR', 'X',\n'UC1', 'X',\
\n'UC2', 'X',\n'UC3', 'X',\n'UC4', 'X',\n'UD1', 'X\
',\n'UD2', 'X',\n'UDP', 'X',\n'UDX', 'X',\n'UFG', \
'X',\n'UFM', 'X',\n'UFP', 'X',\n'UGA', 'X',\n'UIN'\
, 'X',\n'UKP', 'A',\n'UM3', 'X',\n'UMA', 'A',\n'UM\
G', 'X',\n'UMP', 'X',\n'UNA', 'X',\n'UND', 'X',\n'\
UNI', 'X',\n'UNK', 'X',\n'UNN', 'X',\n'UNX', 'X',\\
n'UP5', 'X',\n'UP6', 'X',\n'UPA', 'X',\n'UPF', 'X'\
,\n'UPG', 'X',\n'UPP', 'X',\n'UQ1', 'X',\n'UQ2', '\
X',\n'UQ6', 'X',\n'UR2', 'X',\n'URA', 'X',\n'URE',\
 'X',\n'URF', 'X',\n'URI', 'X',\n'URS', 'X',\n'UTP\
', 'X',\n'UVC', 'X',\n'UVW', 'X',\n'V35', 'X',\n'V\
36', 'X',\n'V4O', 'X',\n'V7O', 'X',\n'VAA', 'V',\n\
'VAC', 'X',\n'VAD', 'V',\n'VAF', 'V',\n'VAG', 'X',\
\n'VAL', 'V',\n'VAN', 'X',\n'VAS', 'X',\n'VAX', 'X\
',\n'VDX', 'X',\n'VDY', 'X',\n'VG1', 'X',\n'VIB', \
'X',\n'VIR', 'X',\n'VIT', 'X',\n'VK3', 'X',\n'VO3'\
, 'X',\n'VO4', 'X',\n'VS1', 'F',\n'VS2', 'F',\n'VS\
3', 'F',\n'VS4', 'F',\n'VXA', 'X',\n'W01', 'X',\n'\
W02', 'X',\n'W03', 'X',\n'W11', 'X',\n'W33', 'X',\\
n'W35', 'X',\n'W42', 'X',\n'W43', 'X',\n'W54', 'X'\
,\n'W56', 'X',\n'W59', 'X',\n'W71', 'X',\n'W84', '\
X',\n'W8R', 'X',\n'W91', 'X',\n'WAY', 'X',\n'WCC',\
 'X',\n'WO2', 'X',\n'WO4', 'X',\n'WRB', 'X',\n'WRR\
', 'X',\n'WRS', 'X',\n'WW7', 'X',\n'X2F', 'X',\n'X\
7O', 'X',\n'XAA', 'X',\n'XAN', 'X',\n'XAO', 'X',\n\
'XBB', 'X',\n'XBP', 'X',\n'XDN', 'X',\n'XDP', 'X',\
\n'XIF', 'X',\n'XIM', 'X',\n'XK2', 'X',\n'XL1', 'X\
',\n'XLS', 'X',\n'XMP', 'X',\n'XN1', 'X',\n'XN2', \
'X',\n'XN3', 'X',\n'XUL', 'X',\n'XV6', 'X',\n'XYD'\
, 'X',\n'XYH', 'X',\n'XYL', 'X',\n'XYP', 'X',\n'XY\
S', 'X',\n'YOF', 'Y',\n'YRR', 'X',\n'YT3', 'X',\n'\
YZ9', 'X',\n'Z34', 'G',\n'Z5A', 'X',\n'ZAF', 'X',\\
n'ZAP', 'X',\n'ZEB', 'X',\n'ZEN', 'X',\n'ZES', 'X'\
,\n'ZID', 'X',\n'ZMR', 'X',\n'ZN3', 'X',\n'ZNH', '\
X',\n'ZNO', 'X',\n'ZO3', 'X',\n'ZPR', 'P',\n'ZRA',\
 'A',\n'ZST', 'X',\n'ZYA', 'A',\n\n\n'ASN','N');\n\
} \n\n\nsub file2head\n      {\n	my $file = shift;\
\n	my $size = shift;\n	my $f= new FileHandle;\n	my\
 $line;\n	open ($f,$file);\n	read ($f,$line, $size\
);\n	close ($f);\n	return $line;\n      }\nsub fil\
e2tail\n      {\n	my $file = shift;\n	my $size = s\
hift;\n	my $f= new FileHandle;\n	my $line;\n	\n	op\
en ($f,$file);\n	seek ($f,$size*-1, 2);\n	read ($f\
,$line, $size);\n	close ($f);\n	return $line;\n   \
   }\n\n\nsub vtmpnam\n      {\n	my $r=rand(100000\
);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n\
	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f\
);\n	return $f;\n      }\n\nsub myexit\n  {\n    m\
y $code=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){re\
turn;}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ###\
 ONLY BARE EXIT\n    exit ($code);\n  }\nsub set_e\
rror_lock\n    {\n      my $name = shift;\n      m\
y $pid=$$;\n\n      \n      &lock4tc ($$,\"LERROR\\
", \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\\
n      return;\n    }\nsub set_lock\n  {\n    my $\
pid=shift;\n    my $msg= shift;\n    my $p=getppid\
();\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$\
msg\\n\");\n  }\nsub unset_lock\n   {\n     \n    \
my $pid=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRE\
LEASE\",\"\");\n  }\nsub shift_lock\n  {\n    my $\
from=shift;\n    my $to=shift;\n    my $from_type=\
shift;\n    my $to_type=shift;\n    my $action=shi\
ft;\n    my $msg;\n    \n    if (!&lock4tc($from, \
$from_type, \"LCHECK\", \"\")){return 0;}\n    $ms\
g=&lock4tc ($from, $from_type, \"LREAD\", \"\");\n\
    &lock4tc ($from, $from_type,\"LRELEASE\", $msg\
);\n    &lock4tc ($to, $to_type, $action, $msg);\n\
    return;\n  }\nsub isshellpid\n  {\n    my $p=s\
hift;\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\"\
)){return 0;}\n    else\n      {\n	my $c=lock4tc($\
p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){re\
turn 1;}\n      }\n    return 0;\n  }\nsub isrootp\
id\n  {\n    if(lock4tc (getppid(), \"LLOCK\", \"L\
CHECK\")){return 0;}\n    else {return 1;}\n  }\ns\
ub lock4tc\n	{\n	  my ($pid,$type,$action,$value)=\
@_;\n	  my $fname;\n	  my $host=hostname;\n	  \n	 \
 if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$\
host.lock4tcoffee\";}\n	  elsif ( $type eq \"LERRO\
R\"){ $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\\
";}\n	  elsif ( $type eq \"LWARNING\"){ $fname=\"$\
LOCKDIR/.$pid.$host.warning4tcoffee\";}\n	  \n	  i\
f ($debug_lock)\n	    {\n	      print STDERR \"\\n\
\\t---lock4tc(tcg): $action => $fname =>$value (RD\
: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($action eq \
\"LCHECK\") {return -e $fname;}\n	  elsif ($action\
 eq \"LREAD\"){return file2string($fname);}\n	  el\
sif ($action eq \"LSET\") {return string2file ($va\
lue, $fname, \">>\");}\n	  elsif ($action eq \"LRE\
SET\") {return string2file ($value, $fname, \">\")\
;}\n	  elsif ($action eq \"LRELEASE\") \n	    {\n	\
      if ( $debug_lock)\n		{\n		  my $g=new FileHa\
ndle;\n		  open ($g, \">>$fname\");\n		  print $g \
\"\\nDestroyed by $$\\n\";\n		  close ($g);\n		  s\
afe_system (\"mv $fname $fname.old\");\n		}\n	    \
  else\n		{\n		  unlink ($fname);\n		}\n	    }\n	 \
 return \"\";\n	}\n	\nsub file2string\n	{\n	  my $\
file=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\n\
	  open ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\
\n	  close ($f);\n	  return $r;\n	}\nsub string2fi\
le \n    {\n    my ($s,$file,$mode)=@_;\n    my $f\
=new FileHandle;\n    \n    open ($f, \"$mode$file\
\");\n    print $f  \"$s\";\n    close ($f);\n  }\\
n\nBEGIN\n    {\n      srand;\n    \n      $SIG{'S\
IGUP'}='signal_cleanup';\n      $SIG{'SIGINT'}='si\
gnal_cleanup';\n      $SIG{'SIGQUIT'}='signal_clea\
nup';\n      $SIG{'SIGILL'}='signal_cleanup';\n   \
   $SIG{'SIGTRAP'}='signal_cleanup';\n      $SIG{'\
SIGABRT'}='signal_cleanup';\n      $SIG{'SIGEMT'}=\
'signal_cleanup';\n      $SIG{'SIGFPE'}='signal_cl\
eanup';\n      \n      $SIG{'SIGKILL'}='signal_cle\
anup';\n      $SIG{'SIGPIPE'}='signal_cleanup';\n \
     $SIG{'SIGSTOP'}='signal_cleanup';\n      $SIG\
{'SIGTTIN'}='signal_cleanup';\n      $SIG{'SIGXFSZ\
'}='signal_cleanup';\n      $SIG{'SIGINFO'}='signa\
l_cleanup';\n      \n      $SIG{'SIGBUS'}='signal_\
cleanup';\n      $SIG{'SIGALRM'}='signal_cleanup';\
\n      $SIG{'SIGTSTP'}='signal_cleanup';\n      $\
SIG{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGV\
TALRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='\
signal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signa\
l_cleanup';\n      $SIG{'SIGTERM'}='signal_cleanup\
';\n      $SIG{'SIGCONT'}='signal_cleanup';\n     \
 $SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGP\
ROF'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='si\
gnal_cleanup';\n\n      $SIG{'SIGSYS'}='signal_cle\
anup';\n      $SIG{'SIGURG'}='signal_cleanup';\n  \
    $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{\
'SIGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINCH\
'}='signal_cleanup';\n      \n      $SIG{'INT'}='s\
ignal_cleanup';\n      $SIG{'TERM'}='signal_cleanu\
p';\n      $SIG{'KILL'}='signal_cleanup';\n      $\
SIG{'QUIT'}='signal_cleanup';\n      \n      our $\
debug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n\
      \n      \n      foreach my $a (@ARGV){$CL.=\\
" $a\";}\n      if ( $debug_lock ){print STDERR \"\
\\n\\n\\n********** START PG: $PROGRAM ***********\
**\\n\";}\n      if ( $debug_lock ){print STDERR \\
"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ ***\
**********\\n\";}\n      if ( $debug_lock ){print \
STDERR \"\\n --- $$ -- $CL\\n\";}\n      \n	     \\
n      \n      \n    }\nsub flush_error\n  {\n    \
my $msg=shift;\n    return add_error ($EXIT_FAILUR\
E,$$, $$,getppid(), $msg, $CL);\n  }\nsub add_erro\
r \n  {\n    my $code=shift;\n    my $rpid=shift;\\
n    my $pid=shift;\n    my $ppid=shift;\n    my $\
type=shift;\n    my $com=shift;\n    \n    $ERROR_\
DONE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\\
"$pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"L\
ERROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n    \
lock4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK:\
 $ppid -> $pid\\n\");\n   \n    return $code;\n  }\
\nsub add_warning \n  {\n    my $rpid=shift;\n    \
my $pid =shift;\n    my $command=shift;\n    my $m\
sg=\"$$ -- WARNING: $command\\n\";\n    print STDE\
RR \"$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSE\
T\", $msg);\n  }\n\nsub signal_cleanup\n  {\n    p\
rint dtderr \"\\n**** $$ (tcg) was killed\\n\";\n \
   &cleanup;\n    exit ($EXIT_FAILURE);\n  }\nsub \
clean_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $\
dir){return ;}\n    elsif (!($dir=~/tmp/)){return \
;}#safety check 1\n    elsif (($dir=~/\\*/)){retur\
n ;}#safety check 2\n    else\n      {\n	`rm -rf $\
dir`;\n      }\n    return;\n  }\nsub cleanup\n  {\
\n    #print stderr \"\\n----tc: $$ Kills $PIDCHIL\
D\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p\
=getppid();\n    $CLEAN_EXIT_STARTED=1;\n    \n   \
 \n    \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\\
", \"\"))\n      {\n	my $ppid=getppid();\n	if (!$E\
RROR_DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"\
LSET\", \"$$ -- STACK: $p -> $$\\n\");\n	    &lock\
4tc($$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\"\
);\n	  }\n      }\n    my $warning=&lock4tc($$, \"\
LWARNING\", \"LREAD\", \"\");\n    my $error=&lock\
4tc($$,  \"LERROR\", \"LREAD\", \"\");\n    #relea\
se error and warning lock if root\n    \n    if (i\
srootpid() && ($warning || $error) )\n      {\n	\n\
	print STDERR \"**************** Summary *********\
****\\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"L\
ERROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNIN\
G\",\"RELEASE\",\"\");\n      } \n    \n    \n    \
foreach my $f (@TMPFILE_LIST)\n      {\n	if (-e $f\
){unlink ($f);} \n      }\n    foreach my $d (@TMP\
DIR_LIST)\n      {\n	clean_dir ($d);\n      }\n   \
 #No More Lock Release\n    #&lock4tc($$,\"LLOCK\"\
,\"LRELEASE\",\"\"); #release lock \n\n    if ( $d\
ebug_lock ){print STDERR \"\\n\\n\\n********** END\
 PG: $PROGRAM ($$) *************\\n\";}\n    if ( \
$debug_lock ){print STDERR \"\\n\\n\\n**********(t\
cg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  }\
\nEND \n  {\n    \n    &cleanup();\n  }\n   \n\nsu\
b safe_system \n{\n  my $com=shift;\n  my $ntry=sh\
ift;\n  my $ctry=shift;\n  my $pid;\n  my $status;\
\n  my $ppid=getppid();\n  if ($com eq \"\"){retur\
n 1;}\n  \n  \n\n  if (($pid = fork ()) < 0){retur\
n (-1);}\n  if ($pid == 0)\n    {\n      set_lock(\
$$, \" -SHELL- $com (tcg)\");\n      exec ($com);\\
n    }\n  else\n    {\n      lock4tc ($$, \"LLOCK\\
", \"LSET\", \"$pid\\n\");#update parent\n      $P\
IDCHILD=$pid;\n    }\n  if ($debug_lock){printf ST\
DERR \"\\n\\t .... safe_system (fasta_seq2hmm)  p:\
 $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,WT\
ERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"\
LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILUR\
E || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\\
n    {\n      if ($ntry && $ctry <$ntry)\n	{\n	  a\
dd_warning ($$,$$,\"$com failed [retry: $ctry]\");\
\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\"\
);\n	  return safe_system ($com, $ntry, ++$ctry);\\
n	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift_\
lock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"\
))\n	    {\n	      add_warning ($$,$$,\"$com faile\
d\");\n	    }\n	  else\n	    {\n	      lock4tc ($p\
id, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  r\
eturn $?;}\n      else\n	{\n	  if (!shift_lock ($p\
id,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n\
	      myexit(add_error ($EXIT_FAILURE,$$,$pid,get\
ppid(), \"UNSPECIFIED system\", $com));\n	    }\n	\
}\n    }\n  return $?;\n}\n\nsub check_configurati\
on \n    {\n      my @l=@_;\n      my $v;\n      f\
oreach my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAI\
L\")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n	\
	add_warning($$,$$,\"Could Not Use EMAIL\");\n		my\
exit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"EM\
AIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p e\
q \"INTERNET\")\n	    {\n	      if ( !&check_inter\
net_connection())\n		{\n		  myexit(add_error ($EXI\
T_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\\
n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n	\
      if (!&pg_is_installed (\"wget\") && !&pg_is_\
installed (\"curl\"))\n		{\n		  myexit(add_error (\
$EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:w\
get\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_\
installed ($p)))\n	    {\n	      myexit(add_error \
($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:\
$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n    \
}\nsub pg_is_installed\n  {\n    my @ml=@_;\n    m\
y $r, $p, $m;\n    my $supported=0;\n    \n    my \
$p=shift (@ml);\n    if ($p=~/::/)\n      {\n	if (\
safe_system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){r\
eturn 1;}\n	else {return 0;}\n      }\n    else\n \
     {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\
\"){return 0;}\n	else {return 1;}\n      }\n  }\n\\
n\n\nsub check_internet_connection\n  {\n    my $i\
nternet;\n    my $tmp;\n    &check_configuration (\
 \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n \
   if     (&pg_is_installed    (\"wget\")){`wget w\
ww.google.com -O$tmp >/dev/null 2>/dev/null`;}\n  \
  elsif  (&pg_is_installed    (\"curl\")){`curl ww\
w.google.com -o$tmp >/dev/null 2>/dev/null`;}\n   \
 \n    if ( !-e $tmp || -s $tmp < 10){$internet=0;\
}\n    else {$internet=1;}\n    if (-e $tmp){unlin\
k $tmp;}\n\n    return $internet;\n  }\nsub check_\
pg_is_installed\n  {\n    my @ml=@_;\n    my $r=&p\
g_is_installed (@ml);\n    if (!$r && $p=~/::/)\n \
     {\n	print STDERR \"\\nYou Must Install the pe\
rl package $p on your system.\\nRUN:\\n\\tsudo per\
l -MCPAN -e 'install $pg'\\n\";\n      }\n    elsi\
f (!$r)\n      {\n	myexit(flush_error(\"\\nProgram\
 $p Supported but Not Installed on your system\"))\
;\n      }\n    else\n      {\n	return 1;\n      }\
\n  }\n\n\nsub remote_is_pdb_name_deprecated\n{\n \
   my $in=@_[0];\n    my ($ref_file, $pdb);\n    m\
y ($value,$value1,$value2);\n    my $max=2;\n    \\
n    \n    \n    $ref_file=\"$cache/pdb_entry_type\
.txt\";\n    \n    if ( $in=~/[^\\w\\d\\:\\_]/){re\
turn 0;}\n    elsif ($no_remote_pdb_dir==1)\n     \
 {\n	my $pdbdir=$ENV{'PDB_DIR'};\n	\n	my $r1=\"$pd\
bdir/derived_data/pdb_entry_type.txt\";\n	my $r2=$\
ref_file;\n	if    (-e $r1){$ref_file=$r1;}\n	elsif\
 (-e $r2){$ref_file=$r2;}\n	else\n	  {\n	    my $p\
=substr ($in,0, 4);\n	    add_warning ($$, $$, \"C\
annot find pdb_entry_type.txt;  $p is assumed to b\
e valid; add ftp://ftp.wwpdb.org/pub/pdb/derived_d\
ata/pdb_entry_type.txt in $cache to check name sta\
tus\");\n	  }\n      }\n    elsif ( !-e $ref_file \
|| (-M $ref_file)>$max || -z $ref_file)\n      {\n\
	&url2file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_d\
ata/pdb_entry_type.txt\", $ref_file);\n      }\n  \
  $pdb=substr ($in,0, 4);\n    chomp(($value1=`gre\
p -c $pdb $ref_file`));\n    $pdb=lc($pdb);\n    c\
homp(($value2=`grep -c $pdb $ref_file`));\n    $va\
lue=($value1 || $value2)?1:0;\n    $value=($value>\
0)?1:0;\n    \n    return $value;\n  }\n","use Cwd\
;\nuse Env;\nuse File::Path;\nuse FileHandle;\nuse\
 strict;\n\n\nour (%MODE, %PG, %ENV_SET, %SUPPORTE\
D_OS);\n\nour $VERSION=\"0\"; #_#UPDATE_VERSION\n\\
nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\nour $\
INTERNET=0;\n\nour $CP=\"cp \"; #was causing a cra\
sh on MacOSX\nour $SILENT=\">/dev/null 2>/dev/null\
\";\nour $WEB_BASE=\"http://www.tcoffee.org\";\nou\
r $TCLINKDB_ADDRESS=\"$WEB_BASE/Resources/tclinkdb\
.txt\";\nour $OS=get_os();\nour $ROOT=&get_root();\
\nour $CD=cwd();\nour $CDIR=$CD;\nour $HOME=$ENV{'\
HOME'};\n\nour $OSNAME=$ENV{'OSNAME'};\nour $OSARC\
H=$ENV{'OSARCH'};\nour $REPO_ROOT=\"\";\n\nour $TC\
DIR;\nour $TCCACHE;\nour $TCTMP;\nour $TCM;\nour $\
TCMETHODS;\nour $TCPLUGINS;\nour $PLUGINS_DIR=\"\"\
;\nour $INSTALL_DIR=\"\";\nour $email;\nour $recom\
pile;\n\nour $CXX=\"g++\";\nour $CXXFLAGS=\"\";\n\\
nour $CPP=\"g++\";\nour $CPPFLAGS=\"\";\n\nour $CC\
=\"gcc\";\nour $CFLAGS=$ENV{'CFLAGS'};\n\nour $FC=\
\"f77\";\nour $FFLAGS=\"\";\n\nmy $install=\"all\"\
;\nmy $default_update_action=\"no_update\";\nmy @r\
equired_applications=(\"wget_OR_curl\");\nmy @smod\
e=(\"all\", \"clean\", \"install\");\n\n&initializ\
e_PG();\nmy $cl=join( \" \", @ARGV);\nif ($#ARGV==\
-1 || ($cl=~/-h/) ||($cl=~/-H/) )\n  {\n     print\
 \"\\n!!!!!!! ./install  t_coffee             --> \
installs t_coffee only\";\n     print \"\\n!!!!!!!\
 ./install  all                  --> installs all \
the modes [mcoffee, expresso, psicoffee,rcoffee..]\
\";\n     print \"\\n!!!!!!! ./install  [mcoffee|r\
coffee|..] --> installs the specified mode\";\n   \
  print \"\\n!!!!!!! ./install  -h                \
   --> print usage\\n\\n\";\n     if ( $#ARGV==-1)\
{exit ($EXIT_FAILURE);}\n   }\n     \nif (($cl=~/-\
h/) ||($cl=~/-H/) )\n  {\n    my $m;\n    print \"\
\\n\\n!!!!!!! advanced mode\\n\";\n    foreach $m \
((keys (%MODE)),@smode)\n      {\n	print \"!!!!!!!\
       ./install $m\\n\";\n      }\n    \n    prin\
t \"!!!!!!! ./install [target:package|mode|] [-upd\
ate|-force|-exec=dir|-dis=dir|-root|-tclinkdb=file\
|-] [CC=|FCC=|CXX=|CFLAGS=|CXXFLAGS=]\\n\";\n    p\
rint \"!!!!!!! ./install clean    [removes all exe\
cutables]\\n\";\n    print \"!!!!!!! ./install [op\
tional:target] -update               [updates pack\
age already installed]\\n\";\n    print \"!!!!!!! \
./install [optional:target] -recompile            \
[forces the recompilation of T-Coffee]\\n\";\n\n  \
  print \"!!!!!!! ./install [optional:target] -for\
ce                [Forces recompilation over every\
thing]\\n\";\n    \n    print \"!!!!!!! ./install \
[optional:target] -root                 [You are r\
unning as root]\\n\";\n    print \"!!!!!!! ./insta\
ll [optional:target] -exec=/foo/bar/       [addres\
s for the T-Coffee executable]\\n\";\n    print \"\
!!!!!!! ./install [optional:target] -dis=/foo/bar/\
        [Address where distributions should be sto\
red]\\n\";\n    print \"!!!!!!! ./install [optiona\
l:target] -tclinkdb=foo|update  [file containing a\
ll the packages to be installed]\\n\";\n    print \
\"!!!!!!! ./install [optional:target] -clean      \
          [clean everything]\\n\";\n    print \"!!\
!!!!! ./install [optional:target] -plugins        \
      [plugins directory]\\n\";\n    print \"!!!!!\
!! ./install [optional:target] -tcdir=/foor/bar   \
   [base path where T-Coffee will be installed - d\
efault ~/.t_coffee]\\n\";\n    print \"!!!!!!! ./i\
nstall [optional:target] -repo=/path/to/repo   [bi\
naries repository root directory]\\n\";\n    print\
 \"!!!!!!! ./install [optional:target] -email=<you\
r email>   [needed for remote BLAST]\\n\";\n    pr\
int \"!!!!!!! ./install [optional:target] -proxy=<\
proxy>   [may be needed to access remote services]\
\\n\";\n    \n    print \"!!!!!!! mode:\";\n    fo\
reach $m (keys(%MODE)){print \"$m \";}\n    print \
\"\\n\";\n    print \"!!!!!!! Packages:\";\n    fo\
reach $m (keys (%PG)){print \"$m \";}\n    print \\
"\\n\";\n    \n    print \"\\n\\n\";\n    exit ($E\
XIT_FAILURE);\n  }\n\n\n\nmy (@argl)=($cl=~/(\\S+=\
[^=]+)\\s\\w+=/g);\npush (@argl, ($cl=~/(\\S+=[^=]\
+\\S)\\s*$/g));\n\nforeach $a (@argl)\n  {\n    if\
 ( ($cl=~/CXX=(.*)/)){$CXX=$1;}\n    if ( ($cl=~/-\
CC=(.*)/    )){$CC=$1;}\n    if ( ($cl=~/-FC=(.*)/\
    )){$FC=$1;}\n    if ( ($cl=~/-CFLAGS=(.*)/)){$\
CFLAGS=$1;}\n    if ( ($cl=~/-CXXFLAGS=(.*)/)){$CX\
XFLAGS=$1;}\n  }\nour ($ROOT_INSTALL, $NO_QUESTION\
, $default_update_action,$BINARIES_ONLY,$force, $d\
efault_update_action, $INSTALL_DIR, $PLUGINS_DIR, \
$DISTRIBUTIONS,$tclinkdb, $proxy, $clean);\nif ( (\
$cl=~/-root/)){$ROOT_INSTALL=1;}\nif ( ($cl=~/-no_\
question/)){$NO_QUESTION=1;}\nif ( ($cl=~/-update/\
)){$default_update_action=\"update\";}\nif ( ($cl=\
~/-recompile/)){$recompile=1;}\n\n\n$BINARIES_ONLY\
=1;\n\nif ( ($cl=~/-nobinaries/)){$BINARIES_ONLY=0\
;}\nif ( ($cl=~/-force/)){$force=1;$default_update\
_action=\"update\"}\nif ( ($cl=~/-exec=\\s*(\\S+)/\
)){$INSTALL_DIR=$1;}\nif ( ($cl=~/-plugins=\\s*(\\\
S+)/)){$PLUGINS_DIR=$1;}\nif ( ($cl=~/-dis=\\s*(\\\
S+)/)){$DISTRIBUTIONS=$1;}\n\nif ( ($cl=~/-tclinkd\
b=\\s*(\\S+)/)){$tclinkdb=$1;}\nif ( ($cl=~/-proxy\
=\\s*(\\S+)/)){$proxy=$1;}\nif ( ($cl=~/-clean/)){\
$clean=1;}\nif ( ($cl=~/-repo=\\s*(\\S+)/)){ $REPO\
_ROOT=$1; }\nif ( ($cl=~/-tcdir=\\s*(\\S+)/)){ $TC\
DIR=$1; }\n\nif ( ($cl=~/-email=\\s*(\\S+)/)){$ema\
il=$1;}\n\n\nif ($tclinkdb){&update_tclinkdb ($tcl\
inkdb);}\n\n\nif( $REPO_ROOT ne \"\" ) {\n	if( $OS\
NAME eq \"\" ) { print \"You have specified the re\
pository folder but the required \\\"OSNAME\\\" en\
viroment variable is missing. \\n\"; exit 1; } \n	\
if( $OSARCH eq \"\" ) { print \"You have specified\
 the repository folder but the required \\\"OSARCH\
\\\" enviroment variable is missing. \\n\"; exit 1\
; } \n}\n\n\nif(!$TCDIR) { $TCDIR=\"$HOME/.t_coffe\
e\"; }\n&add_dir ($TCDIR);\n&add_dir ($TCCACHE=\"$\
TCDIR/cache\");\n&add_dir ($TCTMP=\"$CDIR/tmp\");\\
n&add_dir ($TCM=\"$TCDIR/mcoffee\");\n&add_dir ($T\
CMETHODS=\"$TCDIR/methods\");\n&add_dir ($TCPLUGIN\
S=\"$TCDIR/plugins/$OS\");\n\n\nour $BASE=\"$CD/bi\
n\";\nour $BIN=\"$BASE/cache/binaries/$OS\";\nour \
$DOWNLOAD_DIR=\"$BASE/cache/download\";\nour $DOWN\
LOAD_FILE=\"$DOWNLOAD_DIR/files\";\nour $TMP=\"$BA\
SE/cache/tmp\";\n\n&add_dir($BASE);\n&add_dir($BIN\
);\n&add_dir($DOWNLOAD_DIR);\n&add_dir($DOWNLOAD_F\
ILE);\nif (!$DISTRIBUTIONS){$DISTRIBUTIONS=\"$DOWN\
LOAD_DIR/distributions\";}\n&add_dir ($DISTRIBUTIO\
NS);\n&add_dir ($TMP);\n\n\nif    (!$PLUGINS_DIR &\
& !$ROOT_INSTALL){$PLUGINS_DIR=$TCPLUGINS;}\nelsif\
 (!$PLUGINS_DIR &&  $ROOT_INSTALL){$PLUGINS_DIR=\"\
/usr/local/bin/\";}\n\nif    (!$INSTALL_DIR && !$R\
OOT_INSTALL){$INSTALL_DIR=\"$TCDIR/bin/$OS\";mkpat\
h ($INSTALL_DIR);}\nelsif (!$INSTALL_DIR &&  $ROOT\
_INSTALL){$INSTALL_DIR=\"/usr/local/bin/\";}\n\nif\
 (-d \"mcoffee\"){`cp mcoffee/* $TCM`;}\n\n\nour $\
ENV_FILE=\"$TCDIR/.t_coffee_env\";\nunlink ($ENV_F\
ILE);\n&add2env_file ($ENV_FILE,\"EMAIL_4_TCOFFEE\\
", $email);\n&add2env_file ($ENV_FILE,\"http_proxy\
_4_TCOFFEE\", $proxy);\n&env_file2putenv ($ENV_FIL\
E);\n&set_proxy($proxy);\n\n\n\nmy ($target, $p, $\
r);\n$target=$p;\n\nforeach $p (  ((keys (%PG)),(k\
eys(%MODE)),(@smode)) )\n  {\n    if ($ARGV[0] eq \
$p && $target eq \"\"){$target=$p;}\n  }\nif ($tar\
get eq \"\"){exit ($EXIT_FAILURE);}\n\n\nforeach $\
r (@required_applications)\n  {\n    my @app_list;\
\n    my $i;\n    $i=0;\n    \n    @app_list=split\
 (/_OR_/, $r);\n    foreach my $pg (@app_list)\n  \
    {\n	$i+=&pg_is_installed ($pg);\n      }\n    \
if ($i==0)\n      {\n      print \"One of the foll\
owing packages must be installed to proceed: \";\n\
      foreach my $pg (@app_list)\n	{\n	  print (\"\
$pg \");\n	}\n      die;\n    }\n  }\n\n\n\n\n\n\n\
&sign_license_ni();\n\n\n$PG{C}{compiler}=get_C_co\
mpiler($CC);\n$PG{Fortran}{compiler}=get_F_compile\
r($FC);\n$PG{CXX}{compiler}=$PG{CPP}{compiler}=$PG\
{GPP}{compiler}=get_CXX_compiler($CXX);\nif ($CXXF\
LAGS){$PG{CPP}{options}=$PG{GPP}{options}=$PG{CXX}\
{options}=$CXXFLAGS;}\nif ($CFLAGS ne \"\" ){$PG{C\
}{options}=$CFLAGS;}\nforeach my $c (keys(%PG))\n \
 {\n    my $arguments;\n    if ($PG{$c}{compiler})\
\n      {\n	$arguments=\"$PG{$c}{compiler_flag}=$P\
G{$c}{compiler} \";\n	if ($PG{$c}{options})\n	  {\\
n	    $arguments.=\"$PG{$c}{options_flag}='\" . $P\
G{$c}{options} . \"' \";\n	  }\n	$PG{$c}{arguments\
}=$arguments;\n      }\n  }\n\nif ($PG{$target}){$\
PG{$target}{install}=1;}\nelse\n  {\n    foreach m\
y $pg (keys(%PG))\n      {\n	if ( $target eq \"all\
\" || ($PG{$pg}{mode}=~/$target/))\n	  {\n	    $PG\
{$pg} {install}=1;\n	  }\n      }\n  }\n\nforeach \
my $pg (keys(%PG))\n  {\n    if (!$PG{$pg}{update_\
action}){$PG{$pg}{update_action}=$default_update_a\
ction;}\n    elsif ($PG{$pg}{update_action} eq \"n\
ever\"){$PG{$pg}{install}=0;}\n    if ( $force && \
$PG{$pg}{install})\n      {\n	`rm $BIN/$pg $BIN/$p\
g.exe $SILENT`;\n      }\n    if ($PG{$pg}{update_\
action} eq \"update\" && $PG{$pg}{install}){$PG{$p\
g}{update}=1;}\n  }\n\nif (($target=~/clean/))\n  \
{\n    print \"------- cleaning executables -----\\
\n\";\n    `rm bin/* $SILENT`;\n    exit ($EXIT_SU\
CCESS);\n  }\n\nif ( !$PG{$target}){print \"------\
- Installing T-Coffee Modes\\n\";}\n\nforeach my $\
m (keys(%MODE))\n  {\n    if ( $target eq \"all\" \
|| $target eq $m)\n      {\n	print \"\\n------- Th\
e installer will now install the $m components $MO\
DE{$m}{description}\\n\";\n	foreach my $pg (keys(%\
PG))\n	  {\n	    if ( $PG{$pg}{mode} =~/$m/ && $PG\
{$pg}{install})\n	      {\n		if ($PG{$pg}{touched}\
){print \"------- $PG{$pg}{dname}: already process\
ed\\n\";}\n		else {$PG{$pg}{success}=&install_pg($\
pg);$PG{$pg}{touched}=1;}\n	      }\n	  }\n      }\
\n  }\n\nif ( $PG{$target}){print \"------- Instal\
ling Individual Package\\n\";}\nforeach my $pg (ke\
ys (%PG))\n  {\n    \n    if ( $PG{$pg}{install} &\
& !$PG{$pg}{touched})\n      {\n	print \"\\n------\
- Install $pg\\n\";\n	$PG{$pg}{success}=&install_p\
g($pg);$PG{$pg}{touched}=1;\n      }\n  }\nprint \\
"------- Finishing The installation\\n\";\nmy $fin\
al_report=&install ($INSTALL_DIR);\n\nprint \"\\n\\
";\nprint \"**************************************\
*******************************\\n\";\nprint \"***\
*****              INSTALLATION SUMMARY          *\
****************\\n\";\nprint \"******************\
**************************************************\
*\\n\";\nprint \"------- SUMMARY package Installat\
ion:\\n\";\nprint \"-------   Executable Installed\
 in: $PLUGINS_DIR\\n\";\n\nforeach my $pg (keys(%P\
G))\n  {\n    if ( $PG{$pg}{install})\n      {\n	m\
y $bin_status=($PG{$pg}{from_binary} && $PG{$pg}{s\
uccess})?\"[from binary]\":\"\";\n	if     ( $PG{$p\
g}{new} && !$PG{$pg}{old})                     {pr\
int \"*------        $PG{$pg}{dname}: installed $b\
in_status\\n\"; $PG{$pg}{status}=1;}\n	elsif  ( $P\
G{$pg}{new} &&  $PG{$pg}{old})                    \
 {print \"*------        $PG{$pg}{dname}: updated \
$bin_status\\n\"  ; $PG{$pg}{status}=1;} \n	elsif \
 (!$PG{$pg}{new} &&  $PG{$pg}{old} && !$PG{$pg}{up\
date}){print \"*------        $PG{$pg}{dname}: pre\
vious\\n\" ; $PG{$pg}{status}=1;}\n	elsif  (!$PG{$\
pg}{new} &&  $PG{$pg}{old} &&  $PG{$pg}{update}){p\
rint \"*------        $PG{$pg}{dname}: failed upda\
te (previous installation available)\\n\";$PG{$pg}\
{status}=0;}\n	else                               \
                           {print \"*------       \
 $PG{$pg}{dname}: failed installation\\n\";$PG{$pg\
}{status}=0;}\n      }\n  }\nmy $failure;\n\nif ( \
!$PG{$target}){print \"*------ SUMMARY mode Instal\
lation:\\n\";}\nforeach my $m (keys(%MODE))\n  {\n\
  \n    if ( $target eq \"all\" || $target eq $m)\\
n      {\n	my $succesful=1;\n	foreach my $pg (keys\
(%PG))\n	  {\n	    if (($PG{$pg}{mode}=~/$m/) && $\
PG{$pg}{install} && $PG{$pg}{status}==0)\n	      {\
\n		$succesful=0;\n		print \"*!!!!!!       $PG{$pg\
}{dname}: Missing\\n\";\n	      }\n	  }\n	if ( $su\
ccesful)\n	  {\n	    $MODE{$m}{status}=1;\n	    pr\
int \"*------       MODE $MODE{$m}{dname} SUCCESSF\
ULLY installed\\n\";\n	  }\n	else\n	  {\n	    $fai\
lure++;\n	    $MODE{$m}{status}=0;\n	    print \"*\
!!!!!!       MODE $MODE{$m}{dname} UNSUCCESSFULLY \
installed\\n\";\n	  }\n      }\n  }\n\n    \n     \
 \nif ($clean==1 && ($BASE=~/install4tcoffee/) ){p\
rint \"*------ Clean Installation Directory: $BASE\
\\n\";`rm -rf $BASE`;}\nforeach my $pg (keys(%PG))\
{if ($PG{$pg}{install} && $PG{$pg}{status}==0){exi\
t ($EXIT_FAILURE);}}\n\nif ($failure)\n  {\n    pr\
int \"********************************************\
*************************\\n\";\n    print \"*****\
***     SOME PACKAGES FAILED TO INSTALL        ***\
**************\\n\";\n    print \"****************\
**************************************************\
***\\n\";\n    print \"\\nSome of the reported fai\
lures may be due to connectivity problems\";\n    \
print \"\\nRerun the installation and the installe\
r will specifically try to install the missing pac\
kages\";\n    print \"\\nIf this Fails, go to the \
original website and install the package manually\\
";\n  }\n\nprint \"*******************************\
**************************************\\n\";\nprin\
t \"********              FINALIZE YOUR INSTALLATI\
ON    *****************\\n\";\nprint \"***********\
**************************************************\
********\\n\";\nprint \"------- Your third party e\
xecutables are in:\\n\"; \nprint \"-------       $\
PLUGINS_DIR:\\n\";\nprint \"------- Your t_coffee \
exccutable is in\\n\";\nprint \"-------       $INS\
TALL_DIR:\\n\";\nprint \"------- In order to make \
your installation permanent add these two lines\\n\
\";\nprint \"export PATH=$INSTALL_DIR:\\$PATH\\n\"\
;\nprint \"export PLUGINS_4_TCOFFEE=$PLUGINS_DIR:\\
\n\";\nif ($OS eq \"linux\")\n  {\n    print \"---\
----       to the file: $HOME/.bashrc\\n\";\n  }\n\
else \n  {\n    print \"-------       to the file:\
 $HOME/.profile\\ OR $HOME/.basrc\";\n  }\nexit ($\
EXIT_SUCCESS);  \n  \nsub get_CXX_compiler\n  {\n \
   my $c=@_[0];\n    my (@clist)=(\"g++\");\n    \\
n    return get_compil ($c, @clist);\n }\nsub get_\
C_compiler\n  {\n    my $c=@_[0];\n    my (@clist)\
=(\"gcc\", \"cc\", \"icc\");\n    \n    return get\
_compil ($c, @clist);\n }\n\nsub get_F_compiler\n \
 {\n    my ($c)=@_[0];\n    my @clist=(\"f77\", \"\
g77\",\"g95\", \"gfortran\", \"ifort\");\n    retu\
rn get_compil ($c, @clist);\n  } \n       \nsub ge\
t_compil\n  {\n    my ($fav,@clist)=(@_);\n    \n \
   #return the first compiler found installed in t\
he system. Check first the favorite\n    foreach m\
y $c ($fav,@clist)\n      {\n	if  (&pg_is_installe\
d ($c)){return $c;}\n      }\n    return \"\";\n  \
}\nsub exit_if_pg_not_installed\n  {\n    my (@arg\
)=(@_);\n    \n    foreach my $p (@arg)\n      {\n\
	if ( !&pg_is_installed ($p))\n	  {\n	    print \"\
!!!!!!!! The $p utility must be installed for this\
 installation to proceed [FATAL]\\n\";\n	    die;\\
n	  }\n      }\n    return 1;\n  }\nsub set_proxy\\
n  {\n    my ($proxy)=(@_);\n    my (@list,$p);\n \
   \n    @list= (\"HTTP_proxy\", \"http_proxy\", \\
"HTTP_PROXY\", \"ALL_proxy\", \"all_proxy\",\"HTTP\
_proxy_4_TCOFFEE\",\"http_proxy_4_TCOFFEE\");\n   \
 \n    if (!$proxy)\n      {\n	foreach my $p (@lis\
t)\n	  {\n	    if ( ($ENV_SET{$p}) || $ENV{$p}){$p\
roxy=$ENV{$p};}\n	  }\n      }\n    foreach my $p(\
@list){$ENV{$p}=$proxy;}\n  }\n	\nsub check_intern\
et_connection\n  {\n    my $internet;\n    \n    i\
f ( -e \"x\"){unlink (\"x\");}\n    if     (&pg_is\
_installed    (\"wget\")){`wget www.google.com -Ox\
 >/dev/null 2>/dev/null`;}\n    elsif  (&pg_is_ins\
talled    (\"curl\")){`curl www.google.com -ox >/d\
ev/null 2>/dev/null`;}\n    else\n      {\n	printf\
 stderr \"\\nERROR: No pg for remote file fetching\
 [wget or curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE\
);\n      }\n    \n    if ( !-e \"x\" || -s \"x\" \
< 10){$internet=0;}\n    else {$internet=1;}\n    \
if (-e \"x\"){unlink \"x\";}\n    return $internet\
;\n  }\nsub url2file\n  {\n    my ($cmd, $file,$wg\
et_arg, $curl_arg)=(@_);\n    my ($exit,$flag, $pg\
, $arg);\n\n    \n    \n    if ($INTERNET || check\
_internet_connection ()){$INTERNET=1;}\n    else\n\
      {\n	print STDERR \"ERROR: No Internet Connec\
tion [FATAL:install.pl]\\n\";\n	exit ($EXIT_FAILUR\
E);\n      }\n\n    $cmd=~s/http\\:/https\\:/g;\n \
   \n    if     (&pg_is_installed    (\"wget\")){$\
pg=\"wget\"; $flag=\"-O\";$arg=\"--tries=2 --conne\
ct-timeout=10 --no-check-certificate $wget_arg\";}\
\n    elsif  (&pg_is_installed    (\"curl\")){$pg=\
\"curl\"; $flag=\"-f -o\";$arg=$curl_arg;}\n    el\
se\n      {\n	printf stderr \"\\nERROR: No pg for \
remote file fetching [wget or curl][FATAL]\\n\";\n\
	exit ($EXIT_FAILURE);\n      }\n    \n    \n    i\
f (-e $file){unlink($file);}\n    $exit=system \"$\
pg $cmd $flag$file $arg\";\n    return $exit;\n  }\
\n\nsub pg_is_installed\n  {\n    my ($p, $dir)=(@\
_);\n    my ($r,$m, $ret);\n    my ($supported, $l\
anguage, $compil);\n    \n  \n    if ( $PG{$p})\n \
     {\n	$language=$PG{$p}{language2};\n	$compil=$\
PG{$language}{compiler};\n      }\n    \n    if ( \
$compil eq \"CPAN\")\n      {\n	if ( system (\"per\
l -M$p -e 1\")==$EXIT_SUCCESS){$ret=1;}\n	else {$r\
et=0;}\n      }\n    elsif ($dir)\n      {\n	if (-\
e \"$dir/$p\" || -e \"$dir/$p\\.exe\"){$ret=1;}\n	\
else {$ret=0;}\n      }\n    elsif (-e \"$PLUGINS_\
DIR/$p\" || -e \"$PLUGINS_DIR/$p.exe\"){$ret=1;}\n\
    else\n      {\n	$r=`which $p 2>/dev/null`;\n	i\
f ($r eq \"\"){$ret=0;}\n	else {$ret=1;}\n      }\\
n   \n    return $ret;\n  }\nsub install\n  {\n   \
 my ($new_bin)=(@_);\n    my ($copied, $report);\n\
\n    \n    if (!$ROOT_INSTALL)\n      {\n	`$CP $B\
IN/* $PLUGINS_DIR`;\n	if (-e \"$BIN/t_coffee\")\n	\
  {\n	    `$CP $BIN/t_coffee $INSTALL_DIR`;\n	    \
  unlink(\"$PLUGINS_DIR/t_coffee\");\n	  }\n	$copi\
ed=1;\n      }\n    else\n      {\n	$copied=&root_\
run (\"You must be root to finalize the installati\
on\", \"$CP $BIN/* $PLUGINS_DIR $SILENT\");\n	if (\
-e \"$BIN/t_coffee\")\n	  {\n	    &root_run (\"You\
 must be root to finalize the installation\", \"$C\
P $BIN/t_coffee $INSTALL_DIR\");\n	    &root_run (\
\"You must be root to finalize the installation\",\
 \"rm  $PLUGINS_DIR/t_coffee\");\n	  }\n      }\n \
   \n     \n  if ( !$copied)\n    {\n      $report\
=\"*!!!!!! Installation unsuccesful. The executabl\
es have been left in $BASE/bin\\n\";\n    }\n  els\
if ( $copied && $ROOT)\n    {\n      $report=\"*--\
---- Installation succesful. Your executables have\
 been copied in $new_bin and are on your PATH\\n\"\
;\n    }\n  elsif ( $copied && !$ROOT)\n    {\n   \
   $report= \"*!!!!!! T-Coffee has been installed \
in $INSTALL_DIR\\n\";\n      $report= \"*!!!!!! T-\
Coffee and associated packages have been copied in\
: $PLUGINS_DIR\\n\";\n      $report.=\"*!!!!!! Thi\
s T-Coffee location is NOT on your PATH sytem vari\
able\\n\";\n      if ( $OS eq \"linux\")\n	{\n	  $\
report.=\"*!!!!!! You can do so by adding the foll\
owing line in your ~/.bashrc file:\\n\";\n	}\n    \
  else\n	{\n	  $report.=\"*!!!!!! You can do so by\
 adding the following line in your ~/.profile file\
:\\n\";\n	}\n      $report.=\"*!!!!!! export PATH=\
$INSTALL_DIR:\\$PATH\\n\";\n    }\n  return $repor\
t;\n}\n\nsub sign_license_ni\n  {\n    my $F=new F\
ileHandle;\n    open ($F, \"license.txt\");\n    w\
hile (<$F>)\n      {\n	print \"$_\";\n      }\n   \
 close ($F);\n    \n    return;\n  }\n\nsub instal\
l_pg\n  {\n    my ($pg)=(@_);\n    my ($report, $p\
revious, $language, $compiler, $return);\n    \n  \
  if (!$PG{$pg}{install}){return 1;}\n    \n    $p\
revious=&pg_is_installed ($pg);\n    \n    if ($PG\
{$pg}{update_action} eq \"no_update\" && $previous\
)\n      {\n	$PG{$pg}{old}=1;\n	$PG{$pg}{new}=0;\n\
	$return=1;\n      }\n    else\n      {\n	$PG{$pg}\
{old}=$previous;\n	\n	if ($PG{$pg} {language2} eq \
\"Perl\"){&install_perl_package ($pg);}\n	elsif ($\
BINARIES_ONLY && &install_binary_package ($pg)){$P\
G{$pg}{from_binary}=1;}\n	elsif (&install_source_p\
ackage ($pg)){;}\n	else \n	  {\n	    \n	    if (!&\
supported_os($OS))\n	      {\n		print \"!!!!!!!! $\
pg compilation failed, binary unsupported for $OS\\
\n\"; \n	      }\n	    elsif (!($PG{$pg}{from_bina\
ry}=&install_binary_package ($pg)))\n	      {\n		p\
rint \"!!!!!!!! $pg compilation and  binary instal\
lation failed\\n\";\n	      }\n	  }\n	$PG{$pg}{new\
}=$return=&pg_is_installed ($pg,$BIN);\n      }\n\\
n    \n    return $return;\n  }\nsub install_perl_\
package\n  {\n    my ($pg)=(@_);\n    my ($report,\
 $language, $compiler);\n    \n    $language=$PG{$\
pg} {language2};\n    $compiler=$PG{$language}{com\
piler};\n    \n    if (!&pg_is_installed ($pg))\n \
     {\n	if ( $OS eq \"windows\"){`perl -M$compile\
r -e 'install $pg'`;}\n	elsif ( $ROOT eq \"sudo\")\
{system (\"sudo perl -M$compiler -e 'install $pg'\\
");}\n	else {system (\"su root -c perl -M$compiler\
 -e 'install $pg'\");}\n      }\n    return &pg_is\
_installed ($pg);\n  }\n\n\n\nsub install_source_p\
ackage\n  {\n    my ($pg)=(@_);\n    my ($report, \
$download, $arguments, $language, $address, $name,\
 $ext, $main_dir, $distrib);\n    my $wget_tmp=\"$\
TMP/wget.tmp\";\n    my (@fl);\n    if ( $default_\
update_action ne \"update\" && (-e \"$BIN/$pg\" ||\
 -e \"$BIN/$pg.exe\" )  ){return 1;}\n    \n    #\\
n    # check if the module exists in the repositor\
y cache \n    #\n	if( repo_load($pg) ) {\n		return\
 1;\n	}\n    \n    if ($pg eq \"t_coffee\")  {retu\
rn   &install_t_coffee_source ($pg);}\n    elsif (\
$pg eq \"TMalign\"){return   &install_TMalign ($pg\
);}\n    \n    chdir $DISTRIBUTIONS;\n    \n    $d\
ownload=$PG{$pg}{source};\n    \n    if (($downloa\
d =~/tgz/))\n      {\n	($address,$name,$ext)=($dow\
nload=~/(.+\\/)([^\\/]+)(\\.tgz).*/);\n      }\n  \
  elsif (($download=~/tar\\.gz/))\n      {\n	($add\
ress,$name,$ext)=($download=~/(.+\\/)([^\\/]+)(\\.\
tar\\.gz).*/);\n      }\n    elsif (($download=~/t\
ar/))\n      {\n	($address,$name,$ext)=($download=\
~/(.+\\/)([^\\/]+)(\\.tar).*/);\n      }\n    else\
\n      {\n	($address,$name)=($download=~/(.+\\/)(\
[^\\/]+)/);\n	$ext=\"\";\n      }\n    $distrib=\"\
$name$ext\";\n    \n    if ( !-d $pg){mkdir $pg;}\\
n    chdir $pg;\n   \n    #get the distribution if\
 available\n    if ( -e \"$DOWNLOAD_DIR/$distrib\"\
)\n      {\n	`$CP $DOWNLOAD_DIR/$distrib .`;\n    \
  }\n    #UNTAR and Prepare everything\n    if (!-\
e \"$name.tar\" && !-e \"$name\")\n      {\n	&chec\
k_rm ($wget_tmp);\n	print \"\\n------- Downloading\
/Installing $pg\\n\";\n	\n	if (!-e $distrib && &ur\
l2file (\"$download\", \"$wget_tmp\")==$EXIT_SUCCE\
SS)\n	  {\n	    \n	    `mv $wget_tmp $distrib`;\n	\
    `$CP $distrib $DOWNLOAD_DIR/`;\n	  }\n\n	if (!\
-e $distrib)\n	  {\n	    print \"!!!!!!! Download \
of $pg distribution failed\\n\";\n	    print \"!!!\
!!!! Check Address: $PG{$pg}{source}\\n\";\n	    r\
eturn 0;\n	  }\n	print \"\\n------- unzipping/unta\
ring $name\\n\";\n	if (($ext =~/z/))\n	  { \n	    \
&flush_command (\"gunzip -f $name$ext\");\n	    \n\
	  }\n	if (($ext =~/tar/) || ($ext =~/tgz/))\n	  {\
\n	    &flush_command(\"tar -xvf $name.tar\");\n	 \
 }\n      }\n    #Guess and enter the distribution\
 directory\n    @fl=ls($p);\n    foreach my $f (@f\
l)\n      {\n	if (-d $f)\n	  {\n	    $main_dir=$f;\
\n	  }\n      }\n    if (-d $main_dir)\n	  \n     \
 {\n	chdir $main_dir;}\n    else\n      {\n	print \
\"Error: $main_dir does not exist\";\n      }\n   \
 print \"\\n------- Compiling/Installing $pg\\n\";\
\n    `make clean $SILENT`;\n    \n    \n    #\n  \
  # SAP module\n    #\n    if ($pg eq \"sap\")\n  \
    {\n	if (-e \"./configure\")\n	  {\n	    #new s\
ap distribution\n	    \n	    &flush_command (\"./c\
onfigure\");\n	    &flush_command (\"make clean\")\
;\n	    &flush_command (\"make\");\n	    &check_cp\
 (\"./src/$pg\", \"$BIN\");\n	    repo_store(\"./s\
rc/$pg\");\n	  }\n	else\n	  {\n	    #old style dis\
tribution\n	    `rm *.o sap  sap.exe ./util/aa/*.o\
  ./util/wt/.o $SILENT`;\n	    &flush_command (\"m\
ake $arguments sap\");\n	    &check_cp ($pg, \"$BI\
N\");\n	    repo_store($pg);\n	  }\n      }\n    \\
n    #\n    # CLUSTALW2 module\n    #\n    elsif (\
$pg eq \"clustalw2\")\n      {\n	&flush_command(\"\
./configure\");\n	&flush_command(\"make $arguments\
\");\n	&check_cp (\"./src/$pg\", \"$BIN\");\n	repo\
_store(\"./src/$pg\");\n      }\n\n    #\n    # CL\
USTAL-OMEGA module\n    #\n    elsif ($pg eq \"clu\
stalo\")\n      {\n	&flush_command(\"./configure\"\
);\n	&flush_command(\"make $arguments\");\n	&check\
_cp (\"./src/$pg\", \"$BIN\");\n	repo_store(\"./sr\
c/$pg\");\n      }\n\n    #\n    # STRIKE module\n\
    #\n    elsif ($pg eq \"strike\")\n      {\n	&f\
lush_command(\"make $arguments\");\n	&check_cp (\"\
./bin/$pg\", \"$BIN\");\n	repo_store(\"./bin/$pg\"\
);\n      }\n    \n    #\n    # FSA module\n    # \
\n    elsif ($pg eq \"fsa\")\n      {\n	&flush_com\
mand(\"./configure --prefix=$BIN\");\n	&flush_comm\
and(\"make $arguments\");\n	&flush_command (\"make\
 install\");\n\n	repo_store(\"fsa\", \"$BIN/bin\")\
;\n	`mv $BIN/bin/* $BIN`;\n	`rmdir $BIN/bin`;\n   \
   }\n    \n    #\n    # CLUSTALW module\n    #\n \
   elsif ($pg eq \"clustalw\")\n      {\n	&flush_c\
ommand(\"make $arguments clustalw\");\n	`$CP $pg $\
BIN $SILENT`;\n	repo_store($pg);\n      }\n    \n \
   #\n    # MAFFT module\n    #\n    elsif ($pg eq\
 \"mafft\")\n      {\n	my $base=cwd();\n	my $c;\n	\
\n	#compile core\n	mkpath (\"./mafft/bin\");\n	mkp\
ath (\"./mafft/lib\");\n	chdir \"$base/core\";\n	`\
make clean $SILENT`;\n	&flush_command (\"make $arg\
uments\");\n	&flush_command (\"make install LIBDIR\
=../mafft/lib BINDIR=../mafft/bin\");\n	\n	#compil\
e extension\n	chdir \"$base/extensions\";\n	`make \
clean $SILENT`;\n	&flush_command (\"make $argument\
s\");\n	&flush_command (\"make install LIBDIR=../m\
afft/lib BINDIR=../mafft/bin\");\n	\n	#put everyth\
ing in mafft and copy the compiled stuff in bin\n	\
chdir \"$base\";\n	if ($ROOT_INSTALL)\n	  {\n	    \
&root_run (\"You Must be Root to Install MAFFT\\n\\
", \"mkdir /usr/local/mafft/;$CP mafft/lib/* /usr/\
local/mafft;$CP mafft/lib/mafft* /usr/local/bin ;$\
CP mafft/bin/mafft /usr/local/bin/; \");\n	  }\n	e\
lse\n	  {\n	    `$CP mafft/lib/*  $BIN`;\n	    `$C\
P mafft/bin/mafft  $BIN`;\n	  }\n	`tar -cvf mafft.\
tar mafft`;\n	`gzip mafft.tar`;\n	`mv mafft.tar.gz\
 $BIN`;\n	\n	repo_store(\"mafft/bin/mafft\", \"maf\
ft/lib/\", \"$BIN/mafft.tar.gz\");\n      }\n     \
 \n    #\n    # DIALIGN-TX module\n    #\n    elsi\
f ( $pg eq \"dialign-tx\" )\n      {\n	my $f;\n	my\
 $base=cwd();\n\n	chdir \"./source\";\n	if ($OS eq\
 \"macosx\"){&flush_command (\"cp makefile.MAC_OS \
makefile\");}\n\n	&flush_command (\" make CPPFLAGS\
='-O3 -funroll-loops' all\");\n	\n	chdir \"..\";\n\
	&check_cp (\"./source/$pg\", \"$BIN\");\n	repo_st\
ore(\"./source/$pg\");\n      }\n      \n    #\n  \
  # DIALIGN-T module \n    # (is the same as diali\
gn-tx, but it is mantained for backward name compa\
tibility with tcoffee)\n    #\n    elsif ( $pg eq \
\"dialign-t\" )\n      {\n	my $f;\n	my $base=cwd()\
;\n\n	chdir \"./source\";\n	if ($OS eq \"macosx\")\
{&flush_command (\"cp makefile.MAC_OS makefile\");\
}\n\n	&flush_command (\" make CPPFLAGS='-O3 -funro\
ll-loops' all\");\n	\n	chdir \"..\";\n	&check_cp (\
\"./source/dialign-tx\", \"$BIN/dialign-t\");\n	re\
po_store(\"$BIN/dialign-t\");	\n      }      \n   \
   \n    #\n    # POA module\n    #\n    elsif ($p\
g eq \"poa\")\n      {\n	&flush_command (\"make $a\
rguments poa\");\n	&check_cp (\"$pg\", \"$BIN\");\\
n	repo_store(\"$pg\");\n      }\n     \n     \n   \
 #\n    # PROBCONS module\n    #\n    elsif ( $pg \
eq \"probcons\")\n      {\n	&add_C_libraries(\"./P\
robabilisticModel.h\", \"list\", \"cstring\");\n	\\
n	`rm *.exe $SILENT`;\n	&flush_command (\"make $ar\
guments probcons\");\n	&check_cp(\"$pg\", \"$BIN/$\
pg\");\n	repo_store(\"$pg\");\n      }\n      \n  \
  #\n    # PROBCONS RNA module\n    #\n    elsif (\
 $pg eq \"probconsRNA\")\n      {\n	&add_C_librari\
es(\"./ProbabilisticModel.h\", \"list\", \"cstring\
\");\n	&add_C_libraries(\"./Main.cc\", \"iomanip\"\
, \"cstring\",\"climits\");\n	`rm *.exe $SILENT`;\\
n	&flush_command (\"make $arguments probcons\");\n\
	&check_cp(\"probcons\", \"$BIN/$pg\");\n	repo_sto\
re(\"$BIN/$pg\");\n      }\n\n	#\n	# MUSCLE module\
\n	#\n    elsif (  $pg eq \"muscle\")\n      {	\n	\
`rm *.o muscle muscle.exe $SILENT`;\n	if ($OS eq \\
"macosx\" || $OS eq \"linux\")\n	  {\n	    &replac\
e_line_in_file (\"./Makefile\", \"LDLIBS = -lm -st\
atic\",  \"LDLIBS = -lm\");\n	  }\n	elsif ($OS eq \
\"windows\")\n	  {\n	    &replace_line_in_file (\"\
./intmath.cpp\",  \"double log2e\",      \"double \
cedric_log\");\n	    &replace_line_in_file (\"./in\
tmath.cpp\",  \"double log2\",       \"double log_\
notuse\");\n	    &replace_line_in_file (\"./intmat\
h.cpp\",  \"double cedric_log\", \"double log2e\")\
;\n	  }\n	&flush_command (\"make $arguments all\")\
;\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_store(\"$\
pg\");	\n      }\n      \n     #\n     # MUS4 modu\
le\n     #\n     elsif (  $pg eq \"mus4\")\n      \
{\n	`rm *.o muscle muscle.exe $SILENT`;\n	&flush_c\
ommand (\"./mk\");\n	&check_cp(\"$pg\", \"$BIN\");\
\n	repo_store(\"$pg\");	\n      }\n      \n    #\n\
    # PCMA module\n    #\n    elsif ( $pg eq \"pcm\
a\")\n      {\n	if ($OS eq \"macosx\")\n	  {\n	   \
 &replace_line_in_file (\"./alcomp2.c\", \"malloc.\
h\",  \"\");\n	  }\n	&flush_command (\"make $argum\
ents pcma\");\n	&check_cp(\"$pg\", \"$BIN\");\n	re\
po_store(\"$pg\");	\n      }\n      \n    #\n    #\
 KALIGN module\n    #\n    elsif ($pg eq \"kalign\\
")\n      {\n	&flush_command (\"./configure\");\n	\
&flush_command(\"make $arguments\");\n	&check_cp (\
\"$pg\",$BIN);\n	repo_store(\"$pg\");	\n      }\n \
     \n    #\n    # AMAP module\n    #\n    elsif \
( $pg eq \"amap\")\n      {\n	&add_C_libraries(\".\
/Amap.cc\", \"iomanip\", \"cstring\",\"climits\");\
	\n	`make clean $SILENT`;\n	&flush_command (\"make\
 $arguments all\");\n	&check_cp (\"$pg\", $BIN);\n\
	repo_store(\"$pg\");	\n      }\n      \n    #\n  \
  # PRODA module\n    #\n    elsif ( $pg eq \"prod\
a\")\n      {\n	`sed -i '' 's/int errno = 0;/int e\
rrno; errno = 0;/' Main.cc`;\n	&add_C_libraries(\"\
AlignedFragment.h\", \"vector\", \"iostream\", \"c\
string\",\"cstdlib\");\n	&add_C_libraries(\"Main.c\
c\", \"vector\", \"climits\");	\n	&add_C_libraries\
(\"Sequence.cc\", \"stdlib.h\", \"cstdio\");	\n	&f\
lush_command (\"make $arguments all\");\n	&check_c\
p (\"$pg\", $BIN);\n	repo_store(\"$pg\");	\n      \
}\n      \n    #\n    # PRANK module\n    #\n    e\
lsif ( $pg eq \"prank\")\n      {\n	&flush_command\
 (\"make $arguments all\");\n	&check_cp (\"$pg\", \
$BIN);\n	repo_store(\"$pg\");	\n      }\n      \n \
   #\n    # !!!! MUSTANG module\n    #\n     elsif\
 ( $pg eq \"mustang\")\n      {\n	&flush_command (\
\"rm ./bin/*\");\n	&flush_command (\"make $argumen\
ts all\");\n\n	if ( $OS=~/windows/){&flush_command\
(\"cp ./bin/* $BIN/mustang.exe\");}\n	else {&flush\
_command(\"cp ./bin/* $BIN/mustang\");}\n	\n	repo_\
store(\"$BIN/mustang\");\n      }\n\n	#\n	# RNAplf\
old module\n	#\n    elsif ( $pg eq \"RNAplfold\")\\
n      {\n	&flush_command(\"./configure\");\n	&flu\
sh_command (\"make $arguments all\");\n	&check_cp(\
\"./Progs/RNAplfold\", \"$BIN\");\n	&check_cp(\"./\
Progs/RNAalifold\", \"$BIN\");\n	&check_cp(\"./Pro\
gs/RNAfold\", \"$BIN\");\n	\n	repo_store(\"./Progs\
/RNAplfold\", \"./Progs/RNAalifold\", \"./Progs/RN\
Afold\");\n      }\n      \n    #\n    # !!! RETRE\
E module\n    #\n    elsif ( $pg eq \"retree\")\n \
     {\n	chdir \"src\";\n	&flush_command (\"cp Mak\
efile.unx Makefile\");\n	&flush_command (\"make $a\
rguments all\");\n	&flush_command (\"make put\");\\
n	system \"cp ../exe/* $BIN\";\n	\n	repo_store(\"r\
etree\", \"../exe\");\n      }\n	\n    chdir $CDIR\
;\n    return &pg_is_installed ($pg, $BIN);\n  }\n\
\nsub install_t_coffee_source\n  {\n    my ($pg)=(\
@_);\n    my ($report,$cflags, $arguments, $langua\
ge, $compiler) ;\n\n    #1-Install T-Coffee\n    c\
hdir \"t_coffee_source\";\n    &flush_command (\"m\
ake clean\");\n    print \"\\n------- Compiling T-\
Coffee\\n\";\n    $language=$PG{$pg} {language2};\\
n    $arguments=$PG{$language}{arguments};\n    \n\
    if ( $CC ne \"\")\n      {\n	print \"make -i $\
arguments t_coffee \\n\";\n	&flush_command (\"make\
 -i $arguments t_coffee\");\n      }\n    &check_c\
p ($pg, $BIN);\n    \n    chdir $CDIR;\n    return\
 &pg_is_installed ($pg, $BIN);\n  }\nsub install_T\
Malign\n  {\n    my ($pg)=(@_);\n    my $report;\n\
    chdir \"t_coffee_source\";\n    print \"\\n---\
---- Compiling TMalign\\n\";\n    `rm TMalign TMal\
ign.exe $SILENT`;\n    if ( $FC ne \"\"){&flush_co\
mmand (\"make -i $PG{Fortran}{arguments} TMalign\"\
);}\n    &check_cp ($pg, $BIN);\n    repo_store($p\
g);\n\n    if ( !-e \"$BIN/$pg\" && pg_has_binary_\
distrib ($pg))\n      {\n	print \"!!!!!!! Compilat\
ion of $pg impossible. Will try to install from bi\
nary\\n\";\n	return &install_binary_package ($pg);\
\n      }\n    chdir $CDIR;\n    return &pg_is_ins\
talled ($pg, $BIN);\n  }\n\nsub pg_has_binary_dist\
rib\n  {\n    my ($pg)=(@_);\n    if ($PG{$pg}{win\
dows}){return 1;}\n    elsif ($PG{$pg}{osx}){retur\
n 1;}\n    elsif ($PG{$pg}{macosx}){return 1;}\n\n\
    elsif ($PG{$pg}{linux}){return 1;}\n    return\
 0;\n  }\nsub install_binary_package\n  {\n    my \
($pg)=(@_);\n    my ($base,$report,$name, $downloa\
d, $arguments, $language, $dir);\n    my $isdir;\n\
    &input_os();\n    \n    #\n    # - paolodt - C\
heck if the module exists in the repository cache \
\n    #\n	if( repo_load($pg) ) {\n	    $PG{$pg}{fr\
om_binary}=1;\n		return 1;\n	}\n    # - paolodt - \
end \n    \n    if (!&supported_os($OS)){return 0;\
}\n    if ( $PG{$pg}{binary}){$name=$PG{$pg}{binar\
y};}\n    else {$name=$pg;}\n    if ($name eq \"t_\
coffee\")\n      {\n	#check if local bin is there\\
n	if (-e \"./bin/$OS/t_coffee\")\n	  {\n	    print\
 \"\\n------- Installing  T-Coffee from Pre-Compil\
ed/Pre-Downloaded $OS binary\\n\";\n	    print \"\\
\n------- If you want to trigger a fresh compilati\
on use -recompile\\n\";\n	    &check_cp (\"./bin/$\
OS/t_coffee\", $BIN);\n	    return &pg_is_installe\
d ($pg, $BIN);\n	  }\n	#try to get precompiled bin\
ary -- available from MAC is distribution from MAC\
\n	else\n	  {\n	    $download=\"$WEB_BASE/Packages\
/Binaries/tcoffee/$OS/$name.$VERSION\";\n	  }\n   \
   }\n    else\n      {\n	$download=\"$WEB_BASE/Pa\
ckages/Binaries/plugins/$OS/$name\";\n      }\n   \
 \n    $base=cwd();\n    chdir $TMP;\n    \n    if\
 (!-e $name)\n      {\n	`rm x $SILENT`;\n	if ( url\
2file(\"$download\",\"x\")==$EXIT_SUCCESS)\n	  {\n\
	    `mv x $name`;\n	  }\n      }\n    \n    if (!\
-e $name)\n      {\n	print \"!!!!!!! $PG{$pg}{dnam\
e}: Download of $pg binary failed\\n\";\n	print \"\
!!!!!!! $PG{$pg}{dname}: Check Address: $download\\
\n\";\n	chdir $base;\n	return 0;\n      }\n    pri\
nt \"\\n------- Installing $pg\\n\";\n    \n    if\
 ($name =~/tar\\.gz/)\n      {\n	`gunzip  -f $name\
`;\n	`tar -xvf $pg.tar`;\n	chdir $pg;\n	`chmod u+x\
 *`;\n 	`mv * $BIN`;\n	#if (!($pg=~/\\*/)){`rm -rf\
 $pg`;}\n      }\n    else\n      {\n	&check_cp (\\
"$pg\", \"$BIN\");\n	`chmod u+x $BIN/$pg`; \n	unli\
nk ($pg);\n      }\n    chdir $base;\n    $PG{$pg}\
{from_binary}=1;\n\n    return &pg_is_installed ($\
pg, $BIN);\n  }\n\n	\nsub add_dir\n  {\n    my $di\
r=@_[0];\n    \n    if (!-e $dir && !-d $dir)\n   \
   {\n	my @l;\n	umask (0000);\n	@l=mkpath ($dir,{m\
ode => 0777});\n	\n      }\n    else\n      {\n	re\
turn 0;\n      }\n  }\nsub check_rm \n  {\n    my \
($file)=(@_);\n    \n    if ( -e $file)\n      {\n\
	return unlink($file);\n      }\n    return 0;\n  \
}\nsub check_cp\n  {\n    my ($from, $to)=(@_);\n \
   if ( !-e $from && -e \"$from\\.exe\"){$from=\"$\
from\\.exe\";}\n    if ( !-e $from){return 0;}\n  \
      \n    `$CP $from $to`;\n    return 1;\n  }\n\
\nsub repo_store \n{\n   # check that all required\
 data are available\n   if( $REPO_ROOT eq \"\" ) {\
 return; }\n\n\n    # extract the package name fro\
m the specified path\n    my $pg =`basename $_[0]`\
;\n    chomp($pg);\n	\n    my $VER = $PG{$pg}{vers\
ion};\n    my $CACHE = \"$REPO_ROOT/$pg/$VER/$OSNA\
ME-$OSARCH\"; \n    \n    print \"-------- Storing\
 package: \\\"$pg\\\" to path: $CACHE\\n\";\n    \\
n    # clean the cache path if exists and create i\
t again\n    `rm -rf $CACHE`;\n    `mkdir -p $CACH\
E`;\n    \n 	for my $path (@_) {\n\n	    # check i\
f it is a single file \n	 	if( -f $path ) {\n	    \
	`cp $path $CACHE`;\n		}\n		# .. or a directory, i\
n this case copy all the content \n		elsif( -d $pa\
th ) {\n			opendir(IMD, $path);\n			my @thefiles= \
readdir(IMD);\n			closedir(IMD);\n			\n			for my $\
_file (@thefiles) {\n				if( $_file ne \".\" && $_\
file ne \"..\") {\n	    			`cp $path/$_file $CACHE\
`;\n				}\n			}\n		} \n	}	   \n    \n	\n}   \n\nsu\
b repo_load \n{\n    my ($pg)=(@_);\n\n    #Bypass\
 the Repository Cache\n    return 0;\n    # check \
that all required data are available\n    if( $REP\
O_ROOT eq \"\" ) { return 0; }\n\n    my $VER = $P\
G{$pg}{version};\n    my $CACHE = \"$REPO_ROOT/$pg\
/$VER/$OSNAME-$OSARCH\"; \n    if( !-e \"$CACHE/$p\
g\" ) {\n   	 	print \"-------- Module \\\"$pg\\\"\
 NOT found on repository cache.\\n\";\n    	return\
 0;\n    }\n    \n    print \"-------- Module \\\"\
$pg\\\" found on repository cache. Using copy on p\
ath: $CACHE\\n\";\n    `cp $CACHE/* $BIN`;\n    re\
turn 1;\n}\n\nsub check_file_list_exists \n  {\n  \
  my ($base, @flist)=(@_);\n    my $f;\n\n    fore\
ach $f (@flist)\n      {\n	if ( !-e \"$base/$f\"){\
return 0;}\n      }\n    return 1;\n  }\nsub ls\n \
 {\n    my $f=@_[0];\n    my @fl;\n    chomp(@fl=`\
ls -1 $f`);\n    return @fl;\n  }\nsub flush_comma\
nd\n  {\n    my $command=@_[0];\n    my $F=new Fil\
eHandle;\n    open ($F, \"$command|\");\n    while\
 (<$F>){print \"    --- $_\";}\n    close ($F);\n \
 }    \n\nsub input_installation_directory\n  {\n \
   my $dir=@_[0];\n    my $new;\n    \n    print \\
"------- The current installation directory is: [$\
dir]\\n\";\n    print \"??????? Return to keep the\
 default or new value:\";\n   \n    if ($NO_QUESTI\
ON==0)\n      {\n	chomp ($new=<stdin>);\n	while ( \
$new ne \"\" && !input_yes (\"You have entered $ne\
w. Is this correct? ([y]/n):\"))\n	  {\n	    print\
 \"???????New installation directory:\";\n	    cho\
mp ($new=<stdin>);\n	  }\n	$dir=($new eq \"\")?$di\
r:$new;\n	$dir=~s/\\/$//;\n      }\n    \n    if (\
 -d $dir){return $dir;}\n    elsif (&root_run (\"Y\
ou must be root to create $dir\",\"mkdir $dir\")==\
$EXIT_SUCCESS){return $dir;}\n    else\n      {\n	\
print \"!!!!!!! $dir could not be created\\n\";\n	\
if ( $NO_QUESTION)\n	  {\n	    return \"\";\n	  }\\
n	elsif ( &input_yes (\"??????? Do you want to pro\
vide a new directory([y]/n)?:\"))\n	  {\n	    retu\
rn input_installation_directory ($dir);\n	  }\n	el\
se\n	  {\n	    return \"\";\n	  }\n      }\n    \n\
  }\nsub input_yes\n  {\n    my $question =@_[0];\\
n    my $answer;\n\n    if ($NO_QUESTION==1){retur\
n 1;}\n    \n    if ($question eq \"\"){$question=\
\"??????? Do you wish to proceed ([y]/n)?:\";}\n  \
  print $question;\n    chomp($answer=lc(<STDIN>))\
;\n    if (($answer=~/^y/) || $answer eq \"\"){ret\
urn 1;}\n    elsif ( ($answer=~/^n/)){return 0;}\n\
    else\n      {\n	return input_yes($question);\n\
      }\n  }\nsub root_run\n  {\n    my ($txt, $cm\
d)=(@_);\n    \n    if ( system ($cmd)==$EXIT_SUCC\
ESS){return $EXIT_SUCCESS;}\n    else \n      {\n	\
print \"------- $txt\\n\";\n	if ( $ROOT eq \"sudo\\
"){return system (\"sudo $cmd\");}\n	else {return \
system (\"su root -c \\\"$cmd\\\"\");}\n      }\n \
 }\nsub get_root\n  {\n    if (&pg_is_installed (\\
"sudo\")){return \"sudo\";}\n    else {return \"su\
\";}\n  }\n\nsub get_os\n  {\n    my $raw_os=`unam\
e`;\n    my $os;\n\n    $raw_os=lc ($raw_os);\n   \
 \n    if ($raw_os =~/cygwin/){$os=\"windows\";}\n\
    elsif ($raw_os =~/linux/){$os=\"linux\";}\n   \
 elsif ($raw_os =~/osx/){$os=\"macosx\";}\n    els\
if ($raw_os =~/darwin/){$os=\"macosx\";}\n    else\
\n      {\n	$os=$raw_os;\n      }\n    return $os;\
\n  }\nsub input_os\n  {\n    my $answer;\n    if \
($OS) {return $OS;}\n    \n    print \"??????? whi\
ch os do you use: [w]indows, [l]inux, [m]acosx:?\"\
;\n    $answer=lc(<STDIN>);\n\n    if (($answer=~/\
^m/)){$OS=\"macosx\";}\n    elsif ( ($answer=~/^w/\
)){$OS=\"windows\";}\n    elsif ( ($answer=~/^linu\
x/)){$OS=\"linux\";}\n    \n    else\n      {\n	re\
turn &input_os();\n      }\n    return $OS;\n  }\n\
\nsub supported_os\n  {\n    my ($os)=(@_[0]);\n  \
  return $SUPPORTED_OS{$os};\n  }\n\nsub add2env_f\
ile\n  {\n    my ($env, $var, $value)=(@_);\n    m\
y $F = new FileHandle;\n    my $t;\n    if (!$valu\
e){return;}\n    #make sure new variables do not g\
et duplicated\n    if ( -e $env)\n      {\n	open (\
$F, \"$env\");\n	while (<$F>)\n	  {\n	    my $line\
=$_;\n	    if (!($line=~/$var/)){$t.=$line;}\n	  }\
\n	close ($F);\n      }\n    $t.=\"$var=$value\\n\\
";\n    open ($F, \">$env\");\n    print $F \"$t\"\
;\n    $ENV{$var}=$value;\n    close ($F);\n  }   \
 \n    \n\n\nsub update_tclinkdb \n  {\n    my $fi\
le =@_[0];\n    my $name;\n    my $F=new FileHandl\
e;\n    my ($download, $address, $name, $l, $db);\\
n    \n    if ( $file eq \"update\"){$file=$TCLINK\
DB_ADDRESS;}\n    \n    if ( $file =~/http:\\/\\//\
 || $file =~/ftp:\\/\\//)\n      {\n	($address, $n\
ame)=($download=~/(.*)\\/([^\\/]+)$/);\n	`rm x $SI\
LENT`;\n	if (&url2file ($file,\"x\")==$EXIT_SUCCES\
S)\n	  {\n	    print \"------- Susscessful upload \
of $name\";\n	    `mv x $name`;\n	    $file=$name;\
\n	  }\n      }\n    open ($F, \"$file\");\n    wh\
ile (<$F>)\n      {\n	my $l=$_;\n	if (($l =~/^\\/\\
\//) || ($db=~/^#/)){;}\n	elsif ( !($l =~/\\w/)){;\
}\n	else\n	  {\n	    my @v=split (/\\s+/, $l);\n	 \
   if ( $l=~/^MODE/)\n	      {\n		$MODE{$v[1]}{$v[\
2]}=$v[3];\n	      }\n	    elsif ($l=~/^PG/)\n	   \
   {\n		$PG{$v[1]}{$v[2]}=$v[3];\n	      }\n	  }\n\
      }\n    close ($F);\n    &post_process_PG();\\
n    return;\n  }\n\n\n\nsub initialize_PG\n  {\n\\
n$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\"TCOFFEE\";\n$P\
G{\"t_coffee\"}{\"type\"}=\"sequence_multiple_alig\
ner\";\n$PG{\"t_coffee\"}{\"ADDRESS\"}=\"http://ww\
w.tcoffee.org\";\n$PG{\"t_coffee\"}{\"language\"}=\
\"C++\";\n$PG{\"t_coffee\"}{\"language2\"}=\"CXX\"\
;\n$PG{\"t_coffee\"}{\"source\"}=\"http://www.tcof\
fee.org/Packages/sources/tcoffee/stable/T-COFFEE_d\
istribution.tar.gz\";\n$PG{\"t_coffee\"}{\"update_\
action\"}=\"always\";\n$PG{\"t_coffee\"}{\"binary\\
"}=\"t_coffee\";\n$PG{\"t_coffee\"}{\"mode\"}=\"tc\
offee,mcoffee,rcoffee,expresso,3dcoffee\";\n$PG{\"\
clustalo\"}{\"4_TCOFFEE\"}=\"CLUSTALO\";\n$PG{\"cl\
ustalo\"}{\"type\"}=\"sequence_multiple_aligner\";\
\n$PG{\"clustalo\"}{\"ADDRESS\"}=\"http://www.clus\
tal.org/omega/\";\n$PG{\"clustalo\"}{\"language\"}\
=\"C++\";\n$PG{\"clustalo\"}{\"language2\"}=\"C++\\
";\n$PG{\"clustalo\"}{\"source\"}=\"http://www.clu\
stal.org/omega/clustal-omega-1.2.4.tar.gz\";\n$PG{\
\"clustalo\"}{\"mode\"}=\"mcoffee\";\n$PG{\"clusta\
lo\"}{\"binary\"}=\"clustalo\";\n$PG{\"clustalo\"}\
{\"version\"}=\"1.2.4\";\n$PG{\"strike\"}{\"4_TCOF\
FEE\"}=\"STRIKE\";\n$PG{\"strike\"}{\"type\"}=\"se\
quence_alignment_scoring\";\n$PG{\"strike\"}{\"ADD\
RESS\"}=\"http://www.tcoffee.org/Projects/strike/i\
ndex.html\";\n$PG{\"strike\"}{\"language\"}=\"C++\\
";\n$PG{\"strike\"}{\"language2\"}=\"CXX\";\n$PG{\\
"strike\"}{\"source\"}=\"http://www.tcoffee.org/Pr\
ojects/strike/strike_v1.2.tar.bz2\";\n$PG{\"strike\
\"}{\"mode\"}=\"tcoffee,expresso\";\n$PG{\"strike\\
"}{\"version\"}=\"1.2\";\n$PG{\"strike\"}{\"binary\
\"}=\"strike\";\n$PG{\"clustalw2\"}{\"4_TCOFFEE\"}\
=\"CLUSTALW2\";\n$PG{\"clustalw2\"}{\"type\"}=\"se\
quence_multiple_aligner\";\n$PG{\"clustalw2\"}{\"A\
DDRESS\"}=\"http://www.clustal.org\";\n$PG{\"clust\
alw2\"}{\"language\"}=\"C++\";\n$PG{\"clustalw2\"}\
{\"language2\"}=\"CXX\";\n$PG{\"clustalw2\"}{\"sou\
rce\"}=\"http://www.clustal.org/download/2.0.10/cl\
ustalw-2.0.10-src.tar.gz\";\n$PG{\"clustalw2\"}{\"\
mode\"}=\"mcoffee,rcoffee\";\n$PG{\"clustalw2\"}{\\
"binary\"}=\"clustalw2\";\n$PG{\"clustalw2\"}{\"ve\
rsion\"}=\"2.0.10\";\n$PG{\"clustalw\"}{\"4_TCOFFE\
E\"}=\"CLUSTALW\";\n$PG{\"clustalw\"}{\"type\"}=\"\
sequence_multiple_aligner\";\n$PG{\"clustalw\"}{\"\
ADDRESS\"}=\"http://www.clustal.org\";\n$PG{\"clus\
talw\"}{\"language\"}=\"C\";\n$PG{\"clustalw\"}{\"\
language2\"}=\"C\";\n$PG{\"clustalw\"}{\"source\"}\
=\"http://www.clustal.org/download/1.X/ftp-igbmc.u\
-strasbg.fr/pub/ClustalW/clustalw1.82.UNIX.tar.gz\\
";\n$PG{\"clustalw\"}{\"mode\"}=\"mcoffee,rcoffee\\
";\n$PG{\"clustalw\"}{\"version\"}=\"1.82\";\n$PG{\
\"clustalw\"}{\"binary\"}=\"clustalw\";\n$PG{\"dia\
lign-t\"}{\"4_TCOFFEE\"}=\"DIALIGNT\";\n$PG{\"dial\
ign-t\"}{\"type\"}=\"sequence_multiple_aligner\";\\
n$PG{\"dialign-t\"}{\"ADDRESS\"}=\"http://dialign-\
tx.gobics.de/\";\n$PG{\"dialign-t\"}{\"DIR\"}=\"/u\
sr/share/dialign-tx/\";\n$PG{\"dialign-t\"}{\"lang\
uage\"}=\"C\";\n$PG{\"dialign-t\"}{\"language2\"}=\
\"C\";\n$PG{\"dialign-t\"}{\"source\"}=\"http://di\
align-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz\";\n$PG\
{\"dialign-t\"}{\"mode\"}=\"mcoffee\";\n$PG{\"dial\
ign-t\"}{\"binary\"}=\"dialign-t\";\n$PG{\"dialign\
-t\"}{\"version\"}=\"1.0.2\";\n$PG{\"dialign-tx\"}\
{\"4_TCOFFEE\"}=\"DIALIGNTX\";\n$PG{\"dialign-tx\"\
}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"\
dialign-tx\"}{\"ADDRESS\"}=\"http://dialign-tx.gob\
ics.de/\";\n$PG{\"dialign-tx\"}{\"DIR\"}=\"/usr/sh\
are/dialign-tx/\";\n$PG{\"dialign-tx\"}{\"language\
\"}=\"C\";\n$PG{\"dialign-tx\"}{\"language2\"}=\"C\
\";\n$PG{\"dialign-tx\"}{\"source\"}=\"http://dial\
ign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz\";\n$PG{\\
"dialign-tx\"}{\"mode\"}=\"mcoffee\";\n$PG{\"diali\
gn-tx\"}{\"binary\"}=\"dialign-tx\";\n$PG{\"dialig\
n-tx\"}{\"version\"}=\"1.0.2\";\n$PG{\"poa\"}{\"4_\
TCOFFEE\"}=\"POA\";\n$PG{\"poa\"}{\"type\"}=\"sequ\
ence_multiple_aligner\";\n$PG{\"poa\"}{\"ADDRESS\"\
}=\"http://www.bioinformatics.ucla.edu/poa/\";\n$P\
G{\"poa\"}{\"language\"}=\"C\";\n$PG{\"poa\"}{\"la\
nguage2\"}=\"C\";\n$PG{\"poa\"}{\"source\"}=\"http\
://downloads.sourceforge.net/poamsa/poaV2.tar.gz\"\
;\n$PG{\"poa\"}{\"DIR\"}=\"/usr/share/\";\n$PG{\"p\
oa\"}{\"FILE1\"}=\"blosum80.mat\";\n$PG{\"poa\"}{\\
"mode\"}=\"mcoffee\";\n$PG{\"poa\"}{\"binary\"}=\"\
poa\";\n$PG{\"poa\"}{\"version\"}=\"2.0\";\n$PG{\"\
probcons\"}{\"4_TCOFFEE\"}=\"PROBCONS\";\n$PG{\"pr\
obcons\"}{\"type\"}=\"sequence_multiple_aligner\";\
\n$PG{\"probcons\"}{\"ADDRESS\"}=\"http://probcons\
.stanford.edu/\";\n$PG{\"probcons\"}{\"language2\"\
}=\"CXX\";\n$PG{\"probcons\"}{\"language\"}=\"C++\\
";\n$PG{\"probcons\"}{\"source\"}=\"http://probcon\
s.stanford.edu/probcons_v1_12.tar.gz\";\n$PG{\"pro\
bcons\"}{\"mode\"}=\"mcoffee\";\n$PG{\"probcons\"}\
{\"binary\"}=\"probcons\";\n$PG{\"probcons\"}{\"ve\
rsion\"}=\"1.12\";\n$PG{\"msaprobs\"}{\"4_TCOFFEE\\
"}=\"MSAPROBS\";\n$PG{\"msaprobs\"}{\"type\"}=\"se\
quence_multiple_aligner\";\n$PG{\"msaprobs\"}{\"AD\
DRESS\"}=\"http://msaprobs.sourceforge.net/homepag\
e.htm#latest\";\n$PG{\"msaprobs\"}{\"language2\"}=\
\"CXX\";\n$PG{\"msaprobs\"}{\"language\"}=\"C++\";\
\n$PG{\"msaprobs\"}{\"source\"}=\"https://sourcefo\
rge.net/projects/msaprobs/files/MSAProbs-MPI/MSAPr\
obs-MPI_rel1.0.5.tar.gz\";\n$PG{\"msaprobs\"}{\"mo\
de\"}=\"mcoffee\";\n$PG{\"msaprobs\"}{\"binary\"}=\
\"msaprobs\";\n$PG{\"msaprobs\"}{\"version\"}=\"1.\
05\";\n$PG{\"msaprobs\"}{\"update_action\"}=\"neve\
r\";\n$PG{\"upp\"}{\"4_TCOFFEE\"}=\"UPP\";\n$PG{\"\
upp\"}{\"type\"}=\"sequence_multiple_aligner\";\n$\
PG{\"upp\"}{\"ADDRESS\"}=\"http://www.cs.utexas.ed\
u/users/phylo/software/upp/\";\n$PG{\"upp\"}{\"lan\
guage2\"}=\"CXX\";\n$PG{\"upp\"}{\"language\"}=\"C\
++\";\n$PG{\"upp\"}{\"source\"}=\"https://github.c\
om/smirarab/pasta/archive/upp.zip\";\n$PG{\"upp\"}\
{\"mode\"}=\"mcoffee\";\n$PG{\"upp\"}{\"binary\"}=\
\"upp\";\n$PG{\"upp\"}{\"version\"}=\"1\";\n$PG{\"\
upp\"}{\"update_action\"}=\"never\";\n$PG{\"famsa\\
"}{\"4_TCOFFEE\"}=\"FAMSA\";\n$PG{\"famsa\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"famsa\"\
}{\"ADDRESS\"}=\"https://github.com/refresh-bio/FA\
MSA\";\n$PG{\"famsa\"}{\"language\"}=\"C++\";\n$PG\
{\"famsa\"}{\"language\"}=\"C++\";\n$PG{\"famsa\"}\
{\"source\"}=\"https://github.com/refresh-bio/FAMS\
A.git\";\n$PG{\"famsa\"}{\"mode\"}=\"mcoffee,rcoff\
ee\";\n$PG{\"famsa\"}{\"binary\"}=\"famsa\";\n$PG{\
\"famsa\"}{\"version\"}=\"1.1\";\n$PG{\"mafft\"}{\\
"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"mafft\"}{\"type\"}\
=\"sequence_multiple_aligner\";\n$PG{\"mafft\"}{\"\
ADDRESS\"}=\"http://align.bmr.kyushu-u.ac.jp/mafft\
/online/server/\";\n$PG{\"mafft\"}{\"language\"}=\\
"C\";\n$PG{\"mafft\"}{\"language\"}=\"C\";\n$PG{\"\
mafft\"}{\"source\"}=\"http://mafft.cbrc.jp/alignm\
ent/software/mafft-7.310-with-extensions-src.tgz\"\
;\n$PG{\"mafft\"}{\"mode\"}=\"mcoffee,rcoffee\";\n\
$PG{\"mafft\"}{\"binary\"}=\"mafft.tar.gz\";\n$PG{\
\"mafft\"}{\"version\"}=\"7.310\";\n$PG{\"msa\"}{\\
"4_TCOFFEE\"}=\"MSA\";\n$PG{\"msa\"}{\"type\"}=\"s\
equence_multiple_aligner\";\n$PG{\"msa\"}{\"ADDRES\
S\"}=\"https://www.ncbi.nlm.nih.gov/CBBresearch/Sc\
haffer/msa.html\";\n$PG{\"msa\"}{\"language\"}=\"C\
\";\n$PG{\"msa\"}{\"language\"}=\"C\";\n$PG{\"msa\\
"}{\"source\"}=\"ftp://ftp.ncbi.nih.gov/pub/msa/ms\
a.tar.Z\";\n$PG{\"msa\"}{\"mode\"}=\"mcoffee\";\n$\
PG{\"msa\"}{\"binary\"}=\"msa.pl\";\n$PG{\"msa\"}{\
\"version\"}=\"1.0\";\n$PG{\"msa\"}{\"update_actio\
n\"}=\"never\";\n$PG{\"dca\"}{\"4_TCOFFEE\"}=\"DCA\
\";\n$PG{\"dca\"}{\"type\"}=\"sequence_multiple_al\
igner\";\n$PG{\"dca\"}{\"ADDRESS\"}=\"https://bibi\
serv2.cebitec.uni-bielefeld.de/dca\";\n$PG{\"dca\"\
}{\"language\"}=\"C\";\n$PG{\"dca\"}{\"language\"}\
=\"C\";\n$PG{\"dca\"}{\"source\"}=\"https://bibise\
rv2.cebitec.uni-bielefeld.de/applications/dca/reso\
urces/downloads/dca-1.1-src.tar.gz\";\n$PG{\"dca\"\
}{\"mode\"}=\"mcoffee\";\n$PG{\"dca\"}{\"binary\"}\
=\"dca.pl\";\n$PG{\"dca\"}{\"version\"}=\"1.1\";\n\
$PG{\"dca\"}{\"update_action\"}=\"never\";\n$PG{\"\
muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\";\n$PG{\"muscle\
\"}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\
\"muscle\"}{\"ADDRESS\"}=\"http://www.drive5.com/m\
uscle/\";\n$PG{\"muscle\"}{\"language\"}=\"C++\";\\
n$PG{\"muscle\"}{\"language2\"}=\"GPP\";\n$PG{\"mu\
scle\"}{\"source\"}=\"http://www.drive5.com/muscle\
/downloads3.7/muscle3.7_src.tar.gz\";\n$PG{\"muscl\
e\"}{\"windows\"}=\"http://www.drive5.com/muscle/d\
ownloads3.7/muscle3.7_win32.zip\";\n$PG{\"muscle\"\
}{\"linux\"}=\"http://www.drive5.com/muscle/downlo\
ads3.7/muscle3.7_linux_ia32.tar.gz\";\n$PG{\"muscl\
e\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"muscle\\
"}{\"version\"}=\"3.7\";\n$PG{\"pcma\"}{\"4_TCOFFE\
E\"}=\"PCMA\";\n$PG{\"pcma\"}{\"type\"}=\"sequence\
_multiple_aligner\";\n$PG{\"pcma\"}{\"ADDRESS\"}=\\
"http://prodata.swmed.edu/pcma/pcma.php\";\n$PG{\"\
pcma\"}{\"language\"}=\"C\";\n$PG{\"pcma\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"pcma\"}{\"source\"}=\"http:\
//prodata.swmed.edu/download/pub/PCMA/pcma.tar.gz\\
";\n$PG{\"pcma\"}{\"mode\"}=\"mcoffee\";\n$PG{\"pc\
ma\"}{\"version\"}=\"1.0\";\n$PG{\"kalign\"}{\"4_T\
COFFEE\"}=\"KALIGN\";\n$PG{\"kalign\"}{\"type\"}=\\
"sequence_multiple_aligner\";\n$PG{\"kalign\"}{\"A\
DDRESS\"}=\"http://msa.cgb.ki.se\";\n$PG{\"kalign\\
"}{\"language\"}=\"C\";\n$PG{\"kalign\"}{\"languag\
e2\"}=\"C\";\n$PG{\"kalign\"}{\"source\"}=\"http:/\
/msa.cgb.ki.se/downloads/kalign/current.tar.gz\";\\
n$PG{\"kalign\"}{\"mode\"}=\"mcoffee\";\n$PG{\"kal\
ign\"}{\"version\"}=\"1.0\";\n$PG{\"amap\"}{\"4_TC\
OFFEE\"}=\"AMAP\";\n$PG{\"amap\"}{\"type\"}=\"sequ\
ence_multiple_aligner\";\n$PG{\"amap\"}{\"ADDRESS\\
"}=\"http://bio.math.berkeley.edu/amap/\";\n$PG{\"\
amap\"}{\"language\"}=\"C++\";\n$PG{\"amap\"}{\"la\
nguage2\"}=\"CXX\";\n$PG{\"amap\"}{\"source\"}=\"h\
ttps://github.com/mes5k/amap-align/archive/amap.zi\
p\";\n$PG{\"amap\"}{\"mode\"}=\"mcoffee\";\n$PG{\"\
amap\"}{\"version\"}=\"2.0\";\n$PG{\"amap\"}{\"upd\
ate_action\"}=\"never\";\n$PG{\"proda\"}{\"4_TCOFF\
EE\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\"}=\"seque\
nce_multiple_aligner\";\n$PG{\"proda\"}{\"ADDRESS\\
"}=\"http://proda.stanford.edu\";\n$PG{\"proda\"}{\
\"language\"}=\"C++\";\n$PG{\"proda\"}{\"language2\
\"}=\"CXX\";\n$PG{\"proda\"}{\"source\"}=\"http://\
proda.stanford.edu/proda_1_0.tar.gz\";\n$PG{\"prod\
a\"}{\"mode\"}=\"mcoffee\";\n$PG{\"proda\"}{\"vers\
ion\"}=\"1.0\";\n$PG{\"prank\"}{\"4_TCOFFEE\"}=\"P\
RANK\";\n$PG{\"prank\"}{\"type\"}=\"sequence_multi\
ple_aligner\";\n$PG{\"prank\"}{\"ADDRESS\"}=\"http\
://www.ebi.ac.uk/goldman-srv/prank/\";\n$PG{\"pran\
k\"}{\"language\"}=\"C++\";\n$PG{\"prank\"}{\"lang\
uage2\"}=\"CXX\";\n$PG{\"prank\"}{\"source\"}=\"ht\
tp://www.ebi.ac.uk/goldman-srv/prank/src/prank/pra\
nk.src.100802.tgz\";\n$PG{\"prank\"}{\"mode\"}=\"m\
coffee\";\n$PG{\"prank\"}{\"version\"}=\"100303\";\
\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\n$PG{\"sap\\
"}{\"type\"}=\"structure_pairwise_aligner\";\n$PG{\
\"sap\"}{\"ADDRESS\"}=\"https://mathbio.crick.ac.u\
k/wiki/Software#SAP\";\n$PG{\"sap\"}{\"language\"}\
=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"C\";\n$PG{\\
"sap\"}{\"source\"}=\"https://github.com/jkleinj/S\
AP/archive/v.1.1.3.tar.gz\";\n$PG{\"sap\"}{\"mode\\
"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}{\"version\\
"}=\"1.1.3\";\n$PG{\"sap\"}{\"binary\"}=\"sap\";\n\
$PG{\"TMalign\"}{\"4_TCOFFEE\"}=\"TMALIGN\";\n$PG{\
\"TMalign\"}{\"type\"}=\"structure_pairwise_aligne\
r\";\n$PG{\"TMalign\"}{\"ADDRESS\"}=\"http://zhang\
lab.ccmb.med.umich.edu/TM-align/TMalign.f\";\n$PG{\
\"TMalign\"}{\"language\"}=\"Fortran\";\n$PG{\"TMa\
lign\"}{\"language2\"}=\"Fortran\";\n$PG{\"TMalign\
\"}{\"source\"}=\"http://zhanglab.ccmb.med.umich.e\
du/TM-align/TMalign.f\";\n$PG{\"TMalign\"}{\"linux\
\"}=\"http://zhanglab.ccmb.med.umich.edu/TM-align/\
TMalign_32.gz\";\n$PG{\"TMalign\"}{\"mode\"}=\"exp\
resso,3dcoffee\";\n$PG{\"TMalign\"}{\"version\"}=\\
"2013.05.11\";\n$PG{\"mustang\"}{\"4_TCOFFEE\"}=\"\
MUSTANG\";\n$PG{\"mustang\"}{\"type\"}=\"structure\
_pairwise_aligner\";\n$PG{\"mustang\"}{\"ADDRESS\"\
}=\"http://lcb.infotech.monash.edu.au/mustang/\";\\
n$PG{\"mustang\"}{\"language\"}=\"C++\";\n$PG{\"mu\
stang\"}{\"language2\"}=\"CXX\";\n$PG{\"mustang\"}\
{\"source\"}=\"http://lcb.infotech.monash.edu.au/m\
ustang/mustang_v3.2.3.tgz\";\n$PG{\"mustang\"}{\"m\
ode\"}=\"expresso,3dcoffee\";\n$PG{\"mustang\"}{\"\
version\"}=\"3.2.3\";\n$PG{\"lsqman\"}{\"4_TCOFFEE\
\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\"}=\"struc\
ture_pairwise_aligner\";\n$PG{\"lsqman\"}{\"ADDRES\
S\"}=\"empty\";\n$PG{\"lsqman\"}{\"language\"}=\"e\
mpty\";\n$PG{\"lsqman\"}{\"language2\"}=\"empty\";\
\n$PG{\"lsqman\"}{\"source\"}=\"empty\";\n$PG{\"ls\
qman\"}{\"update_action\"}=\"never\";\n$PG{\"lsqma\
n\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"align\
_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n$PG{\"align\
_pdb\"}{\"type\"}=\"structure_pairwise_aligner\";\\
n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empty\";\n$PG{\\
"align_pdb\"}{\"language\"}=\"empty\";\n$PG{\"alig\
n_pdb\"}{\"language2\"}=\"empty\";\n$PG{\"align_pd\
b\"}{\"source\"}=\"empty\";\n$PG{\"align_pdb\"}{\"\
update_action\"}=\"never\";\n$PG{\"align_pdb\"}{\"\
mode\"}=\"expresso,3dcoffee\";\n$PG{\"fugueali\"}{\
\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fugueali\"}{\"typ\
e\"}=\"structure_pairwise_aligner\";\n$PG{\"fuguea\
li\"}{\"ADDRESS\"}=\"http://mizuguchilab.org/fugue\
/\";\n$PG{\"fugueali\"}{\"language\"}=\"empty\";\n\
$PG{\"fugueali\"}{\"language2\"}=\"empty\";\n$PG{\\
"fugueali\"}{\"source\"}=\"empty\";\n$PG{\"fugueal\
i\"}{\"update_action\"}=\"never\";\n$PG{\"fugueali\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"dalili\
te.pl\"}{\"4_TCOFFEE\"}=\"DALILITEc\";\n$PG{\"dali\
lite.pl\"}{\"type\"}=\"structure_pairwise_aligner\\
";\n$PG{\"dalilite.pl\"}{\"ADDRESS\"}=\"built_in\"\
;\n$PG{\"dalilite.pl\"}{\"ADDRESS2\"}=\"http://www\
.ebi.ac.uk/Tools/webservices/services/dalilite\";\\
n$PG{\"dalilite.pl\"}{\"language\"}=\"Perl\";\n$PG\
{\"dalilite.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"\
dalilite.pl\"}{\"source\"}=\"empty\";\n$PG{\"dalil\
ite.pl\"}{\"update_action\"}=\"never\";\n$PG{\"dal\
ilite.pl\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\
\"probconsRNA\"}{\"4_TCOFFEE\"}=\"PROBCONSRNA\";\n\
$PG{\"probconsRNA\"}{\"type\"}=\"RNA_multiple_alig\
ner\";\n$PG{\"probconsRNA\"}{\"ADDRESS\"}=\"http:/\
/probcons.stanford.edu/\";\n$PG{\"probconsRNA\"}{\\
"language\"}=\"C++\";\n$PG{\"probconsRNA\"}{\"lang\
uage2\"}=\"CXX\";\n$PG{\"probconsRNA\"}{\"source\"\
}=\"http://probcons.stanford.edu/probconsRNA.tar.g\
z\";\n$PG{\"probconsRNA\"}{\"mode\"}=\"mcoffee,rco\
ffee\";\n$PG{\"probconsRNA\"}{\"version\"}=\"1.0\"\
;\n$PG{\"sfold\"}{\"4_TCOFFEE\"}=\"CONSAN\";\n$PG{\
\"sfold\"}{\"type\"}=\"RNA_pairwise_aligner\";\n$P\
G{\"sfold\"}{\"ADDRESS\"}=\"http://selab.janelia.o\
rg/software/consan/\";\n$PG{\"sfold\"}{\"language\\
"}=\"empty\";\n$PG{\"sfold\"}{\"language2\"}=\"emp\
ty\";\n$PG{\"sfold\"}{\"source\"}=\"empty\";\n$PG{\
\"sfold\"}{\"update_action\"}=\"never\";\n$PG{\"sf\
old\"}{\"mode\"}=\"rcoffee\";\n$PG{\"RNAplfold\"}{\
\"4_TCOFFEE\"}=\"RNAPLFOLD\";\n$PG{\"RNAplfold\"}{\
\"type\"}=\"RNA_secondarystructure_predictor\";\n$\
PG{\"RNAplfold\"}{\"ADDRESS\"}=\"http://www.tbi.un\
ivie.ac.at/RNA/\";\n$PG{\"RNAplfold\"}{\"language\\
"}=\"C\";\n$PG{\"RNAplfold\"}{\"language2\"}=\"C\"\
;\n$PG{\"RNAplfold\"}{\"source\"}=\"http://www.tbi\
.univie.ac.at/RNA/packages/source/ViennaRNA-2.1.9.\
tar.gz\";\n$PG{\"RNAplfold\"}{\"mode\"}=\"rcoffee,\
\";\n$PG{\"RNAplfold\"}{\"binary\"}=\"RNAplfold.ta\
r.gz\";\n$PG{\"RNAplfold\"}{\"version\"}=\"2.1.9\"\
;\n$PG{\"retree\"}{\"4_TCOFFEE\"}=\"PHYLIP\";\n$PG\
{\"retree\"}{\"type\"}=\"Phylogeny\";\n$PG{\"retre\
e\"}{\"ADDRESS\"}=\"http://evolution.gs.washington\
.edu/phylip/\";\n$PG{\"retree\"}{\"language\"}=\"C\
\";\n$PG{\"retree\"}{\"language2\"}=\"C\";\n$PG{\"\
retree\"}{\"source\"}=\"http://www.tcoffee.org/Pac\
kages/mirrors/source/phylip-3.66.tar.gz\";\n$PG{\"\
retree\"}{\"mode\"}=\"trmsd,\";\n$PG{\"retree\"}{\\
"binary\"}=\"retree.tar.gz\";\n$PG{\"retree\"}{\"v\
ersion\"}=\"3.66\";\n$PG{\"hmmtop\"}{\"4_TCOFFEE\"\
}=\"HMMTOP\";\n$PG{\"hmmtop\"}{\"type\"}=\"protein\
_secondarystructure_predictor\";\n$PG{\"hmmtop\"}{\
\"ADDRESS\"}=\"www.enzim.hu/hmmtop/\";\n$PG{\"hmmt\
op\"}{\"language\"}=\"C\";\n$PG{\"hmmtop\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"hmmtop\"}{\"source\"}=\"htt\
p://www.tcoffee.org/Packages/mirrors/hmmtop2.1.tgz\
\";\n$PG{\"hmmtop\"}{\"binary\"}=\"hmmtop\";\n$PG{\
\"hmmtop\"}{\"update_action\"}=\"never\";\n$PG{\"h\
mmtop\"}{\"mode\"}=\"psicoffee\";\n$PG{\"hmmtop\"}\
{\"version\"}=\"2.1\";\n$PG{\"gorIV\"}{\"4_TCOFFEE\
\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"type\"}=\"protein_\
secondarystructure_predictor\";\n$PG{\"gorIV\"}{\"\
ADDRESS\"}=\"http://mig.jouy.inra.fr/logiciels/gor\
IV/\";\n$PG{\"gorIV\"}{\"language\"}=\"C\";\n$PG{\\
"gorIV\"}{\"language2\"}=\"C\";\n$PG{\"gorIV\"}{\"\
source\"}=\"http://www.tcoffee.org/Packages/mirror\
s/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}{\"update_action\
\"}=\"never\";\n$PG{\"gorIV\"}{\"mode\"}=\"tcoffee\
\";\n$PG{\"wublast.pl\"}{\"4_TCOFFEE\"}=\"EBIWUBLA\
STc\";\n$PG{\"wublast.pl\"}{\"type\"}=\"protein_ho\
mology_predictor\";\n$PG{\"wublast.pl\"}{\"ADDRESS\
\"}=\"built_in\";\n$PG{\"wublast.pl\"}{\"ADDRESS2\\
"}=\"http://www.ebi.ac.uk/Tools/webservices/servic\
es/wublast\";\n$PG{\"wublast.pl\"}{\"language\"}=\\
"Perl\";\n$PG{\"wublast.pl\"}{\"language2\"}=\"Per\
l\";\n$PG{\"wublast.pl\"}{\"source\"}=\"empty\";\n\
$PG{\"wublast.pl\"}{\"update_action\"}=\"never\";\\
n$PG{\"wublast.pl\"}{\"mode\"}=\"psicoffee,express\
o,accurate\";\n$PG{\"blastpgp.pl\"}{\"4_TCOFFEE\"}\
=\"EBIBLASTPGPc\";\n$PG{\"blastpgp.pl\"}{\"type\"}\
=\"protein_homology_predictor\";\n$PG{\"blastpgp.p\
l\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"blastpgp.pl\
\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/web\
services/services/blastpgp\";\n$PG{\"blastpgp.pl\"\
}{\"language\"}=\"Perl\";\n$PG{\"blastpgp.pl\"}{\"\
language2\"}=\"Perl\";\n$PG{\"blastpgp.pl\"}{\"sou\
rce\"}=\"empty\";\n$PG{\"blastpgp.pl\"}{\"update_a\
ction\"}=\"never\";\n$PG{\"blastpgp.pl\"}{\"mode\"\
}=\"psicoffee,expresso,accurate\";\n$PG{\"blastall\
\"}{\"4_TCOFFEE\"}=\"blastall\";\n$PG{\"blastall\"\
}{\"type\"}=\"protein_homology_predictor\";\n$PG{\\
"blastall\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.gov\
/blast/executables/LATEST\";\n$PG{\"blastall\"}{\"\
language\"}=\"C\";\n$PG{\"blastall\"}{\"language2\\
"}=\"C\";\n$PG{\"blastall\"}{\"source\"}=\"ftp://f\
tp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0\
/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG{\"blastall\"\
}{\"update_action\"}=\"never\";\n$PG{\"blastall\"}\
{\"mode\"}=\"psicoffee,expresso,3dcoffee\";\n$PG{\\
"legacy_blast.pl\"}{\"4_TCOFFEE\"}=\"NCBIBLAST\";\\
n$PG{\"legacy_blast.pl\"}{\"type\"}=\"protein_homo\
logy_predictor\";\n$PG{\"legacy_blast.pl\"}{\"ADDR\
ESS\"}=\"ftp://ftp.ncbi.nih.gov/blast/executables/\
LATEST\";\n$PG{\"legacy_blast.pl\"}{\"language\"}=\
\"C\";\n$PG{\"legacy_blast.pl\"}{\"language2\"}=\"\
C\";\n$PG{\"legacy_blast.pl\"}{\"source\"}=\"ftp:/\
/ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6\
.0/ncbi-blast-2.6.0+-src.tar.gz\";\n$PG{\"legacy_b\
last.pl\"}{\"update_action\"}=\"never\";\n$PG{\"le\
gacy_blast.pl\"}{\"mode\"}=\"psicoffee,expresso,3d\
coffee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOFFEE\"}=\"SO\
APLITE\";\n$PG{\"SOAP::Lite\"}{\"type\"}=\"library\
\";\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\"http://cpa\
nsearch.perl.org/src/MKUTTER/SOAP-Lite-0.710.08/Ma\
kefile.PL\";\n$PG{\"SOAP::Lite\"}{\"language\"}=\"\
Perl\";\n$PG{\"SOAP::Lite\"}{\"language2\"}=\"Perl\
\";\n$PG{\"SOAP::Lite\"}{\"source\"}=\"empty\";\n$\
PG{\"SOAP::Lite\"}{\"update_action\"}=\"never\";\n\
$PG{\"SOAP::Lite\"}{\"mode\"}=\"none\";\n$PG{\"XML\
::Simple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\";\n$PG{\"X\
ML::Simple\"}{\"type\"}=\"library\";\n$PG{\"XML::S\
imple\"}{\"ADDRESS\"}=\"http://search.cpan.org/~gr\
antm/XML-Simple-2.18/lib/XML/Simple.pm\";\n$PG{\"X\
ML::Simple\"}{\"language\"}=\"Perl\";\n$PG{\"XML::\
Simple\"}{\"language2\"}=\"Perl\";\n$PG{\"XML::Sim\
ple\"}{\"source\"}=\"empty\";\n$PG{\"XML::Simple\"\
}{\"mode\"}=\"psicoffee,expresso,accurate\";\n$PG{\
\"x3dna\"}{\"4_TCOFFEE\"}=\"x3dna\";\n$PG{\"x3dna\\
"}{\"type\"}=\"RNA_secondarystructure_predictor\";\
\n$PG{\"x3dna\"}{\"ADDRESS\"}=\"http://x3dna.bio.c\
olumbia.edu/\";\n$PG{\"x3dna\"}{\"source\"}=\"http\
://www.tcoffee.org/Packages/mirrors/source/x3dna-v\
2.3-linux-64bit.tar.gz\";\n$PG{\"x3dna\"}{\"mode\"\
}=\"saracoffee\";\n$PG{\"x3dna\"}{\"update_action\\
"}=\"never\";\n$PG{\"fsa\"}{\"4_TCOFFEE\"}=\"FSA\"\
;\n$PG{\"fsa\"}{\"type\"}=\"sequence_multiple_alig\
ner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"http://fsa.sou\
rceforge.net/\";\n$PG{\"fsa\"}{\"language\"}=\"C++\
\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\";\n$PG{\"f\
sa\"}{\"source\"}=\"http://sourceforge.net/project\
s/fsa/files/fsa-1.15.3.tar.gz/download/\";\n$PG{\"\
fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fsa\"}{\"vers\
ion\"}=\"1.15.3\";\n$PG{\"fsa\"}{\"update_action\"\
}=\"never\";\n$PG{\"mus4\"}{\"4_TCOFFEE\"}=\"MUS4\\
";\n$PG{\"mus4\"}{\"type\"}=\"sequence_multiple_al\
igner\";\n$PG{\"mus4\"}{\"ADDRESS\"}=\"http://www.\
drive5.com/muscle/\";\n$PG{\"mus4\"}{\"language\"}\
=\"C++\";\n$PG{\"mus4\"}{\"language2\"}=\"GPP\";\n\
$PG{\"mus4\"}{\"source\"}=\"http://www.drive5.com/\
muscle/muscle4.0_src.tar.gz\";\n$PG{\"mus4\"}{\"mo\
de\"}=\"mcoffee,rcoffee\";\n$PG{\"mus4\"}{\"versio\
n\"}=\"4.0\";\n$PG{\"mus4\"}{\"update_action\"}=\"\
never\";\n$MODE{\"tcoffee\"}{\"name\"}=\"tcoffee\"\
;\n$MODE{\"rcoffee\"}{\"name\"}=\"rcoffee\";\n$MOD\
E{\"3dcoffee\"}{\"name\"}=\"3dcoffee\";\n$MODE{\"m\
coffee\"}{\"name\"}=\"mcoffee\";\n$MODE{\"expresso\
\"}{\"name\"}=\"expresso\";\n$MODE{\"trmsd\"}{\"na\
me\"}=\"trmsd\";\n$MODE{\"accurate\"}{\"name\"}=\"\
accurate\";\n$MODE{\"seq_reformat\"}{\"name\"}=\"s\
eq_reformat\";\n\n\n$PG{C}{compiler}=\"gcc\";\n$PG\
{C}{compiler_flag}=\"CC\";\n$PG{C}{options}=\"\";\\
n$PG{C}{options_flag}=\"CFLAGS\";\n$PG{C}{type}=\"\
compiler\";\n\n$PG{\"CXX\"}{compiler}=\"g++\";\n$P\
G{\"CXX\"}{compiler_flag}=\"CXX\";\n$PG{\"CXX\"}{o\
ptions}=\"\";\n$PG{\"CXX\"}{options_flag}=\"CXXFLA\
GS\";\n$PG{CXX}{type}=\"compiler\";\n\n$PG{\"CPP\"\
}{compiler}=\"g++\";\n$PG{\"CPP\"}{compiler_flag}=\
\"CPP\";\n$PG{\"CPP\"}{options}=\"\";\n$PG{\"CPP\"\
}{options_flag}=\"CPPFLAGS\";\n$PG{CPP}{type}=\"co\
mpiler\";\n\n$PG{\"GPP\"}{compiler}=\"g++\";\n$PG{\
\"GPP\"}{compiler_flag}=\"GPP\";\n$PG{\"GPP\"}{opt\
ions}=\"\";\n$PG{\"GPP\"}{options_flag}=\"CFLAGS\"\
;\n$PG{GPP}{type}=\"compiler\";\n\n$PG{Fortran}{co\
mpiler}=\"g77\";\n$PG{Fortran}{compiler_flag}=\"FC\
C\";\n$PG{Fortran}{type}=\"compiler\";\n\n$PG{Perl\
}{compiler}=\"CPAN\";\n$PG{Perl}{type}=\"compiler\\
";\n\n$SUPPORTED_OS{macosx}=\"Macintosh\";\n$SUPPO\
RTED_OS{linux}=\"Linux\";\n$SUPPORTED_OS{windows}=\
\"Cygwin\";\n\n\n\n$MODE{t_coffee}{description}=\"\
 for regular multiple sequence alignments\";\n$MOD\
E{rcoffee} {description}=\" for RNA multiple seque\
nce alignments\";\n\n$MODE{psicoffee} {description\
}=\" for Homology Extended multiple sequence align\
ments\";\n$MODE{expresso}{description}=\" for very\
 accurate structure based multiple sequence alignm\
ents\";\n$MODE{\"3dcoffee\"}{description}=\" for m\
ultiple structure alignments\";\n$MODE{mcoffee} {d\
escription}=\" for combining alternative multiple \
sequence alignment packages\\n------- into a uniqu\
e meta-package. The installer will upload several \
MSA packages and compile them\\n\n\";\n\n\n&post_p\
rocess_PG();\nreturn;\n}\n\nsub post_process_PG\n \
 {\n    my $p;\n    \n    %PG=&name2dname (%PG);\n\
    %MODE=&name2dname(%MODE);\n    foreach $p (key\
s(%PG)){if ( $PG{$p}{type} eq \"compiler\"){$PG{$p\
}{update_action}=\"never\";}}\n    \n  }\n\nsub na\
me2dname\n  {\n    my (%L)=(@_);\n    my ($l, $ml)\
;\n    \n    foreach my $pg (keys(%L))\n      {\n	\
$l=length ($pg);\n	if ( $l>$ml){$ml=$l;}\n      }\\
n    $ml+=1;\n    foreach my $pg (keys(%L))\n     \
 {\n	my $name;\n	$l=$ml-length ($pg);\n	$name=$pg;\
\n	for ( $b=0; $b<$l; $b++)\n	  {\n	    $name .=\"\
 \";\n	  }\n	$L{$pg}{dname}=$name;\n      }\n    r\
eturn %L;\n  }\n\nsub env_file2putenv\n  {\n    my\
 $f=@_[0];\n    my $F=new FileHandle;\n    my $n;\\
n    \n    open ($F, \"$f\");\n    while (<$F>)\n \
     {\n	my $line=$_;\n	my($var, $value)=($_=~/(\\\
S+)\\=(\\S*)/);\n	$ENV{$var}=$value;\n	$ENV_SET{$v\
ar}=1;\n	$n++;\n      }\n    close ($F);\n    retu\
rn $n;\n  }\n\nsub replace_line_in_file\n  {\n    \
my ($file, $wordin, $wordout)=@_;\n    my $O=new F\
ileHandle;\n    my $I=new FileHandle;\n    my $l;\\
n    if (!-e $file){return;}\n    \n    system (\"\
mv $file $file.old\");\n    open ($O, \">$file\");\
\n    open ($I, \"$file.old\");\n    while (<$I>)\\
n      {\n	$l=$_;\n	if (!($l=~/$wordin/)){print $O\
 \"$l\";}\n	elsif ( $wordout ne \"\"){$l=~s/$wordi\
n/$wordout/g;print $O \"$l\";}\n      }\n    close\
 ($O);\n    close ($I);\n    return;\n  }\n\nsub a\
dd_C_libraries\n  {\n   my ($file,$first,@list)=@_\
;\n   \n    my $O=new FileHandle;\n    my $I=new F\
ileHandle;\n    my ($l,$anchor);\n    if (!-e $fil\
e){return;}\n   \n    $anchor=\"#include <$first>\\
";\n	 \n    system (\"mv $file $file.old\");\n    \
open ($O, \">$file\");\n    open ($I, \"$file.old\\
");\n    while (<$I>)\n      {\n	$l=$_;\n	print $O\
 \"$l\";\n	if (!($l=~/$anchor/))\n	   {\n	    \n	 \
   foreach my $lib (@list)\n	       {\n           \
       print $O \"#include <$lib>\\n\";\n	       }\
\n           }\n      }\n    close ($O);\n    clos\
e ($I);\n    return;\n    }\n","use Env;\nuse Cwd;\
\n@suffix=(\"tmp\", \"temp\", \"cache\", \"t_coffe\
e\", \"core\", \"tcoffee\");\n\nif ($#ARGV==-1)\n \
 {\n    print \"clean_cache.pl -file <file to add \
in -dir> -dir=<dir> -size=<value in Mb>\\n0: unlim\
ited -1 always.\\nWill only clean directories matc\
hing:[\";\n    foreach $k(@suffix){print \"*$k* \"\
;}\n    print \"]\\n\";\n    exit (EXIT_FAILURE);\\
n  }\n\n$cl=join (\" \",@ARGV);\nif (($cl=~/\\-no_\
action/))\n  {\n    exit (EXIT_SUCCESS);\n  }\n\ni\
f (($cl=~/\\-debug/))\n  {\n    $DEBUG=1;\n  }\nel\
se\n  {\n    $DEBUG=0;\n  }\n\nif (($cl=~/\\-dir=(\
\\S+)/))\n  {\n    $dir=$1;\n  }\nelse\n  {\n    $\
dir=\"./\";\n  }\n\nif ($cl=~/\\-file=(\\S+)/)\n  \
{\n    $file=$1;\n  }\nelse\n  {\n    $file=0;\n  \
}\n\nif ($cl=~/\\-size=(\\S+)/)\n  {\n    $max_siz\
e=$1;\n  }\nelse\n  {\n    $max_size=0;#unlimited\\
n  }\nif ($cl=~/\\-force/)\n  {\n    $force=1;\n  \
}\nelse\n  {\n    $force=0;\n  }\n\nif ($cl=~/\\-a\
ge=(\\S+)/)\n  {\n    $max_age=$1;\n  }\nelse\n  {\
\n    $max_age=0;#unlimited\n  }\n\n$max_size*=100\
0000;\nif ( ! -d $dir)\n  {\n    print STDERR \"\\\
nCannot process $dir: does not exist \\n\";\n    e\
xit (EXIT_FAILURE);\n  }\n\nif ( !($dir=~/^\\//))\\
n  {\n    $base=cwd();\n    $dir=\"$base/$dir\";\n\
  }\n\n$proceed=0;\nforeach $s (@suffix)\n  {\n   \
 \n    if (($dir=~/$s/)){$proceed=1;}\n    $s=uc (\
$s);\n    if (($dir=~/$s/)){$proceed=1;}\n  }\nif \
( $proceed==0)\n  {\n    print STDERR \"Clean_cach\
e.pl can only clean directories whose absolute pat\
h name contains the following strings:\";\n    for\
each $w (@suffix) {print STDERR \"$w \";$w=lc($w);\
 print STDERR \"$w \";}\n    print STDERR \"\\nCan\
not process $dir\\n\";\n    exit (EXIT_FAILURE);\n\
  }\n\n$name_file=\"$dir/name_file.txt\";\n$size_f\
ile=\"$dir/size_file.txt\";\nif ( $force){&create_\
ref_file ($dir,$name_file,$size_file);}\nif ($file\
){&add_file ($dir, $name_file, $size_file, $file);\
}\n&clean_dir ($dir, $name_file, $size_file, $max_\
size,$max_age);\nexit (EXIT_SUCCESS);\n\nsub clean\
_dir \n  {\n    my ($dir, $name_file, $size_file, \
$max_size, $max_age)=@_;\n    my ($tot_size, $size\
, $f, $s);\n\n  \n    $tot_size=&get_tot_size ($di\
r, $name_file, $size_file);\n\n    if ( $tot_size<\
=$max_size){return ;}\n    else {$max_size/=2;}\n \
   \n    #recreate the name file in case some temp\
rary files have not been properly registered\n    \
&create_ref_file ($dir, $name_file, $size_file, $m\
ax_age);\n  \n    $new_name_file=&vtmpnam();\n    \
open (R, \"$name_file\");\n    open (W, \">$new_na\
me_file\");\n    while (<R>)\n      {\n	my $line=$\
_;\n	\n	($f, $s)=($line=~/(\\S+) (\\S+)/);\n	if ( \
!($f=~/\\S/)){next;}\n	\n	elsif ($max_size && $tot\
_size>=$max_size && !($f=~/name_file/))\n	  {\n	  \
  remove ( \"$dir/$f\");\n	    $tot_size-=$s;\n	  \
}\n	elsif ( $max_age && -M(\"$dir/$f\")>=$max_age)\
\n	  {\n	    remove ( \"$dir/$f\");\n	    $tot_siz\
e-=$s;\n	  }\n	else\n	  {\n	    print W \"$f $s\\n\
\";\n	  }\n      }\n    close (R);\n    close (W);\
\n    open (F, \">$size_file\");\n    print F \"$t\
ot_size\";\n    if ( -e $new_name_file){`mv $new_n\
ame_file $name_file`;}\n    close (F);\n  }\nsub g\
et_tot_size\n  {\n    my ($dir, $name_file, $size_\
file)=@_;\n    my $size;\n    \n    if ( !-d $dir)\
{return 0;}\n    if ( !-e $name_file)\n      {\n	\\
n	&create_ref_file ($dir, $name_file, $size_file);\
\n      }\n    open (F, \"$size_file\");\n    $siz\
e=<F>;\n    close (F);\n    chomp ($size);\n    re\
turn $size;\n  }\nsub size \n  {\n    my $f=@_[0];\
\n\n    if ( !-d $f){return -s($f);}\n    else {re\
turn &dir2size($f);}\n  }\nsub dir2size\n  {\n    \
my $d=@_[0];\n    my ($s, $f);\n    \n    if ( !-d\
 $d) {return 0;}\n    \n    foreach $f (&dir2list \
($d))\n      {\n	if ( -d $f){$s+=&dir2size (\"$d/$\
f\");}\n	else {$s+= -s \"$dir/$f\";}\n      }\n   \
 return $s;\n  }\n\nsub remove \n  {\n    my $file\
=@_[0];\n    my ($f);\n    \n    debug_print( \"--\
- $file ---\\n\");\n    if (($file eq \".\") || ($\
file eq \"..\") || ($file=~/\\*/)){return EXIT_FAI\
LURE;}\n    elsif ( !-d $file)\n      {\n	debug_pr\
int (\"unlink $file\\n\");\n	if (-e $file){unlink \
($file);}\n      }\n    elsif ( -d $file)\n      {\
\n	debug_print (\"++++++++ $file +++++++\\n\");\n	\
foreach $f (&dir2list($file))\n	  {\n	    &remove \
(\"$file/$f\");\n	  }\n	debug_print (\"rmdir $file\
\\n\");\n	rmdir $file;\n      }\n    else\n      {\
\n	debug_print (\"????????? $file ????????\\n\");\\
n      }\n    return EXIT_SUCCESS;\n  }\n\nsub dir\
2list\n  {\n    my $dir=@_[0];\n    my (@list1, @l\
ist2,@list3, $l);\n\n    opendir (DIR,$dir);\n    \
@list1=readdir (DIR);\n    closedir (DIR);\n    \n\
    foreach $l (@list1)\n      {\n	if ( $l ne \".\\
" && $l ne \"..\"){@list2=(@list2, $l);}\n      }\\
n    @list3 = sort { (-M \"$dir/$list2[$b]\") <=> \
(-M \"$dir/$list2[$a]\")} @list2;\n    return @lis\
t3;\n    \n  }\n\nsub debug_print\n  {\n    \n    \
if ($DEBUG==1){print @_;}\n    \n  }\nsub create_r\
ef_file\n  {\n    my ($dir,$name_file,$size_file)=\
@_;\n    my ($f, $s, $tot_size, @l);\n    \n    if\
 ( !-d $dir){return;}\n    \n    @l=&dir2list ($di\
r);\n    open (F, \">$name_file\");\n    foreach $\
f (@l)\n      {\n	$s=&size(\"$dir/$f\");\n	$tot_si\
ze+=$s;\n	print F \"$f $s\\n\";\n      }\n    &mye\
cho ($tot_size, \">$size_file\");\n    close (F);\\
n  }\nsub add_file \n  {\n    my ($dir,$name_file,\
$size_file,$file)=@_;\n    my ($s, $tot_size);\n  \
  \n    if ( !-d $dir)   {return;}\n    if ( !-e \\
"$dir/$file\" ) {return;}\n    if ( !-e $name_file\
){&create_ref_file ($dir,$name_file,$size_file);}\\
n					    \n    $s=&size(\"$dir/$file\");\n    ope\
n (F, \">>$name_file\");\n    print F \"$file\\n\"\
;\n    close (F);\n\n    $tot_size=&get_tot_size (\
$dir,$name_file,$size_file);\n    $tot_size+=$s;\n\
    &myecho ($tot_size, \">$size_file\");\n    \n \
 }\n	\nsub myecho\n  {\n    my ($string, $file)=@_\
;\n    open (ECHO, $file) || die;\n    print ECHO \
\"$string\";\n    close (ECHO);\n  }\n    \n		\n	\\
nsub vtmpnam\n  {\n    my $tmp_file_name;\n    $tm\
p_name_counter++;\n    $tmp_file_name=\"tmp_file_f\
or_clean_cache_pdb$$.$tmp_name_counter\";\n    $tm\
p_file_list[$ntmp_file++]=$tmp_file_name;\n    if \
( -e $tmp_file_name) {return &vtmpnam ();}\n    el\
se {return $tmp_file_name;}\n  }\n","\nmy $address\
=\"http://www.tcoffee.org/Data/Datasets/NatureProt\
ocolsDataset.tar.gz\";\nmy $out=\"NatureProtocolsD\
ataset.tar.gz\";\n&url2file ($address,$out);\n\nif\
 ( -e $out)\n  {\n    \n    system (\"gunzip Natur\
eProtocolsDataset.tar.gz\");\n    system (\"tar -x\
vf NatureProtocolsDataset.tar\");\n  	system (\"rm\
 -rf NatureProtocolsDataset.tar\");  \n    print \\
"Your Data Set is in the Folder 'NatureProtocolsDa\
taset'\\n\";\n  }\nelse \n  {\n    print \"Could n\
ot Download Dataset --- Web site may be down -- Tr\
y again later\\n\";\n  }\n\n\n\n\nsub url2file\n{\\
n    my ($address, $out, $wget_arg, $curl_arg)=(@_\
);\n    my ($pg, $flag, $r, $arg, $count);\n    \n\
    if (!$CONFIGURATION){&check_configuration (\"w\
get\", \"INTERNET\", \"gzip\");$CONFIGURATION=1;}\\
n    \n    if (&pg_is_installed (\"wget\"))   {$pg\
=\"wget\"; $flag=\"-O\";$arg=$wget_arg;}\n    elsi\
f (&pg_is_installed (\"curl\")){$pg=\"curl\"; $fla\
g=\"-o\";$arg=$curl_arg;}\n    return system (\"$p\
g $address $flag $out>/dev/null 2>/dev/null\");\n\\
n}\n\nsub pg_is_installed\n  {\n    my @ml=@_;\n  \
  my $r, $p, $m;\n    my $supported=0;\n    \n    \
my $p=shift (@ml);\n    if ($p=~/::/)\n      {\n	i\
f (system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){ret\
urn 1;}\n	else {return 0;}\n      }\n    else\n   \
   {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"\
){return 0;}\n	else {return 1;}\n      }\n  }\nsub\
 check_configuration \n    {\n      my @l=@_;\n   \
   my $v;\n      foreach my $p (@l)\n	{\n	  \n	  i\
f   ( $p eq \"EMAIL\")\n	    { \n	      if ( !($EM\
AIL=~/@/))\n		{\n		  exit (EXIT_FAILURE);\n		}\n	 \
   }\n	  elsif( $p eq \"INTERNET\")\n	    {\n	    \
  if ( !&check_internet_connection())\n		{\n		  ex\
it (EXIT_FAILURE);\n		}\n	    }\n	  elsif( $p eq \\
"wget\")\n	    {\n	      if (!&pg_is_installed (\"\
wget\") && !&pg_is_installed (\"curl\"))\n		{\n		 \
 exit (EXIT_FAILURE);\n		}\n	    }\n	  elsif( !(&p\
g_is_installed ($p)))\n	    {\n	      exit (EXIT_F\
AILURE);\n	    }\n	}\n      return 1;\n    }\nsub \
check_internet_connection\n  {\n    my $internet;\\
n    my $tmp;\n    &check_configuration ( \"wget\"\
); \n    \n    $tmp=&vtmpnam ();\n    \n    if    \
 (&pg_is_installed    (\"wget\")){`wget www.google\
.com -O$tmp >/dev/null 2>/dev/null`;}\n    elsif  \
(&pg_is_installed    (\"curl\")){`curl www.google.\
com -o$tmp >/dev/null 2>/dev/null`;}\n    \n    if\
 ( !-e $tmp || -s $tmp < 10){$internet=0;}\n    el\
se {$internet=1;}\n    if (-e $tmp){unlink $tmp;}\\
n\n    return $internet;\n  }\n\nsub vtmpnam\n    \
  {\n	my $r=rand(100000);\n	my $f=\"file.$r.$$\";\\
n	while (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\n	\
push (@TMPFILE_LIST, $f);\n	return $f;\n      }\n\\
n","\n$t_coffee=\"t_coffee\";\n\nforeach $value ( \
@ARGV)\n  {\n    $seq_file=$seq_file.\" \".$value;\
\n  }\n\n$name=$ARGV[0];\n$name=~s/\\.[^\\.]*$//;\\
n$lib_name=\"$name.mocca_lib\";\n$type=`t_coffee $\
seq_file -get_type -quiet`;\nchop ($type);\n\nif (\
 $type eq \"PROTEIN\"){$lib_mode=\"lalign_rs_s_pai\
r -lalign_n_top 20\";}\nelsif ( $type eq\"DNA\"){$\
lib_mode=\"lalign_rs_s_dna_pair -lalign_n_top 40\"\
;}\n\nif ( !(-e $lib_name))\n  {\n	  \n  $command=\
\"$t_coffee -mocca -seq_weight=no -cosmetic_penalt\
y=0 -mocca_interactive -in $lib_mode -out_lib $lib\
_name -infile $seq_file\";\n  \n  }\nelsif ( (-e $\
lib_name))\n  {\n  $command=\"$t_coffee -mocca -se\
q_weight=no -cosmetic_penalty=0 -mocca_interactive\
 -in $lib_name -infile $seq_file\";\n  \n  }\n\nsy\
stem ($command);\n\nexit;\n\n","my $WSDL = 'http:/\
/www.ebi.ac.uk/Tools/webservices/wsdl/WSDaliLite.w\
sdl';\n\nuse SOAP::Lite;\nuse Data::Dumper;\nuse G\
etopt::Long qw(:config no_ignore_case bundling);\n\
use File::Basename;\n\nmy $checkInterval = 5;\n\nm\
y %params=(\n	    'async' => '1', # Use async mode\
 and simulate sync mode in client\n	    );\nGetOpt\
ions(\n    'pdb1=s'     => \\$params{'sequence1'},\
\n    'chainid1=s' => \\$params{'chainid1'},\n    \
'pdb2=s'     => \\$params{'sequence2'},\n    'chai\
nid2=s' => \\$params{'chainid2'},\n    \"help|h\"	\
 => \\$help, # Usage info\n    \"async|a\"	 => \\$\
async, # Asynchronous submission\n    \"polljob\"	\
 => \\$polljob, # Get results\n    \"status\"	 => \
\\$status, # Get status\n    \"jobid|j=s\"  => \\$\
jobid, # JobId\n    \"email|S=s\"  => \\$params{em\
ail}, # E-mail address\n    \"trace\"      => \\$t\
race, # SOAP messages\n    \"sequence=s\" => \\$se\
quence, # Input PDB\n    );\n\nmy $scriptName = ba\
sename($0, ());\nif($help) {\n    &usage();\n    e\
xit(0);\n}\n\nif($trace) {\n    print \"Tracing ac\
tive\\n\";\n    SOAP::Lite->import(+trace => 'debu\
g');\n}\n\nmy $soap = SOAP::Lite\n    ->service($W\
SDL)\n    ->on_fault(sub {\n        my $soap = shi\
ft;\n        my $res = shift;\n        # Throw an \
exception for all faults\n        if(ref($res) eq \
'') {\n            die($res);\n        } else {\n \
           die($res->faultstring);\n        }\n   \
     return new SOAP::SOM;\n    }\n               \
);\n\nif( !($polljob || $status) &&\n    !( define\
d($params{'sequence1'}) && defined($params{'sequen\
ce2'}) )\n    ) {\n    print STDERR 'Error: bad op\
tion combination', \"\\n\";\n    &usage();\n    ex\
it(1);\n}\nelsif($polljob && defined($jobid)) {\n \
   print \"Getting results for job $jobid\\n\";\n \
   getResults($jobid);\n}\nelsif($status && define\
d($jobid)) {\n    print STDERR \"Getting status fo\
r job $jobid\\n\";\n    my $result = $soap->checkS\
tatus($jobid);\n    print STDOUT \"$result\", \"\\\
n\";\n    if($result eq 'DONE') {\n	print STDERR \\
"To get results: $scriptName --polljob --jobid $jo\
bid\\n\";\n    }\n}\nelse {\n    if(-f $params{'se\
quence1'}) {\n	$params{'sequence1'} = read_file($p\
arams{'sequence1'});\n    }\n    if(-f $params{'se\
quence2'}) {\n	$params{'sequence2'} = read_file($p\
arams{'sequence2'});\n    }\n\n    my $jobid;\n   \
 my $paramsData = SOAP::Data->name('params')->type\
(map=>\\%params);\n    # For SOAP::Lite 0.60 and e\
arlier parameters are passed directly\n    if($SOA\
P::Lite::VERSION eq '0.60' || $SOAP::Lite::VERSION\
 =~ /0\\.[1-5]/) {\n        $jobid = $soap->runDal\
iLite($paramsData);\n    }\n    # For SOAP::Lite 0\
.69 and later parameter handling is different, so \
pass\n    # undef's for templated params, and then\
 pass the formatted args.\n    else {\n        $jo\
bid = $soap->runDaliLite(undef,\n				     $paramsD\
ata);\n    }\n\n    if (defined($async)) {\n	print\
 STDOUT $jobid, \"\\n\";\n        print STDERR \"T\
o check status: $scriptName --status --jobid $jobi\
d\\n\";\n    } else { # Synchronous mode\n        \
print STDERR \"JobId: $jobid\\n\";\n        sleep \
1;\n        getResults($jobid);\n    }\n}\n\nsub c\
lientPoll($) {\n    my $jobid = shift;\n    my $re\
sult = 'PENDING';\n    # Check status and wait if \
not finished\n    #print STDERR \"Checking status:\
 $jobid\\n\";\n    while($result eq 'RUNNING' || $\
result eq 'PENDING') {\n        $result = $soap->c\
heckStatus($jobid);\n        print STDERR \"$resul\
t\\n\";\n        if($result eq 'RUNNING' || $resul\
t eq 'PENDING') {\n            # Wait before polli\
ng again.\n            sleep $checkInterval;\n    \
    }\n    }\n}\n\nsub getResults($) {\n    $jobid\
 = shift;\n    # Check status, and wait if not fin\
ished\n    clientPoll($jobid);\n    # Use JobId if\
 output file name is not defined\n    unless(defin\
ed($outfile)) {\n        $outfile=$jobid;\n    }\n\
    # Get list of data types\n    my $resultTypes \
= $soap->getResults($jobid);\n    # Get the data a\
nd write it to a file\n    if(defined($outformat))\
 { # Specified data type\n        my $selResultTyp\
e;\n        foreach my $resultType (@$resultTypes)\
 {\n            if($resultType->{type} eq $outform\
at) {\n                $selResultType = $resultTyp\
e;\n            }\n        }\n        $res=$soap->\
poll($jobid, $selResultType->{type});\n        wri\
te_file($outfile.'.'.$selResultType->{ext}, $res);\
\n    } else { # Data types available\n        # W\
rite a file for each output type\n        for my $\
resultType (@$resultTypes){\n            #print \"\
Getting $resultType->{type}\\n\";\n            $re\
s=$soap->poll($jobid, $resultType->{type});\n     \
       write_file($outfile.'.'.$resultType->{ext},\
 $res);\n        }\n    }\n}\n\nsub read_file($) {\
\n    my $filename = shift;\n    open(FILE, $filen\
ame);\n    my $content;\n    my $buffer;\n    whil\
e(sysread(FILE, $buffer, 1024)) {\n	$content.= $bu\
ffer;\n    }\n    close(FILE);\n    return $conten\
t;\n}\n\nsub write_file($$) {\n    my ($tmp,$entit\
y) = @_;\n    print STDERR \"Creating result file:\
 \".$tmp.\"\\n\";\n    unless(open (FILE, \">$tmp\\
")) {\n	return 0;\n    }\n    syswrite(FILE, $enti\
ty);\n    close (FILE);\n    return 1;\n}\n\nsub u\
sage {\n    print STDERR <<EOF\nDaliLite\n========\
\n\nPairwise comparison of protein structures\n\n[\
Required]\n\n  --pdb1                : str  : PDB \
ID for structure 1\n  --pdb2                : str \
 : PDB ID for structure 2\n\n[Optional]\n\n  --cha\
in1              : str  : Chain identifer in struc\
ture 1\n  --chain2              : str  : Chain ide\
ntifer in structure 2\n\n[General]\n\n  -h, --help\
            :      : prints this help text\n  -S, \
--email           : str  : user email address\n  -\
a, --async           :      : asynchronous submiss\
ion\n      --status          :      : poll for the\
 status of a job\n      --polljob         :      :\
 poll for the results of a job\n  -j, --jobid     \
      : str  : jobid for an asynchronous job\n  -O\
, --outfile         : str  : file name for results\
 (default is jobid)\n      --trace	        :      \
: show SOAP messages being interchanged \n\nSynchr\
onous job:\n\n  The results/errors are returned as\
 soon as the job is finished.\n  Usage: $scriptNam\
e --email <your\\@email> [options] pdbFile [--outf\
ile string]\n  Returns: saves the results to disk\\
n\nAsynchronous job:\n\n  Use this if you want to \
retrieve the results at a later time. The results \
\n  are stored for up to 24 hours. \n  The asynchr\
onous submission mode is recommended when users ar\
e submitting \n  batch jobs or large database sear\
ches	\n  Usage: $scriptName --email <your\\@email>\
 --async [options] pdbFile\n  Returns: jobid\n\n  \
Use the jobid to query for the status of the job. \
\n  Usage: $scriptName --status --jobid <jobId>\n \
 Returns: string indicating the status of the job:\
\n    DONE - job has finished\n    RUNNING - job i\
s running\n    NOT_FOUND - job cannot be found\n  \
  ERROR - the jobs has encountered an error\n\n  W\
hen done, use the jobid to retrieve the status of \
the job. \n  Usage: $scriptName --polljob --jobid \
<jobId> [--outfile string]\n\n[Help]\n\n  For more\
 detailed help information refer to\n  http://www.\
ebi.ac.uk/DaliLite/\nEOF\n;\n}\n","my $WSDL = 'htt\
p://www.ebi.ac.uk/Tools/webservices/wsdl/WSWUBlast\
.wsdl';\n\nuse strict;\nuse SOAP::Lite;\nuse Getop\
t::Long qw(:config no_ignore_case bundling);\nuse \
File::Basename;\n\nmy $checkInterval = 15;\n\nmy $\
numOpts = scalar(@ARGV);\nmy ($outfile, $outformat\
, $help, $async, $polljob, $status, $ids, $jobid, \
$trace, $sequence);\nmy %params= ( # Defaults\n	  \
    'async' => 1, # Force into async mode\n	      \
'exp' => 10.0, # E-value threshold\n	      'numal'\
 => 50, # Maximum number of alignments\n	      'sc\
ores' => 100, # Maximum number of scores\n        \
    );\nGetOptions( # Map the options into variabl\
es\n    \"program|p=s\"     => \\$params{program},\
 # BLAST program\n    \"database|D=s\"    => \\$pa\
rams{database}, # Search database\n    \"matrix|m=\
s\"      => \\$params{matrix}, # Scoring matrix\n \
   \"exp|E=f\"         => \\$params{exp}, # E-valu\
e threshold\n    \"echofilter|e\"    => \\$params{\
echofilter}, # Display filtered sequence\n    \"fi\
lter|f=s\"      => \\$params{filter}, # Low comple\
xity filter name\n    \"alignments|b=i\"  => \\$pa\
rams{numal}, # Number of alignments\n    \"scores|\
s=i\"      => \\$params{scores}, # Number of score\
s\n    \"sensitivity|S=s\" => \\$params{sensitivit\
y}, # Search sensitivity\n    \"sort|t=s\"	      =\
> \\$params{sort}, # Sort hits by...\n    \"stats|\
T=s\"       => \\$params{stats}, # Scoring statist\
ic to use\n    \"strand|d=s\"      => \\$params{st\
rand}, # Strand to use in DNA vs. DNA search\n    \
\"topcombon|c=i\"   => \\$params{topcombon}, # Con\
sistent sets of HSPs\n    \"outfile=s\"       => \\
\$outfile, # Output file\n    \"outformat|o=s\"   \
=> \\$outformat, # Output format\n    \"help|h\"	 \
     => \\$help, # Usage info\n    \"async|a\"	   \
   => \\$async, # Asynchronous mode\n    \"polljob\
\"	      => \\$polljob, # Get results\n    \"statu\
s\"	      => \\$status, # Get job status\n    \"id\
s\"             => \\$ids, # Get ids from result\n\
    \"jobid|j=s\"       => \\$jobid, # JobId\n    \
\"email=s\"         => \\$params{email}, # E-mail \
address\n    \"trace\"           => \\$trace, # SO\
AP trace\n    \"sequence=s\"      => \\$sequence, \
# Query sequence\n    );\n\nmy $scriptName = basen\
ame($0, ());\nif($help || $numOpts == 0) {\n    &u\
sage();\n    exit(0);\n}\n\nif($trace){\n    print\
 STDERR \"Tracing active\\n\";\n    SOAP::Lite->im\
port(+trace => 'debug');\n}\n\nmy $soap = SOAP::Li\
te\n    ->service($WSDL)\n    ->proxy('http://loca\
lhost/',\n    #proxy => ['http' => 'http://your.pr\
oxy.server/'], # HTTP proxy\n    timeout => 600, #\
 HTTP connection timeout\n    )\n    ->on_fault(su\
b { # SOAP fault handler\n        my $soap = shift\
;\n        my $res = shift;\n        # Throw an ex\
ception for all faults\n        if(ref($res) eq ''\
) {\n            die($res);\n        } else {\n   \
         die($res->faultstring);\n        }\n     \
   return new SOAP::SOM;\n    }\n               );\
\n\nif( !($polljob || $status || $ids) &&\n    !( \
defined($ARGV[0]) || defined($sequence) )\n    ) {\
\n    print STDERR 'Error: bad option combination'\
, \"\\n\";\n    &usage();\n    exit(1);\n}\nelsif(\
$polljob && defined($jobid)) {\n    print \"Gettin\
g results for job $jobid\\n\";\n    getResults($jo\
bid);\n}\nelsif($status && defined($jobid)) {\n   \
 print STDERR \"Getting status for job $jobid\\n\"\
;\n    my $result = $soap->checkStatus($jobid);\n \
   print STDOUT \"$result\\n\";\n    if($result eq\
 'DONE') {\n	print STDERR \"To get results: $scrip\
tName --polljob --jobid $jobid\\n\";\n    }\n}  \n\
elsif($ids && defined($jobid)) {\n    print STDERR\
 \"Getting ids from job $jobid\\n\";\n    getIds($\
jobid);\n}\nelse {\n    # Prepare input data\n    \
my $content;\n    my (@contents) = ();\n    if(-f \
$ARGV[0] || $ARGV[0] eq '-') {	\n	$content={type=>\
'sequence',content=>read_file($ARGV[0])};	\n    }\\
n    if($sequence) {	\n	if(-f $sequence || $sequen\
ce eq '-') {	\n	    $content={type=>'sequence',con\
tent=>read_file($ARGV[0])};	\n	} else {\n	    $con\
tent={type=>'sequence',content=>$sequence};\n	}\n \
   }\n    push @contents, $content;\n\n    # Submi\
t the job\n    my $paramsData = SOAP::Data->name('\
params')->type(map=>\\%params);\n    my $contentDa\
ta = SOAP::Data->name('content')->value(\\@content\
s);\n    # For SOAP::Lite 0.60 and earlier paramet\
ers are passed directly\n    if($SOAP::Lite::VERSI\
ON eq '0.60' || $SOAP::Lite::VERSION =~ /0\\.[1-5]\
/) {\n        $jobid = $soap->runWUBlast($paramsDa\
ta, $contentData);\n    }\n    # For SOAP::Lite 0.\
69 and later parameter handling is different, so p\
ass\n    # undef's for templated params, and then \
pass the formatted args.\n    else {\n        $job\
id = $soap->runWUBlast(undef, undef,\n				   $para\
msData, $contentData);\n    }\n\n    # Asynchronou\
s mode: output jobid and exit.\n    if (defined($a\
sync)) {\n	print STDOUT $jobid, \"\\n\";\n        \
print STDERR \"To check status: $scriptName --stat\
us --jobid $jobid\\n\";\n    }\n    # Synchronous \
mode: try to get results\n    else {\n        prin\
t STDERR \"JobId: $jobid\\n\";\n        sleep 1;\n\
        getResults($jobid);\n    }\n}\n\nsub getId\
s($) {\n    my $jobid = shift;\n    my $results = \
$soap->getIds($jobid);\n    for my $result (@$resu\
lts){\n	print \"$result\\n\";\n    }\n}\n\nsub cli\
entPoll($) {\n    my $jobid = shift;\n    my $resu\
lt = 'PENDING';\n    # Check status and wait if no\
t finished\n    while($result eq 'RUNNING' || $res\
ult eq 'PENDING') {\n        $result = $soap->chec\
kStatus($jobid);\n        print STDERR \"$result\\\
n\";\n        if($result eq 'RUNNING' || $result e\
q 'PENDING') {\n            # Wait before polling \
again.\n            sleep $checkInterval;\n       \
 }\n    }\n}\n\nsub getResults($) {\n    my $jobid\
 = shift;\n    my $res;\n    # Check status, and w\
ait if not finished\n    clientPoll($jobid);\n    \
# Use JobId if output file name is not defined\n  \
  unless(defined($outfile)) {\n        $outfile=$j\
obid;\n    }\n    # Get list of data types\n    my\
 $resultTypes = $soap->getResults($jobid);\n    # \
Get the data and write it to a file\n    if(define\
d($outformat)) { # Specified data type\n	if($outfo\
rmat eq 'xml') {$outformat = 'toolxml';}\n	if($out\
format eq 'txt') {$outformat = 'tooloutput';}\n   \
     my $selResultType;\n        foreach my $resul\
tType (@$resultTypes) {\n            if($resultTyp\
e->{type} eq $outformat) {\n                $selRe\
sultType = $resultType;\n            }\n        }\\
n        $res=$soap->poll($jobid, $selResultType->\
{type});\n	if($outfile eq '-') {\n	     write_file\
($outfile, $res);\n	} else {\n	    write_file($out\
file.'.'.$selResultType->{ext}, $res);\n	}\n    } \
else { # Data types available\n        # Write a f\
ile for each output type\n        for my $resultTy\
pe (@$resultTypes){\n            #print STDERR \"G\
etting $resultType->{type}\\n\";\n            $res\
=$soap->poll($jobid, $resultType->{type});\n	    i\
f($outfile eq '-') {\n		write_file($outfile, $res)\
;\n	    } else {\n		write_file($outfile.'.'.$resul\
tType->{ext}, $res);\n	    }\n        }\n    }\n}\\
n\nsub read_file($) {\n    my $filename = shift;\n\
    my ($content, $buffer);\n    if($filename eq '\
-') {\n	while(sysread(STDIN, $buffer, 1024)) {\n	 \
   $content .= $buffer;\n	}\n    }\n    else { # F\
ile\n	open(FILE, $filename) or die \"Error: unable\
 to open input file\";\n	while(sysread(FILE, $buff\
er, 1024)) {\n	    $content .= $buffer;\n	}\n	clos\
e(FILE);\n    }\n    return $content;\n}\n\nsub wr\
ite_file($$) {\n    my ($filename, $data) = @_;\n \
   print STDERR 'Creating result file: ' . $filena\
me . \"\\n\";\n    if($filename eq '-') {\n	print \
STDOUT $data;\n    }\n    else {\n	open(FILE, \">$\
filename\") or die \"Error: unable to open output \
file\";\n	syswrite(FILE, $data);\n	close(FILE);\n \
   }\n}\n\nsub usage {\n    print STDERR <<EOF\nWU\
-BLAST\n========\n\nRapid sequence database search\
 programs utilizing the BLAST algorithm.\n   \n[Re\
quired]\n\n      --email       : str  : user email\
 address \n  -p, --program	    : str  : BLAST prog\
ram to use: blastn, blastp, blastx, \n            \
                 tblastn or tblastx\n  -D, --datab\
ase    : str  : database to search\n  seqFile     \
      : file : query sequence data file (\"-\" for\
 STDIN)\n\n[Optional]\n\n  -m, --matrix	    : str \
 : scoring matrix\n  -E, --exp	    : real : 0<E<= \
1000. Statistical significance threshold\n        \
                     for reporting database sequen\
ce matches.\n  -e, --echofilter  :      : display \
the filtered query sequence in the output\n  -f, -\
-filter	    : str  : activates filtering of the qu\
ery sequence\n  -b, --alignments  : int  : number \
of alignments to be reported\n  -s, --scores	    :\
 int  : number of scores to be reported\n  -S, --s\
ensitivity : str  :\n  -t, --sort	    : str  :\n  \
-T, --stats       : str  :\n  -d, --strand      : \
str  : DNA strand to search with in DNA vs. DNA se\
arches \n  -c, --topcombon   :      :\n\n[General]\
	\n\n  -h, --help       :      : prints this help \
text\n  -a, --async      :      : forces to make a\
n asynchronous query\n      --status     :      : \
poll for the status of a job\n      --polljob    :\
      : poll for the results of a job\n  -j, --job\
id      : str  : jobid that was returned when an a\
synchronous job \n                            was \
submitted.\n  -O, --outfile    : str  : name of th\
e file results should be written to \n            \
                (default is based on the jobid; \"\
-\" for STDOUT)\n  -o, --outformat  : str  : txt o\
r xml output (no file is written)\n      --trace	 \
  :      : show SOAP messages being interchanged \\
n\nSynchronous job:\n\n  The results/errors are re\
turned as soon as the job is finished.\n  Usage: $\
scriptName --email <your\\@email> [options...] seq\
File\n  Returns: saves the results to disk\n\nAsyn\
chronous job:\n\n  Use this if you want to retriev\
e the results at a later time. The results \n  are\
 stored for up to 24 hours. \n  The asynchronous s\
ubmission mode is recommended when users are submi\
tting \n  batch jobs or large database searches	\n\
  Usage: $scriptName --async --email <your\\@email\
> [options...] seqFile\n  Returns : jobid\n\n  Use\
 the jobid to query for the status of the job. \n \
 Usage: $scriptName --status --jobid <jobId>\n  Re\
turns : string indicating the status of the job:\n\
    DONE - job has finished\n    RUNNING - job is \
running\n    NOT_FOUND - job cannot be found\n    \
ERROR - the jobs has encountered an error\n\n  Whe\
n done, use the jobid to retrieve the status of th\
e job. \n  Usage: $scriptName --polljob --jobid <j\
obId> [--outfile string]\n  Returns: saves the res\
ults to disk\n\n[Help]\n\nFor more detailed help i\
nformation refer to \nhttp://www.ebi.ac.uk/blast2/\
WU-Blast2_Help_frame.html\n \nEOF\n;\n}\n","\nmy $\
WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsd\
l/WSBlastpgp.wsdl';\n\nuse SOAP::Lite;\nuse Getopt\
::Long qw(:config no_ignore_case bundling);\nuse F\
ile::Basename;\n\nmy $checkInterval = 15;\n\nmy %p\
arams=(\n	    'async' => '1', # Use async mode and\
 simulate sync mode in client\n	    );\nGetOptions\
(\n    \"mode=s\"           => \\$params{mode}, # \
Search mode: PSI-Blast or PHI-Blast\n    \"databas\
e|d=s\"     => \\$params{database}, # Database to \
search\n    \"matrix|M=s\"       => \\$params{matr\
ix},# Scoring maxtrix\n    \"exp|e=f\"          =>\
 \\$params{exp}, # E-value\n    \"expmulti|h=f\"  \
   => \\$params{expmulti}, # E-value\n    \"filter\
|F=s\"       => \\$params{filter}, # Low complexit\
y filter\n    \"dropoff|X=i\"      => \\$params{dr\
opoff}, # Dropoff score\n    \"finaldropoff|Z=i\" \
=> \\$params{finaldropoff}, # Final dropoff score\\
n    \"scores|v=i\"       => \\$params{scores}, # \
Max number of scores\n    \"align=i\"          => \
\\$params{align}, # Alignment view\n    \"startreg\
ion|S=i\"  => \\$params{startregion}, # Start of r\
egion in query\n    \"endregion|H=i\"    => \\$par\
ams{endregion}, # End of region in query\n    \"ma\
xpasses|j=i\"    => \\$params{maxpasses}, # Number\
 of PSI iterations\n    \"opengap|G=i\"      => \\\
$params{opengap}, # Gap open penalty\n    \"extend\
gap|E=i\"    => \\$params{extendgap}, # Gap extens\
ion penalty\n    \"pattern=s\"        => \\$params\
{pattern}, # PHI-BLAST pattern\n    \"usagemode|p=\
s\"    => \\$params{usagemode}, # PHI-BLAST progra\
m\n    \"appxml=s\"         => \\$params{appxml}, \
# Application XML\n    \"sequence=s\"       => \\$\
sequence, # Query sequence\n    \"help\"	       =>\
 \\$help, # Usage info\n    \"polljob\"	       => \
\\$polljob, # Get results\n    \"status\"	       =\
> \\$status, # Get status\n    \"ids\"      	     \
  => \\$ids, # Get ids from result\n    \"jobid=s\\
"          => \\$jobid, # JobId\n    \"outfile=s\"\
        => \\$outfile, # Output filename\n    \"ou\
tformat|o=s\"    => \\$outformat, # Output file fo\
rmat\n    \"async|a\"	       => \\$async, # Async \
submission\n    \"email=s\"          => \\$params{\
email}, # User e-mail address\n    \"trace\"      \
      => \\$trace, # Show SOAP messages\n    );\n\\
nmy $scriptName = basename($0, ());\nif($help) {\n\
    &usage();\n    exit(0);\n}\n\nif ($trace){\n  \
  print \"Tracing active\\n\";\n    SOAP::Lite->im\
port(+trace => 'debug');\n}\n\nmy $soap = SOAP::Li\
te\n    ->service($WSDL)\n    ->on_fault(sub {\n  \
      my $soap = shift;\n        my $res = shift;\\
n        # Throw an exception for all faults\n    \
    if(ref($res) eq '') {\n            die($res);\\
n        } else {\n            die($res->faultstri\
ng);\n        }\n        return new SOAP::SOM;\n  \
  }\n               );\n\nif( !($polljob || $statu\
s || $ids) &&\n    !( (defined($ARGV[0]) && -f $AR\
GV[0]) || defined($sequence) )\n    ) {\n    print\
 STDERR 'Error: bad option combination', \"\\n\";\\
n    &usage();\n    exit(1);\n}\nelsif($polljob &&\
 defined($jobid)) {\n    print \"Getting results f\
or job $jobid\\n\";\n    getResults($jobid);\n}\ne\
lsif($status && defined($jobid)) {\n    print STDE\
RR \"Getting status for job $jobid\\n\";\n    my $\
result = $soap->checkStatus($jobid);\n    print ST\
DOUT $result, \"\\n\";\n    if($result eq 'DONE') \
{\n	print STDERR \"To get results: $scriptName --p\
olljob --jobid $jobid\\n\";\n    }\n}  \nelsif($id\
s && defined($jobid)) {\n    print STDERR \"Gettin\
g ids from job $jobid\\n\";\n    getIds($jobid);\n\
}\nelse {\n    if(-f $ARGV[0]) {	\n	$content={type\
=>'sequence', content=>read_file($ARGV[0])};	\n   \
 }\n    if($sequence) {	\n	if(-f $sequence) {\n	  \
  $content={type=>'sequence', content=>read_file($\
sequence)};	\n	} else {\n	    $content={type=>'seq\
uence', content=>$sequence};\n	}\n    }\n    push \
@content, $content;\n\n    my $jobid;\n    my $par\
amsData = SOAP::Data->name('params')->type(map=>\\\
%params);\n    my $contentData = SOAP::Data->name(\
'content')->value(\\@content);\n    # For SOAP::Li\
te 0.60 and earlier parameters are passed directly\
\n    if($SOAP::Lite::VERSION eq '0.60' || $SOAP::\
Lite::VERSION =~ /0\\.[1-5]/) {\n        $jobid = \
$soap->runBlastpgp($paramsData, $contentData);\n  \
  }\n    # For SOAP::Lite 0.69 and later parameter\
 handling is different, so pass\n    # undef's for\
 templated params, and then pass the formatted arg\
s.\n    else {\n        $jobid = $soap->runBlastpg\
p(undef, undef,\n				    $paramsData, $contentData\
);\n    }\n\n    if (defined($async)) {\n	print ST\
DOUT $jobid, \"\\n\";\n        print STDERR \"To c\
heck status: $scriptName --status --jobid $jobid\\\
n\";\n    } else { # Synchronous mode\n        pri\
nt STDERR \"JobId: $jobid\\n\";\n        sleep 1;\\
n        getResults($jobid);\n    }\n}\n\nsub getI\
ds($) {\n    $jobid = shift;\n    my $results = $s\
oap->getIds($jobid);\n    for $result (@$results){\
\n	print \"$result\\n\";\n    }\n}\n\nsub clientPo\
ll($) {\n    my $jobid = shift;\n    my $result = \
'PENDING';\n    # Check status and wait if not fin\
ished\n    #print STDERR \"Checking status: $jobid\
\\n\";\n    while($result eq 'RUNNING' || $result \
eq 'PENDING') {\n        $result = $soap->checkSta\
tus($jobid);\n        print STDERR \"$result\\n\";\
\n        if($result eq 'RUNNING' || $result eq 'P\
ENDING') {\n            # Wait before polling agai\
n.\n            sleep $checkInterval;\n        }\n\
    }\n}\n\nsub getResults($) {\n    $jobid = shif\
t;\n    # Check status, and wait if not finished\n\
    clientPoll($jobid);\n    # Use JobId if output\
 file name is not defined\n    unless(defined($out\
file)) {\n        $outfile=$jobid;\n    }\n    # G\
et list of data types\n    my $resultTypes = $soap\
->getResults($jobid);\n    # Get the data and writ\
e it to a file\n    if(defined($outformat)) { # Sp\
ecified data type\n        my $selResultType;\n   \
     foreach my $resultType (@$resultTypes) {\n   \
         if($resultType->{type} eq $outformat) {\n\
                $selResultType = $resultType;\n   \
         }\n        }\n        $res=$soap->poll($j\
obid, $selResultType->{type});\n        write_file\
($outfile.'.'.$selResultType->{ext}, $res);\n    }\
 else { # Data types available\n        # Write a \
file for each output type\n        for my $resultT\
ype (@$resultTypes){\n            #print \"Getting\
 $resultType->{type}\\n\";\n            $res=$soap\
->poll($jobid, $resultType->{type});\n            \
write_file($outfile.'.'.$resultType->{ext}, $res);\
\n        }\n    }\n}\n\nsub read_file($) {\n    m\
y $filename = shift;\n    open(FILE, $filename);\n\
    my $content;\n    my $buffer;\n    while(sysre\
ad(FILE, $buffer, 1024)) {\n	$content.= $buffer;\n\
    }\n    close(FILE);  \n    return $content;\n}\
\n\nsub write_file($$) {\n    my ($tmp,$entity) = \
@_;\n    print STDERR \"Creating result file: \".$\
tmp.\"\\n\";\n    unless(open (FILE, \">$tmp\")) {\
\n	return 0;\n    }\n    syswrite(FILE, $entity);\\
n    close (FILE);\n    return 1;\n}\n\nsub usage \
{\n    print STDERR <<EOF\nBlastpgp\n========\n   \
\nThe blastpgp program implements the PSI-BLAST an\
d PHI-BLAST variations\nof NCBI BLAST.\n\nFor more\
 detailed help information refer to\nhttp://www.eb\
i.ac.uk/blastpgp/blastpsi_help_frame.html\n \nBlas\
tpgp specific options:\n\n[Required]\n\n      --mo\
de            : str  : search mode to use: PSI-Bla\
st or PHI-Blast\n  -d, --database        : str  : \
protein database to search\n  seqFile             \
  : file : query sequence\n\n[Optional]\n\n  -M, -\
-matrix          : str  : scoring matrix\n  -e, --\
exp             : real : Expectation value\n  -h, \
--expmulti        : real : threshold (multipass mo\
del)\n  -F, --filter          : str  : filter quer\
y sequence with SEG [T,F]\n  -m, --align          \
 : int  : alignment view option:\n                \
                 0 - pairwise, 1 - M/S identities,\
\n                                 2 - M/S non-ide\
ntities, 3 - Flat identities,\n                   \
              4 - Flat non-identities\n  -G, --ope\
ngap         : int  : cost to open a gap\n  -E, --\
extendgap       : int  : cost to extend a gap\n  -\
g, --gapalign        : str  : Gapped [T,F]\n  -v, \
--scores          : int  : number of scores to be \
reported\n  -j, --maxpasses       : int  : number \
of iterations\n  -X, --dropoff         : int  : Dr\
opoff score\n  -Z, --finaldropoff    : int  : Drop\
off for final alignment\n  -S, --startregion     :\
 int  : Start of required region in query\n  -H, -\
-endregion       : int  : End of required region i\
n query\n  -k, --pattern         : str  : Hit File\
 (PHI-BLAST only)\n  -p, --usagemode       : str  \
: Program option (PHI-BLAST only):\n              \
                   blastpgp, patseedp, seedp\n\n[G\
eneral]\n\n      --help            :      : prints\
 this help text\n  -a, --async           :      : \
forces to make an asynchronous query\n      --stat\
us          :      : poll for the status of a job\\
n      --polljob         :      : poll for the res\
ults of a job\n      --jobid           : str  : jo\
bid of an asynchronous job\n      --ids           \
  :      : get hit identifiers for result \n  -O, \
--outfile         : str  : name of the file result\
s should be written to\n                          \
       (default is based on the jobid)\n  -o, --ou\
tformat       : str  : txt or xml output (no file \
is written)\n      --trace           :      : show\
 SOAP messages being interchanged\n\nSynchronous j\
ob:\n\n  The results/errors are returned as soon a\
s the job is finished.\n  Usage: blastpgp.pl --ema\
il <your@email> [options...] seqfile\n  Returns: s\
aves the results to disk\n\nAsynchronous job:\n\n \
 Use this if you want to retrieve the results at a\
 later time. The results\n  are stored for up to 2\
4 hours.\n  The asynchronous submission mode is re\
commended when users are submitting\n  batch jobs \
or large database searches\n  Usage: blastpgp.pl -\
-email <your@email> --async [options...] seqFile\n\
  Returns: jobid\n\n  Use the jobid to query for t\
he status of the job.\n  Usage: blastpgp.pl --stat\
us --jobid <jobId>\n  Returns: string indicating t\
he status of the job\n    DONE - job has finished\\
n    RUNNING - job is running\n    NOT_FOUND - job\
 cannot be found\n    ERROR - the jobs has encount\
ered an error\n\n  When done, use the jobid to ret\
rieve the results of the job.\n  Usage: blastpgp.p\
l --polljob --jobid <jobId> [--outfile <fileName>]\
\n  Returns: saves the results to disk\nEOF\n;\n}\\
n","\n=head1 NAME\n\nncbiblast.pl\n\n=head1 DESCRI\
PTION\n\nNCBI Blast (REST) web service Perl client\
 using L<LWP>.\n\nTested with:\n\n=over\n\n=item *\
\nL<LWP> 6.35, L<XML::Simple> 2.25 and Perl 5.22.0\
 (MacOS 10.13.6)\n\n=back\n\nFor further informati\
on see:\n\n=over\n\n=item *\nL<https://www.ebi.ac.\
uk/Tools/webservices/>\n\n=back\n\n=head1 LICENSE\\
n\nCopyright 2012-2018 EMBL - European Bioinformat\
ics Institute\n\nLicensed under the Apache License\
, Version 2.0 (the \"License\");\nyou may not use \
this file except in compliance with the License.\n\
You may obtain a copy of the License at\n\n    htt\
p://www.apache.org/licenses/LICENSE-2.0\n\nUnless \
required by applicable law or agreed to in writing\
, software\ndistributed under the License is distr\
ibuted on an \"AS IS\" BASIS,\nWITHOUT WARRANTIES \
OR CONDITIONS OF ANY KIND, either express or impli\
ed.\nSee the License for the specific language gov\
erning permissions and\nlimitations under the Lice\
nse.\n\nPerl Client Automatically generated with:\\
nhttps://github.com/ebi-wp/webservice-clients-gene\
rator\n\n=cut\n\nuse strict;\nuse warnings;\n\nuse\
 English;\nuse LWP;\nuse XML::Simple;\nuse Getopt:\
:Long qw(:config no_ignore_case bundling);\nuse Fi\
le::Basename;\nuse Data::Dumper;\nuse Time::HiRes \
qw(usleep);\n\nmy $baseUrl = 'https://www.ebi.ac.u\
k/Tools/services/rest/ncbiblast';\nmy $version = '\
2019-07-03 16:26';\n\nmy $checkInterval = 3;\n\nmy\
 $maxErrorStatusCount = 3;\n\nmy $outputLevel = 1;\
\n\nmy $numOpts = scalar(@ARGV);\nmy %params = (\n\
    'debugLevel' => 0,\n    'maxJobs'    => 1\n);\\
n\nGetOptions(\n    # Tool specific options\n    '\
program=s'       => \\$params{'program'},        #\
 The BLAST program to be used for the Sequence Sim\
ilarity Search.\n    'task=s'          => \\$param\
s{'task'},           # Task option (only selectabl\
e for blastn)\n    'matrix=s'        => \\$params{\
'matrix'},         # (Protein searches) The substi\
tution matrix used for scoring alignments when sea\
rching the database.\n    'alignments=i'    => \\$\
params{'alignments'},     # Maximum number of matc\
h alignments reported in the result output.\n    '\
scores=i'        => \\$params{'scores'},         #\
 Maximum number of match score summaries reported \
in the result output.\n    'exp=s'           => \\\
$params{'exp'},            # Limits the number of \
scores and alignments reported based on the expect\
ation value. This is the maximum number of times t\
he match is expected to occur by chance.\n    'dro\
poff=i'       => \\$params{'dropoff'},        # Th\
e amount a score can drop before gapped extension \
of word hits is halted\n    'match_scores=s'  => \\
\$params{'match_scores'},   # (Nucleotide searches\
) The match score is the bonus to the alignment sc\
ore when matching the same base. The mismatch is t\
he penalty when failing to match.\n    'gapopen=i'\
       => \\$params{'gapopen'},        # Penalty t\
aken away from the score when a gap is created in \
sequence. Increasing the gap openning penalty will\
 decrease the number of gaps in the final alignmen\
t.\n    'gapext=i'        => \\$params{'gapext'}, \
        # Penalty taken away from the score for ea\
ch base or residue in the gap. Increasing the gap \
extension penalty favors short gaps in the final a\
lignment, conversly decreasing the gap extension p\
enalty favors long gaps in the final alignment.\n \
   'filter=s'        => \\$params{'filter'},      \
   # Filter regions of low sequence complexity. Th\
is can avoid issues with low complexity sequences \
where matches are found due to composition rather \
than meaningful sequence similarity. However in so\
me cases filtering also masks regions of interest \
and so should be used with caution.\n    'seqrange\
=s'      => \\$params{'seqrange'},       # Specify\
 a range or section of the input sequence to use i\
n the search. Example: Specifying '34-89' in an in\
put sequence of total length 100, will tell BLAST \
to only use residues 34 to 89, inclusive.\n    'ga\
palign'        => \\$params{'gapalign'},       # T\
his is a true/false setting that tells the program\
 the perform optimised alignments within regions i\
nvolving gaps. If set to true, the program will pe\
rform an alignment using gaps. Otherwise, if it is\
 set to false, it will report only individual HSP \
where two sequence match each other, and thus will\
 not produce alignments with gaps.\n    'wordsize=\
i'      => \\$params{'wordsize'},       # Word siz\
e for wordfinder algorithm\n    'compstats=s'     \
=> \\$params{'compstats'},      # Use composition-\
based statistics.\n    'align=i'         => \\$par\
ams{'align'},          # Formating for the alignme\
nts\n    'transltable=i'   => \\$params{'transltab\
le'},    # Query Genetic code to use in translatio\
n\n    'stype=s'         => \\$params{'stype'},   \
       # Indicates if the sequence is protein or D\
NA/RNA.\n    'sequence=s'      => \\$params{'seque\
nce'},       # The query sequence can be entered d\
irectly into this form. The sequence can be in GCG\
, FASTA, EMBL (Nucleotide only), GenBank, PIR, NBR\
F, PHYLIP or UniProtKB/Swiss-Prot (Protein only) f\
ormat. A partially formatted sequence is not accep\
ted. Adding a return to the end of the sequence ma\
y help certain applications understand the input. \
Note that directly using data from word processors\
 may yield unpredictable results as hidden/control\
 characters may be present.\n    'database=s'     \
 => \\$params{'database'},       # Database\n    #\
 Generic options\n    'email=s'         => \\$para\
ms{'email'},          # User e-mail address\n    '\
title=s'         => \\$params{'title'},          #\
 Job title\n    'outfile=s'       => \\$params{'ou\
tfile'},        # Output file name\n    'outformat\
=s'     => \\$params{'outformat'},      # Output f\
ile type\n    'jobid=s'         => \\$params{'jobi\
d'},          # JobId\n    'help|h'          => \\\
$params{'help'},           # Usage help\n    'asyn\
cjob'        => \\$params{'asyncjob'},       # Asy\
nchronous submission\n    'polljob'         => \\$\
params{'polljob'},        # Get results\n    'poll\
Freq=f'      => \\$params{'pollFreq'},       # Pol\
l Frequency\n    'resultTypes'     => \\$params{'r\
esultTypes'},    # Get result types\n    'status' \
         => \\$params{'status'},         # Get sta\
tus\n    'params'          => \\$params{'params'},\
         # List input parameters\n    'paramDetail\
=s'   => \\$params{'paramDetail'},    # Get detail\
s for parameter\n    'multifasta'      => \\$param\
s{'multifasta'},     # Multiple fasta input\n    '\
useSeqId'        => \\$params{'useSeqId'},       #\
 Seq Id file name\n    'maxJobs=i'       => \\$par\
ams{'maxJobs'},        # Max. parallel jobs\n\n   \
 'verbose'         => \\$params{'verbose'},       \
 # Increase output level\n    'version'         =>\
 \\$params{'version'},        # Prints out the ver\
sion of the Client and exit.\n    'quiet'         \
  => \\$params{'quiet'},          # Decrease outpu\
t level\n    'debugLevel=i'    => \\$params{'debug\
Level'},     # Debugging level\n    'baseUrl=s'   \
    => \\$baseUrl,                  # Base URL for\
 service.\n);\nif ($params{'verbose'}) {$outputLev\
el++}\nif ($params{'quiet'}) {$outputLevel--}\nif \
($params{'pollFreq'}) {$checkInterval = $params{'p\
ollFreq'} * 1000 * 1000}\nif ($params{'baseUrl'}) \
{$baseUrl = $params{'baseUrl'}}\n\n&print_debug_me\
ssage('MAIN', 'LWP::VERSION: ' . $LWP::VERSION,\n \
   1);\n\n&print_debug_message('MAIN', \"params:\\\
n\" . Dumper(\\%params), 11);\n\nmy $ua;\n\nmy $sc\
riptName = basename($0, ());\n\nif ($params{'help'\
} || $numOpts == 0) {\n    &usage();\n    exit(0);\
\n}\n\n&print_debug_message('MAIN', 'baseUrl: ' . \
$baseUrl, 1);\nif (\n    !(\n        $params{'poll\
job'}\n            || $params{'resultTypes'}\n    \
        || $params{'status'}\n            || $para\
ms{'params'}\n            || $params{'paramDetail'\
}\n            || $params{'version'}\n    )\n     \
   && !(defined($ARGV[0]) || defined($params{'sequ\
ence'}))\n) {\n\n    # Bad argument combination, s\
o print error message and usage\n    print STDERR \
'Error: bad option combination', \"\\n\";\n    &us\
age();\n    exit(1);\n}\nelsif ($params{'params'})\
 {\n    &print_tool_params();\n}\n\nelsif ($params\
{'paramDetail'}) {\n    &print_param_details($para\
ms{'paramDetail'});\n}\n\nelsif ($params{'version'\
}) {\n  print STDOUT 'Revision: ' . $version, \"\\\
n\";\n  exit(1);\n}\n\nelsif ($params{'status'} &&\
 defined($params{'jobid'})) {\n    &print_job_stat\
us($params{'jobid'});\n}\n\nelsif ($params{'result\
Types'} && defined($params{'jobid'})) {\n    &prin\
t_result_types($params{'jobid'});\n}\n\nelsif ($pa\
rams{'polljob'} && defined($params{'jobid'})) {\n \
   &get_results($params{'jobid'});\n}\n\nelse {\n \
   # Multiple input sequence mode, assume fasta fo\
rmat.\n    if (defined($params{'multifasta'}) && $\
params{'multifasta'}) {\n        &multi_submit_job\
();\n    }\n\n    # Entry identifier list file.\n \
   elsif ((defined($params{'sequence'}) && $params\
{'sequence'} =~ m/^\\@/)\n        || (defined($ARG\
V[0]) && $ARGV[0] =~ m/^\\@/)) {\n        my $list\
_filename = $params{'sequence'} || $ARGV[0];\n    \
    $list_filename =~ s/^\\@//;\n        &list_fil\
e_submit_job($list_filename);\n    }\n    # Defaul\
t: single sequence/identifier.\n    else {\n      \
  # Warn for invalid batch only option use.\n     \
   if (defined($params{'useSeqId'}) && $params{'us\
eSeqId'}) {\n            print STDERR \"Warning: -\
-useSeqId option ignored.\\n\";\n            delet\
e $params{'useSeqId'};\n        }\n        if (def\
ined($params{'maxJobs'}) && $params{'maxJobs'} > 1\
) {\n            print STDERR \"Warning: --maxJobs\
 option ignored.\\n\";\n            $params{'maxJo\
bs'} = 1;\n        }\n        # Load the sequence \
data and submit.\n        &submit_job(&load_data()\
);\n    }\n}\n\n\n\n=head1 FUNCTIONS\n\n=cut\n\n\n\
=head2 rest_user_agent()\n\nGet a LWP UserAgent to\
 use to perform REST requests.\n\n  my $ua = &rest\
_user_agent();\n\n=cut\n\nsub rest_user_agent() {\\
n    print_debug_message('rest_user_agent', 'Begin\
', 21);\n    # Create an LWP UserAgent for making \
HTTP calls.\n    my $ua = LWP::UserAgent->new();\n\
    # Set 'User-Agent' HTTP header to identifiy th\
e client.\n    my $revisionNumber = 0;\n    $revis\
ionNumber = \"Revision: \" . $version;\n    $ua->a\
gent(\"EBI-Sample-Client/$revisionNumber ($scriptN\
ame; $OSNAME) \" . $ua->agent());\n    # Configure\
 HTTP proxy support from environment.\n    $ua->en\
v_proxy;\n    print_debug_message('rest_user_agent\
', 'End', 21);\n    return $ua;\n}\n\n=head2 rest_\
error()\n\nCheck a REST response for an error cond\
ition. An error is mapped to a die.\n\n  &rest_err\
or($response, $content_data);\n\n=cut\n\nsub rest_\
error() {\n    print_debug_message('rest_error', '\
Begin', 21);\n    my $response = shift;\n    my $c\
ontentdata;\n    if (scalar(@_) > 0) {\n        $c\
ontentdata = shift;\n    }\n    if (!defined($cont\
entdata) || $contentdata eq '') {\n        $conten\
tdata = $response->content();\n    }\n    # Check \
for HTTP error codes\n    if ($response->is_error)\
 {\n        my $error_message = '';\n        # HTM\
L response.\n        if ($contentdata =~ m/<h1>([^\
<]+)<\\/h1>/) {\n            $error_message = $1;\\
n        }\n        #  XML response.\n        elsi\
f ($contentdata =~ m/<description>([^<]+)<\\/descr\
iption>/) {\n            $error_message = $1;\n   \
     }\n        # die 'http status: ' . $response-\
>code . ' ' . $response->message . '  ' . $error_m\
essage;\n    }\n    print_debug_message('rest_erro\
r', 'End', 21);\n}\n\n=head2 rest_request()\n\nPer\
form a REST request (HTTP GET).\n\n  my $response_\
str = &rest_request($url);\n\n=cut\n\nsub rest_req\
uest {\n    print_debug_message('rest_request', 'B\
egin', 11);\n    my $requestUrl = shift;\n    prin\
t_debug_message('rest_request', 'URL: ' . $request\
Url, 11);\n\n    # Get an LWP UserAgent.\n    $ua \
= &rest_user_agent() unless defined($ua);\n    # A\
vailable HTTP compression methods.\n    my $can_ac\
cept;\n    eval {\n        $can_accept = HTTP::Mes\
sage::decodable();\n    };\n    $can_accept = '' u\
nless defined($can_accept);\n    # Perform the req\
uest\n    my $response = $ua->get($requestUrl,\n  \
      'Accept-Encoding' => $can_accept, # HTTP com\
pression.\n    );\n    print_debug_message('rest_r\
equest', 'HTTP status: ' . $response->code,\n     \
   11);\n    print_debug_message('rest_request',\n\
        'response length: ' . length($response->co\
ntent()), 11);\n    print_debug_message('rest_requ\
est',\n        'request:' . \"\\n\" . $response->r\
equest()->as_string(), 32);\n    print_debug_messa\
ge('rest_request',\n        'response: ' . \"\\n\"\
 . $response->as_string(), 32);\n    # Unpack poss\
ibly compressed response.\n    my $retVal;\n    if\
 (defined($can_accept) && $can_accept ne '') {\n  \
      $retVal = $response->decoded_content();\n   \
 }\n    # If unable to decode use orginal content.\
\n    $retVal = $response->content() unless define\
d($retVal);\n    # Check for an error.\n    &rest_\
error($response, $retVal);\n    print_debug_messag\
e('rest_request', 'retVal: ' . $retVal, 12);\n    \
print_debug_message('rest_request', 'End', 11);\n\\
n    # Return the response data\n    return $retVa\
l;\n}\n\n=head2 rest_get_parameters()\n\nGet list \
of tool parameter names.\n\n  my (@param_list) = &\
rest_get_parameters();\n\n=cut\n\nsub rest_get_par\
ameters {\n    print_debug_message('rest_get_param\
eters', 'Begin', 1);\n    my $url = $baseUrl . '/p\
arameters/';\n    my $param_list_xml_str = rest_re\
quest($url);\n    my $param_list_xml = XMLin($para\
m_list_xml_str);\n    my (@param_list) = @{$param_\
list_xml->{'id'}};\n    print_debug_message('rest_\
get_parameters', 'End', 1);\n    return(@param_lis\
t);\n}\n\n=head2 rest_get_parameter_details()\n\nG\
et details of a tool parameter.\n\n  my $paramDeta\
il = &rest_get_parameter_details($param_name);\n\n\
=cut\n\nsub rest_get_parameter_details {\n    prin\
t_debug_message('rest_get_parameter_details', 'Beg\
in', 1);\n    my $parameterId = shift;\n    print_\
debug_message('rest_get_parameter_details',\n     \
   'parameterId: ' . $parameterId, 1);\n    my $ur\
l = $baseUrl . '/parameterdetails/' . $parameterId\
;\n    my $param_detail_xml_str = rest_request($ur\
l);\n    my $param_detail_xml = XMLin($param_detai\
l_xml_str);\n    print_debug_message('rest_get_par\
ameter_details', 'End', 1);\n    return($param_det\
ail_xml);\n}\n\n=head2 rest_run()\n\nSubmit a job.\
\n\n  my $job_id = &rest_run($email, $title, \\%pa\
rams );\n\n=cut\n\nsub rest_run {\n    print_debug\
_message('rest_run', 'Begin', 1);\n    my $email =\
 shift;\n    my $title = shift;\n    my $params = \
shift;\n    $email = '' if (!$email);\n    print_d\
ebug_message('rest_run', 'email: ' . $email, 1);\n\
    if (defined($title)) {\n        print_debug_me\
ssage('rest_run', 'title: ' . $title, 1);\n    }\n\
    print_debug_message('rest_run', 'params: ' . D\
umper($params), 1);\n\n    # Get an LWP UserAgent.\
\n    $ua = &rest_user_agent() unless defined($ua)\
;\n\n    # Clean up parameters\n    my (%tmp_param\
s) = %{$params};\n    $tmp_params{'email'} = $emai\
l;\n    $tmp_params{'title'} = $title;\n    foreac\
h my $param_name (keys(%tmp_params)) {\n        if\
 (!defined($tmp_params{$param_name})) {\n         \
   delete $tmp_params{$param_name};\n        }\n  \
  }\n\n    # Submit the job as a POST\n    my $url\
 = $baseUrl . '/run';\n    my $response = $ua->pos\
t($url, \\%tmp_params);\n    print_debug_message('\
rest_run', 'HTTP status: ' . $response->code, 11);\
\n    print_debug_message('rest_run',\n        're\
quest:' . \"\\n\" . $response->request()->as_strin\
g(), 11);\n    print_debug_message('rest_run',\n  \
      'response: ' . length($response->as_string()\
) . \"\\n\" . $response->as_string(), 11);\n\n    \
# Check for an error.\n    &rest_error($response);\
\n\n    # The job id is returned\n    my $job_id =\
 $response->content();\n    print_debug_message('r\
est_run', 'End', 1);\n    return $job_id;\n}\n\n=h\
ead2 rest_get_status()\n\nCheck the status of a jo\
b.\n\n  my $status = &rest_get_status($job_id);\n\\
n=cut\n\nsub rest_get_status {\n    print_debug_me\
ssage('rest_get_status', 'Begin', 1);\n    my $job\
_id = shift;\n    print_debug_message('rest_get_st\
atus', 'jobid: ' . $job_id, 2);\n    my $status_st\
r = 'UNKNOWN';\n    my $url = $baseUrl . '/status/\
' . $job_id;\n    $status_str = &rest_request($url\
);\n    print_debug_message('rest_get_status', 'st\
atus_str: ' . $status_str, 2);\n    print_debug_me\
ssage('rest_get_status', 'End', 1);\n    return $s\
tatus_str;\n}\n\n=head2 rest_get_result_types()\n\\
nGet list of result types for finished job.\n\n  m\
y (@result_types) = &rest_get_result_types($job_id\
);\n\n=cut\n\nsub rest_get_result_types {\n    pri\
nt_debug_message('rest_get_result_types', 'Begin',\
 1);\n    my $job_id = shift;\n    print_debug_mes\
sage('rest_get_result_types', 'jobid: ' . $job_id,\
 2);\n    my (@resultTypes);\n    my $url = $baseU\
rl . '/resulttypes/' . $job_id;\n    my $result_ty\
pe_list_xml_str = &rest_request($url);\n    my $re\
sult_type_list_xml = XMLin($result_type_list_xml_s\
tr);\n    (@resultTypes) = @{$result_type_list_xml\
->{'type'}};\n    print_debug_message('rest_get_re\
sult_types',\n        scalar(@resultTypes) . ' res\
ult types', 2);\n    print_debug_message('rest_get\
_result_types', 'End', 1);\n    return(@resultType\
s);\n}\n\n=head2 rest_get_result()\n\nGet result d\
ata of a specified type for a finished job.\n\n  m\
y $result = rest_get_result($job_id, $result_type)\
;\n\n=cut\n\nsub rest_get_result {\n    print_debu\
g_message('rest_get_result', 'Begin', 1);\n    my \
$job_id = shift;\n    my $type = shift;\n    print\
_debug_message('rest_get_result', 'jobid: ' . $job\
_id, 1);\n    print_debug_message('rest_get_result\
', 'type: ' . $type, 1);\n    my $url = $baseUrl .\
 '/result/' . $job_id . '/' . $type;\n    my $resu\
lt = &rest_request($url);\n    print_debug_message\
('rest_get_result', length($result) . ' characters\
',\n        1);\n    print_debug_message('rest_get\
_result', 'End', 1);\n    return $result;\n}\n\n\n\
=head2 print_debug_message()\n\nPrint debug messag\
e at specified debug level.\n\n  &print_debug_mess\
age($method_name, $message, $level);\n\n=cut\n\nsu\
b print_debug_message {\n    my $function_name = s\
hift;\n    my $message = shift;\n    my $level = s\
hift;\n    if ($level <= $params{'debugLevel'}) {\\
n        print STDERR '[', $function_name, '()] ',\
 $message, \"\\n\";\n    }\n}\n\n=head2 print_tool\
_params()\n\nPrint list of tool parameters.\n\n  &\
print_tool_params();\n\n=cut\n\nsub print_tool_par\
ams {\n    print_debug_message('print_tool_params'\
, 'Begin', 1);\n    my (@param_list) = &rest_get_p\
arameters();\n    foreach my $param (sort (@param_\
list)) {\n        print $param, \"\\n\";\n    }\n \
   print_debug_message('print_tool_params', 'End',\
 1);\n}\n\n=head2 print_param_details()\n\nPrint d\
etails of a tool parameter.\n\n  &print_param_deta\
ils($param_name);\n\n=cut\n\nsub print_param_detai\
ls {\n    print_debug_message('print_param_details\
', 'Begin', 1);\n    my $paramName = shift;\n    p\
rint_debug_message('print_param_details', 'paramNa\
me: ' . $paramName, 2);\n    my $paramDetail = &re\
st_get_parameter_details($paramName);\n    print $\
paramDetail->{'name'}, \"\\t\", $paramDetail->{'ty\
pe'}, \"\\n\";\n    print $paramDetail->{'descript\
ion'}, \"\\n\";\n    if (defined($paramDetail->{'v\
alues'}->{'value'})) {\n        if (ref($paramDeta\
il->{'values'}->{'value'}) eq 'ARRAY') {\n        \
    foreach my $value (@{$paramDetail->{'values'}-\
>{'value'}}) {\n                &print_param_value\
($value);\n            }\n        }\n        else \
{\n            &print_param_value($paramDetail->{'\
values'}->{'value'});\n        }\n    }\n    print\
_debug_message('print_param_details', 'End', 1);\n\
}\n\n=head2 print_param_value()\n\nPrint details o\
f a tool parameter value.\n\n  &print_param_detail\
s($param_value);\n\nUsed by print_param_details() \
to handle both singluar and array values.\n\n=cut\\
n\nsub print_param_value {\n    my $value = shift;\
\n    print $value->{'value'};\n    if ($value->{'\
defaultValue'} eq 'true') {\n        print \"\\t\"\
, 'default';\n    }\n    print \"\\n\";\n    print\
 \"\\t\", $value->{'label'}, \"\\n\";\n    if (def\
ined($value->{'properties'})) {\n        foreach\n\
        my $key (sort ( keys(%{$value->{'propertie\
s'}{'property'}}) )) {\n            if (ref($value\
->{'properties'}{'property'}{$key}) eq 'HASH'\n   \
             && defined($value->{'properties'}{'pr\
operty'}{$key}{'value'})\n            ) {\n       \
         print \"\\t\", $key, \"\\t\",\n          \
          $value->{'properties'}{'property'}{$key}\
{'value'}, \"\\n\";\n            }\n            el\
se {\n                print \"\\t\", $value->{'pro\
perties'}{'property'}{'key'},\n                   \
 \"\\t\", $value->{'properties'}{'property'}{'valu\
e'}, \"\\n\";\n                last;\n            \
}\n        }\n    }\n}\n\n=head2 print_job_status(\
)\n\nPrint status of a job.\n\n  &print_job_status\
($job_id);\n\n=cut\n\nsub print_job_status {\n    \
print_debug_message('print_job_status', 'Begin', 1\
);\n    my $jobid = shift;\n    print_debug_messag\
e('print_job_status', 'jobid: ' . $jobid, 1);\n   \
 if ($outputLevel > 0) {\n        print STDERR 'Ge\
tting status for job ', $jobid, \"\\n\";\n    }\n \
   my $result = &rest_get_status($jobid);\n    pri\
nt \"$result\\n\";\n    if ($result eq 'FINISHED' \
&& $outputLevel > 0) {\n        print STDERR \"To \
get results: perl $scriptName --polljob --jobid \"\
 . $jobid\n            . \"\\n\";\n    }\n    prin\
t_debug_message('print_job_status', 'End', 1);\n}\\
n\n=head2 print_result_types()\n\nPrint available \
result types for a job.\n\n  &print_result_types($\
job_id);\n\n=cut\n\nsub print_result_types {\n    \
print_debug_message('result_types', 'Begin', 1);\n\
    my $jobid = shift;\n    print_debug_message('r\
esult_types', 'jobid: ' . $jobid, 1);\n    if ($ou\
tputLevel > 0) {\n        print STDERR 'Getting re\
sult types for job ', $jobid, \"\\n\";\n    }\n   \
 my $status = &rest_get_status($jobid);\n    if ($\
status eq 'PENDING' || $status eq 'RUNNING') {\n  \
      print STDERR 'Error: Job status is ', $statu\
s,\n            '. To get result types the job mus\
t be finished.', \"\\n\";\n    }\n    else {\n    \
    my (@resultTypes) = &rest_get_result_types($jo\
bid);\n        if ($outputLevel > 0) {\n          \
  print STDOUT 'Available result types:', \"\\n\";\
\n        }\n        foreach my $resultType (@resu\
ltTypes) {\n            print STDOUT $resultType->\
{'identifier'}, \"\\n\";\n            if (defined(\
$resultType->{'label'})) {\n                print \
STDOUT \"\\t\", $resultType->{'label'}, \"\\n\";\n\
            }\n            if (defined($resultType\
->{'description'})) {\n                print STDOU\
T \"\\t\", $resultType->{'description'}, \"\\n\";\\
n            }\n            if (defined($resultTyp\
e->{'mediaType'})) {\n                print STDOUT\
 \"\\t\", $resultType->{'mediaType'}, \"\\n\";\n  \
          }\n            if (defined($resultType->\
{'fileSuffix'})) {\n                print STDOUT \\
"\\t\", $resultType->{'fileSuffix'}, \"\\n\";\n   \
         }\n        }\n        if ($status eq 'FIN\
ISHED' && $outputLevel > 0) {\n            print S\
TDERR \"\\n\", 'To get results:', \"\\n\",\n      \
          \"  perl $scriptName --polljob --jobid \\
" . $params{'jobid'} . \"\\n\",\n                \\
"  perl $scriptName --polljob --outformat <type> -\
-jobid \"\n                    . $params{'jobid'} \
. \"\\n\";\n        }\n    }\n    print_debug_mess\
age('result_types', 'End', 1);\n}\n\n=head2 submit\
_job()\n\nSubmit a job to the service.\n\n  &submi\
t_job($seq);\n\n=cut\n\nsub submit_job {\n    prin\
t_debug_message('submit_job', 'Begin', 1);\n\n    \
# Set input sequence\n    $params{'sequence'} = sh\
ift;\n    my $seq_id = shift;\n\n    # Load parame\
ters\n    &load_params();\n\n    # Submit the job\\
n    my $jobid = &rest_run($params{'email'}, $para\
ms{'title'}, \\%params);\n\n    # Asynchronous sub\
mission.\n    if (defined($params{'asyncjob'})) {\\
n        print STDOUT $jobid, \"\\n\";\n        if\
 ($outputLevel > 0) {\n            print STDERR\n \
               \"To check status: perl $scriptName\
 --status --jobid $jobid\\n\";\n        }\n    }\n\
\n    # Simulate synchronous submission serial mod\
e.\n    else {\n        if ($outputLevel > 0) {\n \
           print STDERR \"JobId: $jobid\\n\";\n   \
     } else {\n            print STDERR \"$jobid\\\
n\";\n        }\n        usleep($checkInterval);\n\
        # Get results.\n        &get_results($jobi\
d, $seq_id);\n\n    }\n    print_debug_message('su\
bmit_job', 'End', 1);\n    return $jobid;\n}\n=hea\
d2 multi_submit_job()\n\nSubmit multiple jobs assu\
ming input is a collection of fasta formatted sequ\
ences.\n\n  &multi_submit_job();\n\n=cut\n\nsub mu\
lti_submit_job {\n    print_debug_message('multi_s\
ubmit_job', 'Begin', 1);\n    my (@filename_list) \
= ();\n\n    # Query sequence\n    if (defined($AR\
GV[0])) {                  # Bare option\n        \
if (-f $ARGV[0] || $ARGV[0] eq '-') { # File\n    \
        push(@filename_list, $ARGV[0]);\n        }\
\n        else {\n            warn 'Warning: Input\
 file \"' . $ARGV[0] . '\" does not exist';\n     \
   }\n    }\n    if ($params{'sequence'}) {       \
                               # Via --sequence\n \
       if (-f $params{'sequence'} || $params{'sequ\
ence'} eq '-') { # File\n            push(@filenam\
e_list, $params{'sequence'});\n        }\n        \
else {\n            warn 'Warning: Input file \"'\\
n                . $params{'sequence'}\n          \
      . '\" does not exist';\n        }\n    }\n\n\
    # Job identifier tracking for parallel executi\
on.\n    my @jobid_list = ();\n    my $job_number \
= 0;\n    $/ = '>';\n    foreach my $filename (@fi\
lename_list) {\n        my $INFILE;\n        if ($\
filename eq '-') { # STDIN.\n            open($INF\
ILE, '<-')\n                or die 'Error: unable \
to STDIN (' . $! . ')';\n        }\n        else {\
 # File.\n            open($INFILE, '<', $filename\
)\n                or die 'Error: unable to open f\
ile '\n                . $filename . ' ('\n       \
         . $! . ')';\n        }\n        while (<$\
INFILE>) {\n            my $seq = $_;\n           \
 $seq =~ s/>$//;\n            if ($seq =~ m/(\\S+)\
/) {\n                my $seq_id = $1;\n          \
      print STDERR \"Submitting job for: $seq_id\\\
n\"\n                    if ($outputLevel > 0);\n \
               $seq = '>' . $seq;\n               \
 &print_debug_message('multi_submit_job', $seq, 11\
);\n                $job_number++;\n              \
  my $job_id = &submit_job($seq, $seq_id);\n\n    \
            my $job_info_str = sprintf('%s %d %d',\
 $job_id, 0, $job_number);\n\n                push\
(@jobid_list, $job_info_str);\n            }\n\n  \
          # Parallel mode, wait for job(s) to fini\
sh to free slots.\n            while ($params{'max\
Jobs'} > 1\n                && scalar(@jobid_list)\
 >= $params{'maxJobs'}) {\n                &_job_l\
ist_poll(\\@jobid_list);\n                print_de\
bug_message('multi_submit_job',\n                 \
   'Remaining jobs: ' . scalar(@jobid_list), 1);\n\
            }\n        }\n        close $INFILE;\n\
    }\n\n    # Parallel mode, wait for remaining j\
obs to finish.\n    while ($params{'maxJobs'} > 1 \
&& scalar(@jobid_list) > 0) {\n        &_job_list_\
poll(\\@jobid_list);\n        print_debug_message(\
'multi_submit_job',\n            'Remaining jobs: \
' . scalar(@jobid_list), 1);\n    }\n    print_deb\
ug_message('multi_submit_job', 'End', 1);\n}\n\n\n\
=head2 _job_list_poll()\n\nPoll the status of a li\
st of jobs and fetch results for finished jobs.\n\\
n  while(scalar(@jobid_list) > 0) {\n    &_job_lis\
t_poll(\\@jobid_list);\n  }\n\n=cut\n\nsub _job_li\
st_poll {\n    print_debug_message('_job_list_poll\
', 'Begin', 1);\n    my $jobid_list = shift;\n    \
print_debug_message('_job_list_poll', 'Num jobs: '\
 . scalar(@$jobid_list),\n        11);\n\n    # Lo\
op though job Id list polling job status.\n    for\
 (my $jobNum = (scalar(@$jobid_list) - 1); $jobNum\
 > -1; $jobNum--) {\n        my ($jobid, $seq_id, \
$error_count, $job_number) =\n            split(/\\
\s+/, $jobid_list->[$jobNum]);\n        print_debu\
g_message('_job_list_poll', 'jobNum: ' . $jobNum, \
12);\n        print_debug_message('_job_list_poll'\
,\n            'Job info: ' . $jobid_list->[$jobNu\
m], 12);\n\n        # Get job status.\n        my \
$job_status = &rest_get_status($jobid);\n        p\
rint_debug_message('_job_list_poll', 'Status: ' . \
$job_status, 12);\n\n        # Fetch results and r\
emove finished/failed jobs from list.\n        if \
(\n            !(\n                $job_status eq \
'RUNNING'\n                    || $job_status eq '\
PENDING'\n                    || ($job_status eq '\
ERROR'\n                    && $error_count < $max\
ErrorStatusCount)\n            )\n        ) {\n   \
         if ($job_status eq 'ERROR' || $job_status\
 eq 'FAILED') {\n                print STDERR\n   \
                 \"Warning: job $jobid failed for \
sequence $job_number: $seq_id\\n\";\n            }\
\n			# Duplicated getting results.\n            #&\
get_results($jobid, $seq_id);\n            splice(\
@$jobid_list, $jobNum, 1);\n        }\n        els\
e {\n\n            # Update error count, increment\
 for new error or clear old errors.\n            i\
f ($job_status eq 'ERROR') {\n                $err\
or_count++;\n            }\n            elsif ($er\
ror_count > 0) {\n                $error_count--;\\
n            }\n\n            # Update job trackin\
g info.\n            my $job_info_str = sprintf('%\
s %s %d %d',\n                $jobid, $seq_id, $er\
ror_count, $job_number);\n            $jobid_list-\
>[$jobNum] = $job_info_str;\n        }\n    }\n   \
 print_debug_message('_job_list_poll', 'Num jobs: \
' . scalar(@$jobid_list),\n        11);\n    print\
_debug_message('_job_list_poll', 'End', 1);\n}\n\n\
=head2 list_file_submit_job()\n\nSubmit multiple j\
obs using a file containing a list of entry identi\
fiers as\ninput.\n\n  &list_file_submit_job($list_\
filename)\n\n=cut\n\nsub list_file_submit_job {\n \
   print_debug_message('list_file_submit_job', 'Be\
gin', 1);\n    my $filename = shift;\n\n    # Open\
 the file of identifiers.\n    my $LISTFILE;\n    \
if ($filename eq '-') { # STDIN.\n        open($LI\
STFILE, '<-')\n            or die 'Error: unable t\
o STDIN (' . $! . ')';\n    }\n    else { # File.\\
n        open($LISTFILE, '<', $filename)\n        \
    or die 'Error: unable to open file ' . $filena\
me . ' (' . $! . ')';\n    }\n\n    # Job identifi\
er tracking for parallel execution.\n    my @jobid\
_list = ();\n    my $job_number = 0;\n\n    # Iter\
ate over identifiers, submitting each job\n    whi\
le (<$LISTFILE>) {\n        my $line = $_;\n      \
  chomp($line);\n        if ($line ne '') {\n     \
       &print_debug_message('list_file_submit_job'\
, 'line: ' . $line, 2);\n            if ($line =~ \
m/\\w:\\w/) {\n                # Check this is an \
identifier\n                my $seq_id = $line;\n \
               print STDERR \"Submitting job for: \
$seq_id\\n\"\n                    if ($outputLevel\
 > 0);\n                $job_number++;\n          \
      my $job_id = &submit_job($seq_id, $seq_id);\\
n                my $job_info_str =\n             \
       sprintf('%s %s %d %d', $job_id, $seq_id, 0,\
 $job_number);\n                push(@jobid_list, \
$job_info_str);\n            }\n            else {\
\n                print STDERR\n                  \
  \"Warning: line \\\"$line\\\" is not recognised \
as an identifier\\n\";\n            }\n\n         \
   # Parallel mode, wait for job(s) to finish to f\
ree slots.\n            while ($params{'maxJobs'} \
> 1\n                && scalar(@jobid_list) >= $pa\
rams{'maxJobs'}) {\n                &_job_list_pol\
l(\\@jobid_list);\n                print_debug_mes\
sage('list_file_submit_job',\n                    \
'Remaining jobs: ' . scalar(@jobid_list), 1);\n   \
         }\n        }\n    }\n    close $LISTFILE;\
\n\n    # Parallel mode, wait for remaining jobs t\
o finish.\n    while ($params{'maxJobs'} > 1 && sc\
alar(@jobid_list) > 0) {\n        &_job_list_poll(\
\\@jobid_list);\n        print_debug_message('list\
_file_submit_job',\n            'Remaining jobs: '\
 . scalar(@jobid_list), 1);\n    }\n    print_debu\
g_message('list_file_submit_job', 'End', 1);\n}\n\\
n\n=head2 load_data()\n\nLoad sequence data from f\
ile or option specified on the command-line.\n\n  \
&load_data();\n\n=cut\n\nsub load_data {\n    prin\
t_debug_message('load_data', 'Begin', 1);\n    my \
$retSeq;\n\n    # Query sequence\n    if (defined(\
$ARGV[0])) {                  # Bare option\n     \
   if (-f $ARGV[0] || $ARGV[0] eq '-') { # File\n \
           $retSeq = &read_file($ARGV[0]);\n      \
  }\n        else { # DB:ID or sequence\n         \
   $retSeq = $ARGV[0];\n        }\n    }\n    if (\
$params{'sequence'}) {                            \
          # Via --sequence\n        if (-f $params\
{'sequence'} || $params{'sequence'} eq '-') { # Fi\
le\n            $retSeq = &read_file($params{'sequ\
ence'});\n        }\n        else { # DB:ID or seq\
uence\n            $retSeq = $params{'sequence'};\\
n        }\n    }\n    print_debug_message('load_d\
ata', 'End', 1);\n    return $retSeq;\n}\n\n=head2\
 load_params()\n\nLoad job parameters from command\
-line options.\n\n  &load_params();\n\n=cut\n\nsub\
 load_params {\n    print_debug_message('load_para\
ms', 'Begin', 1);\n\n    # Pass default values and\
 fix bools (without default value)\n    if ($param\
s{'stype'} eq 'protein') {\n        if (!$params{'\
task'}) {\n            $params{'task'} = 'blastp'\\
n        }\n    }\n    if ($params{'stype'} eq 'nu\
cleotide') {\n        if (!$params{'task'}) {\n   \
         $params{'task'} = 'blastn'\n        }\n  \
  }\n    if ($params{'stype'} eq 'vector') {\n    \
    if (!$params{'task'}) {\n            $params{'\
task'} = 'blastn'\n        }\n    }\n\n    if ($pa\
rams{'stype'} eq 'protein') {\n        if (!$param\
s{'matrix'}) {\n            $params{'matrix'} = 'B\
LOSUM62'\n        }\n    }\n    if ($params{'stype\
'} eq 'nucleotide') {\n        if (!$params{'matri\
x'}) {\n            $params{'matrix'} = 'NONE'\n  \
      }\n    }\n    if ($params{'stype'} eq 'vecto\
r') {\n        if (!$params{'matrix'}) {\n        \
    $params{'matrix'} = 'NONE'\n        }\n    }\n\
\n    if (!$params{'alignments'}) {\n        $para\
ms{'alignments'} = 50\n    }\n\n    if (!$params{'\
scores'}) {\n        $params{'scores'} = 50\n    }\
\n\n    if (!$params{'exp'}) {\n        $params{'e\
xp'} = '10'\n    }\n\n    if (!$params{'dropoff'})\
 {\n        $params{'dropoff'} = 0\n    }\n\n    i\
f ($params{'stype'} eq 'nucleotide') {\n        if\
 (!$params{'match_scores'}) {\n            $params\
{'match_scores'} = '1,-3'\n        }\n    }\n    i\
f ($params{'stype'} eq 'vector') {\n        if (!$\
params{'match_scores'}) {\n            $params{'ma\
tch_scores'} = '1,-3'\n        }\n    }\n\n    if \
(!$params{'gapopen'}) {\n        $params{'gapopen'\
} = -1\n    }\n\n    if (!$params{'gapext'}) {\n  \
      $params{'gapext'} = -1\n    }\n\n    if ($pa\
rams{'stype'} eq 'protein') {\n        if (!$param\
s{'filter'}) {\n            $params{'filter'} = 'F\
'\n        }\n    }\n    if ($params{'stype'} eq '\
nucleotide') {\n        if (!$params{'filter'}) {\\
n            $params{'filter'} = 'T'\n        }\n \
   }\n    if ($params{'stype'} eq 'vector') {\n   \
     if (!$params{'filter'}) {\n            $param\
s{'filter'} = 'T'\n        }\n    }\n\n    if (!$p\
arams{'gapalign'}) {\n        $params{'gapalign'} \
= 'true'\n    }\n\n    if (!$params{'compstats'}) \
{\n        $params{'compstats'} = 'F'\n    }\n\n  \
  if (!$params{'align'}) {\n        $params{'align\
'} = 0\n    }\n\n    if (!$params{'transltable'}) \
{\n        $params{'transltable'} = 1\n    }\n\n  \
  print_debug_message('load_params', 'End', 1);\n}\
\n\n=head2 client_poll()\n\nClient-side job pollin\
g.\n\n  &client_poll($job_id);\n\n=cut\n\nsub clie\
nt_poll {\n    print_debug_message('client_poll', \
'Begin', 1);\n    my $jobid = shift;\n    my $stat\
us = 'PENDING';\n\n    # Check status and wait if \
not finished. Terminate if three attempts get \"ER\
ROR\".\n    my $errorCount = 0;\n    while ($statu\
s eq 'RUNNING'\n        || $status eq 'PENDING'\n \
       || ($status eq 'ERROR' && $errorCount < 2))\
 {\n        $status = rest_get_status($jobid);\n  \
      print STDERR \"$status\\n\" if ($outputLevel\
 > 0);\n        if ($status eq 'ERROR') {\n       \
     $errorCount++;\n        }\n        elsif ($er\
rorCount > 0) {\n            $errorCount--;\n     \
   }\n        if ($status eq 'RUNNING'\n          \
  || $status eq 'PENDING'\n            || $status \
eq 'ERROR') {\n\n            # Wait before polling\
 again.\n            usleep($checkInterval);\n    \
    }\n    }\n    print_debug_message('client_poll\
', 'End', 1);\n    return $status;\n}\n\n=head2 ge\
t_results()\n\nGet the results for a job identifie\
r.\n\n  &get_results($job_id);\n\n=cut\n\nsub get_\
results {\n    print_debug_message('get_results', \
'Begin', 1);\n    my $jobid = shift;\n    print_de\
bug_message('get_results', 'jobid: ' . $jobid, 1);\
\n    my $seq_id = shift;\n    print_debug_message\
('get_results', 'seq_id: ' . $seq_id, 1) if ($seq_\
id);\n\n    my $output_basename = $jobid;\n\n    #\
 Verbose\n    if ($outputLevel > 1) {\n        pri\
nt 'Getting results for job ', $jobid, \"\\n\";\n \
   }\n\n    # Check status, and wait if not finish\
ed\n    client_poll($jobid);\n\n    # Default outp\
ut file names use JobId, however the name can be s\
pecified...\n    if (defined($params{'outfile'})) \
{\n        $output_basename = $params{'outfile'};\\
n    }\n    # Or use sequence identifer.\n    elsi\
f (defined($params{'useSeqId'} && defined($seq_id)\
 && $seq_id ne '')) {\n        $output_basename = \
$seq_id;\n\n        # Make safe to use as a file n\
ame.\n        $output_basename =~ s/\\W/_/g;\n    \
}\n\n    # Use JobId if output file name is not de\
fined\n    else {\n        unless (defined($params\
{'outfile'})) {\n            #$params{'outfile'} =\
 $jobid;\n            $output_basename = $jobid;\n\
        }\n    }\n\n    # Get list of data types\n\
    my (@resultTypes) = rest_get_result_types($job\
id);\n\n\n    # Get the data and write it to a fil\
e\n    if (defined($params{'outformat'})) {\n     \
   # Specified data type\n        # check to see i\
f there are multiple formats (comma separated)\n  \
      my $sep = \",\";\n        my (@multResultTyp\
es);\n        if ($params{'outformat'} =~ /$sep/) \
{\n            @multResultTypes = split(',', $para\
ms{'outformat'});\n        }\n        else {\n    \
        $multResultTypes[0] = $params{'outformat'}\
;\n        }\n        # check if the provided form\
ats are recognised\n        foreach my $inputType \
(@multResultTypes) {\n            my $expectation \
= 0;\n            foreach my $resultType (@resultT\
ypes) {\n                if ($resultType->{'identi\
fier'} eq $inputType && $expectation eq 0) {\n    \
                $expectation = 1;\n               \
 }\n            }\n            if ($expectation ne\
 1) {\n                die 'Error: unknown result \
format \"' . $inputType . '\"';\n            }\n  \
      }\n        # if so get the files\n        my\
 $selResultType;\n        foreach my $resultType (\
@resultTypes) {\n            if (grep {$_ eq $resu\
ltType->{'identifier'}} @multResultTypes) {\n     \
           $selResultType = $resultType;\n        \
        my $result = rest_get_result($jobid, $selR\
esultType->{'identifier'});\n                if (d\
efined($params{'outfile'}) && $params{'outfile'} e\
q '-') {\n                    write_file($params{'\
outfile'}, $result);\n                }\n         \
       else {\n                    write_file(\n  \
                      $output_basename . '.'\n    \
                        . $selResultType->{'identi\
fier'} . '.'\n                            . $selRe\
sultType->{'fileSuffix'},\n                       \
 $result\n                    );\n                \
}\n            }\n        }\n    }\n    else { # D\
ata types available\n        # Write a file for ea\
ch output type\n        for my $resultType (@resul\
tTypes) {\n            if ($outputLevel > 1) {\n  \
              print STDERR 'Getting ', $resultType\
->{'identifier'}, \"\\n\";\n            }\n       \
     my $result = rest_get_result($jobid, $resultT\
ype->{'identifier'});\n            if (defined($pa\
rams{'outfile'}) && $params{'outfile'} eq '-') {\n\
                write_file($params{'outfile'}, $re\
sult);\n            }\n            else {\n       \
         write_file(\n                    $output_\
basename . '.'\n                        . $resultT\
ype->{'identifier'} . '.'\n                       \
 . $resultType->{'fileSuffix'},\n                 \
   $result\n                );\n            }\n   \
     }\n    }\n    print_debug_message('get_result\
s', 'End', 1);\n}\n\n=head2 read_file()\n\nRead a \
file into a scalar. The special filename '-' can b\
e used to read from\nstandard input (STDIN).\n\n  \
my $data = &read_file($filename);\n\n=cut\n\nsub r\
ead_file {\n    print_debug_message('read_file', '\
Begin', 1);\n    my $filename = shift;\n    print_\
debug_message('read_file', 'filename: ' . $filenam\
e, 2);\n    my ($content, $buffer);\n    if ($file\
name eq '-') {\n        while (sysread(STDIN, $buf\
fer, 1024)) {\n            $content .= $buffer;\n \
       }\n    }\n    else {\n        # File\n     \
   open(my $FILE, '<', $filename)\n            or \
die \"Error: unable to open input file $filename (\
$!)\";\n        while (sysread($FILE, $buffer, 102\
4)) {\n            $content .= $buffer;\n        }\
\n        close($FILE);\n    }\n    print_debug_me\
ssage('read_file', 'End', 1);\n    return $content\
;\n}\n\n=head2 write_file()\n\nWrite data to a fil\
e. The special filename '-' can be used to write t\
o\nstandard output (STDOUT).\n\n  &write_file($fil\
ename, $data);\n\n=cut\n\nsub write_file {\n    pr\
int_debug_message('write_file', 'Begin', 1);\n    \
my ($filename, $data) = @_;\n    print_debug_messa\
ge('write_file', 'filename: ' . $filename, 2);\n  \
  if ($outputLevel > 0) {\n        print STDERR 'C\
reating result file: ' . $filename . \"\\n\";\n   \
 }\n    if ($filename eq '-') {\n        print STD\
OUT $data;\n    }\n    else {\n        open(my $FI\
LE, '>', $filename)\n            or die \"Error: u\
nable to open output file $filename ($!)\";\n     \
   syswrite($FILE, $data);\n        close($FILE);\\
n    }\n    print_debug_message('write_file', 'End\
', 1);\n}\n\n=head2 usage()\n\nPrint program usage\
 message.\n\n  &usage();\n\n=cut\n\nsub usage {\n \
   print STDERR <<EOF\nEMBL-EBI NCBI Blast Perl Cl\
ient:\n\nSequence similarity search with NCBI Blas\
t.\n\n[Required (for job submission)]\n  --email  \
             E-mail address.\n  --program         \
    The BLAST program to be used for the Sequence \
Similarity\n                        Search.\n  --s\
type               Indicates if the sequence is pr\
otein or DNA/RNA.\n  --sequence            The que\
ry sequence can be entered directly into this form\
.\n                        The sequence can be in \
GCG, FASTA, EMBL (Nucleotide only),\n             \
           GenBank, PIR, NBRF, PHYLIP or UniProtKB\
/Swiss-Prot (Protein\n                        only\
) format. A partially formatted sequence is not\n \
                       accepted. Adding a return t\
o the end of the sequence may\n                   \
     help certain applications understand the inpu\
t. Note that\n                        directly usi\
ng data from word processors may yield\n          \
              unpredictable results as hidden/cont\
rol characters may be\n                        pre\
sent.\n  --database            Database.\n\n[Optio\
nal]\n  --task                Task option (only se\
lectable for blastn).\n  --matrix              (Pr\
otein searches) The substitution matrix used for s\
coring\n                        alignments when se\
arching the database.\n  --alignments          Max\
imum number of match alignments reported in the re\
sult\n                        output.\n  --scores \
             Maximum number of match score summari\
es reported in the\n                        result\
 output.\n  --exp                 Limits the numbe\
r of scores and alignments reported based on\n    \
                    the expectation value. This is\
 the maximum number of times\n                    \
    the match is expected to occur by chance.\n  -\
-dropoff             The amount a score can drop b\
efore gapped extension of word\n                  \
      hits is halted.\n  --match_scores        (Nu\
cleotide searches) The match score is the bonus to\
 the\n                        alignment score when\
 matching the same base. The mismatch is\n        \
                the penalty when failing to match.\
\n  --gapopen             Penalty taken away from \
the score when a gap is created in\n              \
          sequence. Increasing the gap openning pe\
nalty will decrease\n                        the n\
umber of gaps in the final alignment.\n  --gapext \
             Penalty taken away from the score for\
 each base or residue\n                        in \
the gap. Increasing the gap extension penalty favo\
rs\n                        short gaps in the fina\
l alignment, conversly decreasing the\n           \
             gap extension penalty favors long gap\
s in the final\n                        alignment.\
\n  --filter              Filter regions of low se\
quence complexity. This can avoid\n               \
         issues with low complexity sequences wher\
e matches are found\n                        due t\
o composition rather than meaningful sequence\n   \
                     similarity. However in some c\
ases filtering also masks\n                       \
 regions of interest and so should be used with ca\
ution.\n  --seqrange            Specify a range or\
 section of the input sequence to use in\n        \
                the search. Example: Specifying '3\
4-89' in an input sequence\n                      \
  of total length 100, will tell BLAST to only use\
 residues 34\n                        to 89, inclu\
sive.\n  --gapalign            This is a true/fals\
e setting that tells the program the\n            \
            perform optimised alignments within re\
gions involving gaps.\n                        If \
set to true, the program will perform an alignment\
 using\n                        gaps. Otherwise, i\
f it is set to false, it will report only\n       \
                 individual HSP where two sequence\
 match each other, and thus\n                     \
   will not produce alignments with gaps.\n  --wor\
dsize            Word size for wordfinder algorith\
m.\n  --compstats           Use composition-based \
statistics.\n  --align               Formating for\
 the alignments.\n  --transltable         Query Ge\
netic code to use in translation.\n\n[General]\n  \
-h, --help            Show this help message and e\
xit.\n  --asyncjob            Forces to make an as\
ynchronous query.\n  --title               Title f\
or job.\n  --status              Get job status.\n\
  --resultTypes         Get available result types\
 for job.\n  --polljob             Poll for the st\
atus of a job.\n  --pollFreq            Poll frequ\
ency in seconds (default 3s).\n  --jobid          \
     JobId that was returned when an asynchronous \
job was submitted.\n  --outfile             File n\
ame for results (default is JobId; for STDOUT).\n \
 --multifasta          Treat input as a set of fas\
ta formatted sequences.\n  --useSeqId            U\
se sequence identifiers for output filenames.\n   \
                     Only available in multi-fasta\
 and multi-identifier modes.\n  --maxJobs         \
    Maximum number of concurrent jobs. Only\n     \
                   available in multifasta or list\
 file modes.\n  --outformat           Result forma\
t(s) to retrieve. It accepts comma-separated value\
s.\n  --params              List input parameters.\
\n  --paramDetail         Display details for inpu\
t parameter.\n  --quiet               Decrease out\
put.\n  --verbose             Increase output.\n  \
--version             Prints out the version of th\
e Client and exit.\n  --baseUrl             Base U\
RL. Defaults to:\n                        https://\
www.ebi.ac.uk/Tools/services/rest/ncbiblast\n\nSyn\
chronous job:\n  The results/errors are returned a\
s soon as the job is finished.\n  Usage: perl $scr\
iptName --email <your\\@email.com> [options...] <S\
eqFile|SeqID(s)>\n  Returns: results as an attachm\
ent\n\nAsynchronous job:\n  Use this if you want t\
o retrieve the results at a later time. The result\
s\n  are stored for up to 24 hours.\n  Usage: perl\
 $scriptName --asyncjob --email <your\\@email.com>\
 [options...] <SeqFile|SeqID(s)>\n  Returns: jobid\
\n\nCheck status of Asynchronous job:\n  Usage: pe\
rl $scriptName --status --jobid <jobId>\n\nRetriev\
e job data:\n  Use the jobid to query for the stat\
us of the job. If the job is finished,\n  it also \
returns the results/errors.\n  Usage: perl $script\
Name --polljob --jobid <jobId> [--outfile string]\\
n  Returns: string indicating the status of the jo\
b and if applicable, results\n  as an attachment.\\
n\nFurther information:\n  https://www.ebi.ac.uk/T\
ools/webservices and\n    https://github.com/ebi-w\
p/webservice-clients\n\nSupport/Feedback:\n  https\
://www.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FEEDBA\
CK/SUPPORT\n\nPlease contact us at L<https://www.e\
bi.ac.uk/support/> if you have any\nfeedback, sugg\
estions or issues with the service or this client.\
\n\n=cut\n","\n=head1 NAME\n\nwublast_lwp.pl\n\n=h\
ead1 DESCRIPTION\n\nWU-BLAST (REST) web service Pe\
rl client using L<LWP>.\n\nTested with:\n\n=over\n\
\n=item *\nL<LWP> 5.79, L<XML::Simple> 2.12 and Pe\
rl 5.8.3\n\n=item *\nL<LWP> 5.808, L<XML::Simple> \
2.18 and Perl 5.8.8 (Ubuntu 8.04 LTS)\n\n=item *\n\
L<LWP> 5.834, L<XML::Simple> 2.18 and Perl 5.10.1 \
(Ubuntu 10.04 LTS)\n\n=item *\nL<LWP> 6.03, L<XML:\
:Simple> 2.18 and Perl 5.14.2 (Ubuntu 12.04 LTS)\n\
\n=back\n\nFor further information see:\n\n=over\n\
\n=item *\nL<http://www.ebi.ac.uk/Tools/webservice\
s/services/sss/wu_blast_rest>\n\n=item *\nL<http:/\
/www.ebi.ac.uk/Tools/webservices/tutorials/perl>\n\
\n=back\n\n=head1 LICENSE\n\nCopyright 2012-2013 E\
MBL - European Bioinformatics Institute\n\nLicense\
d under the Apache License, Version 2.0 (the \"Lic\
ense\");\nyou may not use this file except in comp\
liance with the License.\nYou may obtain a copy of\
 the License at\n\n    http://www.apache.org/licen\
ses/LICENSE-2.0\n\nUnless required by applicable l\
aw or agreed to in writing, software\ndistributed \
under the License is distributed on an \"AS IS\" B\
ASIS,\nWITHOUT WARRANTIES OR CONDITIONS OF ANY KIN\
D, either express or implied.\nSee the License for\
 the specific language governing permissions and\n\
limitations under the License.\n\n=head1 VERSION\n\
\n$Id: wublast_lwp.pl 2560 2013-03-20 12:56:31Z hp\
m $\n\n=cut\n\nuse strict;\nuse warnings;\n\nuse E\
nglish;\nuse LWP;\nuse XML::Simple;\nuse Getopt::L\
ong qw(:config no_ignore_case bundling);\nuse File\
::Basename;\nuse Data::Dumper;\n\nmy $baseUrl = 'h\
ttp://www.ebi.ac.uk/Tools/services/rest/wublast';\\
n\nmy $checkInterval = 3;\n\nmy $outputLevel = 1;\\
n\nmy $numOpts = scalar(@ARGV);\nmy %params = ( 'd\
ebugLevel' => 0 );\n\nmy %tool_params = ();\nGetOp\
tions(\n\n	# Tool specific options\n	'program|p=s'\
     => \\$tool_params{'program'},      # BLAST pr\
ogram\n	'database|D=s'    => \\$params{'database'}\
,     # Search database\n	'matrix|m=s'      => \\$\
tool_params{'matrix'},       # Scoring matrix\n	'e\
xp|E=f'         => \\$tool_params{'exp'},         \
 # E-value threshold\n	'viewfilter|e'    => \\$too\
l_params{'viewfilter'},   # Display filtered seque\
nce\n	'filter|f=s'      => \\$tool_params{'filter'\
},       # Low complexity filter name\n	'alignment\
s|n=i'  => \\$tool_params{'alignments'},   # Numbe\
r of alignments\n	'scores|s=i'      => \\$tool_par\
ams{'scores'},       # Number of scores\n	'sensiti\
vity|S=s' => \\$tool_params{'sensitivity'},  # Sea\
rch sensitivity\n	'sort|t=s'        => \\$tool_par\
ams{'sort'},         # Sort hits by...\n	'stats|T=\
s'       => \\$tool_params{'stats'},        # Scor\
ing statistic to use\n	'strand|d=s'      => \\$too\
l_params{'strand'},       # Strand to use\n	'topco\
mbon|c=i'   => \\$tool_params{'topcombon'},    # C\
onsistent sets of HSPs\n	'align|A=i'       => \\$t\
ool_params{'align'},   # Pairwise alignment format\
\n	'stype=s' => \\$tool_params{'stype'},    # Sequ\
ence type 'protein' or 'dna'\n	'sequence=s' => \\$\
params{'sequence'},         # Query sequence file \
or DB:ID\n	'multifasta' => \\$params{'multifasta'}\
,       # Multiple fasta input\n\n	# Compatability\
 options, old command-line.\n	'echofilter|e'    =>\
 \\$params{'echofilter'},   # Display filtered seq\
uence\n	'b=i'  => \\$params{'numal'},        # Num\
ber of alignments\n	'appxml=s'        => \\$params\
{'appxml'},       # Application XML\n\n	# Generic \
options\n	'email=s'       => \\$params{'email'},  \
        # User e-mail address\n	'title=s'       =>\
 \\$params{'title'},          # Job title\n	'outfi\
le=s'     => \\$params{'outfile'},        # Output\
 file name\n	'outformat=s'   => \\$params{'outform\
at'},      # Output file type\n	'jobid=s'       =>\
 \\$params{'jobid'},          # JobId\n	'help|h'  \
      => \\$params{'help'},           # Usage help\
\n	'async'         => \\$params{'async'},         \
 # Asynchronous submission\n	'polljob'       => \\\
$params{'polljob'},        # Get results\n	'result\
Types'   => \\$params{'resultTypes'},    # Get res\
ult types\n	'status'        => \\$params{'status'}\
,         # Get status\n	'params'        => \\$par\
ams{'params'},         # List input parameters\n	'\
paramDetail=s' => \\$params{'paramDetail'},    # G\
et details for parameter\n	'quiet'         => \\$p\
arams{'quiet'},          # Decrease output level\n\
	'verbose'       => \\$params{'verbose'},        #\
 Increase output level\n	'debugLevel=i'  => \\$par\
ams{'debugLevel'},     # Debug output level\n	'bas\
eUrl=s'     => \\$baseUrl,                  # Base\
 URL for service.\n);\nif ( $params{'verbose'} ) {\
 $outputLevel++ }\nif ( $params{'quiet'} )  { $out\
putLevel-- }\n\n&print_debug_message( 'MAIN', 'LWP\
::VERSION: ' . $LWP::VERSION,\n	1 );\n\n&print_deb\
ug_message( 'MAIN', \"params:\\n\" . Dumper( \\%pa\
rams ),           11 );\n&print_debug_message( 'MA\
IN', \"tool_params:\\n\" . Dumper( \\%tool_params \
), 11 );\n\nmy $ua;\n\nmy $scriptName = basename( \
$0, () );\n\nif ( $params{'help'} || $numOpts == 0\
 ) {\n	&usage();\n	exit(0);\n}\n\n&print_debug_mes\
sage( 'MAIN', 'baseUrl: ' . $baseUrl, 1 );\n\nif (\
\n	!(\n		   $params{'polljob'}\n		|| $params{'resu\
ltTypes'}\n		|| $params{'status'}\n		|| $params{'p\
arams'}\n		|| $params{'paramDetail'}\n	)\n	&& !( d\
efined( $ARGV[0] ) || defined( $params{'sequence'}\
 ) )\n  )\n{\n\n	# Bad argument combination, so pr\
int error message and usage\n	print STDERR 'Error:\
 bad option combination', \"\\n\";\n	&usage();\n	e\
xit(1);\n}\n\nelsif ( $params{'params'} ) {\n	&pri\
nt_tool_params();\n}\n\nelsif ( $params{'paramDeta\
il'} ) {\n	&print_param_details( $params{'paramDet\
ail'} );\n}\n\nelsif ( $params{'status'} && define\
d( $params{'jobid'} ) ) {\n	&print_job_status( $pa\
rams{'jobid'} );\n}\n\nelsif ( $params{'resultType\
s'} && defined( $params{'jobid'} ) ) {\n	&print_re\
sult_types( $params{'jobid'} );\n}\n\nelsif ( $par\
ams{'polljob'} && defined( $params{'jobid'} ) ) {\\
n	&get_results( $params{'jobid'} );\n}\n\nelse {\n\
\n	# Multiple input sequence mode, assume fasta fo\
rmat.\n	if ( $params{'multifasta'} ) {\n		&multi_s\
ubmit_job();\n	}\n\n	# Entry identifier list file.\
\n	elsif (( defined( $params{'sequence'} ) && $par\
ams{'sequence'} =~ m/^\\@/ )\n		|| ( defined( $ARG\
V[0] ) && $ARGV[0] =~ m/^\\@/ ) )\n	{\n		my $list_\
filename = $params{'sequence'} || $ARGV[0];\n		$li\
st_filename =~ s/^\\@//;\n		&list_file_submit_job(\
$list_filename);\n	}\n\n	# Default: single sequenc\
e/identifier.\n	else {\n\n		# Load the sequence da\
ta and submit.\n		&submit_job( &load_data() );\n	}\
\n}\n\n=head1 FUNCTIONS\n\n=cut\n\n\n=head2 rest_u\
ser_agent()\n\nGet a LWP UserAgent to use to perfo\
rm REST requests.\n\n  my $ua = &rest_user_agent()\
;\n\n=cut\n\nsub rest_user_agent() {\n	print_debug\
_message( 'rest_user_agent', 'Begin', 21 );\n	# Cr\
eate an LWP UserAgent for making HTTP calls.\n	my \
$ua = LWP::UserAgent->new();\n	# Set 'User-Agent' \
HTTP header to identifiy the client.\n	'$Revision:\
 2560 $' =~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-C\
lient/$1 ($scriptName; $OSNAME) \" . $ua->agent())\
;\n	# Configure HTTP proxy support from environmen\
t.\n	$ua->env_proxy;\n	print_debug_message( 'rest_\
user_agent', 'End', 21 );\n	return $ua;\n}\n\n=hea\
d2 rest_error()\n\nCheck a REST response for an er\
ror condition. An error is mapped to a die.\n\n  &\
rest_error($response, $content_data);\n\n=cut\n\ns\
ub rest_error() {\n	print_debug_message( 'rest_err\
or', 'Begin', 21 );\n	my $response = shift;\n	my $\
contentdata;\n	if(scalar(@_) > 0) {\n		$contentdat\
a = shift;\n	}\n	if(!defined($contentdata) || $con\
tentdata eq '') {\n		$contentdata = $response->con\
tent();\n	}\n	# Check for HTTP error codes\n	if ( \
$response->is_error ) {\n		my $error_message = '';\
\n		# HTML response.\n		if(	$contentdata =~ m/<h1>\
([^<]+)<\\/h1>/ ) {\n			$error_message = $1;\n		}\\
n		#  XML response.\n		elsif($contentdata =~ m/<de\
scription>([^<]+)<\\/description>/) {\n			$error_m\
essage = $1;\n		}\n		die 'http status: ' . $respon\
se->code . ' ' . $response->message . '  ' . $erro\
r_message;\n	}\n	print_debug_message( 'rest_error'\
, 'End', 21 );\n}\n\n=head2 rest_request()\n\nPerf\
orm a REST request (HTTP GET).\n\n  my $response_s\
tr = &rest_request($url);\n\n=cut\n\nsub rest_requ\
est {\n	print_debug_message( 'rest_request', 'Begi\
n', 11 );\n	my $requestUrl = shift;\n	print_debug_\
message( 'rest_request', 'URL: ' . $requestUrl, 11\
 );\n\n	# Get an LWP UserAgent.\n	$ua = &rest_user\
_agent() unless defined($ua);\n	# Available HTTP c\
ompression methods.\n	my $can_accept;\n	eval {\n	 \
   $can_accept = HTTP::Message::decodable();\n	};\\
n	$can_accept = '' unless defined($can_accept);\n	\
# Perform the request\n	my $response = $ua->get($r\
equestUrl,\n		'Accept-Encoding' => $can_accept, # \
HTTP compression.\n	);\n	print_debug_message( 'res\
t_request', 'HTTP status: ' . $response->code,\n		\
11 );\n	print_debug_message( 'rest_request',\n		'r\
esponse length: ' . length($response->content()), \
11 );\n	print_debug_message( 'rest_request',\n		'r\
equest:' .\"\\n\" . $response->request()->as_strin\
g(), 32 );\n	print_debug_message( 'rest_request',\\
n		'response: ' . \"\\n\" . $response->as_string()\
, 32 );\n	# Unpack possibly compressed response.\n\
	my $retVal;\n	if ( defined($can_accept) && $can_a\
ccept ne '') {\n	    $retVal = $response->decoded_\
content();\n	}\n	# If unable to decode use orginal\
 content.\n	$retVal = $response->content() unless \
defined($retVal);\n	# Check for an error.\n	&rest_\
error($response, $retVal);\n	print_debug_message( \
'rest_request', 'retVal: ' . $retVal, 12 );\n	prin\
t_debug_message( 'rest_request', 'End', 11 );\n\n	\
# Return the response data\n	return $retVal;\n}\n\\
n=head2 rest_get_parameters()\n\nGet list of tool \
parameter names.\n\n  my (@param_list) = &rest_get\
_parameters();\n\n=cut\n\nsub rest_get_parameters \
{\n	print_debug_message( 'rest_get_parameters', 'B\
egin', 1 );\n	my $url                = $baseUrl . \
'/parameters/';\n	my $param_list_xml_str = rest_re\
quest($url);\n	my $param_list_xml     = XMLin($par\
am_list_xml_str);\n	my (@param_list)       = @{ $p\
aram_list_xml->{'id'} };\n	print_debug_message( 'r\
est_get_parameters', 'End', 1 );\n	return (@param_\
list);\n}\n\n=head2 rest_get_parameter_details()\n\
\nGet details of a tool parameter.\n\n  my $paramD\
etail = &rest_get_parameter_details($param_name);\\
n\n=cut\n\nsub rest_get_parameter_details {\n	prin\
t_debug_message( 'rest_get_parameter_details', 'Be\
gin', 1 );\n	my $parameterId = shift;\n	print_debu\
g_message( 'rest_get_parameter_details',\n		'param\
eterId: ' . $parameterId, 1 );\n	my $url          \
        = $baseUrl . '/parameterdetails/' . $param\
eterId;\n	my $param_detail_xml_str = rest_request(\
$url);\n	my $param_detail_xml     = XMLin($param_d\
etail_xml_str);\n	print_debug_message( 'rest_get_p\
arameter_details', 'End', 1 );\n	return ($param_de\
tail_xml);\n}\n\n=head2 rest_run()\n\nSubmit a job\
.\n\n  my $job_id = &rest_run($email, $title, \\%p\
arams );\n\n=cut\n\nsub rest_run {\n	print_debug_m\
essage( 'rest_run', 'Begin', 1 );\n	my $email  = s\
hift;\n	my $title  = shift;\n	my $params = shift;\\
n	print_debug_message( 'rest_run', 'email: ' . $em\
ail, 1 );\n	if ( defined($title) ) {\n		print_debu\
g_message( 'rest_run', 'title: ' . $title, 1 );\n	\
}\n	print_debug_message( 'rest_run', 'params: ' . \
Dumper($params), 1 );\n\n	# Get an LWP UserAgent.\\
n	$ua = &rest_user_agent() unless defined($ua);\n\\
n	# Clean up parameters\n	my (%tmp_params) = %{$pa\
rams};\n	$tmp_params{'email'} = $email;\n	$tmp_par\
ams{'title'} = $title;\n	foreach my $param_name ( \
keys(%tmp_params) ) {\n		if ( !defined( $tmp_param\
s{$param_name} ) ) {\n			delete $tmp_params{$param\
_name};\n		}\n	}\n\n	# Submit the job as a POST\n	\
my $url = $baseUrl . '/run';\n	my $response = $ua-\
>post( $url, \\%tmp_params );\n	print_debug_messag\
e( 'rest_run', 'HTTP status: ' . $response->code, \
11 );\n	print_debug_message( 'rest_run',\n		'reque\
st:' .\"\\n\" . $response->request()->as_string(),\
 11 );\n	print_debug_message( 'rest_run',\n		'resp\
onse: ' . length($response->as_string()) . \"\\n\"\
 . $response->as_string(), 11 );\n\n	# Check for a\
n error.\n	&rest_error($response);\n\n	# The job i\
d is returned\n	my $job_id = $response->content();\
\n	print_debug_message( 'rest_run', 'End', 1 );\n	\
return $job_id;\n}\n\n=head2 rest_get_status()\n\n\
Check the status of a job.\n\n  my $status = &rest\
_get_status($job_id);\n\n=cut\n\nsub rest_get_stat\
us {\n	print_debug_message( 'rest_get_status', 'Be\
gin', 1 );\n	my $job_id = shift;\n	print_debug_mes\
sage( 'rest_get_status', 'jobid: ' . $job_id, 2 );\
\n	my $status_str = 'UNKNOWN';\n	my $url        = \
$baseUrl . '/status/' . $job_id;\n	$status_str = &\
rest_request($url);\n	print_debug_message( 'rest_g\
et_status', 'status_str: ' . $status_str, 2 );\n	p\
rint_debug_message( 'rest_get_status', 'End', 1 );\
\n	return $status_str;\n}\n\n=head2 rest_get_resul\
t_types()\n\nGet list of result types for finished\
 job.\n\n  my (@result_types) = &rest_get_result_t\
ypes($job_id);\n\n=cut\n\nsub rest_get_result_type\
s {\n	print_debug_message( 'rest_get_result_types'\
, 'Begin', 1 );\n	my $job_id = shift;\n	print_debu\
g_message( 'rest_get_result_types', 'jobid: ' . $j\
ob_id, 2 );\n	my (@resultTypes);\n	my $url        \
              = $baseUrl . '/resulttypes/' . $job_\
id;\n	my $result_type_list_xml_str = &rest_request\
($url);\n	my $result_type_list_xml     = XMLin($re\
sult_type_list_xml_str);\n	(@resultTypes) = @{ $re\
sult_type_list_xml->{'type'} };\n	print_debug_mess\
age( 'rest_get_result_types',\n		scalar(@resultTyp\
es) . ' result types', 2 );\n	print_debug_message(\
 'rest_get_result_types', 'End', 1 );\n	return (@r\
esultTypes);\n}\n\n=head2 rest_get_result()\n\nGet\
 result data of a specified type for a finished jo\
b.\n\n  my $result = rest_get_result($job_id, $res\
ult_type);\n\n=cut\n\nsub rest_get_result {\n	prin\
t_debug_message( 'rest_get_result', 'Begin', 1 );\\
n	my $job_id = shift;\n	my $type   = shift;\n	prin\
t_debug_message( 'rest_get_result', 'jobid: ' . $j\
ob_id, 1 );\n	print_debug_message( 'rest_get_resul\
t', 'type: ' . $type,    1 );\n	my $url    = $base\
Url . '/result/' . $job_id . '/' . $type;\n	my $re\
sult = &rest_request($url);\n	print_debug_message(\
 'rest_get_result', length($result) . ' characters\
',\n		1 );\n	print_debug_message( 'rest_get_result\
', 'End', 1 );\n	return $result;\n}\n\n\n=head2 pr\
int_debug_message()\n\nPrint debug message at spec\
ified debug level.\n\n  &print_debug_message($meth\
od_name, $message, $level);\n\n=cut\n\nsub print_d\
ebug_message {\n	my $function_name = shift;\n	my $\
message       = shift;\n	my $level         = shift\
;\n	if ( $level <= $params{'debugLevel'} ) {\n		pr\
int STDERR '[', $function_name, '()] ', $message, \
\"\\n\";\n	}\n}\n\n=head2 print_tool_params()\n\nP\
rint list of tool parameters.\n\n  &print_tool_par\
ams();\n\n=cut\n\nsub print_tool_params {\n	print_\
debug_message( 'print_tool_params', 'Begin', 1 );\\
n	my (@param_list) = &rest_get_parameters();\n	for\
each my $param ( sort(@param_list) ) {\n		print $p\
aram, \"\\n\";\n	}\n	print_debug_message( 'print_t\
ool_params', 'End', 1 );\n}\n\n=head2 print_param_\
details()\n\nPrint details of a tool parameter.\n\\
n  &print_param_details($param_name);\n\n=cut\n\ns\
ub print_param_details {\n	print_debug_message( 'p\
rint_param_details', 'Begin', 1 );\n	my $paramName\
 = shift;\n	print_debug_message( 'print_param_deta\
ils', 'paramName: ' . $paramName, 2 );\n	my $param\
Detail = &rest_get_parameter_details($paramName);\\
n	print $paramDetail->{'name'}, \"\\t\", $paramDet\
ail->{'type'}, \"\\n\";\n	print $paramDetail->{'de\
scription'}, \"\\n\";\n	if(defined($paramDetail->{\
'values'}->{'value'})) {\n		if(ref($paramDetail->{\
'values'}->{'value'}) eq 'ARRAY') {\n			foreach my\
 $value ( @{ $paramDetail->{'values'}->{'value'} }\
 ) {\n				&print_param_value($value);\n			}\n		}\n\
		else {\n				&print_param_value($paramDetail->{'v\
alues'}->{'value'});\n		}\n	}\n	print_debug_messag\
e( 'print_param_details', 'End', 1 );\n}\n\n=head2\
 print_param_value()\n\nPrint details of a tool pa\
rameter value.\n\n  &print_param_details($param_va\
lue);\n\nUsed by print_param_details() to handle b\
oth singluar and array values.\n\n=cut\n\nsub prin\
t_param_value {\n	my $value = shift;\n	print $valu\
e->{'value'};\n	if ( $value->{'defaultValue'} eq '\
true' ) {\n		print \"\\t\", 'default';\n	}\n	print\
 \"\\n\";\n	print \"\\t\", $value->{'label'}, \"\\\
n\";\n	if ( defined( $value->{'properties'} ) ) {\\
n		foreach\n		  my $key ( sort( keys( %{ $value->{\
'properties'}{'property'} } ) ) )\n		{\n			if ( re\
f( $value->{'properties'}{'property'}{$key} ) eq '\
HASH'\n				&& defined( $value->{'properties'}{'pro\
perty'}{$key}{'value'} )\n			  )\n			{\n				print \
\"\\t\", $key, \"\\t\",\n				  $value->{'propertie\
s'}{'property'}{$key}{'value'}, \"\\n\";\n			}\n		\
	else {\n				print \"\\t\", $value->{'properties'}\
{'property'}{'key'},\n				  \"\\t\", $value->{'pro\
perties'}{'property'}{'value'}, \"\\n\";\n				last\
;\n			}\n		}\n	}\n}\n\n=head2 print_job_status()\n\
\nPrint status of a job.\n\n  &print_job_status($j\
ob_id);\n\n=cut\n\nsub print_job_status {\n	print_\
debug_message( 'print_job_status', 'Begin', 1 );\n\
	my $jobid = shift;\n	print_debug_message( 'print_\
job_status', 'jobid: ' . $jobid, 1 );\n	if ( $outp\
utLevel > 0 ) {\n		print STDERR 'Getting status fo\
r job ', $jobid, \"\\n\";\n	}\n	my $result = &rest\
_get_status($jobid);\n	print \"$result\\n\";\n	if \
( $result eq 'FINISHED' && $outputLevel > 0 ) {\n	\
	print STDERR \"To get results: $scriptName --poll\
job --jobid \" . $jobid\n		  . \"\\n\";\n	}\n	prin\
t_debug_message( 'print_job_status', 'End', 1 );\n\
}\n\n=head2 print_result_types()\n\nPrint availabl\
e result types for a job.\n\n  &print_result_types\
($job_id);\n\n=cut\n\nsub print_result_types {\n	p\
rint_debug_message( 'result_types', 'Begin', 1 );\\
n	my $jobid = shift;\n	print_debug_message( 'resul\
t_types', 'jobid: ' . $jobid, 1 );\n	if ( $outputL\
evel > 0 ) {\n		print STDERR 'Getting result types\
 for job ', $jobid, \"\\n\";\n	}\n	my $status = &r\
est_get_status($jobid);\n	if ( $status eq 'PENDING\
' || $status eq 'RUNNING' ) {\n		print STDERR 'Err\
or: Job status is ', $status,\n		  '. To get resul\
t types the job must be finished.', \"\\n\";\n	}\n\
	else {\n		my (@resultTypes) = &rest_get_result_ty\
pes($jobid);\n		if ( $outputLevel > 0 ) {\n			prin\
t STDOUT 'Available result types:', \"\\n\";\n		}\\
n		foreach my $resultType (@resultTypes) {\n			pri\
nt STDOUT $resultType->{'identifier'}, \"\\n\";\n	\
		if ( defined( $resultType->{'label'} ) ) {\n				\
print STDOUT \"\\t\", $resultType->{'label'}, \"\\\
n\";\n			}\n			if ( defined( $resultType->{'descri\
ption'} ) ) {\n				print STDOUT \"\\t\", $resultTy\
pe->{'description'}, \"\\n\";\n			}\n			if ( defin\
ed( $resultType->{'mediaType'} ) ) {\n				print ST\
DOUT \"\\t\", $resultType->{'mediaType'}, \"\\n\";\
\n			}\n			if ( defined( $resultType->{'fileSuffix\
'} ) ) {\n				print STDOUT \"\\t\", $resultType->{\
'fileSuffix'}, \"\\n\";\n			}\n		}\n		if ( $status\
 eq 'FINISHED' && $outputLevel > 0 ) {\n			print S\
TDERR \"\\n\", 'To get results:', \"\\n\",\n			  \\
"  $scriptName --polljob --jobid \" . $params{'job\
id'} . \"\\n\",\n			  \"  $scriptName --polljob --\
outformat <type> --jobid \"\n			  . $params{'jobid\
'} . \"\\n\";\n		}\n	}\n	print_debug_message( 'res\
ult_types', 'End', 1 );\n}\n\n=head2 submit_job()\\
n\nSubmit a job to the service.\n\n  &submit_job($\
seq);\n\n=cut\n\nsub submit_job {\n	print_debug_me\
ssage( 'submit_job', 'Begin', 1 );\n\n	# Set input\
 sequence\n	$tool_params{'sequence'} = shift;\n\n	\
# Load parameters\n	&load_params();\n\n	# Submit t\
he job\n	my $jobid = &rest_run( $params{'email'}, \
$params{'title'}, \\%tool_params );\n\n	# Simulate\
 sync/async mode\n	if ( defined( $params{'async'} \
) ) {\n		print STDOUT $jobid, \"\\n\";\n		if ( $ou\
tputLevel > 0 ) {\n			print STDERR\n			  \"To chec\
k status: $scriptName --status --jobid $jobid\\n\"\
;\n		}\n	}\n	else {\n		if ( $outputLevel > 0 ) {\n\
			print STDERR \"JobId: $jobid\\n\";\n		}\n		slee\
p 1;\n		&get_results($jobid);\n	}\n	print_debug_me\
ssage( 'submit_job', 'End', 1 );\n}\n\n=head2 mult\
i_submit_job()\n\nSubmit multiple jobs assuming in\
put is a collection of fasta formatted sequences.\\
n\n  &multi_submit_job();\n\n=cut\n\nsub multi_sub\
mit_job {\n	print_debug_message( 'multi_submit_job\
', 'Begin', 1 );\n	my $jobIdForFilename = 1;\n	$jo\
bIdForFilename = 0 if ( defined( $params{'outfile'\
} ) );\n	my (@filename_list) = ();\n\n	# Query seq\
uence\n	if ( defined( $ARGV[0] ) ) {    # Bare opt\
ion\n		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    \
# File\n			push( @filename_list, $ARGV[0] );\n		}\\
n		else {\n			warn 'Warning: Input file \"' . $ARG\
V[0] . '\" does not exist'\n		}\n	}\n	if ( $params\
{'sequence'} ) {                   # Via --sequenc\
e\n		if ( -f $params{'sequence'} || $params{'seque\
nce'} eq '-' ) {    # File\n			push( @filename_lis\
t, $params{'sequence'} );\n		}\n		else {\n			warn \
'Warning: Input file \"' . $params{'sequence'} . '\
\" does not exist'\n		}\n	}\n\n	$/ = '>';\n	foreac\
h my $filename (@filename_list) {\n		my $INFILE;\n\
		if($filename eq '-') { # STDIN.\n			open( $INFIL\
E, '<-' )\n			  or die 'Error: unable to STDIN (' \
. $! . ')';\n		} else { # File.\n			open( $INFILE,\
 '<', $filename )\n			  or die 'Error: unable to o\
pen file ' . $filename . ' (' . $! . ')';\n		}\n		\
while (<$INFILE>) {\n			my $seq = $_;\n			$seq =~ \
s/>$//;\n			if ( $seq =~ m/(\\S+)/ ) {\n				print \
STDERR \"Submitting job for: $1\\n\"\n				  if ( $\
outputLevel > 0 );\n				$seq = '>' . $seq;\n				&p\
rint_debug_message( 'multi_submit_job', $seq, 11 )\
;\n				&submit_job($seq);\n				$params{'outfile'} \
= undef if ( $jobIdForFilename == 1 );\n			}\n		}\\
n		close $INFILE;\n	}\n	print_debug_message( 'mult\
i_submit_job', 'End', 1 );\n}\n\n=head2 list_file_\
submit_job()\n\nSubmit multiple jobs using a file \
containing a list of entry identifiers as \ninput.\
\n\n  &list_file_submit_job($list_filename)\n\n=cu\
t\n\nsub list_file_submit_job {\n	print_debug_mess\
age( 'list_file_submit_job', 'Begin', 11 );\n	my $\
filename         = shift;\n	my $jobIdForFilename =\
 1;\n	$jobIdForFilename = 0 if ( defined( $params{\
'outfile'} ) );\n\n	# Iterate over identifiers, su\
bmitting each job\n	my $LISTFILE;\n	if($filename e\
q '-') { # STDIN.\n		open( $LISTFILE, '<-' )\n		  \
or die 'Error: unable to STDIN (' . $! . ')';\n	} \
else { # File.\n		open( $LISTFILE, '<', $filename \
)\n		  or die 'Error: unable to open file ' . $fil\
ename . ' (' . $! . ')';\n	}\n	while (<$LISTFILE>)\
 {\n		my $line = $_;\n		chomp($line);\n		if ( $lin\
e ne '' ) {\n			&print_debug_message( 'list_file_s\
ubmit_job', 'line: ' . $line, 2 );\n			if ( $line \
=~ m/\\w:\\w/ ) {    # Check this is an identifier\
\n				print STDERR \"Submitting job for: $line\\n\\
"\n				  if ( $outputLevel > 0 );\n				&submit_job\
($line);\n			}\n			else {\n				print STDERR\n\"War\
ning: line \\\"$line\\\" is not recognised as an i\
dentifier\\n\";\n			}\n		}\n		$params{'outfile'} =\
 undef if ( $jobIdForFilename == 1 );\n	}\n	close \
$LISTFILE;\n	print_debug_message( 'list_file_submi\
t_job', 'End', 11 );\n}\n\n=head2 load_data()\n\nL\
oad sequence data from file or option specified on\
 the command-line.\n\n  &load_data();\n\n=cut\n\ns\
ub load_data {\n	print_debug_message( 'load_data',\
 'Begin', 1 );\n	my $retSeq;\n\n	# Query sequence\\
n	if ( defined( $ARGV[0] ) ) {    # Bare option\n	\
	if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File\
\n			$retSeq = &read_file( $ARGV[0] );\n		}\n		els\
e {                                     # DB:ID or\
 sequence\n			$retSeq = $ARGV[0];\n		}\n	}\n	if ( \
$params{'sequence'} ) {                   # Via --\
sequence\n		if ( -f $params{'sequence'} || $params\
{'sequence'} eq '-' ) {    # File\n			$retSeq = &r\
ead_file( $params{'sequence'} );\n		}\n		else {   \
 # DB:ID or sequence\n			$retSeq = $params{'sequen\
ce'};\n		}\n	}\n	print_debug_message( 'load_data',\
 'End', 1 );\n	return $retSeq;\n}\n\n=head2 load_p\
arams()\n\nLoad job parameters from command-line o\
ptions.\n\n  &load_params();\n\n=cut\n\nsub load_p\
arams {\n	print_debug_message( 'load_params', 'Beg\
in', 1 );\n\n	# Database(s) to search\n	my (@dbLis\
t) = split /[ ,]/, $params{'database'};\n	$tool_pa\
rams{'database'} = \\@dbList;\n\n	# Compatability \
options, old command-line.\n	if(!$tool_params{'vie\
wfilter'} && $params{'echofilter'}) {\n		$tool_par\
ams{'viewfilter'} = 'true';\n	}\n	if(!$tool_params\
{'alignments'} && $params{'numal'}) {\n		$tool_par\
ams{'alignments'} = $params{'numal'};\n	}\n	# TODO\
: set alignment format option to get NCBI BLAST XM\
L.\n	if($params{'appxml'}) {\n		$tool_params{'alig\
n'} = '';\n	}\n\n	print_debug_message( 'load_param\
s', 'End', 1 );\n}\n\n=head2 client_poll()\n\nClie\
nt-side job polling.\n\n  &client_poll($job_id);\n\
\n=cut\n\nsub client_poll {\n	print_debug_message(\
 'client_poll', 'Begin', 1 );\n	my $jobid  = shift\
;\n	my $status = 'PENDING';\n\n	my $errorCount = 0\
;\n	while ($status eq 'RUNNING'\n		|| $status eq '\
PENDING'\n		|| ( $status eq 'ERROR' && $errorCount\
 < 2 ) )\n	{\n		$status = rest_get_status($jobid);\
\n		print STDERR \"$status\\n\" if ( $outputLevel \
> 0 );\n		if ( $status eq 'ERROR' ) {\n			$errorCo\
unt++;\n		}\n		elsif ( $errorCount > 0 ) {\n			$er\
rorCount--;\n		}\n		if (   $status eq 'RUNNING'\n	\
		|| $status eq 'PENDING'\n			|| $status eq 'ERROR\
' )\n		{\n\n			# Wait before polling again.\n			sl\
eep $checkInterval;\n		}\n	}\n	print_debug_message\
( 'client_poll', 'End', 1 );\n	return $status;\n}\\
n\n=head2 get_results()\n\nGet the results for a j\
ob identifier.\n\n  &get_results($job_id);\n\n=cut\
\n\nsub get_results {\n	print_debug_message( 'get_\
results', 'Begin', 1 );\n	my $jobid = shift;\n	pri\
nt_debug_message( 'get_results', 'jobid: ' . $jobi\
d, 1 );\n\n	# Verbose\n	if ( $outputLevel > 1 ) {\\
n		print 'Getting results for job ', $jobid, \"\\n\
\";\n	}\n\n	# Check status, and wait if not finish\
ed\n	client_poll($jobid);\n\n	# Use JobId if outpu\
t file name is not defined\n	unless ( defined( $pa\
rams{'outfile'} ) ) {\n		$params{'outfile'} = $job\
id;\n	}\n\n	# Get list of data types\n	my (@result\
Types) = rest_get_result_types($jobid);\n\n	# Get \
the data and write it to a file\n	if ( defined( $p\
arams{'outformat'} ) ) {    # Specified data type\\
n		my $selResultType;\n		foreach my $resultType (@\
resultTypes) {\n			if ( $resultType->{'identifier'\
} eq $params{'outformat'} ) {\n				$selResultType \
= $resultType;\n			}\n		}\n		if ( defined($selResu\
ltType) ) {\n			my $result =\n			  rest_get_result\
( $jobid, $selResultType->{'identifier'} );\n			if\
 ( $params{'outfile'} eq '-' ) {\n				write_file( \
$params{'outfile'}, $result );\n			}\n			else {\n	\
			write_file(\n					$params{'outfile'} . '.'\n			\
		  . $selResultType->{'identifier'} . '.'\n					 \
 . $selResultType->{'fileSuffix'},\n					$result\n\
				);\n			}\n		}\n		else {\n			die 'Error: unknow\
n result format \"' . $params{'outformat'} . '\"';\
\n		}\n	}\n	else {    # Data types available\n		  \
    # Write a file for each output type\n		for my \
$resultType (@resultTypes) {\n			if ( $outputLevel\
 > 1 ) {\n				print STDERR 'Getting ', $resultType\
->{'identifier'}, \"\\n\";\n			}\n			my $result = \
rest_get_result( $jobid, $resultType->{'identifier\
'} );\n			if ( $params{'outfile'} eq '-' ) {\n				\
write_file( $params{'outfile'}, $result );\n			}\n\
			else {\n				write_file(\n					$params{'outfile'\
} . '.'\n					  . $resultType->{'identifier'} . '.\
'\n					  . $resultType->{'fileSuffix'},\n					$re\
sult\n				);\n			}\n		}\n	}\n	print_debug_message(\
 'get_results', 'End', 1 );\n}\n\n=head2 read_file\
()\n\nRead a file into a scalar. The special filen\
ame '-' can be used to read from \nstandard input \
(STDIN).\n\n  my $data = &read_file($filename);\n\\
n=cut\n\nsub read_file {\n	print_debug_message( 'r\
ead_file', 'Begin', 1 );\n	my $filename = shift;\n\
	print_debug_message( 'read_file', 'filename: ' . \
$filename, 2 );\n	my ( $content, $buffer );\n	if (\
 $filename eq '-' ) {\n		while ( sysread( STDIN, $\
buffer, 1024 ) ) {\n			$content .= $buffer;\n		}\n\
	}\n	else {    # File\n		open( my $FILE, '<', $fil\
ename )\n		  or die \"Error: unable to open input \
file $filename ($!)\";\n		while ( sysread( $FILE, \
$buffer, 1024 ) ) {\n			$content .= $buffer;\n		}\\
n		close($FILE);\n	}\n	print_debug_message( 'read_\
file', 'End', 1 );\n	return $content;\n}\n\n=head2\
 write_file()\n\nWrite data to a file. The special\
 filename '-' can be used to write to \nstandard o\
utput (STDOUT).\n\n  &write_file($filename, $data)\
;\n\n=cut\n\nsub write_file {\n	print_debug_messag\
e( 'write_file', 'Begin', 1 );\n	my ( $filename, $\
data ) = @_;\n	print_debug_message( 'write_file', \
'filename: ' . $filename, 2 );\n	if ( $outputLevel\
 > 0 ) {\n		print STDERR 'Creating result file: ' \
. $filename . \"\\n\";\n	}\n	if ( $filename eq '-'\
 ) {\n		print STDOUT $data;\n	}\n	else {\n		open( \
my $FILE, '>', $filename )\n		  or die \"Error: un\
able to open output file $filename ($!)\";\n		sysw\
rite( $FILE, $data );\n		close($FILE);\n	}\n	print\
_debug_message( 'write_file', 'End', 1 );\n}\n\n=h\
ead2 usage()\n\nPrint program usage message.\n\n  \
&usage();\n\n=cut\n\nsub usage {\n	print STDERR <<\
EOF\nWU-BLAST\n========\n   \nRapid sequence datab\
ase search programs utilizing the BLAST algorithm\\
n    \n[Required]\n\n  -p, --program      : str  :\
 BLAST program to use, see --paramDetail program\n\
  -D, --database     : str  : database(s) to searc\
h, space separated. See\n                         \
     --paramDetail database\n      --stype        \
: str  : query sequence type, see --paramDetail st\
ype\n  seqFile            : file : query sequence \
(\"-\" for STDIN, \\@filename for\n               \
               identifier list file)\n\n[Optional]\
\n\n  -m, --matrix       : str  : scoring matrix, \
see --paramDetail matrix\n  -e, --exp          : r\
eal : 0<E<= 1000. Statistical significance thresho\
ld \n                              for reporting d\
atabase sequence matches.\n  -e, --viewfilter   : \
     : display the filtered query sequence\n  -f, \
--filter       : str  : filter the query sequence \
for low complexity \n                             \
 regions, see --paramDetail filter\n  -A, --align \
       : int  : pairwise alignment format, see --p\
aramDetail align\n  -s, --scores       : int  : nu\
mber of scores to be reported\n  -b, --alignments \
  : int  : number of alignments to report\n  -S, -\
-sensitivity  : str  : sensitivity of the search, \
\n                              see --paramDetail \
sensitivity\n  -t, --sort	     : str  : sort order\
 for hits, see --paramDetail sort\n  -T, --stats  \
      : str  : statistical model, see --paramDetai\
l stats\n  -d, --strand       : str  : DNA strand \
to search with,\n                              see\
 --paramDetail strand\n  -c, --topcombon    : str \
 : consistent sets of HSPs\n      --multifasta   :\
      : treat input as a set of fasta formatted se\
quences\n\n[General]\n\n  -h, --help         :    \
  : prints this help text\n      --async        : \
     : forces to make an asynchronous query\n     \
 --email        : str  : e-mail address\n      --t\
itle        : str  : title for job\n      --status\
       :      : get job status\n      --resultType\
s  :      : get available result types for job\n  \
    --polljob      :      : poll for the status of\
 a job\n      --jobid        : str  : jobid that w\
as returned when an asynchronous job \n           \
                   was submitted.\n      --outfile\
      : str  : file name for results (default is j\
obid;\n                              \"-\" for STD\
OUT)\n      --outformat    : str  : result format \
to retrieve\n      --params       :      : list in\
put parameters\n      --paramDetail  : str  : disp\
lay details for input parameter\n      --quiet    \
    :      : decrease output\n      --verbose     \
 :      : increase output\n   \nSynchronous job:\n\
\n  The results/errors are returned as soon as the\
 job is finished.\n  Usage: $scriptName --email <y\
our\\@email> [options...] seqFile\n  Returns: resu\
lts as an attachment\n\nAsynchronous job:\n\n  Use\
 this if you want to retrieve the results at a lat\
er time. The results \n  are stored for up to 24 h\
ours. 	\n  Usage: $scriptName --async --email <you\
r\\@email> [options...] seqFile\n  Returns: jobid\\
n\n  Use the jobid to query for the status of the \
job. If the job is finished, \n  it also returns t\
he results/errors.\n  Usage: $scriptName --polljob\
 --jobid <jobId> [--outfile string]\n  Returns: st\
ring indicating the status of the job and if appli\
cable, results \n  as an attachment.\n\nFurther in\
formation:\n\n  http://www.ebi.ac.uk/Tools/webserv\
ices/services/sss/wu_blast_rest\n  http://www.ebi.\
ac.uk/Tools/webservices/tutorials/perl\n\nSupport/\
Feedback:\n\n  http://www.ebi.ac.uk/support/\nEOF\\
n}\n\n=head1 FEEDBACK/SUPPORT\n\nPlease contact us\
 at L<http://www.ebi.ac.uk/support/> if you have a\
ny \nfeedback, suggestions or issues with the serv\
ice or this client.\n\n=cut\n","\n\n\nmy $PROBTRES\
H = 0.3;# base pairs below this prob threshold wil\
l be ignored\nmy $WEIGHT = 100.0; # float!!\nmy $N\
UCALPH = \"ACGTUNRYMKSWHBVD\";\nuse vars qw($NUCAL\
PH $WEIGHT);\n\nmy $myname = basename($0);\n\nuse \
strict;\nuse warnings;\n\nuse File::Basename;\nuse\
 Getopt::Long;\nuse File::Glob ':glob';\nuse File:\
:Spec;\nuse File::Temp qw/ tempfile tempdir /;\n\n\
\n\n\nsub tcoffeelib_header($;$)\n{\n    my ($nseq\
, $fd) = @_;\n    if (! defined($fd)) {\n        $\
fd = *STDOUT;\n    }\n    printf $fd \"! TC_LIB_FO\
RMAT_01\\n\";\n    printf $fd \"%d\\n\", $nseq;\n}\
\n\n\nsub tcoffeelib_header_addseq($$;$)\n{\n    m\
y ($id, $seq, $fd) = @_;\n    if (! defined($fd)) \
{\n        $fd = *STDOUT;\n    }\n    printf $fd \\
"%s %d %s\\n\", $id, length($seq), $seq;\n}\n\n\ns\
ub tcoffeelib_comment($;$)\n{\n    my ($comment, $\
fd) = @_;\n    if (! defined($fd)) {\n        $fd \
= *STDOUT;\n    }\n    printf $fd \"!\" . $comment\
 . \"\\n\";\n}\n\n\nsub tcoffeelib_struct($$$;$)\n\
{\n    my ($nseq, $len, $bpm, $fd) = @_;\n\n    if\
 (! defined($fd)) {\n        $fd = *STDOUT;\n    }\
\n\n    # output basepair indices with fixed weigh\
t\n    printf $fd \"#%d %d\\n\", $nseq, $nseq;\n  \
  # output basepairs (only once) and with unit-off\
set\n    for (my $i=0; $i<$len; $i++) {\n        f\
or (my $j=$i+1; $j<$len; $j++) {\n            if (\
! defined($bpm->[$i][$j])) {\n                prin\
t STDERR \"ERROR: \\$bpm->[$i][$j] undefined\\n\";\
\n            }\n            if ($bpm->[$i][$j]>0)\
 {\n                print $fd $i+1;\n             \
   print $fd \" \";\n                print $fd $j+\
1;\n                print $fd \" \" . $bpm->[$i][$\
j] . \"\\n\";\n            }\n        }\n    }\n}\\
n\n\nsub tcoffeelib_footer(;$)\n{\n    my ($fd) = \
@_;\n    if (! defined($fd)) {\n        $fd = *STD\
OUT;\n    }\n    print $fd \"! SEQ_1_TO_N\\n\";\n}\
\n\n\n    \nsub plfold($$$)\n{    \n    my ($id, $\
seq, $probtresh) = @_;\n    my (@struct);# return\\
n    my ($templ, $fhtmp, $fnametmp, $cmd, $ctr, $w\
indow_size);\n    our $ntemp++;\n    \n    $templ \
= $myname . \".pid-\" . $$ .$ntemp .\".XXXXXX\";\n\
    ($fhtmp, $fnametmp) = tempfile($templ, UNLINK \
=> 1); \n    print $fhtmp \">$id\\n$seq\\n\";\n\n \
   # --- init basepair array\n    #\n    for (my $\
i=0; $i<length($seq); $i++) {\n        for (my $j=\
$i+1; $j<length($seq); $j++) {\n            $struc\
t[$i][$j]=0;\n        }\n    }\n\n\n    # --- call\
 rnaplfold and drop a readme\n    #\n    $window_s\
ize=(length($seq)<70)?length($seq):70;\n    $cmd =\
 \"RNAplfold -W $window_size < $fnametmp >/dev/nul\
l\";\n    system($cmd);\n    \n    if ($? != 0) {\\
n        printf STDERR \"ERROR: RNAplfold ($cmd) e\
xited with error status %d\\n\", $? >> 8;\n       \
 return;\n    }\n    #unlink($fnametmp);\n    my $\
fps = sprintf(\"%s_dp.ps\", $id); # check long nam\
e\n    \n    if (! -s $fps) {\n      {\n\n	$fps = \
sprintf(\"%s_dp.ps\", substr($id,0,12)); # check s\
hort name\n 	if (! -s $fps)\n	  {\n	    die(\"coul\
dn't find expected file $fps\\n\");\n	    return;\\
n	  }\n      }\n    }\n\n    \n    # --- read base\
 pairs from created postscript\n    #\n    open(FH\
, $fps);\n    while (my $line = <FH>) {\n        m\
y ($nti, $ntj, $prob);\n        chomp($line);     \
   \n        # line: bp bp sqrt-prob ubox\n       \
 my @match = ($line =~ m/^([0-9]+) +([0-9]+) +([0-\
9\\.]+) +ubox$/);\n        if (scalar(@match)) {\n\
            $nti=$1;\n            $ntj=$2;\n      \
      $prob=$3*$3;# prob stored as square root\n\n\
            if ($prob>$probtresh) {\n             \
   #printf STDERR \"\\$struct[$nti][$ntj] sqrtprob\
=$3 prob=$prob > $probtresh\\n\";\n               \
 $struct[$nti-1][$ntj-1] = $WEIGHT\n            }\\
n            # store with zero-offset\n        }\n\
    }\n    close(FH);\n\n    # remove or gzi posts\
cript\n    #\n    unlink($fps);\n    #\n    # or g\
zip\n    #$cmd = \"gzip -qf $fps\";\n    #system($\
cmd);\n    #if ($? != 0) {\n    #    printf STDERR\
 \"ERROR: gzip ($cmd) exited with error status %d\\
\n\", $? >> 8;\n    #}\n\n    return \\@struct;\n}\
\n\n\n\n\n\nsub rnaseqfmt($)\n{\n    my ($seq) = @\
_;\n    # remove gaps\n    $seq =~ s/-//g;\n    # \
uppercase RNA\n    $seq = uc($seq);\n    # T -> U\\
n    $seq =~ s/T/U/g;\n    # check for invalid cha\
raters\n    $_ = $seq;\n    s/[^$NUCALPH]//g;\n   \
 return $_;\n}\n\n\n\n\nsub usage(;$)\n{    \n    \
my ($errmsg) = @_;\n    if ($errmsg) {\n        pr\
int STDERR \"ERROR: $errmsg\\n\";\n    }\n    prin\
t STDERR << \"EOF\";\n$myname:\n Creates a T-Coffe\
e RNA structure library from RNAplfold prediction.\
\n See FIXME:citation\nUsage:\n $myname -in seq_fi\
le -out tcoffee_lib\nEOF\n    exit(1);\n}\n\nsub r\
ead_fasta_seq \n  {\n    my $f=$_[0];\n    my %hse\
q;\n    my (@seq, @com, @name);\n    my ($a, $s,$n\
seq);\n\n    open (F, $f);\n    while (<F>)\n     \
 {\n	$s.=$_;\n      }\n    close (F);\n\n    \n   \
 @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    @seq\
 =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>(\\S\
*)(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1;\n  \\
n    for ($a=0; $a<$nseq; $a++)\n      {\n	my $n=$\
name[$a];\n	my $s;\n	$hseq{$n}{name}=$n;\n	$s=$seq\
[$a];$s=~s/\\s//g;\n	\n	$hseq{$n}{seq}=$s;\n	$hseq\
{$n}{com}=$com[$a];\n      }\n    return %hseq;\n \
 }\n\n\n\n\n\n\n\nmy $fmsq = \"\";\nmy $flib = \"\\
";\nmy %OPTS;\nmy %seq;\nmy ($id, $nseq, $i);\nmy \
@nl;\n\nGetOptions(\"in=s\" => \\$fmsq, \"out=s\" \
=> \\$flib);\n\nif (! -s $fmsq) {\n    usage(\"emp\
ty or non-existant file \\\"$fmsq\\\"\")\n}\nif (l\
ength($flib)==0) {\n    usage(\"empty out-filename\
\")\n}\n\n\n\n\n\n\n%seq=read_fasta_seq($fmsq);\n\\
n\n@nl=keys(%seq);\n\n$nseq=$#nl+1;\nopen FD_LIB, \
\">$flib\" or die \"can't open $flib!\";\ntcoffeel\
ib_header($nseq, *FD_LIB);\nforeach $id (keys (%se\
q))\n  {\n    my ($seq, $fmtseq);\n    \n    $seq \
= $seq{$id}{seq};\n    \n    $fmtseq = rnaseqfmt($\
seq);# check here, formatting for folding importan\
t later\n    if (length($seq)!=length($fmtseq)) {\\
n        print STDERR \"ERROR: invalid sequence $i\
d is not an RNA sequence. read seq is: $seq\\n\";\\
n        exit\n      }\n   \n    tcoffeelib_header\
_addseq($id, uc($seq), *FD_LIB);\n  }\ntcoffeelib_\
comment(\"generated by $myname on \" . localtime()\
, *FD_LIB);\n\n\n\n$i=0;\nforeach $id (keys (%seq)\
)\n  {\n    my ($cleanid, $seq, $bpm);\n    $seq=$\
seq{$id}{seq};\n    $cleanid = $id;\n    $cleanid \
=~ s,[/ ],_,g;# needed for rnaplfold\n    $seq = r\
naseqfmt($seq);\n    \n    $bpm = plfold($cleanid,\
 rnaseqfmt($seq), $PROBTRESH);       \n    \n    t\
coffeelib_struct($i+1, length($seq), $bpm, *FD_LIB\
);\n    $i++;\n}\n\n\ntcoffeelib_footer(*FD_LIB);\\
nclose FD_LIB;\nexit (0);\n\n","\n\n\n\n\n$cmd=joi\
n ' ', @ARGV;\nif ($cmd=~/-infile=(\\S+)/){ $seqfi\
le=$1;}\nif ($cmd=~/-outfile=(\\S+)/){ $libfile=$1\
;}\n\n\n\n%s=read_fasta_seq ($seqfile);\n\nopen (F\
, \">$libfile\");\nforeach $name (keys (%s))\n  {\\
n    my $tclib=\"$name.RNAplfold_tclib\";\n    pri\
nt (F \">$name _F_ $tclib\\n\");\n    seq2RNAplfol\
d2tclib ($name, $s{$name}{seq}, $tclib);\n  }\nclo\
se (F);\nexit (EXIT_SUCCESS);\n\nsub seq2RNAplfold\
2tclib\n  {\n    my ($name, $seq, $tclib)=@_;\n   \
 my ($tmp);\n    $n++;\n    $tmp=\"tmp4seq2RNAplfo\
ld_tclib.$$.$n.pep\";\n    open (RF, \">$tmp\");\n\
    print (RF \">$name\\n$seq\\n\");\n    close (R\
F);\n    \n    system \"t_coffee -other_pg RNAplfo\
ld2tclib.pl -in=$tmp -out=$tclib\";\n    \n    unl\
ink ($tmp);\n    return $tclib;\n  }\n    \n    \n\
sub read_fasta_seq \n  {\n    my $f=@_[0];\n    my\
 %hseq;\n    my (@seq, @com, @name);\n    my ($a, \
$s,$nseq);\n\n    open (F, $f);\n    while (<F>)\n\
      {\n	$s.=$_;\n      }\n    close (F);\n\n    \
\n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n   \
 @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/\
>\\S*(.*)\\n([^>]*)/g);\n\n    \n    $nseq=$#name+\
1;\n    \n    for ($a=0; $a<$nseq; $a++)\n      {\\
n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$\
n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n    \
  }\n    return %hseq;\n  }\n","use Getopt::Long;\\
nuse File::Path;\nuse Env;\nuse FileHandle;\nuse C\
wd;\nuse Sys::Hostname;\nour $PIDCHILD;\nour $ERRO\
R_DONE;\nour @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\\
nour $EXIT_SUCCESS=0;\n\nour $REFDIR=getcwd;\nour \
$EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\n\nour $PRO\
GRAM=\"tc_generic_method.pl\";\nour $CL=$PROGRAM;\\
n\nour $CLEAN_EXIT_STARTED;\nour $debug_lock=$ENV{\
\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TC\
OFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour\
 $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERR\
ORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock (\
$$);\nif (isshellpid(getppid())){lock4tc(getppid()\
, \"LLOCK\", \"LSET\", \"$$\\n\");}\n      \nour $\
print;\nmy ($fmsq1, $fmsq2, $output, $outfile, $ar\
ch, $psv, $hmmtop_home, $trim, $cov, $sample, $mod\
e, $gor_home, $gor_seq, $gor_obs);\n\nGetOptions(\\
"-in=s\" => \\$fmsq1,\"-output=s\" =>\\$output ,\"\
-out=s\" => \\$outfile, \"-arch=s\" => \\$arch,\"-\
psv=s\" => \\$psv, \"-hmmtop_home=s\", \\$hmmtop_h\
ome,\"-trim=s\" =>\\$trim ,\"-print=s\" =>\\$print\
,\"-cov=s\" =>\\$cov , \"-sample=s\" =>\\$sample, \
\"-mode=s\" =>\\$mode, \"-gor_home=s\"=>\\$gor_hom\
e, \"-gor_seq=s\"=>\\$gor_seq,\"-gor_obs=s\"=>\\$g\
or_obs);\n\n\nif (!$mode){$mode = \"hmmtop\"}\nels\
if ($mode eq \"hmmtop\"){;}\nelsif ($mode eq \"gor\
\"){;}\nelse {myexit(flush_error (\"-mode=$mode is\
 unknown\"));}\n\n\nour $HOME=$ENV{\"HOME\"};\nour\
 $MCOFFEE=($ENV{\"MCOFFEE_4_TCOFFEE\"})?$ENV{\"MCO\
FFEE_4_TCOFFEE\"}:\"$HOME/.t_coffee/mcoffee\";\n\n\
if ($mode eq \"hmmtop\")\n  {\n    \n    check_con\
figuration (\"hmmtop\");\n    if (-e $arch){$ENV{'\
HMMTOP_ARCH'}=$arch;}\n    elsif (-e $ENV{HMMTOP_A\
RCH}){$arch=$ENV{HMMTOP_ARCH};}\n    elsif (-e \"$\
MCOFFEE/hmmtop.arch\"){$arch=$ENV{'HMMTOP_ARCH'}=\\
"$MCOFFEE/hmmtop.arch\";}\n    elsif (-e \"$hmmtop\
_home/hmmtop.arch\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$\
hmmtop_home/hmmtop.arch\";}\n    else {myexit(flus\
h_error ( \"Could not find ARCH file for hmmtop\")\
);}\n    \n    \n    if (-e $psv){$ENV{'HMMTOP_PSV\
'}=$psv;}\n    elsif (-e $ENV{HMMTOP_PSV}){$psv=$E\
NV{HMMTOP_PSV};}\n    elsif (-e \"$MCOFFEE/hmmtop.\
psv\"){$psv=$ENV{'HMMTOP_PSV'}=\"$MCOFFEE/hmmtop.p\
sv\";}\n    elsif (-e \"$hmmtop_home/hmmtop.psv\")\
{$psv=$ENV{'HMMTOP_PSV'}=\"$hmmtop_home/hmmtop.psv\
\";}\n    else {myexit(flush_error ( \"Could not f\
ind PSV file for hmmtop\"));}\n\n  }\nelsif ($mode\
 eq \"gor\")\n  {\n    our $GOR_SEQ;\n    our $GOR\
_OBS;\n    \n    check_configuration (\"gorIV\");\\
n    if (-e $gor_seq){$GOR_SEQ=$gor_seq;}\n    els\
if (-e $ENV{GOR_SEQ}){$GOR_SEQ=$ENV{GOR_SEQ};}\n  \
  elsif (-e \"$MCOFFEE/New_KS.267.seq\"){$GOR_SEQ=\
\"$MCOFFEE/New_KS.267.seq\";}\n    elsif (-e \"$go\
r_home/New_KS.267.seq\"){$GOR_SEQ=\"$gor_home/New_\
KS.267.seq\";}\n    else {myexit(flush_error ( \"C\
ould not find SEQ file for gor\"));}\n\n    if (-e\
 $gor_obs){$GOR_OBS=$gor_obs;}\n    elsif (-e $ENV\
{GOR_OBS}){$GOR_OBS=$ENV{GOR_OBS};}\n    elsif (-e\
 \"$MCOFFEE/New_KS.267.obs\"){$GOR_OBS=\"$MCOFFEE/\
New_KS.267.obs\";}\n    elsif (-e \"$gor_home/New_\
KS.267.obs\"){$GOR_OBS=\"$gor_home/New_KS.267.obs\\
";}\n    else {myexit(flush_error ( \"Could not fi\
nd OBS file for gor\"));}\n  }\n\n\nif ( ! -e $fms\
q1){myexit(flush_error (\"Could Not Read Input fil\
e $fmsq1\"));}\n\n\nmy $fmsq2=vtmpnam();\nmy $fmsq\
3=vtmpnam();\nmy $tmpfile=vtmpnam();\nmy $predfile\
=vtmpnam();\n\nif ($trim){$trim_action=\" +trim _a\
ln_%%$trim\\_K1 \";}\nif ($cov) {$cov_action= \" +\
sim_filter _aln_c$cov \";}\n&safe_system(\"t_coffe\
e -other_pg seq_reformat -in $fmsq1 -action +conve\
rt 'BOUJXZ-' $cov_action $trim_action -output fast\
a_aln -out $fmsq2\");\nmy (%pred, %seq, %predA);\n\
\n\n%seq=read_fasta_seq($fmsq2);\n%seq=fasta2sampl\
e(\\%seq, $sample);\n\nif (1==2 &&$mode eq \"hmmto\
p\" && $output eq \"cons\")\n  {\n    fasta2hmmtop\
_cons($outfile,\\%seq);\n  }\nelse\n  {\n   \n    \
%pred=fasta2pred(\\%seq, $mode);\n    %predA=pred2\
aln (\\%pred, \\%seq);\n    \n    \n    if (!$outp\
ut || $output eq \"prediction\"){output_fasta_seq \
(\\%predA, $outfile);}\n    elsif ($output eq \"co\
lor_html\"){pred2color (\\%pred,\\%seq, $outfile);\
}\n    elsif ($output eq \"cons\"){pred2cons($outf\
ile,\\%predA);}\n    else {flush_error (\"$output \
is an unknown output mode\");}\n  }\n\nsub fasta2s\
ample\n  {\n    my $SR=shift;\n    my $it=shift;\n\
    my %S=%$SR;\n    \n    my $seq=index2seq_name \
(\\%S, 1);\n    my $l=length($S{$seq}{seq});\n    \
my @sl=keys(%S);\n    my $nseq=$#sl+1;\n    my $in\
dex=$nseq;\n  \n    if (!$sample) {return %S;}\n  \
  for (my $a=0; $a<$it; $a++)\n      {\n	my $newse\
q=\"\";\n	my $nname=\"$seq\\_sampled_$index\";\n	f\
or (my $p=0; $p<$l; $p++)\n	  {\n	    my $i=int(ra\
nd($nseq));\n	    \n	    my $name = $sl[$i];\n	   \
 my $seq=$S{$name}{seq};\n	    my $r=substr ($seq,\
 $p, 1);\n	    $newseq.=$r;\n	  }\n	$S{$nname}{nam\
e}=$nname;\n	$S{$nname}{seq}=$newseq;\n	$S{$nname}\
{com}=\"sampled\";\n	$S{$nname}{index}=++$index;\n\
      }\n    return %S;\n  }\n	      \nsub fasta2p\
red\n  {\n    my $s=shift;\n    my $mode=shift;\n\\
n    if ( $mode eq \"hmmtop\"){return fasta2hmmtop\
_pred($s);}\n    elsif ($mode eq \"gor\"){return f\
asta2gor_pred ($s);}\n  }\nsub fasta2hmmtop_cons\n\
  {\n    my $outfile=shift;\n    my $SR=shift;\n  \
  \n    my $o = new FileHandle;\n    my $i = new F\
ileHandle;\n    my $tmp_in =vtmpnam();\n    my $tm\
p_out=vtmpnam();\n    my %seq=%$SR;\n    my %pred;\
\n    my $N=keys(%seq);\n    \n    output_fasta_se\
q (\\%seq,$tmp_in, \"seq\");\n    `hmmtop -pi=mpre\
d -if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\\
n    open ($o, \">$outfile\");\n    open ($i, \"$t\
mp_out\");\n    while (<$i>)\n      {\n	my $l=$_;\\
n	if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\\
">$2 NSEQ: $N\\n\";print $o \"$line\";}\n	elsif ( \
($l=~/.*pred(.*)/))  {my $line=\"$1\\n\";print $o \
\"$line\";}\n      }\n    close ($o);\n    close (\
$i);\n    return read_fasta_seq($tmp);\n  }\nsub f\
asta2hmmtop_pred\n  {\n    my $SR=shift;\n    my $\
o = new FileHandle;\n    my $i = new FileHandle;\n\
    my $tmp    =vtmpnam();\n    my $tmp_in =vtmpna\
m();\n    my $tmp_out=vtmpnam();\n    my %seq=%$SR\
;\n    my %pred;\n    \n\n    output_fasta_seq (\\\
%seq,$tmp_in, \"seq\");\n\n    \n    `hmmtop -if=$\
tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\n    \n\
\n    \n    \n    open ($o, \">$tmp\");\n    open \
($i, \"$tmp_out\");\n    while (<$i>)\n      {\n	m\
y $l=$_;\n	if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){m\
y $line=\">$2\\n\";print $o \"$line\";}\n	elsif ( \
($l=~/.*pred(.*)/))  {my $line=\"$1\\n\";print $o \
\"$line\";}\n      }\n    close ($o);\n    close (\
$i);\n    return read_fasta_seq($tmp);\n  }\n    \\
n	\n	\n	    \n	\n	\n\n	\nsub fasta2gor_pred\n  {\n\
    my $SR=shift;\n    my $o = new FileHandle;\n  \
  my $i = new FileHandle;\n    my $tmp    =vtmpnam\
();\n    my $tmp_in =vtmpnam();\n    my $tmp_out=v\
tmpnam();\n    my %seq=%$SR;\n    my %pred;\n    \\
n\n    output_fasta_seq (\\%seq,$tmp_in, \"seq\");\
\n    `gorIV -prd $tmp_in -seq $GOR_SEQ -obs $GOR_\
OBS >$tmp_out`;\n    open ($o, \">$tmp\");\n    op\
en ($i, \"$tmp_out\");\n    while (<$i>)\n      {\\
n	my $l=$_;\n\n	\n	if ( $l=~/>/){print $o \"$l\";}\
\n	elsif ( $l=~/Predicted Sec. Struct./){$l=~s/Pre\
dicted Sec. Struct\\.//;print $o \"$l\";}\n      }\
\n    close ($o);\n    close ($i);\n    return rea\
d_fasta_seq($tmp);\n  }\n			\n			     \nsub index2\
seq_name\n  {\n    \n    my $SR=shift;\n    my $in\
dex=shift;\n    \n    \n    my %S=%$SR;\n    \n   \
 foreach my $s (%S)\n      {\n	if ( $S{$s}{index}=\
=$index){return $s;}\n      }\n    return \"\";\n \
 }\n\nsub pred2cons\n  {\n    my $outfile=shift;\n\
    my $predR=shift;\n    my $seq=shift;\n    my %\
P=%$predR;\n    my %C;\n    my ($s,@r,$nseq);\n   \
 my $f= new FileHandle;\n\n    open ($f, \">$outfi\
le\");\n\n    if (!$seq){$seq=index2seq_name(\\%P,\
1);}\n    foreach my $s (keys(%P))\n      {\n	$nse\
q++;\n	$string= $P{$s}{seq};\n	$string = uc $strin\
g;\n	my @r=split (//,$string);\n	for (my $a=0; $a<\
=$#r; $a++)\n	  {\n	    if (($r[$a]=~/[OHICE]/)){$\
C{$a}{$r[$a]}++;}\n	  }\n      }\n    @l=keys(%C);\
\n    \n    \n    $s=$P{$seq}{seq};\n    print $f \
\">$seq pred based on $nseq\\n\";\n    @r=split (/\
/,$s);\n    \n    for (my $x=0; $x<=$#r; $x++)\n  \
    {\n	if ($r[$x] ne \"-\")\n	  {\n	    my $h=$C{\
$x}{H};\n	    my $i=$C{$x}{I};\n	    my $o=$C{$x}{\
O};\n	    my $c=$C{$x}{C};\n	    my $e=$C{$x}{E};\\
n	    my $l=$i+$o;\n	    \n	    if ($h>=$i && $h>=\
$o && $h>=$c && $h>=$e){$r[$x]='H';}\n	    elsif (\
$i>=$o && $i>=$c && $i>=$e){$r[$x]='I';}\n	    els\
if ($o>=$c && $o>=$e){$r[$x]='O';}\n	    elsif ($c\
>=$e){$r[$x]='C';}\n	    else {$r[$x]='E';}\n	  }\\
n      }\n    $j=join ('', @r);\n    print $f \"$j\
\\n\";\n    close ($f);\n    return $j;\n  }\n\nsu\
b pred2aln\n  {\n    my $PR=shift;\n    my $AR=shi\
ft;\n    \n    my $f=new FileHandle;\n    my %P=%$\
PR;\n    my %A=%$AR;\n    my %PA;\n    my $tmp=vtm\
pnam();\n    my $f= new FileHandle;\n    \n    ope\
n ($f, \">$tmp\");\n    foreach my $s (sort{$A{$a}\
{index}<=>$A{$b}{index}}(keys (%A)))\n      {\n	my\
 (@list, $seq, @plist, @pseq, $L, $PL, $c, $w);\n	\
my $seq;\n	my $seq=$A{$s}{seq};\n	my $pred=$P{$s}{\
seq};\n	$seq=pred2alnS($P{$s}{seq},$A{$s}{seq});\n\
	print $f \">$s\\n$seq\\n\";\n      }\n    close (\
$f);\n    return read_fasta_seq ($tmp);\n  }\nsub \
pred2alnS\n  {\n    my $pred=shift;\n    my $aln= \
shift;\n    my ($j,$a,$b);\n    my @P=split (//, $\
pred);\n    my @A=split (//, $aln);\n    for ($a=$\
b=0;$a<=$#A; $a++)\n      {\n	if ($A[$a] ne \"-\")\
{$A[$a]=$P[$b++];}\n      }\n    if ($b!= ($#P+1))\
{add_warning (\"Could not thread sequence: $b $#P\\
");}\n    \n    $j= join ('', @A);\n    return $j;\
\n  }\nsub pred2color\n  {\n    my $predP=shift;\n\
    my $alnP=shift;\n    my $out=shift;\n    my $F\
=new FileHandle;\n    my $struc=vtmpnam();\n    my\
 $aln=vtmpnam();\n    \n\n    output_fasta_seq ($a\
lnP, $aln);\n    my %p=%$predP;\n    \n    open ($\
F, \">$struc\");\n    \n    \n    foreach my $s (k\
eys(%p))\n      {\n	\n	print $F \">$s\\n\";\n	my $\
s=uc($p{$s}{seq});\n	\n	$s=~s/[Oo]/0/g;\n	$s=~s/[E\
e]/0/g;\n	\n	$s=~s/[Ii]/5/g;\n	$s=~s/[Cc]/5/g;\n	\\
n	$s=~s/[Hh]/9/g;\n	\n	print $F \"$s\\n\";\n      \
}\n    close ($F);\n    \n    \n    \n    safe_sys\
tem ( \"t_coffee -other_pg seq_reformat -in $aln -\
struc_in $struc -struc_in_f number_fasta -output c\
olor_html -out $out\");\n    return;\n  }\n	  \n  \
  \nsub display_fasta_seq\n  {\n    my $SR=shift;\\
n    my %S=%$SR;\n    \n    foreach my $s (sort{$S\
{$a}{index}<=>$S{$b}{index}}(keys (%S)))\n      {\\
n	print STDERR \">$s\\n$S{$s}{seq}\\n\";\n      }\\
n    close ($f);\n  }\nsub output_fasta_seq\n  {\n\
    my $SR=shift;\n    my $outfile=shift;\n    my \
$mode =shift;\n    my $f= new FileHandle;\n    my \
%S=%$SR;\n    \n    \n    open ($f, \">$outfile\")\
;\n    foreach my $s (sort{$S{$a}{index}<=>$S{$b}{\
index}}(keys (%S)))\n      {\n	my $seq=$S{$s}{seq}\
;\n	if ( $mode eq \"seq\"){$seq=~s/\\-//g;}\n	prin\
t $f \">$s\\n$seq\\n\";\n      }\n    close ($f);\\
n  }\n      \nsub read_fasta_seq \n  {\n    my $f=\
$_[0];\n    my %hseq;\n    my (@seq, @com, @name);\
\n    my ($a, $s,$nseq);\n    my $index;\n    open\
 (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n   \
   }\n    close (F);\n\n    \n    @name=($s=~/>(\\\
S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n(\
[^>]*)/g);\n    @com =($s=~/>.*(.*)\\n([^>]*)/g);\\
n\n\n    $nseq=$#name+1;\n    \n  \n    for ($a=0;\
 $a<$nseq; $a++)\n      {\n	my $n=$name[$a];\n	my \
$s;\n	$hseq{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s/\
/g;\n	$hseq{$n}{index}=++$index;\n	$hseq{$n}{seq}=\
$s;\n	$hseq{$n}{com}=$com[$a];\n      }\n    retur\
n %hseq;\n  }\n\n\nsub file2head\n      {\n	my $fi\
le = shift;\n	my $size = shift;\n	my $f= new FileH\
andle;\n	my $line;\n	open ($f,$file);\n	read ($f,$\
line, $size);\n	close ($f);\n	return $line;\n     \
 }\nsub file2tail\n      {\n	my $file = shift;\n	m\
y $size = shift;\n	my $f= new FileHandle;\n	my $li\
ne;\n	\n	open ($f,$file);\n	seek ($f,$size*-1, 2);\
\n	read ($f,$line, $size);\n	close ($f);\n	return \
$line;\n      }\n\n\nsub vtmpnam\n      {\n	my $r=\
rand(100000);\n	my $f=\"file.$r.$$\";\n	while (-e \
$f)\n	  {\n	    $f=vtmpnam();\n	  }\n	push (@TMPFI\
LE_LIST, $f);\n	return $f;\n      }\n\nsub myexit\\
n  {\n    my $code=@_[0];\n    if ($CLEAN_EXIT_STA\
RTED==1){return;}\n    else {$CLEAN_EXIT_STARTED=1\
;}\n    ### ONLY BARE EXIT\n    exit ($code);\n  }\
\nsub set_error_lock\n    {\n      my $name = shif\
t;\n      my $pid=$$;\n\n      \n      &lock4tc ($\
$,\"LERROR\", \"LSET\", \"$$ -- ERROR: $name $PROG\
RAM\\n\");\n      return;\n    }\nsub set_lock\n  \
{\n    my $pid=shift;\n    my $msg= shift;\n    my\
 $p=getppid();\n    &lock4tc ($pid,\"LLOCK\",\"LRE\
SET\",\"$p$msg\\n\");\n  }\nsub unset_lock\n   {\n\
     \n    my $pid=shift;\n    &lock4tc ($pid,\"LL\
OCK\",\"LRELEASE\",\"\");\n  }\nsub shift_lock\n  \
{\n    my $from=shift;\n    my $to=shift;\n    my \
$from_type=shift;\n    my $to_type=shift;\n    my \
$action=shift;\n    my $msg;\n    \n    if (!&lock\
4tc($from, $from_type, \"LCHECK\", \"\")){return 0\
;}\n    $msg=&lock4tc ($from, $from_type, \"LREAD\\
", \"\");\n    &lock4tc ($from, $from_type,\"LRELE\
ASE\", $msg);\n    &lock4tc ($to, $to_type, $actio\
n, $msg);\n    return;\n  }\nsub isshellpid\n  {\n\
    my $p=shift;\n    if (!lock4tc ($p, \"LLOCK\",\
 \"LCHECK\")){return 0;}\n    else\n      {\n	my $\
c=lock4tc($p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-\
SHELL-/){return 1;}\n      }\n    return 0;\n  }\n\
sub isrootpid\n  {\n    if(lock4tc (getppid(), \"L\
LOCK\", \"LCHECK\")){return 0;}\n    else {return \
1;}\n  }\nsub lock4tc\n	{\n	  my ($pid,$type,$acti\
on,$value)=@_;\n	  my $fname;\n	  my $host=hostnam\
e;\n	  \n	  if ($type eq \"LLOCK\"){$fname=\"$LOCK\
DIR/.$pid.$host.lock4tcoffee\";}\n	  elsif ( $type\
 eq \"LERROR\"){ $fname=\"$LOCKDIR/.$pid.$host.err\
or4tcoffee\";}\n	  elsif ( $type eq \"LWARNING\"){\
 $fname=\"$LOCKDIR/.$pid.$host.warning4tcoffee\";}\
\n	  \n	  if ($debug_lock)\n	    {\n	      print S\
TDERR \"\\n\\t---lock4tc(tcg): $action => $fname =\
>$value (RD: $LOCKDIR)\\n\";\n	    }\n\n	  if    (\
$action eq \"LCHECK\") {return -e $fname;}\n	  els\
if ($action eq \"LREAD\"){return file2string($fnam\
e);}\n	  elsif ($action eq \"LSET\") {return strin\
g2file ($value, $fname, \">>\");}\n	  elsif ($acti\
on eq \"LRESET\") {return string2file ($value, $fn\
ame, \">\");}\n	  elsif ($action eq \"LRELEASE\") \
\n	    {\n	      if ( $debug_lock)\n		{\n		  my $g\
=new FileHandle;\n		  open ($g, \">>$fname\");\n		\
  print $g \"\\nDestroyed by $$\\n\";\n		  close (\
$g);\n		  safe_system (\"mv $fname $fname.old\");\\
n		}\n	      else\n		{\n		  unlink ($fname);\n		}\\
n	    }\n	  return \"\";\n	}\n	\nsub file2string\n\
	{\n	  my $file=@_[0];\n	  my $f=new FileHandle;\n\
	  my $r;\n	  open ($f, \"$file\");\n	  while (<$f\
>){$r.=$_;}\n	  close ($f);\n	  return $r;\n	}\nsu\
b string2file \n    {\n    my ($s,$file,$mode)=@_;\
\n    my $f=new FileHandle;\n    \n    open ($f, \\
"$mode$file\");\n    print $f  \"$s\";\n    close \
($f);\n  }\n\nBEGIN\n    {\n      srand;\n    \n  \
    $SIG{'SIGUP'}='signal_cleanup';\n      $SIG{'S\
IGINT'}='signal_cleanup';\n      $SIG{'SIGQUIT'}='\
signal_cleanup';\n      $SIG{'SIGILL'}='signal_cle\
anup';\n      $SIG{'SIGTRAP'}='signal_cleanup';\n \
     $SIG{'SIGABRT'}='signal_cleanup';\n      $SIG\
{'SIGEMT'}='signal_cleanup';\n      $SIG{'SIGFPE'}\
='signal_cleanup';\n      \n      $SIG{'SIGKILL'}=\
'signal_cleanup';\n      $SIG{'SIGPIPE'}='signal_c\
leanup';\n      $SIG{'SIGSTOP'}='signal_cleanup';\\
n      $SIG{'SIGTTIN'}='signal_cleanup';\n      $S\
IG{'SIGXFSZ'}='signal_cleanup';\n      $SIG{'SIGIN\
FO'}='signal_cleanup';\n      \n      $SIG{'SIGBUS\
'}='signal_cleanup';\n      $SIG{'SIGALRM'}='signa\
l_cleanup';\n      $SIG{'SIGTSTP'}='signal_cleanup\
';\n      $SIG{'SIGTTOU'}='signal_cleanup';\n     \
 $SIG{'SIGVTALRM'}='signal_cleanup';\n      $SIG{'\
SIGUSR1'}='signal_cleanup';\n\n\n      $SIG{'SIGSE\
GV'}='signal_cleanup';\n      $SIG{'SIGTERM'}='sig\
nal_cleanup';\n      $SIG{'SIGCONT'}='signal_clean\
up';\n      $SIG{'SIGIO'}='signal_cleanup';\n     \
 $SIG{'SIGPROF'}='signal_cleanup';\n      $SIG{'SI\
GUSR2'}='signal_cleanup';\n\n      $SIG{'SIGSYS'}=\
'signal_cleanup';\n      $SIG{'SIGURG'}='signal_cl\
eanup';\n      $SIG{'SIGCHLD'}='signal_cleanup';\n\
      $SIG{'SIGXCPU'}='signal_cleanup';\n      $SI\
G{'SIGWINCH'}='signal_cleanup';\n      \n      $SI\
G{'INT'}='signal_cleanup';\n      $SIG{'TERM'}='si\
gnal_cleanup';\n      $SIG{'KILL'}='signal_cleanup\
';\n      $SIG{'QUIT'}='signal_cleanup';\n      \n\
      our $debug_lock=$ENV{\"DEBUG_LOCK\"};\n     \
 \n      \n      \n      \n      foreach my $a (@A\
RGV){$CL.=\" $a\";}\n      if ( $debug_lock ){prin\
t STDERR \"\\n\\n\\n********** START PG: $PROGRAM \
*************\\n\";}\n      if ( $debug_lock ){pri\
nt STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOC\
KDIR $$ *************\\n\";}\n      if ( $debug_lo\
ck ){print STDERR \"\\n --- $$ -- $CL\\n\";}\n    \
  \n	     \n      \n      \n    }\nsub flush_error\
\n  {\n    my $msg=shift;\n    return add_error ($\
EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);\n  }\ns\
ub add_error \n  {\n    my $code=shift;\n    my $r\
pid=shift;\n    my $pid=shift;\n    my $ppid=shift\
;\n    my $type=shift;\n    my $com=shift;\n    \n\
    $ERROR_DONE=1;\n    lock4tc ($rpid, \"LERROR\"\
,\"LSET\",\"$pid -- ERROR: $type\\n\");\n    lock4\
tc ($$, \"LERROR\",\"LSET\", \"$pid -- COM: $com\\\
n\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pi\
d -- STACK: $ppid -> $pid\\n\");\n   \n    return \
$code;\n  }\nsub add_warning \n  {\n    my $rpid=s\
hift;\n    my $pid =shift;\n    my $command=shift;\
\n    my $msg=\"$$ -- WARNING: $command\\n\";\n   \
 print STDERR \"$msg\";\n    lock4tc ($$, \"LWARNI\
NG\", \"LSET\", $msg);\n  }\n\nsub signal_cleanup\\
n  {\n    print dtderr \"\\n**** $$ (tcg) was kill\
ed\\n\";\n    &cleanup;\n    exit ($EXIT_FAILURE);\
\n  }\nsub clean_dir\n  {\n    my $dir=@_[0];\n   \
 if ( !-d $dir){return ;}\n    elsif (!($dir=~/tmp\
/)){return ;}#safety check 1\n    elsif (($dir=~/\\
\*/)){return ;}#safety check 2\n    else\n      {\\
n	`rm -rf $dir`;\n      }\n    return;\n  }\nsub c\
leanup\n  {\n    #print stderr \"\\n----tc: $$ Kil\
ls $PIDCHILD\\n\";\n    #kill (SIGTERM,$PIDCHILD);\
\n    my $p=getppid();\n    $CLEAN_EXIT_STARTED=1;\
\n    \n    \n    \n    if (&lock4tc($$,\"LERROR\"\
, \"LCHECK\", \"\"))\n      {\n	my $ppid=getppid()\
;\n	if (!$ERROR_DONE) \n	  {\n	    &lock4tc($$,\"L\
ERROR\", \"LSET\", \"$$ -- STACK: $p -> $$\\n\");\\
n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- CO\
M: $CL\\n\");\n	  }\n      }\n    my $warning=&loc\
k4tc($$, \"LWARNING\", \"LREAD\", \"\");\n    my $\
error=&lock4tc($$,  \"LERROR\", \"LREAD\", \"\");\\
n    #release error and warning lock if root\n    \
\n    if (isrootpid() && ($warning || $error) )\n \
     {\n	\n	print STDERR \"**************** Summar\
y *************\\n$error\\n$warning\\n\";\n\n	&loc\
k4tc($$,\"LERROR\",\"RELEASE\",\"\");\n	&lock4tc($\
$,\"LWARNING\",\"RELEASE\",\"\");\n      } \n    \\
n    \n    foreach my $f (@TMPFILE_LIST)\n      {\\
n	if (-e $f){unlink ($f);} \n      }\n    foreach \
my $d (@TMPDIR_LIST)\n      {\n	clean_dir ($d);\n \
     }\n    #No More Lock Release\n    #&lock4tc($\
$,\"LLOCK\",\"LRELEASE\",\"\"); #release lock \n\n\
    if ( $debug_lock ){print STDERR \"\\n\\n\\n***\
******* END PG: $PROGRAM ($$) *************\\n\";}\
\n    if ( $debug_lock ){print STDERR \"\\n\\n\\n*\
*********(tcg) LOCKDIR: $LOCKDIR $$ *************\\
\n\";}\n  }\nEND \n  {\n    \n    &cleanup();\n  }\
\n   \n\nsub safe_system \n{\n  my $com=shift;\n  \
my $ntry=shift;\n  my $ctry=shift;\n  my $pid;\n  \
my $status;\n  my $ppid=getppid();\n  if ($com eq \
\"\"){return 1;}\n  \n  \n\n  if (($pid = fork ())\
 < 0){return (-1);}\n  if ($pid == 0)\n    {\n    \
  set_lock($$, \" -SHELL- $com (tcg)\");\n      ex\
ec ($com);\n    }\n  else\n    {\n      lock4tc ($\
$, \"LLOCK\", \"LSET\", \"$pid\\n\");#update paren\
t\n      $PIDCHILD=$pid;\n    }\n  if ($debug_lock\
){printf STDERR \"\\n\\t .... safe_system (fasta_s\
eq2hmm)  p: $$ c: $pid COM: $com\\n\";}\n\n  waitp\
id ($pid,WTERMSIG);\n\n  shift_lock ($pid,$$, \"LW\
ARNING\",\"LWARNING\", \"LSET\");\n\n  if ($? == $\
EXIT_FAILURE || lock4tc($pid, \"LERROR\", \"LCHECK\
\", \"\"))\n    {\n      if ($ntry && $ctry <$ntry\
)\n	{\n	  add_warning ($$,$$,\"$com failed [retry:\
 $ctry]\");\n	  lock4tc ($pid, \"LRELEASE\", \"LER\
ROR\", \"\");\n	  return safe_system ($com, $ntry,\
 ++$ctry);\n	}\n      elsif ($ntry == -1)\n	{\n	  \
if (!shift_lock ($pid, $$, \"LERROR\", \"LWARNING\\
", \"LSET\"))\n	    {\n	      add_warning ($$,$$,\\
"$com failed\");\n	    }\n	  else\n	    {\n	      \
lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	\
    }\n	  return $?;}\n      else\n	{\n	  if (!shi\
ft_lock ($pid,$$, \"LERROR\",\"LERROR\", \"LSET\")\
)\n	    {\n	      myexit(add_error ($EXIT_FAILURE,\
$$,$pid,getppid(), \"UNSPECIFIED system\", $com));\
\n	    }\n	}\n    }\n  return $?;\n}\n\nsub check_\
configuration \n    {\n      my @l=@_;\n      my $\
v;\n      foreach my $p (@l)\n	{\n	  \n	  if   ( $\
p eq \"EMAIL\")\n	    { \n	      if ( !($EMAIL=~/@\
/))\n		{\n		add_warning($$,$$,\"Could Not Use EMAI\
L\");\n		myexit(add_error ($EXIT_FAILURE,$$,$$,get\
ppid(),\"EMAIL\",\"$CL\"));\n	      }\n	    }\n	  \
elsif( $p eq \"INTERNET\")\n	    {\n	      if ( !&\
check_internet_connection())\n		{\n		  myexit(add_\
error ($EXIT_FAILURE,$$,$$,getppid(),\"INTERNET\",\
\"$CL\"));\n		}\n	    }\n	  elsif( $p eq \"wget\")\
\n	    {\n	      if (!&pg_is_installed (\"wget\") \
&& !&pg_is_installed (\"curl\"))\n		{\n		  myexit(\
add_error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_\
INSTALLED:wget\",\"$CL\"));\n		}\n	    }\n	  elsif\
( !(&pg_is_installed ($p)))\n	    {\n	      myexit\
(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT\
_INSTALLED:$p\",\"$CL\"));\n	    }\n	}\n      retu\
rn 1;\n    }\nsub pg_is_installed\n  {\n    my @ml\
=@_;\n    my $r, $p, $m;\n    my $supported=0;\n  \
  \n    my $p=shift (@ml);\n    if ($p=~/::/)\n   \
   {\n	if (safe_system (\"perl -M$p -e 1\")==$EXIT\
_SUCCESS){return 1;}\n	else {return 0;}\n      }\n\
    else\n      {\n	$r=`which $p 2>/dev/null`;\n	i\
f ($r eq \"\"){return 0;}\n	else {return 1;}\n    \
  }\n  }\n\n\n\nsub check_internet_connection\n  {\
\n    my $internet;\n    my $tmp;\n    &check_conf\
iguration ( \"wget\"); \n    \n    $tmp=&vtmpnam (\
);\n    \n    if     (&pg_is_installed    (\"wget\\
")){`wget www.google.com -O$tmp >/dev/null 2>/dev/\
null`;}\n    elsif  (&pg_is_installed    (\"curl\"\
)){`curl www.google.com -o$tmp >/dev/null 2>/dev/n\
ull`;}\n    \n    if ( !-e $tmp || -s $tmp < 10){$\
internet=0;}\n    else {$internet=1;}\n    if (-e \
$tmp){unlink $tmp;}\n\n    return $internet;\n  }\\
nsub check_pg_is_installed\n  {\n    my @ml=@_;\n \
   my $r=&pg_is_installed (@ml);\n    if (!$r && $\
p=~/::/)\n      {\n	print STDERR \"\\nYou Must Ins\
tall the perl package $p on your system.\\nRUN:\\n\
\\tsudo perl -MCPAN -e 'install $pg'\\n\";\n      \
}\n    elsif (!$r)\n      {\n	myexit(flush_error(\\
"\\nProgram $p Supported but Not Installed on your\
 system\"));\n      }\n    else\n      {\n	return \
1;\n      }\n  }\n\n\n\n","\n\n\n\n\nmy $FMODEL =\\
"\"; \nmy $TMPDIR = \"/tmp\";\n\n\n\n\nmy $NUCALPH\
 = \"ACGTUNRYMKSWHBVD\";\nmy $PRIMNUCALPH = \"ACGT\
UN\";\nuse vars qw($NUCALPH $PRIMNUCALPH $TMPDIR);\
\n\n\nmy $errmsg;\nuse vars qw($errmsg);\n\n\n\nus\
e Getopt::Long;\nuse Cwd;\nuse File::Basename;\nus\
e File::Temp qw/ tempfile tempdir /;\nuse File::Co\
py;\nuse File::Path;\n\n\n\nsub usage(;$)\n{\n    \
my ($errmsg) = @_;\n    my $myname = basename($0);\
\n\n    if ($errmsg) {\n        print STDERR \"ERR\
OR: $errmsg\\n\";\n    }\n\n    print STDERR << \"\
EOF\";\n    \n$myname: align two sequences by mean\
s of consan\\'s sfold\nUsage:\n $myname -i file -o\
 file -d path\nOptions:\n -i|--in : pairwise input\
 sequence file\n -o|--out: output alignment\n -d|-\
-directory containing data\n\nEOF\n}\n\nsub read_s\
tk_aln \n  {\n    my $f=$_[0];\n    my ($seq, $id)\
;\n    \n    my %hseq;\n\n    open (STK, \"$f\");\\
n    while (<STK>)\n      {\n	if ( /^#/ || /^\\/\\\
// || /^\\s*$/){;}\n	else\n	  {\n	    ($id,$seq)=/\
(\\S+)\\s+(\\S+)/;\n	    $hseq{$id}{'seq'}.=$seq;\\
n	  }\n      }\n    close (STK);\n    return %hseq\
;\n  }\nsub read_fasta_seq \n  {\n    my $f=$_[0];\
\n    my %hseq;\n    my (@seq, @com, @name);\n    \
my ($a, $s,$nseq);\n\n    open (F, $f);\n    while\
 (<F>)\n      {\n	$s.=$_;\n      }\n    close (F);\
\n\n    \n    @name=($s=~/>(.*).*\\n[^>]*/g);\n   \
 \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =\
($s=~/>.*(.*)\\n([^>]*)/g);\n\n    \n    $nseq=$#n\
ame+1;\n    \n    for ($a=0; $a<$nseq; $a++)\n    \
  {\n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$hs\
eq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n\
      }\n    return %hseq;\n  }\n\n\n\nsub sfold_p\
arseoutput($$)\n{\n    my ($frawout, $foutfa) = @_\
;\n    my %haln;\n    my ($fstk, $cmd, $id);\n    \
open FOUTFA, \">$foutfa\";\n    \n    $fstk = $fra\
wout . \".stk\";\n    \n    # first line of raw ou\
t contains info\n    # remaining stuff is stockhol\
m formatted\n    $cmd = \"sed -e '1d' $frawout\";\\
n    system(\"$cmd > $fstk\");\n    if ($? != 0) {\
\n        $errmsg = \"command failed with exit sta\
tus $?.\";\n        $errmsg .=  \"Command was \\\"\
$cmd\\\"\";\n        return -1;\n    }\n\n    # th\
is gives an error message. just ignore it...\n    \
%haln=read_stk_aln ( $fstk);\n    foreach $i (keys\
 (%haln))\n      {\n	my $s;\n	$s=$haln{$i}{'seq'};\
\n	$s =~ s/\\./-/g;\n	print FOUTFA \">$i\\n$s\\n\"\
;\n      }\n    close FOUTFA;\n    return 0;\n}\n\\
n\n\n\nsub sfold_wrapper($$$$)\n{\n    \n    my ($\
fs1, $fs2, $fmodel, $foutfa) = @_;\n    \n\n    my\
 ($cmd, $frawout, $ferrlog, $freadme, $ftimelog, $\
fstk);\n\n    # add  basename($fmsqin) (unknown he\
re!)\n    $frawout = \"sfold.log\";\n    $ferrlog \
= \"sfold.err\";\n    $ftimelog = \"sfold.time\";\\
n    $freadme =  \"sfold.README\";\n    $fstk = \"\
sfold.stk\";\n    \n    # prepare execution...\n  \
  #\n    # ./tmp is essential for dswpalign\n    #\
 otherwise you'll get a segfault\n    mkdir \"./tm\
p\";\n    \n    $cmd = \"sfold -m $fmodel $fs1 $fs\
2\";\n    open(FREADME,\">$freadme\");\n    print \
FREADME \"$cmd\\n\"; \n    close(FREADME);\n\n    \
# and go\n    #\n    system(\"/usr/bin/time -p -o \
$ftimelog $cmd >$frawout 2>$ferrlog\");\n    if ($\
? != 0) {\n        $errmsg = \"command failed with\
 exit status $?\";\n        $errmsg .= \"command w\
as \\\"$cmd\\\". See \" . getcwd . \"\\n\";\n     \
   return -1;\n    }\n\n    return sfold_parseoutp\
ut($frawout, $foutfa);\n}\n\n\n\n\n\n\n\nmy ($help\
, $fmsqin, $fmsaout);\nGetOptions(\"help\"  => \\$\
help,\n           \"in=s\" => \\$fmsqin,\n        \
   \"out=s\" => \\$fmsaout,\n	   \"data=s\" => \\$\
ref_dir);\n\n\n\nif ($help) {\n    usage();\n    e\
xit(0);\n}\nif (! defined($fmsqin)) {\n    usage('\
missing input filename');\n    exit(1);\n}\nif (! \
defined($fmsaout)) {\n    usage('missing output fi\
lename');\n    exit(1);\n\n}\nif (scalar(@ARGV)) {\
\n    usage('Unknown remaining args');\n    exit(1\
);\n}\n\n$FMODEL = \"$ref_dir/mix80.mod\";\nif (! \
-e \"$FMODEL\") {\n    die(\"couldn't find sfold g\
rammar model file. Expected $FMODEL\\n\");\n}\n\n\\
nmy %hseq=read_fasta_seq ($fmsqin);\nmy $id;\n\nfo\
reach $id (keys(%hseq))\n  {\n    push(@seq_array,\
 $hseq{$id});\n  }\n\nif ( scalar(@seq_array) != 2\
 ) {\n    die(\"Need *exactly* two sequences as in\
put (pairwise alignment!).\")\n}\n\n\n\nmy ($sec, \
$min, $hour, $mday, $mon, $year, $wday, $yday, $is\
dst) = localtime(time);\nmy $datei = sprintf(\"%4d\
-%02d-%02d\", $year+1900, $mon+1, $mday);\nmy $tem\
pl = basename($0) . \".\" . $datei . \".pid-\" . $\
$ . \".XXXXXX\";\nmy $wd = tempdir ( $templ, DIR =\
> $TMPDIR);\n\ncopy($fmsqin, \"$wd/\" . basename($\
fmsqin) . \".org\"); # for reproduction\ncopy($FMO\
DEL, \"$wd\");\nmy $fmodel = basename($FMODEL);\nm\
y $orgwd = getcwd;\nchdir $wd;\n\n\n\nmy @sepseqfi\
les;\nforeach $id (keys(%hseq)) {\n    my ($seq, $\
orgseq, $fname, $sout);\n    $seq=$hseq{$id}{'seq'\
};\n    \n    $fname = basename($fmsqin) . \"_$id.\
fa\";\n    # replace funnies in file/id name (e.g.\
 \"/\" \" \" etc)\n    $fname =~ s,[/ ],_,g;\n    \
open (PF, \">$fname\");\n    print (PF \">$id\\n$s\
eq\\n\");\n    close (PF);\n\n    push(@sepseqfile\
s, $fname);\n}\n\nmy ($f1, $f2, $fout);\n$f1 = $se\
pseqfiles[0];\n$f2 = $sepseqfiles[1];\n$fout = $wd\
 . basename($fmsqin) . \".out.fa\";\nif (sfold_wra\
pper($f1, $f2, $fmodel, \"$fout\") != 0) {\n    pr\
intf STDERR \"ERROR: See logs in $wd\\n\";\n    ex\
it(1);\n} else {\n    chdir $orgwd;\n    copy($fou\
t, $fmsaout);\n    rmtree($wd);\n   exit(0);\n}\n"\
,"\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env \
qw(USER);\n\n\n$tmp=clean_cr ($ARGV[0]);\nopen (F,\
 $tmp);\n\nwhile ( <F>)\n  {\n    my $l=$_;\n    i\
f ( $l=~/^# STOCKHOLM/){$stockholm=1;}\n    elsif \
( $stockholm && $l=~/^#/)\n      {\n	$l=~/^#(\\S+)\
\\s+(\\S+)\\s+(\\S*)/g;\n	$l=\"_stockholmhasch_$1\\
\_stockholmspace_$2 $3\\n\";\n      }\n    $file.=\
$l;\n  }\nclose (F);\nunlink($tmp);\n$file1=$file;\
\n\n$file=~s/\\#/_hash_symbol_/g;\n$file=~s/\\@/_a\
robase_symbol_/g;\n\n\n$file=~s/\\n[\\.:*\\s]+\\n/\
\\n\\n/g;\n\n$file=~s/\\n[ \\t\\r\\f]+(\\b)/\\n\\1\
/g;\n\n\n$file=~s/(\\n\\S+)(\\s+)(\\S)/\\1_blank_\\
\3/g;\n\n$file=~s/[ ]//g;\n$file=~s/_blank_/ /g;\n\
\n\n\n$file =~s/\\n\\s*\\n/#/g;\n\n$file.=\"#\";\n\
$file =~s/\\n/@/g;\n\n\n\n\n@blocks=split /\\#/, $\
file;\nshift (@blocks);\n@s=split /\\@/, $blocks[0\
];\n$nseq=$#s+1;\n\n\n\n$file=join '@', @blocks;\n\
@lines=split /\\@/,$file;\n\n$c=0;\n\nforeach $l (\
@lines)\n  {\n    if (!($l=~/\\S/)){next;}\n    el\
sif ($stockholm && ($l=~/^\\/\\// || $l=~/STOCKHOL\
M/)){next;}#get read of STOCHOLM Terminator\n   \n\
    $l=~/(\\S+)\\s+(\\S*)/g;\n    $n=$1; $s=$2;\n \
   \n    $seq[$c].=$s;\n    $name[$c]=$n;\n    $c+\
+;\n    \n    if ( $c==$nseq){$c=0;}\n    \n  } \n\
\nif ( $c!=0)\n      {\n	print STDERR \"ERROR: $AR\
GV[0] is NOT an MSA in Clustalw format: make sure \
there is no blank line within a block [ERROR]\\n\"\
;\n	exit (EXIT_FAILURE);\n      }\n\nfor ($a=0; $a\
< $nseq; $a++)\n  {\n    $name[$a]=cleanstring ($n\
ame[$a]);\n    $seq[$a]=cleanstring ($seq[$a]);\n \
   $seq[$a]=breakstring($seq[$a], 60);\n    \n    \
$line=\">$name[$a]\\n$seq[$a]\\n\";\n    \n    pri\
nt \"$line\";\n  }\nexit (EXIT_SUCCESS);\n\nsub cl\
eanstring\n  {\n    my $s=@_[0];\n    $s=~s/_hash_\
symbol_/\\#/g;\n    $s=~s/_arobase_symbol_/\\@/g;\\
n    $s=~s/[ \\t]//g;\n    return $s;\n  }\nsub br\
eakstring\n  {\n    my $s=@_[0];\n    my $size=@_[\
1];\n    my @list;\n    my $n,$ns, $symbol;\n    \\
n    @list=split //,$s;\n    $n=0;$ns=\"\";\n    f\
oreach $symbol (@list)\n      {\n	if ( $n==$size)\\
n	  {\n	    $ns.=\"\\n\";\n	    $n=0;\n	  }\n	$ns.\
=$symbol;\n	$n++;\n      }\n    return $ns;\n    }\
\n\nsub clean_cr\n  {\n    my $f=@_[0];\n    my $f\
ile;\n    \n    $tmp=\"f$.$$\";\n    \n    \n    o\
pen (IN, $f);\n    open (OUT, \">$tmp\");\n    \n \
   while ( <IN>)\n      {\n	$file=$_;\n	$file=~s/\\
\r\\n/\\n/g;\n	$file=~s/\\n\\r/\\n/g;\n	$file=~s/\\
\r\\r/\\n/g;\n	$file=~s/\\r/\\n/g;\n	print OUT \"$\
file\";\n      }\n    \n    close (IN);\n    close\
 (OUT);\n    return $tmp;\n  }\n","use strict;\nus\
e FileHandle;\nuse Env qw(HOST);\nuse Env qw(HOME)\
;\nuse Env qw(USER);\n\nmy $format=file2format ($A\
RGV[0]);\n\nif    ($format eq \"clustalw\"){clusta\
lw2fasta($ARGV[0]);}\nelsif ($format eq \"fasta\")\
   {fasta2fasta($ARGV[0]);}\nelsif ($format eq \"m\
sf\")   {msf2fasta($ARGV[0]);}\nelsif ($format eq \
\"phylip\")   {phylip2fasta($ARGV[0]);}\nelsif ($f\
ormat eq \"nameseq\") {display_file ($ARGV[0]);}\n\
 \nexit (0);\n\nsub file2format\n  {\n    my $f=sh\
ift;\n    \n    my $l=file2n_lines($f,2);\n    \n \
   if ( $l=~/^CLUSTAL/){return \"clustalw\";}\n   \
 elsif ($l=~/^SAGA/){return \"clustalw\";}\n    el\
sif ($l=~/^>/){return \"fasta\";}\n    elsif ($l=~\
/^PileUp/){return \"msf\";}\n    elsif ($l=~/\\s+\\
\d+\\s+\\d+\\s/){return \"phylip\";}\n    elsif ($\
l=~/\\#NAMESEQ_01/){return \"nameseq\";}\n    else\
 \n      {\n	print STDERR \"ERROR: $f FILE is NOT \
a supported format [ERROR]\\n\";\n	system (\"cp $f\
 /Users/cnotredame/cedric1.txt\");\n	exit (1);\n  \
    }\n  }\nsub display_file\n    {\n       my $fi\
le=shift;\n       my $F= new FileHandle;\n       o\
pen ($F, $file);\n       while (<$F>){print \"$_\"\
;}\n       close ($F);\n     }\nsub phylip2fasta\n\
    {\n      my $file=shift;\n      my $F= new Fil\
eHandle;\n      my ($seq, $name,$seq);\n      my $\
query_start=-1;\n      my $query_end=-1;\n      my\
 $in_aln=0;\n      my %list;\n      my ($first,$se\
q,$name, $cn, $nseq, $l,%len);\n      \n      open\
 ($F, $file);\n      <$F>;\n      $l=$_;\n      $l\
=~/\\s*(\\d+)\\s*(\\d+)/;\n      $first=1;\n      \
$cn=0;\n      while (<$F>)\n	{\n	  my $l=$_;\n	  i\
f (!($l=~/\\S/))\n	    {\n	      $cn=0;\n	      $f\
irst=0;\n	    }\n	  elsif ($first==1)\n	    {\n	  \
    $l=~/\\s*(\\S+)(.*)/;\n	      my $name=$1;\n	 \
     my $seq=$2;\n	      chomp ($seq);\n	      $se\
q=~s/\\s//g;\n	      $list{$cn}{'name'}=$name;\n	 \
     $list{$cn}{'seq'}.=$seq;\n	      $cn++;\n	   \
   $nseq++;\n	    }\n	  else\n	    {\n	      chomp\
 ($l);\n	      $l=~s/\\s//g;\n	      $list{$cn}{'s\
eq'}.=$l;\n	      $cn++;\n	    }\n	}\n      close \
($F);\n      \n      for (my $a=0; $a<$nseq; $a++)\
\n	{\n	  print \">$list{$a}{'name'}\\n$list{$a}{'s\
eq'}\\n\";\n	}\n    }\n      \nsub msf2fasta\n    \
{\n      my $file=shift;\n      my $F= new FileHan\
dle;\n      my ($seq, $name,$seq);\n      my $quer\
y_start=-1;\n      my $query_end=-1;\n      my $in\
_aln=0;\n      my %list;\n      my ($seq,$name, $n\
, $nseq, $l,%len);\n      \n      open ($F, $file)\
;\n      while (<$F>)\n	{\n	  if ( /\\/\\//){$in_a\
ln=1;}\n	  elsif ( $in_aln && /(\\S+)\\s+(.*)/)\n	\
    {\n	      $name=$1;\n	      $seq=$2;\n	      $\
seq=~s/\\s//g;\n	      $seq=~s/\\~/\\-/g;\n	      \
$seq=~s/\\./\\-/g;\n	      if ( $list{$n}{'name'} \
&& $list{$n}{'name'} ne $name)\n		{\n		  print \"$\
list{$n}{'name'} Vs $name\";\n		  \n		  exit (1);\\
n		}\n	      else\n		{\n		  $list{$n}{'name'}= $na\
me;\n		}\n	      \n	      $list{$n}{'seq'}=$list{$\
n}{'seq'}.$seq;\n	      \n	      $nseq=++$n;\n	   \
   \n	    }\n	  else\n	    {$n=0;}\n	}\n      clos\
e ($F);\n      \n      for (my $a=0; $a<$nseq; $a+\
+)\n	{\n	  my $nl=length ($list{$a}{'name'});\n	  \
my $sl=length ($list{$a}{'seq'});\n	  print \">$li\
st{$a}{'name'}\\n$list{$a}{'seq'}\\n\";\n	}\n    }\
\n    \nsub fasta2fasta\n    {\n      my $file=shi\
ft;\n      my $F= new FileHandle;\n      my ($seq,\
 $name,$n,$l,%len);\n      my $started=0;\n      o\
pen ($F, $file);\n      while (<$F>)\n	{\n	  if ( \
/^>(\\S+)/){$n++;$seq=\"\";$name=$1;}\n	  else\n	 \
   {\n	      $l=$_;\n	      chomp ($l);\n	      \n\
	      $seq.=$l;\n	      $len{$name}=length($seq);\
\n	    }\n	}\n      close ($F);\n      \n      ope\
n ($F, $file);\n      while (<$F>)\n	{\n	  my $l=$\
_;\n	  $l=~s/\\r[\\n]*/\\n/gm;\n	  if ( ($l=~/^>(\\
\S+)(.*)\\n/))\n	    {\n	      my $name=$1;\n	    \
  my $comment=$2;\n	      my $nl=length ($name);\n\
	      my $sl=$len{$name};\n	      if ($comment)\n\
		{\n		  $comment=~s/^\\s+//g;\n		  my $cl=length \
($comment);\n		}\n	      if (!$started){$started=1\
;print \">$name\\n\";}\n	      else {print \"\\n>$\
name\\n\"}\n	    }\n	  else\n	    {\n	      $l=$_;\
\n	      chomp ($l);\n	      $l=~s/\\W//g;\n	     \
 print \"$l\";\n	    }\n	}\n      print \"\\n\";\n\
      close ($F);\n    }\nsub clustalw2fasta\n  {\\
n    my $fname=shift;\n    my ($file1, $file);\n  \
  my (@blocks, @lines,@s, $n,$nseq, $c);\n    my (\
@name, @seq);\n    my $F= new FileHandle;\n    my \
$stockholm;\n   \n    \n    open ($F, $fname);\n  \
  \n    while ( <$F>)\n      {\n	my $l=$_;\n	$l=cl\
ean_cr($l);\n	if ( $l=~/^# STOCKHOLM/){$stockholm=\
1;}\n	elsif ( $stockholm && $l=~/^#/)\n	  {\n	    \
$l=~/^#(\\S+)\\s+(\\S+)\\s+(\\S*)/g;\n	    $l=\"_s\
tockholmhasch_$1\\_stockholmspace_$2 $3\\n\";\n	  \
}\n	$file.=$l;\n      }\n    close ($F);\n        \
\n    #Protect # and @\n    $file=~s/\\#/_hash_sym\
bol_/g;\n    $file=~s/\\@/_arobase_symbol_/g;\n   \
 \n    \n    #Remove annotation\n    $file=~s/\\n[\
\\.:*\\s]+\\n/\\n\\n/g;\n    \n    #Remove White s\
paces before the sequence name\n    $file=~s/\\n[ \
\\t\\r\\f]+(\\b)/\\n\\1/g;\n    \n    \n    #Remov\
e Internal Blanks\n    $file=~s/(\\n\\S+)(\\s+)(\\\
S)/\\1_blank_\\3/g;\n    \n    $file=~s/[ ]//g;\n \
   $file=~s/_blank_/ /g;\n    \n    \n    #Identif\
y Double Blank lines\n    \n    $file =~s/\\n\\s*\\
\n/#/g;\n    \n    $file.=\"#\";\n    $file =~s/\\\
n/@/g;\n    \n    \n    \n    \n    #count nseq\n \
   @blocks=split /\\#/, $file;\n    shift (@blocks\
);\n    @s=split /\\@/, $blocks[0];\n    $nseq=$#s\
+1;\n    \n    #Merge all the sequences and split \
every Nseq\n    \n    \n    $file=join '@', @block\
s;\n    @lines=split /\\@/,$file;\n    \n    $c=0;\
\n    \n    foreach my $l (@lines)\n      {\n	my (\
$n, $s);\n	\n	if (!($l=~/\\S/)){next;}\n	elsif ($s\
tockholm && ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){ne\
xt;}#get read of STOCHOLM Terminator\n	\n	$l=~/(\\\
S+)\\s+(\\S*)/g;\n	$n=$1; $s=$2;\n	\n	$seq[$c].=$s\
;\n	$name[$c]=$n;\n	$c++;\n	\n	if ( $c==$nseq){$c=\
0;}\n	\n      } \n    \n    if ( $c!=0)\n      {\n\
	print STDERR \"ERROR: $fname is NOT an MSA in Clu\
stalw format: make sure there is no blank line wit\
hin a block [ERROR]\\n\";\n	exit (1);\n      }\n  \
  \n    \n    for (my $a=0; $a< $nseq; $a++)\n    \
  {\n	$name[$a]=cleanstring ($name[$a]);\n	$seq[$a\
]=cleanstring ($seq[$a]);\n	print \">$name[$a]\\n$\
seq[$a]\\n\";\n      }\n  }\nsub cleanstring\n    \
{\n      my $s=@_[0];\n      $s=~s/_hash_symbol_/\\
\#/g;\n      $s=~s/_arobase_symbol_/\\@/g;\n      \
$s=~s/[ \\t]//g;\n      return $s;\n    }\n\nsub c\
lean_cr\n  {\n    my $f=shift;\n    $f=~s/\\r\\n/\\
\n/g;\n    $f=~s/\\n\\r/\\n/g;\n    $f=~s/\\r\\r/\\
\n/g;\n    $f=~s/\\r/\\n/g;\n    return $f;\n  }\n\
\nsub file2n_lines\n    {\n      my $file=shift;\n\
      my $nl=shift;\n      my $ret;\n      my $F=n\
ew FileHandle;\n      my $n=0;\n      open ($F, $f\
ile);\n\n      while (<$F>)\n	{\n	  $ret.=$_;\n	  \
$n++;\n	  \n	  if ($n>=$n){close ($F); return $ret\
;}\n	}\n      close ($F);\n      return $ret;\n   \
 }\n","use strict;\nuse FileHandle;\nuse Env qw(HO\
ST);\nuse Env qw(HOME);\nuse Env qw(USER);\nmy %na\
me;\nmy $nseq;\nmy $fasta;\nif ($ARGV[2] eq \"-fas\
ta\"){$fasta=1;}\nmy $F= new FileHandle;\n\nopen (\
$F, $ARGV[1]);\nwhile(<$F>)\n  {\n    my $l=$_;\n \
   if ($l=~/^#/){;}\n    elsif (($l=~/\\d+\\s+\\d+\
\\s+(\\S+)\\s+(\\S+)/))\n      {\n	my $n=$1;\n	$na\
me{$1}++;\n      }\n  }\nclose ($F);\n\nopen ($F, \
$ARGV[0]);\nwhile(<$F>)\n  {\n    my $l=$_;\n    i\
f ($l=~/^#/){;}\n    elsif ($l=~/\\d+\\s+\\d+\\s+(\
\\S+)\\s+(\\S+)/)\n      {\n	my $n=$1;\n	$name{$n}\
++;\n	if ($name{$n}==2){$nseq++;}\n      }\n  }\nc\
lose ($F);\n\nif (!$fasta && $nseq>0)\n  {\n    pr\
int \"#NAMESEQ_01\\n\";\n    print \"# $nseq\\n\";\
\n  }\nopen ($F, $ARGV[0]);\nwhile(<$F>)\n  {\n   \
 my $l=$_;\n    if ($l=~/^#/){;}\n    elsif ($l=~/\
.\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)/)\n      {\n	my \
$n=$1;\n	my $s=$2;\n	if ($name{$n}==2)\n	  {\n	   \
 if ($fasta)\n	      {\n		print \">$n\\n$s\\n\";\n\
	      }\n	    else\n	      {\n		print \"$l\";\n	 \
     }\n	  }\n      }\n  }\nclose ($F);\nexit (0);\
\n\n\n","use Env qw(HOST);\nuse Env qw(HOME);\nuse\
 Env qw(USER);\n\n\n$query_start=-1;\n$query_end=-\
1;\n\nwhile (<>)\n  {\n    if ( /\\/\\//){$in_aln=\
1;}\n    elsif ( $in_aln && /(\\S+)\\s+(.*)/)\n   \
   {\n\n\n	$name=$1;\n	\n\n	$seq=$2;\n	$seq=~s/\\s\
//g;\n        $seq=~s/\\~/\\-/g;\n	$seq=~s/\\./\\-\
/g;\n	if ( $list{$n}{'name'} && $list{$n}{'name'} \
ne $name)\n	  {\n	    print \"$list{$n}{'name'} Vs\
 $name\";\n	    \n	    exit (EXIT_FAILURE);\n	  }\\
n	else\n	  {\n	    $list{$n}{'name'}= $name;\n	  }\
\n\n	$list{$n}{'seq'}=$list{$n}{'seq'}.$seq;\n	\n	\
$nseq=++$n;\n	\n      }\n    else\n      {$n=0;}\n\
  }\n\n\nfor ($a=0; $a<$nseq; $a++)\n  {\n    prin\
t \">$list{$a}{'name'}\\n$list{$a}{'seq'}\\n\";\n \
 }\n      \n","$run_anyway=2;\nmy $msaf=\"msa.in.t\
mp.$$\";\nmy $msaoutf=\"msa.out.tmp.$$\";\nmy $err\
=\"msa.out.err.$$\";\nopen  (F, $ARGV[0]);\nopen  \
(OUT, \">$msaf\");\n$nseq=0;\nwhile (<F>)\n  {\n  \
  $l=$_;\n    if ( $l=~/^>(\\S+)/)\n      {\n	$s=$\
seqname{$nseq++}=$1;\n	print OUT \"$l\";\n	\n     \
 }\n    else \n      {\n	$l=uc($l);\n	print OUT \"\
$l\";\n      }\n  }\n\nclose (F);\nclose(OUT);\n\n\
system (\"msa $msaf > $msaoutf 2>$err\");\nopen (F\
, \"$msaoutf\");\n$read=0;\n$cn=0;\nwhile (<F>)\n \
 {\n    $l=$_;\n    if ($read)\n      {\n	if ($l=~\
/End gaps not penalized/){$read=0;}\n	elsif (!($l=\
~/\\S/))\n	  {\n	    $cn=0;\n	  }\n	else\n	  {\n	 \
   \n	    chomp ($l);\n	    $seqal{$cn++}.=$l;\n	 \
   $tot++;\n	  }\n      }\n    elsif ($l=~/Optimal\
 Multiple Alignment/)\n      {\n	$read=1;\n      }\
\n  }\nclose (F);\n\nif ($tot<1 && $run_anyway==1)\
\n  {\n    print STDERR \"\\nWarning: MSA returned\
 a NULL file -- Use T-Coffee instead\\n\";\n    op\
en (F,$err);\n    while (<F>){print \"$_\";}\n    \
  \n    system (\"t_coffee -seq $msaf -outfile $AR\
GV[1]  -quiet\");\n  }\nelsif ($tot<1 && $run_anyw\
ay==2)\n  {\n    \n    \n    $nseq/=2;\n    $nseq=\
int ($nseq);\n    if ($nseq<2){$nseq=2;}\n    prin\
t \"RUN MSA with NSeq=$nseq\\n\";\n    #print (\"t\
_coffee -dpa -dpa_nseq $nseq -seq $ARGV[0] -dpa_tr\
ee codnd -outfile $ARGV[1] -dpa_method msa_msa\");\
\n    system (\"t_coffee -dpa -dpa_nseq $nseq -seq\
 $ARGV[0] -dpa_tree codnd -outfile $ARGV[1] -dpa_m\
ethod msa_msa>/dev/null\");\n\n  }\nelsif ($tot<1)\
\n  {\n    exit (EXIT_FAILURE);\n  }\nelse\n  {\n \
   open (OUT, \">$ARGV[1]\");\n    for ($a=0; $a<$\
nseq;$a++)\n      {\n	print OUT \">$seqname{$a}\\n\
$seqal{$a}\\n\";\n      }\n    close (OUT);\n  }\n\
\n\n\nunlink ($msaf);\nunlink ($msaoutf);\nunlink \
($err);\n","use strict;\nuse Cwd;\nuse File::Basen\
ame;\nmy $test=0;\n\nmy $tmpdir=\"/tmp/tco/aligner\
s/upp/\";\nmymkdir ($tmpdir);\n\n\n\nif ($ARGV[0] \
eq \"one\")\n  {\n    seq2msa ($ARGV[1], $ARGV[2])\
;\n  }\nelsif ($ARGV[0] eq \"all\")\n  {\n    list\
seq2listmsa ($ARGV[1]);\n  }\n\nsub listseq2listms\
a\n  {\n    my $list=shift;\n    my $cdir = getcwd\
;\n    my $dir=random_string(10);\n    $dir=\"$tmp\
dir/$dir/\";\n    my %h;\n    my $n;\n    mkdir  (\
$dir);\n\n    open (F, $list);\n    while (<F>)\n \
     {\n        my $l=$_;\n\n        chomp($l);\n \
       my @f=split (/\\s+/, $l);\n	if ( -e $f[0])\\
n          {\n            $h{$n}{in}=$f[0];\n     \
       ($h{$n}{name},$h{$n}{path})=fileparse ($f[0\
]);\n            $h{$n}{NFin}= \"$dir/$h{$n}{name}\
.seq\";\n	    \n            $h{$n}{NFout}=\"$dir/$\
h{$n}{name}.aln\";\n\n            $h{$n}{out}=$f[1\
];\n\n            fasta2fastaupp ($h{$n}{in}, $h{$\
n}{NFin});\n            $n++;\n          }\n      \
}\n    close (F);\n    chdir ($dir);\n    \n    if\
 (!$test)\n      {\n	system (\"fbname=\\$(basename\
 `ls *.seq` .seq); \\\n             run_upp.py -s \
\\${fbname}.seq -m amino --cpu 1 -d outdir -o \\${\
fbname}.aln; \\\n             mv outdir/\\${fbname\
}.aln_alignment.fasta \\${fbname}.aln;\");\n      \
}\n    \n    foreach my $n (keys (%h))\n      {\n	\
if ($test)\n	  {\n	    system (\"cp $h{$n}{NFin} $\
h{$n}{NFout}\");\n	    print \"$h{$n}{NFin} $h{$n}\
{NFout} $h{$n}{out}\\n\";\n	  }\n        fastaupp2\
fasta ($h{$n}{NFout},$h{$n}{out});\n      }\n    c\
hdir ($cdir);\n  }\n\nsub seq2msa\n    {\n      my\
 ($in, $out)=@_;\n      my $cdir=getcwd;\n      \n\
      \n      if (!($in=~/\\//)){$in=$cdir.\"/\".$\
in;}\n      if (!($out=~/\\//)){$out=$cdir.\"/\".$\
out;}\n      \n      my $file=random_string(10);\n\
      $file=\"$tmpdir/$file\";\n      open (F, \">\
$file\");\n      print F \"$in $out\\n\";\n      c\
lose (F);\n      \n      return listseq2listmsa ($\
file);\n    }\n	\nsub fasta2fastaupp\n  {\n    my \
($in, $out)=@_;\n    my ($name, $seq, $n);\n    \n\
    if (!-e $in){return;}\n    \n    open (IN, \"$\
in\");\n    open (OUT, \">$out\");\n    local $/ =\
 \"\\n>\";  # read by FASTA record\n    \n    whil\
e (<IN>)\n      {\n	my $l=$_;\n	$l=~s/>//g;\n	$l=\\
">\".$l;\n	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	$l=~s\
/^>*.+\\n//;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	$seq=\
~s/u/x/g;\n	$seq=~s/U/X/g;\n	print OUT \">$name\\n\
$seq\\n\";\n	$n++;\n      }\n    if ($n==2)\n     \
 {\n	print OUT \">fake_seq4upp\\n$seq\\n\";\n     \
 }\n    close (IN);\n    close (OUT);\n    local $\
/=\"\\n\";\n  }\n\nsub fastaupp2fasta\n  {\n    my\
 ($in, $out)=@_;\n    my ($name, $seq, $n);\n    \\
n    if (!-e $in){return;}\n    \n    open (IN, \"\
$in\");\n    open (OUT, \">$out\");\n    local $/ \
= \"\\n>\";  # read by FASTA record\n    \n    whi\
le (<IN>)\n      {\n	my $l=$_;\n	$l=~s/>//g;\n	$l=\
\">\".$l;\n	\n	$l=~/^>(.*)/;\n	$name=$1;\n	\n	$l=~\
s/^>*.+\\n//;\n	$l=~s/\\n//g;\n	$seq=$l;\n	\n	$seq\
=~s/x/u/g;\n	$seq=~s/X/U/g;\n	\n	if (!($name=~/fak\
e_seq4upp/))\n	  {\n	    print OUT \">$name\\n$seq\
\\n\";\n	  }\n      }\n    close (IN);\n    close \
(OUT);\n    local $/=\"\\n\";\n  }\n\nsub random_s\
tring\n    {\n      my $len=shift;\n      my @char\
s = (\"A\"..\"Z\", \"a\"..\"z\");\n      my $strin\
g;\n      $string .= $chars[rand @chars] for 1..$l\
en;\n      return $string;\n    }\n\nsub mymkdir\n\
      {\n	my $d=shift;\n	my $cd='/';\n	\n	foreach \
my $e (split (/\\//, $d))\n	  {\n	    $cd.=\"$e/\"\
;\n	    if ( !-d $cd){mkdir ($cd);}\n	  }\n	return\
;\n      }\n      \n			  \n      \n","use strict;\\
nuse Cwd;\nuse File::Basename;\n\n\nmy $tmpdir=\"/\
tmp/tco/aligners/clustalo/\";\nmymkdir ($tmpdir);\\
n\n\n\nif ($ARGV[0] eq \"one\")\n  {\n    seq2msa \
($ARGV[1], $ARGV[2]);\n  }\nelsif ($ARGV[0] eq \"a\
ll\")\n  {\n    listseq2listmsa ($ARGV[1]);\n  }\n\
\n\n\nsub listseq2listmsa\n  {\n    my $list=shift\
;\n    my $cdir = getcwd;\n    my $dir=random_stri\
ng(10);\n    $dir=\"$tmpdir/$dir/\";\n    my %h;\n\
    my $n;\n    mkdir  ($dir);\n    \n    open (F,\
 $list);\n    while (<F>)\n      {\n	my $l=$_;\n\n\
	chomp($l);\n	my @f=split (/\\s+/, $l);\n	#print \\
"$l: 0:$f[0], 1:$f[1]\\n\";\n	if ( -e $f[0])\n	  {\
\n	    $h{$n}{in}=$f[0];\n	    ($h{$n}{name},$h{$n\
}{path})=fileparse ($f[0]);\n	    $h{$n}{NFin}= \"\
$dir/$h{$n}{name}.seq4nf\";\n	    $h{$n}{NFout}=\"\
$dir/$h{$n}{name}.aln\";\n	    \n	    $h{$n}{out}=\
$f[1];\n	    \n	    translate_fasta_seq (\"uU\", \\
"X\",$h{$n}{in}, $h{$n}{NFin});\n	    $n++;\n	  }\\
n      }\n    close (F);\n    \n    \n    chdir ($\
dir);\n    dump_nf (\"nf\");\n    dump_config ();\\
n   \n    #system (\"nextflow run nf  --name \\'*.\
seq4nf\\' >/dev/null 2>/dev/null\");\n    system (\
\"nextflow run nf  --name \\'*.seq4nf\\'\");\n    \
foreach my $n (keys (%h))\n      {\n	translate_fas\
ta_seq (\"uU\", \"X\",$h{$n}{NFout},$h{$n}{out});\\
n      }\n    chdir ($cdir);\n  }\nsub seq2msa\n  \
  {\n      my ($in, $out)=@_;\n      my $cdir=getc\
wd;\n      \n      \n      if (!($in=~/\\//)){$in=\
$cdir.\"/\".$in;}\n      if (!($out=~/\\//)){$out=\
$cdir.\"/\".$out;}\n      \n      my $file=random_\
string(10);\n      $file=\"$tmpdir/$file\";\n     \
 open (F, \">$file\");\n      print F \"$in $out\\\
n\";\n      close (F);\n      \n      return lists\
eq2listmsa ($file);\n    }\n	\nsub seq2msa_old\n  \
{\n    my ($in, $out)=@_;\n    my $cdir = getcwd;\\
n    my $dir=random_string(10);\n    $dir=\"/tmp/u\
pp.nf4tcoffee/$dir\";\n    my $seq=random_string(1\
0);\n    $seq.=\".fa\";\n    my $aln=$seq;\n    $a\
ln.=\".aln\";\n    \n    mkdir ($dir);\n    transl\
ate_fasta_seq (\"uU\", \"X\",$in, \"$dir/$seq\");\\
n    chdir ($dir);\n    \n    dump_nf (\"nf\");\n \
   dump_config ();\n    print \"IN: $in OUT: $cdir\
/$out\\nDIR: $dir\\nnextflow run nf  --name \\'*.f\
a\\' \\n\";\n    system (\"nextflow run nf  --name\
 \\'*.fa\\' \");\n    print \"$dir/$aln $cdir/$out\
\\n\";\n    translate_fasta_seq (\"xX\", \"U\",$al\
n, \"$cdir/$out\");\n    chdir ($cdir);\n   } \nsu\
b translate_fasta_seq\n  {\n    my ($from, $to, $i\
n, $out)=@_;\n    my $n;\n    my $skip;\n    my $l\
;\n    my $cseq;\n    if (!-e $in){return;}\n    \\
n    open (IN, \"$in\");\n    open (OUT, \">$out\"\
);\n   \n    while (<IN>)\n      {\n	$l=$_;\n	if (\
$l=~\">\"){$n++;$cseq=\"\";}\n	else { $l=~s/[$from\
]/$to/;$cseq.=$l;}\n\n	if ($skip){$skip=0;}\n	elsi\
f ($l=~/>fake_seq/){$skip=1;}\n	else\n	  {\n	    p\
rint OUT \"$l\";\n	  }\n      }\n    if ($n==2 && \
$from eq \"uU\")\n      {\n	print OUT \">fake_seq\\
\n$cseq\\n\";\n      }\n    close (IN);\n    close\
 (OUT);\n  }\n\nsub dump_config\n    {\n      open\
 (F, \">nextflow.config\");\n\n      print F \"doc\
ker.enabled = true\\n\";\n      print F \"process.\
container = \\'cbcrg/benchfam_large_scale\\'\\n\";\
\n      close (F);\n    }\n\nsub dump_nf\n  {\n   \
 my $nff=shift;\n    open (F,\">$nff\");\n    prin\
t F \"#!/usr/bin/env nextflow\\n\";\n    print F \\
"params.base_dir=\\\"./\\\"\\n\";\n    print F \"p\
arams.out_dir=\\\"./\\\"\\n\";\n    print F \"Chan\
nel.fromPath(params.name)\\n\";\n    print F \"\\t\
.map{ tuple(it.baseName, it) }\\n\";\n    \n    pr\
int F \"\\t.set{ file_names_1 }\\n\";\n    print F\
 \"process clustalo_align{\\n\";\n    print F \"\\\
tpublishDir params.out_dir, mode: \\\"copy\\\"\\n\\
";\n    print F \"tag \\\"\\${name}\\\"\";\n    pr\
int F \"\\n\";\n    print F \"\\tinput:\\n\";\n   \
 print F \"\\tset name, file(seq_file) from file_n\
ames_1\\n\";\n    print F \"\\toutput:\\n\";\n    \
print F \"\\tfile \\\"\\${name}.aln\\\"\\n\";\n   \
 print F \"\\n\";\n    print F \" \\\"\\\"\\\"\\n\\
";\n    print F \" clustalo -i \\$seq_file -o \\${\
name}.aln\\n\";\n    print F \"\\\"\\\"\\\"\\n\\n\\
";\n    print F \"}\\n\";\n    close (F);\n  }\n\n\
sub random_string\n    {\n      my $len=shift;\n  \
    my @chars = (\"A\"..\"Z\", \"a\"..\"z\");\n   \
   my $string;\n      $string .= $chars[rand @char\
s] for 1..$len;\n      return $string;\n    }\n\ns\
ub mymkdir\n      {\n	my $d=shift;\n	my $cd='/';\n\
	\n	foreach my $e (split (/\\//, $d))\n	  {\n	    \
$cd.=\"$e/\";\n	    if ( !-d $cd){mkdir ($cd);}\n	\
  }\n	return;\n      }\n      \n			  \n      \n","\
\nmy $msaf=\"msa.in.tmp.$$\";\nmy $msaoutf=\"msa.o\
ut.tmp.$$\";\nmy $cost=\"blosum62.tmp.$$\";\n\nope\
n  (F, $ARGV[0]);\nopen  (OUT, \">$msaf\");\n$nseq\
=0;\nwhile (<F>)\n  {\n    $l=$_;\n    if ( $l=~/^\
>(\\S+)/)\n      {\n	my $simple=\"Seq$nseq\";\n	$s\
=$seqname{$nseq++}=$1;\n	$translate{$simple}=$s;\n\
	\n	print OUT \">$simple\\n\";\n	\n      }\n    el\
se\n      {\n	$l=uc($l);\n	print OUT \"$l\";\n    \
  }\n  }\nclose (F);\nclose(OUT);\n\ndump_blosum (\
$cost);\nsystem (\"dca -c $cost -q $msaf> $msaoutf\
 2>/dev/null\");\nopen (F, \"$msaoutf\");\nopen (O\
UT, \">$ARGV[1]\");\n\n$read=0;\nwhile (<F>)\n  {\\
n    $l=$_;\n    if ($l=~/^>(\\S+)/)\n      {\n	$r\
ead=1;\n	$name=$translate{$1};\n	print OUT \">$nam\
e\\n\";\n      }\n    elsif ($read && ($l=~/\\S/))\
\n      {\n	print OUT \"$l\";\n      }\n    else\n\
      {\n	$read=0;\n      }\n  }\nclose (F);\n\nun\
link ($cost);\nunlink ($msaf);\nunlink ($msaoutf);\
\n\nsub dump_blosum\n  {\n    my $f=shift;\n    op\
en (F, \">$f\");\n\n    print F \"6\\n\";\n    pri\
nt F \"- -   0\\n\";\n    print F \"W W   0\\n\";\\
n    print F \"Y Y   4\\n\";\n    print F \"F F   \
5\\n\";\n    print F \"V V   7\\n\";\n    print F \
\"L L   7\\n\";\n    print F \"I I   7\\n\";\n    \
print F \"M M   6\\n\";\n    print F \"K K   6\\n\\
";\nprint F \"R R   6\\n\";\n    print F \"H H   3\
\\n\";\n    print F \"Q Q   6\\n\";\n    print F \\
"E E   6\\n\";\n    print F \"D D   5\\n\";\n    p\
rint F \"N N   5\\n\";\n    print F \"G G   5\\n\"\
;\n    print F \"A A   7\\n\";\n    print F \"P P \
  4\\n\";\n    print F \"T T   6\\n\";\n    print \
F \"S S   7\\n\";\n    print F \"C C   2\\n\";\n  \
  print F \"- C  10 \\n\";\n    print F \"- S  10\\
\n\";\n    print F \"- T  10 \\n\";\n    print F \\
"- P  10\\n\";\n    print F \"- A  10 \\n\";\n    \
print F \"- G  10\\n\";\n    print F \"- N  10 \\n\
\";\n    print F \"- D  10\\n\";\n    print F \"- \
E  10 \\n\";\n    print F \"- Q  10\\n\";\nprint F\
 \"- H  10 \\n\";\n    print F \"- R  10\\n\";\n  \
  print F \"- K  10 \\n\";\n    print F \"- M  10\\
\n\";\n    print F \"- I  10 \\n\";\n    print F \\
"- L  10\\n\";\n    print F \"- V  10 \\n\";\n    \
print F \"- F  10\\n\";\n    print F \"- Y  10 \\n\
\";\n    print F \"- W  10\\n\";\n    print F \"W \
C  13 \\n\";\n    print F \"W S  14\\n\";\n    pri\
nt F \"W T  13 \\n\";\n    print F \"W P  15\\n\";\
\n    print F \"W A  14 \\n\";\n    print F \"W G \
 13\\n\";\n    print F \"W N  15 \\n\";\n    print\
 F \"W D  15\\n\";\n    print F \"W E  14 \\n\";\n\
    print F \"W Q  13\\n\";\n    print F \"W H  13\
 \\n\";\n    print F \"W R  14\\n\";\n    print F \
\"W K  14 \\n\";\n    print F \"W M  12\\n\";\n   \
 print F \"W I  14 \\n\";\n    print F \"W L  13\\\
n\";\n    print F \"W V  14 \\n\";\n    print F \"\
W F  10\\n\";\n    print F \"W Y   9 \\n\";\n    p\
rint F \"Y C  13\\n\";\n    print F \"Y S  13 \\n\\
";\n    print F \"Y T  13\\n\";\n    print F \"Y P\
  14 \\n\";\n    print F \"Y A  13\\n\";\n    prin\
t F \"Y G  14 \\n\";\n    print F \"Y N  13\\n\";\\
n    print F \"Y D  14 \\n\";\n    print F \"Y E  \
13\\n\";\n    print F \"Y Q  12 \\n\";\n    print \
F \"Y H   9\\n\";\n    print F \"Y R  13 \\n\";\n \
   print F \"Y K  13\\n\";\n    print F \"Y M  12 \
\\n\";\n    print F \"Y I  12\\n\";\n    print F \\
"Y L  12 \\n\";\n    print F \"Y V  12\\n\";\n    \
print F \"Y F   8 \\n\";\n    print F \"F C  13\\n\
\";\nprint F \"F S  13 \\n\";\n    print F \"F T  \
13\\n\";\n    print F \"F P  15 \\n\";\n    print \
F \"F A  13\\n\";\n    print F \"F G  14 \\n\";\n \
   print F \"F N  14\\n\";\n    print F \"F D  14 \
\\n\";\n    print F \"F E  14\\n\";\n    print F \\
"F Q  14 \\n\";\n    print F \"F H  12\\n\";\n    \
print F \"F R  14 \\n\";\n    print F \"F K  14\\n\
\";\n    print F \"F M  11 \\n\";\n    print F \"F\
 I  11\\n\";\n    print F \"F L  11 \\n\";\n    pr\
int F \"F V  12\\n\";\n    print F \"V C  12 \\n\"\
;\n    print F \"V S  13\\n\";\n    print F \"V T \
 11 \\n\";\n    print F \"V P  13\\n\";\n    print\
 F \"V A  11 \\n\";\n    print F \"V G  14\\n\";\n\
    print F \"V N  14 \\n\";\n    print F \"V D  1\
4\\n\";\nprint F \"V E  13 \\n\";\nprint F \"V Q  \
13\\n\";\nprint F \"V H  14 \\n\";\nprint F \"V R \
 14\\n\";\nprint F \"V K  13 \\n\";\nprint F \"V M\
  10\\n\";\nprint F \"V I   8 \\n\";\nprint F \"V \
L  10\\n\";\nprint F \"L C  12 \\n\";\nprint F \"L\
 S  13\\n\";\nprint F \"L T  12 \\n\";\nprint F \"\
L P  14\\n\";\nprint F \"L A  12 \\n\";\nprint F \\
"L G  15\\n\";\nprint F \"L N  14 \\n\";\nprint F \
\"L D  15\\n\";\nprint F \"L E  14 \\n\";\nprint F\
 \"L Q  13\\n\";\nprint F \"L H  14 \\n\";\nprint \
F \"L R  13\\n\";\nprint F \"L K  13 \\n\";\nprint\
 F \"L M   9\\n\";\nprint F \"L I   9 \\n\";\nprin\
t F \"I C  12\\n\";\nprint F \"I S  13 \\n\";\npri\
nt F \"I T  12\\n\";\nprint F \"I P  14 \\n\";\npr\
int F \"I A  12\\n\";\nprint F \"I G  15 \\n\";\np\
rint F \"I N  14\\n\";\nprint F \"I D  14 \\n\";\n\
print F \"I E  14\\n\";\nprint F \"I Q  14 \\n\";\\
nprint F \"I H  14\\n\";\nprint F \"I R  14 \\n\";\
\nprint F \"I K  14\\n\";\nprint F \"I M  10 \\n\"\
;\nprint F \"M C  12\\n\";\nprint F \"M S  12 \\n\\
";\nprint F \"M T  12\\n\";\nprint F \"M P  13 \\n\
\";\nprint F \"M A  12\\n\";\nprint F \"M G  14 \\\
n\";\nprint F \"M N  13\\n\";\nprint F \"M D  14 \\
\n\";\nprint F \"M E  13\\n\";\nprint F \"M Q  11 \
\\n\";\nprint F \"M H  13\\n\";\nprint F \"M R  12\
 \\n\";\nprint F \"M K  12\\n\";\nprint F \"K C  1\
4 \\n\";\nprint F \"K S  11\\n\";\nprint F \"K T  \
12 \\n\";\nprint F \"K P  12\\n\";\nprint F \"K A \
 12 \\n\";\nprint F \"K G  13\\n\";\nprint F \"K N\
  11 \\n\";\nprint F \"K D  12\\n\";\nprint F \"K \
E  10 \\n\";\nprint F \"K Q  10\\n\";\nprint F \"K\
 H  12 \\n\";\nprint F \"K R   9\\n\";\nprint F \"\
R C  14 \\n\";\nprint F \"R S  12\\n\";\nprint F \\
"R T  12 \\n\";\nprint F \"R P  13\\n\";\nprint F \
\"R A  12 \\n\";\nprint F \"R G  13\\n\";\nprint F\
 \"R N  11 \\n\";\nprint F \"R D  13\\n\";\nprint \
F \"R E  11 \\n\";\nprint F \"R Q  10\\n\";\nprint\
 F \"R H  11 \\n\";\nprint F \"H C  14\\n\";\nprin\
t F \"H S  12 \\n\";\nprint F \"H T  13\\n\";\npri\
nt F \"H P  13 \\n\";\nprint F \"H A  13\\n\";\npr\
int F \"H G  13 \\n\";\nprint F \"H N  10\\n\";\np\
rint F \"H D  12 \\n\";\nprint F \"H E  11\\n\";\n\
print F \"H Q  11 \\n\";\nprint F \"Q C  14\\n\";\\
nprint F \"Q S  11 \\n\";\nprint F \"Q T  12\\n\";\
\nprint F \"Q P  12 \\n\";\nprint F \"Q A  12\\n\"\
;\nprint F \"Q G  13 \\n\";\nprint F \"Q N  11\\n\\
";\nprint F \"Q D  11 \\n\";\nprint F \"Q E   9\\n\
\";\nprint F \"E C  15 \\n\";\nprint F \"E S  11\\\
n\";\nprint F \"E T  12 \\n\";\nprint F \"E P  12\\
\n\";\nprint F \"E A  12 \\n\";\nprint F \"E G  13\
\\n\";\nprint F \"E N  11 \\n\";\nprint F \"E D   \
9\\n\";\nprint F \"D C  14 \\n\";\nprint F \"D S  \
11\\n\";\nprint F \"D T  12 \\n\";\nprint F \"D P \
 12\\n\";\nprint F \"D A  13 \\n\";\nprint F \"D G\
  12\\n\";\nprint F \"D N  10 \\n\";\nprint F \"N \
C  14\\n\";\nprint F \"N S  10 \\n\";\nprint F \"N\
 T  11\\n\";\nprint F \"N P  13 \\n\";\nprint F \"\
N A  13\\n\";\nprint F \"N G  11 \\n\";\nprint F \\
"G C  14\\n\";\nprint F \"G S  11 \\n\";\nprint F \
\"G T  13\\n\";\nprint F \"G P  13 \\n\";\nprint F\
 \"G A  11\\n\";\nprint F \"A C  11 \\n\";\nprint \
F \"A S  10\\n\";\nprint F \"A T  11 \\n\";\nprint\
 F \"A P  12\\n\";\nprint F \"P C  14 \\n\";\nprin\
t F \"P S  12\\n\";\nprint F \"P T  12 \\n\";\npri\
nt F \"T C  12\\n\";\nprint F \"T S  10 \\n\";\npr\
int F \"S C  12\\n\";\nclose (F);\n    return;\n  \
}\n    \n","\nuse Env qw(HOST);\nuse Env qw(HOME);\
\nuse Env qw(USER);\n\n                           \
                             \nuse strict;        \
                                     \nuse warning\
s;\nuse diagnostics;\n\nmy $in_hit_list, my $in_al\
n=0, my(%name_list)=(),my (%list)=(),my $n_seq=0; \
my $test=0;\nmy($j)=0, my $n=0, my $nom, my $lg_qu\
ery, my %vu=();\n\nopen (F, \">tmp\");\n\n$/=\"\\n\
\";\nwhile (<>)\n{\n    print F $_;\n    if($_ =~ \
/Query=\\s*(.+?)\\s/i) { $nom=$1;}\n\n    if ( /Se\
quences producing significant alignments/){$in_hit\
_list=1;}\n    \n    if ($_=~ /^pdb\\|/i) { $_=~ s\
/pdb\\|//g; }\n    if ($_=~ /^(1_\\d+)\\s+\\d+/) {\
 $_=~ s/$1/QUERY/;}\n      \n    if ( /^(\\S+).+?\\
\s+[\\d.]+\\s+([\\de.-]+)\\s+$/ && $in_hit_list)	\\
n    {\n	my($id)=$1; # \n	$id=~ s/\\|/_/g; #\n	if \
($id =~ /.+_$/) { chop($id) }; #\n	$name_list{$n_s\
eq++}=$id;\n	$name_list{$n_seq-1}=~ s/.*\\|//g;   \
  \n    }\n  \n    if (/query/i) {$in_aln=1;}\n   \
 if ( /^(\\S+)\\s+(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)\
/ || /^(\\S+)(\\s+)(\\-+)(\\s+)/ && ($in_aln == 1)\
)\n    {\n	my $name=$1;\n	my $start=$2;\n	my $seq=\
$3;\n	my $end=$4;\n		\n	if ($name =~ /QUERY/i) { $\
lg_query=length($seq); }\n\n	unless ($test > $n) #\
m\n	{\n	    my(@seqq)= split('',$seq);\n	    my($g\
ap_missing)= scalar(@seqq);\n	    \n	    while ($g\
ap_missing != $lg_query)  { unshift (@seqq,\"-\");\
 $gap_missing= scalar(@seqq); }\n	    $seq=join(''\
,@seqq);  #m\n	}\n	\n	if ($name =~ /QUERY/i)\n	{\n\
	    $n=0; %vu=(); $j=0;\n	    $list{$n}{'real_nam\
e'}=\"$nom\";\n	}	\n	else\n	{\n	    unless (exists\
 $vu{$name}) { ++$j;}	\n	    $list{$n}{'real_name'\
}=$name_list{$j-1};\n	}\n		\n	$list{$n}{'name'}=$n\
ame;\n\n	$seq=~tr/a-z/A-Z/;\n	$list{$n}{'seq'}=$li\
st{$n}{'seq'};\n	$list{$n}{'seq'}.=$seq;\n\n	$n++;\
\n	$vu{$name}++;\n	$test++;\n   } \n    \n}\n\nmy \
@numero=();\n\nfor (my $a=0; $a<$n; $a++) #m\n{\n \
   my $long=length($list{0}{'seq'});  \n    my $lo\
ng1= length($list{$a}{'seq'});\n  \n    while ($lo\
ng1 ne $long)\n    {\n	$list{$a}{'seq'}.=\"-\";\n	\
$long1= length ($list{$a}{'seq'});\n    } \n \n   \
 push (@numero,\"$list{$a}{'name'} $list{$a}{'real\
_name'}\\n\");\n}\n\nmy %dejavu=();\n\n\nfor (my $\
i=0; $i<=$#numero; $i++)\n{\n    my $s=\">$list{$i\
}{'real_name'}\\n$list{$i}{'seq'}\\n\";\n    my $k\
=0;\n    \n    if (exists $dejavu{$numero[$i]}) {n\
ext;}\n    else\n    {	\n	for ($j=0; $j<$n ; $j++)\
\n	{\n	    if (\"$numero[$i]\" eq \"$numero[$j]\" \
&& $j != $i )\n	    {\n		++$k;\n		$s .=\">$list{$j\
}{'real_name'}\\n$list{$j}{'seq'}\\n\";\n	    }\n	\
}	\n    }\n    \n    if ($k>0) \n    {\n	my $cons;\
\n	open (SOR,\">tempo_aln2cons\"); print SOR $s;  \
close SOR ;\n	open (COM,\"t_coffee -other_pg seq_r\
eformat -in tempo_aln2cons -action +aln2cons +uppe\
r |\") ; \n     	while (<COM>)\n	{	\n	    if (/^>/\
) { $cons =\">$list{$i}{'real_name'}\\n\"; next;}\\
n	    $_=~ s/\\n//g;\n	    $cons .=$_;\n	}\n	close\
 COM; unlink (\"tempo_aln2cons\");\n	print $cons,\\
"\\n\"; print F $cons,\"\\n\";\n    }	\n    else  \
{ print $s;  print F $s; }\n    \n    $dejavu{$num\
ero[$i]}++;\n} #m\n\nexit;\n\n\n\n\n\n\n\n\n\n\n\n\
","use Env;\n\n\n$tmp_dir=\"\";\n$init_dir=\"\";\n\
$program=\"tc_generic_method.pl\";\n\n$blast=@ARGV\
[0];\n\n$name=\"query\";$seq=\"\";\n%p=blast_xml2p\
rofile($name,$seq,100, 0, 0, $blast);\n&output_pro\
file (%p);\n\n\nsub output_profile\n  {\n    my (%\
profile)=(@_);\n    my ($a);\n    for ($a=0; $a<$p\
rofile{n}; $a++)\n      {\n	\n	print \">$profile{$\
a}{name} $profile{$a}{comment}\\n$profile{$a}{seq}\
\\n\";\n      }\n    return;\n  }\nsub file_contai\
ns \n  {\n    my ($file, $tag, $max)=(@_);\n    my\
 ($n);\n    $n=0;\n    \n    if ( !-e $file && ($f\
ile =~/$tag/)) {return 1;}\n    elsif ( !-e $file)\
{return 0;}\n    else \n      {\n	open (FC, \"$fil\
e\");\n	while ( <FC>)\n	  {\n	    if ( ($_=~/$tag/\
))\n	      {\n		close (FC);\n		return 1;\n	      }\
\n	    elsif ($max && $n>$max)\n	      {\n		close \
(FC);\n		return 0;\n	      }\n	    $n++;\n	  }\n  \
    }\n    close (FC);\n    return 0;\n  }\n	    \\
n	  \nsub file2string\n  {\n    my $f=@_[0];\n    \
my $string, $l;\n    open (F,\"$f\");\n    while (\
<F>)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$string.\
=$l;\n      }\n    close (F);\n    $string=~s/\\r\\
\n//g;\n    $string=~s/\\n//g;\n    return $string\
;\n  }\n\n\n\nsub tag2value \n  {\n    \n    my $t\
ag=(@_[0]);\n    my $word=(@_[1]);\n    my $return\
;\n    \n    $tag=~/$word=\"([^\"]+)\"/;\n    $ret\
urn=$1;\n    return $return;\n  }\n      \nsub hit\
_tag2pdbid\n  {\n    my $tag=(@_[0]);\n    my $pdb\
id;\n       \n    $tag=~/id=\"(\\S+)\"/;\n    $pdb\
id=$1;\n    $pdbid=~s/_//;\n    return $pdbid;\n  \
}\nsub id2pdbid \n  {\n    my $id=@_[0];\n  \n    \
if ($id =~/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	$id=\
$1;\n      }\n    $id=~s/[|�_]//g;\n    return $id\
;\n  }\nsub set_blast_type \n  {\n    my $file =@_\
[0];\n    if (&file_contains ($file,\"EBIApplicati\
onResult\",100)){$BLAST_TYPE=\"EBI\";}\n    elsif \
(&file_contains ($file,\"NCBI_BlastOutput\",100)) \
{$BLAST_TYPE=\"NCBI\";}\n    else\n      {\n	$BLAS\
T_TYPE=\"\";\n      }\n    return $BLAST_TYPE;\n  \
}\nsub blast_xml2profile \n  {\n    my ($name,$seq\
,$maxid, $minid, $mincov, $file)=(@_);\n    my (%p\
, $a, $string, $n);\n    \n\n\n    if ($BLAST_TYPE\
 eq \"EBI\" || &file_contains ($file,\"EBIApplicat\
ionResult\",100)){%p=ebi_blast_xml2profile(@_);}\n\
    elsif ($BLAST_TYPE eq \"NCBI\" || &file_contai\
ns ($file,\"NCBI_BlastOutput\",100)){%p=ncbi_blast\
_xml2profile(@_);}\n    else \n      {\n	print \"*\
*********** ERROR: Blast Returned an unknown XML F\
ormat **********************\";\n	die;\n      }\n \
   for ($a=0; $a<$p{n}; $a++)\n      {\n	my $name=\
$p{$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n     \
 }\n    return %p;\n  }\nsub ncbi_blast_xml2profil\
e \n  {\n    my ($name,$seq,$maxid, $minid, $minco\
v, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhi\
ts,@identifyerL);\n    \n    \n    $seq=~s/[^a-zA-\
Z]//g;\n    $L=length ($seq);\n    \n    %hit=&xml\
2tag_list ($string, \"Hit\");\n    \n    \n    for\
 ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\n	my (\
$ldb,$id, $identity, $expectation, $start, $end, $\
coverage, $r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\"\
;\n\n	%ID=&xml2tag_list ($hit{$a}{body}, \"Hit_id\\
");\n	$identifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_\
list ($hit{$a}{body}, \"Hit_def\");\n	$definition=\
$DE{0}{body};\n	\n	%HSP=&xml2tag_list ($hit{$a}{bo\
dy}, \"Hsp\");\n	for ($b=0; $b<$HSP{n}; $b++)\n	  \
{\n	    my (%START,%END,%E,%I,%Q,%M);\n\n	 \n	    \
%START=&xml2tag_list ($HSP{$b}{body}, \"Hsp_query-\
from\");\n	    %HSTART=&xml2tag_list ($HSP{$b}{bod\
y}, \"Hsp_hit-from\");\n	    \n	    %LEN=  &xml2ta\
g_list ($HSP{$b}{body}, \"Hsp_align-len\");\n	    \
%END=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_query-\
to\");\n	    %HEND=  &xml2tag_list ($HSP{$b}{body}\
, \"Hsp_hit-to\");\n	    %E=&xml2tag_list     ($HS\
P{$b}{body}, \"Hsp_evalue\");\n	    %I=&xml2tag_li\
st     ($HSP{$b}{body}, \"Hsp_identity\");\n	    %\
Q=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_qseq\")\
;\n	    %M=&xml2tag_list     ($HSP{$b}{body}, \"Hs\
p_hseq\");\n	    \n	    for ($e=0; $e<$Q{n}; $e++)\
\n\n	      {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{bo\
dy};\n		if ($seq eq\"\"){$seq=$qs;$L=length($seq);\
}\n		\n		$expectation=$E{$e}{body};\n		$identity=(\
$LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*1\
00;\n		$start=$START{$e}{body};\n		$end=$END{$e}{b\
ody};\n		$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND\
{$e}{body};\n	\n		$coverage=(($end-$start)*100)/$L\
;\n\n	\n		if ($identity>$maxid || $identity<$minid\
 || $coverage<$mincov){next;}\n		@lr1=(split (//,$\
qs));\n		@lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		\
for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		fo\
r ($d=0,$c=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$\
c];\n		    if ( $r=~/[A-Za-z]/)\n		      {\n			\n	\
		$p[$nhits][$d + $start-1]=$lr2[$c];\n			$d++;\n	\
	      }\n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nh\
its]=$ms;\n		$QstartL[$nhits]=$start;\n		$HstartL[\
$nhits]=$Hstart;\n		$identityL[$nhits]=$identity;\\
n		$endL[$nhits]=$end;\n		$definitionL[$nhits]=$de\
finition;\n		$identifyerL[$nhits]=$identifyer;\n		\
$comment[$nhits]=\"$ldb|$identifyer [Eval=$expecta\
tion][id=$identity%][start=$Hstart end=$Hend]\";\n\
		$nhits++;\n	      }\n	  }\n      }\n    \n    $p\
rofile{n}=0;\n    $profile{$profile{n}}{name}=$nam\
e;\n    $profile{$profile{n}}{seq}=$seq;\n    $pro\
file {n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\
\n      {\n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$na\
me\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n\
}{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\
\n	$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$\
n}{Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=\
$identityL[$a];\n	$profile{$n}{definition}=$defini\
tionL[$a];\n	$profile{$n}{identifyer}=$identifyerL\
[$a];\n	$profile{$n}{comment}=$comment[$a];\n	for \
($b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	\
      {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      \
}\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\"\
;\n	      }\n	  }\n      }\n    \n    $profile{n}=\
$nhits+1;\n    return %profile;\n  }\nsub ebi_blas\
t_xml2profile \n  {\n    my ($name,$seq,$maxid, $m\
inid, $mincov, $string)=(@_);\n    my ($L,$l, $a,$\
b,$c,$d,$nhits,@identifyerL,$identifyer);\n    \n\\
n    \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($\
seq);\n    %hit=&xml2tag_list ($string, \"hit\");\\
n    \n    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n\
      {\n	my ($ldb,$id, $identity, $expectation, $\
start, $end, $coverage, $r);\n	my (%Q,%M,%E,%I);\n\
	\n	$ldb=&tag2value ($hit{$a}{open}, \"database\")\
;\n	$identifyer=&tag2value ($hit{$a}{open}, \"id\"\
);\n\n	$description=&tag2value ($hit{$a}{open}, \"\
description\");\n	\n	%Q=&xml2tag_list ($hit{$a}{bo\
dy}, \"querySeq\");\n	%M=&xml2tag_list ($hit{$a}{b\
ody}, \"matchSeq\");\n	%E=&xml2tag_list ($hit{$a}{\
body}, \"expectation\");\n	%I=&xml2tag_list ($hit{\
$a}{body}, \"identity\");\n	\n\n	for ($b=0; $b<$Q{\
n}; $b++)\n	  {\n	    \n	    \n	    $qs=$Q{$b}{bod\
y};\n	    $ms=$M{$b}{body};\n	    if ($seq eq\"\")\
{$seq=$qs;$L=length($seq);}\n\n	    $expectation=$\
E{$b}{body};\n	    $identity=$I{$b}{body};\n	    \\
n	    	    \n	    $start=&tag2value ($Q{$b}{open},\
 \"start\");\n	    $end=&tag2value ($Q{$b}{open}, \
\"end\");\n	    $startM=&tag2value ($M{$b}{open}, \
\"start\");\n	    $endM=&tag2value ($M{$b}{open}, \
\"end\");\n	    $coverage=(($end-$start)*100)/$L;\\
n	    \n	   # print \"$id: ID: $identity COV: $cov\
erage [$start $end]\\n\";\n	    \n	    \n	    if (\
$identity>$maxid || $identity<$minid || $coverage<\
$mincov){next;}\n	    # print \"KEEP\\n\";\n\n	   \
 \n	    @lr1=(split (//,$qs));\n	    @lr2=(split (\
//,$ms));\n	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;\
$c++){$p[$nhits][$c]=\"-\";}\n	    for ($d=0,$c=0;\
 $c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n		if ( $r\
=~/[A-Za-z]/)\n		  {\n		    \n		    $p[$nhits][$d \
+ $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	      \
}\n	  \n	    \n	    $identifyerL[$nhits]=$identify\
er;\n	    $comment[$nhits]=\"$ldb|$identifyer [Eva\
l=$expectation][id=$identity%][start=$startM end=$\
endM]\";\n	    $nhits++;\n	  }\n      }\n    \n   \
 $profile{n}=0;\n    $profile{$profile{n}}{name}=$\
name;\n    $profile{$profile{n}}{seq}=$seq;\n    $\
profile {n}++;\n    \n    for ($a=0; $a<$nhits; $a\
++)\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$na\
me\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n\
}{identifyer}=$identifyerL[$a];\n	\n	$profile{$n}{\
comment}=$comment[$a];\n	for ($b=0; $b<$L; $b++)\n\
	  {\n	    if ($p[$a][$b])\n	      {\n		$profile{$\
n}{seq}.=$p[$a][$b];\n	      }\n	    else\n	      \
{\n		$profile{$n}{seq}.=\"-\";\n	      }\n	  }\n  \
    }\n    $profile{n}=$nhits+1;\n    \n    return\
 %profile;\n  }\n\nsub blast_xml2hit_list\n  {\n  \
  my $string=(@_[0]);\n    return &xml2tag_list ($\
string, \"hit\");\n  }\nsub xml2tag_list  \n  {\n \
   my ($string_in,$tag)=@_;\n    my $tag_in, $tag_\
out;\n    my %tag;\n    \n    if (-e $string_in)\n\
      {\n	$string=&file2string ($string_in);\n    \
  }\n    else\n      {\n	$string=$string_in;\n    \
  }\n    $tag_in1=\"<$tag \";\n    $tag_in2=\"<$ta\
g>\";\n    $tag_out=\"/$tag>\";\n    $string=~s/>/\
>##1/g;\n    $string=~s/</##2</g;\n    $string=~s/\
##1/<#/g;\n    $string=~s/##2/#>/g;\n    @l=($stri\
ng=~/(\\<[^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;\
$n=-1;\n  \n \n\n    foreach $t (@l)\n      {\n\n	\
$t=~s/<#//;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_in1/\
 || $t=~/$tag_in2/)\n	  {\n	 \n	    $in=1;\n	    $\
tag{$tag{n}}{open}=$t;\n	    $n++;\n	    \n	  }\n	\
elsif ($t=~/$tag_out/)\n	  {\n	    \n\n	    $tag{$\
tag{n}}{close}=$t;\n	    $tag{n}++;\n	    $in=0;\n\
	  }\n	elsif ($in)\n	  {\n	   \n	    $tag{$tag{n}}\
{body}.=$t;\n	  }\n      }\n  \n    return %tag;\n\
  }\n\n\n\n\n","use Env qw(HOST);\nuse Env qw(HOME\
);\nuse Env qw(USER);\nwhile (<>)\n  {\n    if ( /\
^>(\\S+)/)\n      {\n	if ($list{$1})\n	  {\n	    p\
rint \">$1_$list{$1}\\n\";\n	    $list{$1}++;\n	  \
}\n	else\n	  {\n	    print $_;\n	    $list{$1}=1;\\
n	  }\n      }\n    else\n      {\n	print $_;\n   \
   }\n  }\n      \n","\n\n\nuse Env qw(HOST);\nuse\
 Env qw(HOME);\nuse Env qw(USER);\n\n\nopen (F,$AR\
GV[0]);\nwhile ( <>)\n  {\n    @x=/([^:,;\\)\\(\\s\
]+):[^:,;\\)\\(]*/g;\n    @list=(@list,@x);\n  }\n\
$n=$#list+1;\nforeach $n(@list){print \">$n\\nsequ\
ence\\n\";}\n\n\nclose (F);\n","\nopen (F, $ARGV[0\
]);\n\nwhile ( <F>)\n  {\n    @l=($_=~/(\\S+)/g);\\
n    \n    $name=shift @l;\n    \n    print STDOUT\
 \"\\n>$name\\n\";\n    foreach $e (@l){$e=($e eq \
\"0\")?\"O\":\"I\";print \"$e\";}\n  }\nclose (F);\
\n\n		       \n    \n","use strict;\nuse FileHandl\
e;\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env \
qw(USER);\nmy %name;\nmy $nseq;\nmy $F= new FileHa\
ndle;\nopen ($F, $ARGV[0]);\nwhile(<$F>)\n  {\n   \
 \n    my $l=$_;\n    if ($l=~/^#/){;}\n    elsif \
(($l=~/\\d+\\s+\\d+\\s+(\\S+)\\s+(\\S+)/))\n      \
{\n	my $name=$1;\n	my $seq=$2;\n	print \">$name\\n\
$seq\\n\";\n      }\n  }\nclose ($F);\nexit (0);\n\
\n\n","use Env qw(HOST);\nuse Env qw(HOME);\nuse E\
nv qw(USER);\n\n$tmp=\"$ARGV[0].$$\";\nopen (IN, $\
ARGV[0]);\nopen (OUT, \">$tmp\");\n\nwhile ( <IN>)\
\n  {\n    $file=$_;\n    $file=~s/\\r\\n/\\n/g;\n\
    $file=~s/\\n\\r/\\n/g;\n    $file=~s/\\r\\r/\\\
n/g;\n    $file=~s/\\r/\\n/g;\n    print OUT \"$fi\
le\";\n  }\nclose (IN);\nclose (OUT);\n\nopen (OUT\
, \">$ARGV[0]\");\nopen (IN, \"$tmp\");\n\nwhile (\
 <IN>)\n{\n  print OUT \"$_\";\n}\nclose (IN);\ncl\
ose (OUT);\nunlink ($tmp);\n\n"};
