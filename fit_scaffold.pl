#!/usr/bin/perl

# define the reasonable pdb atom name
%index = (
	"GLY_N", 1, "GLY_CA", 2, "GLY_C", 3, "GLY_O", 4,
	"ALA_N", 5, "ALA_CA", 6, "ALA_C", 7, "ALA_O", 8, "ALA_CB", 9,
	"VAL_N", 10, "VAL_CA", 11, "VAL_C", 12, "VAL_O", 13, "VAL_CB", 14, "VAL_CG1", 15, "VAL_CG2", 16,
	"LEU_N", 17, "LEU_CA", 18, "LEU_C", 19, "LEU_O", 20, "LEU_CB", 21, "LEU_CG", 22, "LEU_CD1", 23, "LEU_CD2", 24,
	"ILE_N", 25, "ILE_CA", 26, "ILE_C", 27, "ILE_O", 28, "ILE_CB", 29, "ILE_CG1", 30, "ILE_CG2", 31, "ILE_CD1", 32,
	"SER_N", 33, "SER_CA", 34, "SER_C", 35, "SER_O", 36, "SER_CB", 37, "SER_OG", 38,
	"THR_N", 39, "THR_CA", 40, "THR_C", 41, "THR_O", 42, "THR_CB", 43, "THR_OG1", 44, "THR_CG2", 45,
	"CYS_N", 46, "CYS_CA", 47, "CYS_C", 48, "CYS_O", 49, "CYS_CB", 50, "CYS_SG", 51,
	"PRO_N", 52, "PRO_CA", 53, "PRO_C", 54, "PRO_O", 55, "PRO_CB", 56, "PRO_CG", 57, "PRO_CD", 58,
	"PHE_N", 59, "PHE_CA", 60, "PHE_C", 61, "PHE_O", 62, "PHE_CB", 63, "PHE_CG", 64, "PHE_CD1", 65, "PHE_CD2", 66,
	"PHE_CE1", 67, "PHE_CE2", 68, "PHE_CZ", 69,
	"TYR_N", 70, "TYR_CA", 71, "TYR_C", 72, "TYR_O", 73, "TYR_CB", 74, "TYR_CG", 75, "TYR_CD1", 76, "TYR_CD2", 77,
	"TYR_CE1", 78, "TYR_CE2", 79, "TYR_CZ", 80, "TYR_OH", 81,
	"TRP_N", 82, "TRP_CA", 83, "TRP_C", 84, "TRP_O", 85, "TRP_CB", 86, "TRP_CG", 87, "TRP_CD1", 88, "TRP_CD2", 89,
	"TRP_NE1", 90, "TRP_CE2", 91, "TRP_CE3", 92, "TRP_CZ2", 93, "TRP_CZ3", 94, "TRP_CH2", 95,
	"HIS_N", 96, "HIS_CA", 97, "HIS_C", 98, "HIS_O", 99, "HIS_CB", 100, "HIS_CG", 101, "HIS_ND1", 102, "HIS_CD2", 103,
	"HIS_CE1", 104, "HIS_NE2", 105,
	"ASP_N", 106, "ASP_CA", 107, "ASP_C", 108, "ASP_O", 109, "ASP_CB", 110, "ASP_CG", 111, "ASP_OD1", 112, "ASP_OD2", 113,
	"ASN_N", 114, "ASN_CA", 115, "ASN_C", 116, "ASN_O", 117, "ASN_CB", 118, "ASN_CG", 119, "ASN_OD1", 120, "ASN_ND2", 121,
	"GLU_N", 122, "GLU_CA", 123, "GLU_C", 124, "GLU_O", 125, "GLU_CB", 126, "GLU_CG", 127, "GLU_CD", 128, "GLU_OE1", 129,
	"GLU_OE2", 130,
	"GLN_N", 131, "GLN_CA", 132, "GLN_C", 133, "GLN_O", 134, "GLN_CB", 135, "GLN_CG", 136, "GLN_CD", 137, "GLN_OE1", 138,
	"GLN_NE2", 139,
	"MET_N", 140, "MET_CA", 141, "MET_C", 142, "MET_O", 143, "MET_CB", 144, "MET_CG", 145, "MET_SD", 146, "MET_CE", 147,
	"LYS_N", 148, "LYS_CA", 149, "LYS_C", 150, "LYS_O", 151, "LYS_CB", 152, "LYS_CG", 153, "LYS_CD", 154, "LYS_CE", 155,
	"LYS_NZ", 156,
	"ARG_N", 157, "ARG_CA", 158, "ARG_C", 159, "ARG_O", 160, "ARG_CB", 161, "ARG_CG", 162, "ARG_CD", 163, "ARG_NE", 164,
	"ARG_CZ", 165, "ARG_NH1", 166, "ARG_NH2", 167);
@list = keys(%index);
$nindex = @list;

# define the sequence in 3-letter and 1-letter formats
@amino3 = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "ASP", "ASN", "GLU", "GLN", "MET", "LYS", "ARG");
@amino1_upper = ("G", "A", "V", "L", "I", "S", "T", "C", "P", "F", "Y", "W", "H", "D", "N", "E", "Q", "M", "K", "R");
@amino1_lower = ("g", "a", "v", "l", "i", "s", "t", "c", "p", "f", "y", "w", "h", "d", "n", "e", "q", "m", "k", "r");
$namino = @amino3;

$Param = "~/bin/Parameters.inp";

if (@ARGV != 3) {
	printf STDERR "fit_scaffold.pl <arg1> <arg2> <arg3>\n";
	printf STDERR "<arg1>: loop epitope file\n";
	printf STDERR "<arg2>: scaffold pdb file\n";
	printf STDERR "<arg3>: scaffolding method\n";
	printf STDERR "1 - click\n";
	printf STDERR "2 - fast\n";
	printf STDERR "3 - mammoth\n";
	printf STDERR "4 - spalign\n";
	printf STDERR "5 - tmalignc\n";
	printf STDERR "6 - tmalignf\n";
	exit;
}

# check and read scaffold
if (not -f $ARGV[0]) {
	printf "error, epitope - $ARGV[0] doesn't exist\n";
	exit;
}

$strlen = length($ARGV[0]);
if (substr($ARGV[0], $strlen - 4, 4) ne ".pdb") {
	printf "error, epitope pdb file - $ARGV[0] doesn't end with .pdb\n";
	exit;
}
$ename = substr($ARGV[0], 0, $strlen - 4);

# read the pdb file
&readpdb($ARGV[0]);
$natm1 = $natm;
$nres1 = $nres;
@anam1 = @anam;
@rnam1 = @rnam;
@chid1 = @chid;
@rseq1 = @rseq;
@aseq1 = @aseq;
@xpdb1 = @xpdb;
@ypdb1 = @ypdb;
@zpdb1 = @zpdb;
@rca1 = @rca;

# check and then read loop epitope
if (not -f $ARGV[1]) {
	printf "error, scaffold protein - $ARGV[1] doesn't exist\n";
	exit;
}

# check scaffold name
$strlen = length($ARGV[1]);
if (substr($ARGV[1], $strlen - 4, 4) ne ".pdb") {
	printf "error, scaffold protein - $ARGV[1] doesn't end with .pdb\n";
	exit;
}
$sname = substr($ARGV[1], 0, $strlen - 4);

# read the pdb file
&readpdb($ARGV[1]);
$natm2 = $natm;
$nres2 = $nres;
@anam2 = @anam;
@rnam2 = @rnam;
@chid2 = @chid;
@rseq2 = @rseq;
@aseq2 = @aseq;
@xpdb2 = @xpdb;
@ypdb2 = @ypdb;
@zpdb2 = @zpdb;
@rca2 = @rca;

# fit the scaffold protein onto epitope
if ($ARGV[2] == 1) {
	system "cp $Param ./";
	system "click $ARGV[0] $ARGV[1] > test.out";

	# 1. get the specific tm-score
	open(OUT, "$ename-$sname.pdb.1.clique");
	@outtxt = <OUT>;
	$outlen = @outtxt;
	close(OUT);
	unlink "test.out", "Parameters.inp";
	unlink "$ename-$sname.1.pdb", "$sname-$ename.1.pdb";
	unlink "$ename-$sname.pdb.1.clique";

	# calculate the pair-wise correspondence
	for ($nmatch = 0, $j = 0; $j < $outlen; $j++) {
		if ($outtxt[$j] =~ /Chain/) {
			for ($k = $j + 1; $k < $outlen; $k++) {
				chomp $outtxt[$k];
				$pair[0][$nmatch] = substr($outtxt[$k], 6, 10) + 0;
				$pair[1][$nmatch] = substr($outtxt[$k], 54, 10) + 0;
				$nmatch++;
			}
		}
	}

	# check connectivity
	for ($nonlinear = 0, $j = 0; $j < $nmatch - 1; $j++) {
		if ($pair[1][$j + 1] < $pair[1][$j]) {
			$nonlinear = 1;
			last;
		}
	}

	# print out warning message
	if ($nonlinear == 1) {
		printf "warning: nonlinear region matched by click\n";
	}
}
elsif ($ARGV[2] == 2) {
	system "fast $ARGV[0] $ARGV[1] > test.out";

	# 1. get the alignment score - rmsd
	open(OUT, "test.out");
	@outtxt = <OUT>;
	$outlen = @outtxt;
	close(OUT);
	unlink "test.out";

	# 2. manage gaps in the alignment
	$alg_seq1 = $alg_corr = $alg_seq2 = "";
	for ($alg_seq2 = "", $j = 0; $j < $outlen; $j++) {
		if ($outtxt[$j] =~ /^ 1:/) {
			chomp $outtxt[$j];
			@arr1 = split /\t+/, $outtxt[$j];
			$narr1 = @arr1;
			$alg_seq1. = $arr1[$narr1 - 1];
		}
		if ($outtxt[$j] =~ /^ 2:/) {
			chomp $outtxt[$j];
			@arr2 = split /\t+/, $outtxt[$j];
			$narr2 = @arr2;
			$alg_seq2. = $arr2[$narr2 - 1];
		}
	}
	$len1 = length($alg_seq1);
	$len2 = length($alg_seq2);

	# check if alignment strings match
	if ($len1 != $len2) {
		printf STDERR "fast alignment error for $pdblist[$i]\n";
		exit;
	}
	$alg_size = length($alg_seq1);

	# search for the first match
	for ($js = 0; $js < $alg_size; $js++) {
		last if (substr($alg_seq1, $js, 1) ne "-" and substr($alg_seq2, $js, 1) ne "-");
	}

	# search for the last match (skip the *)
	for ($je = $alg_size - 2; $je >= 0; $je--) {
		last if (substr($alg_seq1, $je, 1) ne "-" and substr($alg_seq2, $je, 1) ne "-");
	}
	$strlen = $je - $js + 1;

	# find the correspondence section
	for ($iline = 0, $j = $outlen - 1; $j >= 0; $j--) {
		if ($outtxt[$j] =~ /^ 2:/) {
			$iline = $j + 2;
			last;
		}
	}

	# calculate the pair-wise correspondence
	for ($nmatch = 0, $j = $iline; $j < $outlen; $j++) {
		chomp $outtxt[$j];
		next if (length($outtxt[$j]) == 0);
		@arr1 = split / +/, $outtxt[$j];
		$narr1 = @arr1;
		$pair[0][$nmatch] = $arr1[$narr1 - 6];
		$pair[1][$nmatch] = $arr1[$narr1 - 1];
		$nmatch++;
	}
}
elsif ($ARGV[2] == 3) {
	system "mammoth -e $ARGV[0] -p $ARGV[1] > test.out";

	# 1. get the specific tm-score
	open(OUT, "test.out");
	@outtxt = <OUT>;
	$outlen = @outtxt;
	close(OUT);
	unlink "test.out";

	# strucutre alignment-derived sequence alignment
	$alg_seq1 = $alg_corr = $alg_seq2 = "";
	for ($alg_seq2 = "", $j = 0; $j < $outlen; $j++) {
		if ($outtxt[$j] =~ /^Prediction/ and $outtxt[$j + 1] =~ /^Prediction/) {
			chomp $outtxt[$j];
			$alg_seq2 .= substr($outtxt[$j], 11, length($outtxt[$j]) - 11);
		}
		if ($outtxt[$j] =~ /^Experiment/ and $outtxt[$j + 1] =~ /^Experiment/) {
			chomp $outtxt[$j + 1];
			$alg_seq1 .= substr($outtxt[$j+1], 11, length($outtxt[$j + 1]) - 11);
		}
		if ($outtxt[$j] =~ /^Prediction/ and $outtxt[$j + 1] =~ /^Prediction/) {
			chomp $outtxt[$j - 1];
			$alg_corr .= substr($outtxt[$j - 1], 11, length($outtxt[$j - 1]) - 11);
		}
	}

	# skip failed alignment
	next if ($alg_seq1 eq "" and $alg_seq2 eq "" and $alg_corr eq "");

	# calculate the pair-wise correspondence
	$strlen = length($alg_corr);
	$idx_seq1 = $idx_seq2 = 0;
	for ($nmatch = 0, $j = 0; $j < $strlen; $j++) {
		if (substr($alg_corr, $j, 1) eq "*") {
			$pair[0][$nmatch] = $rca1[$idx_seq1];
			$pair[1][$nmatch] = $rca2[$idx_seq2];
			$nmatch++;
		}
		if (substr($alg_seq1, $j, 1) ne "-") {
			$idx_seq1++;
		}
		if (substr($alg_seq2, $j, 1) ne "-") {
			$idx_seq2++;
		}
	}
}
elsif ($ARGV[2] == 4) {
	system "SP-align.gnu $ARGV[0] $ARGV[1] > test.out";

	# 1. get the specific tm-score
	open(OUT, "test.out");
	@outtxt = <OUT>;
	close(OUT);
	unlink "test.out";

	# manage gaps in the alignment
	$alg_seq1 = $outtxt[12];
	chomp $alg_seq1;
	$alg_corr = $outtxt[13];
	chomp $alg_corr;
	$alg_seq2 = $outtxt[14];
	chomp $alg_seq2;
	$alg_size = length($alg_corr);

	# search for the first : or .
	for ($js = 0; $js < $alg_size; $js++) {
		last if (substr($alg_corr, $js, 1) eq ":" or substr($alg_corr, $js, 1) eq ".");
	}

	# search for the last : or .
	for ($je = $alg_size - 1; $je >= 0; $je--) {
		last if (substr($alg_corr, $je, 1) eq ":" or substr($alg_corr, $je, 1) eq ".");
	}
	$strlen = $je - $js + 1;

	# starting & ending residues of seq 1
	$seq1_str = substr($alg_seq1, 0, $js);
	$seq1_end = substr($alg_seq1, $je + 1, length($alg_seq2) - $je - 2);

	# starting & ending residues of seq 2
	$seq2_str = substr($alg_seq2, 0, $js);
	$seq2_end = substr($alg_seq2, $je + 1, length($alg_seq2) - $je - 2);

	# reassign the aligned region
	$alg_seq1 = substr($alg_seq1, $js, $strlen);
	$alg_corr = substr($alg_corr, $js, $strlen);
	$alg_seq2 = substr($alg_seq2, $js, $strlen);

	# 2. calculate the pair-wise correspondence
	for ($nmatch = 0, $j = 0; $j < $strlen; $j++) {
		if (substr($alg_corr, $j, 1) eq ":" or substr($alg_corr, $j, 1) eq ".") {
			$pair[0][$nmatch] = $rca1[$seq1_str + $j - $ngap1];
			$pair[1][$nmatch] = $rca2[$seq2_str + $j - $ngap2];
			$nmatch++;
		}
		elsif (substr($alg_seq1, $j, 1) eq "-") {
			$ngap1++;
		}
		elsif (substr($alg_seq2, $j, 1) eq "-") {
			$ngap2++;
		}
	}
}
elsif ($ARGV[2] == 5) {
	system "TMalignC -A $ARGV[0] -B $ARGV[1] -a T -o TM.sup >& /dev/null";
	open(SUP, "TM.sup");
	@suptxt = <SUP>;
	$suplen = @suptxt;
	close(SUP);
	unlink "TM.sup", "TM.sup_all";
	$nmatch = 0;
	for ($idx = 1, $j = 0; $j < $suplen; $j++) {
		chomp $suptxt[$j];
		if ($suptxt[$j] = ~/^ATOM/) {
			for ($nca = 0, $k = $j; $k < $suplen; $k++) {
				if ($suptxt[$k] = ~/^TER/) {
					$nmatch = $nca;
					$nca = 0;
					$j = $k;
					$idx--;
					last;
				}
				$pair[$idx][$nca] = substr($suptxt[$k], 22, 4) + 0;
				$nca++;
			}
		}
	}
}
elsif ($ARGV[2] == 6) {
	system "TMalignF $ARGV[0] $ARGV[1] -o TM.sup >& /dev/null";
	open(SUP, "TM.sup");
	@suptxt = <SUP>;
	$suplen = @suptxt;
	close(SUP);
	unlink "TM.sup", "TM.sup_all";
	$nmatch = 0;
	for ($idx = 0, $i = 0; $i < $suplen; $i++) {
		chomp $suptxt[$i];
		if ($suptxt[$i] = ~/^ATOM/) {
			for ($nca = 0, $j = $i; $j < $suplen; $j++) {
				if ($suptxt[$j] = ~/^TER/) {
					$nmatch = $nca;
					$nca = 0;
					$i = $j;
					$idx++;
					last;
				};
				$pair[$idx][$nca] = substr($suptxt[$j], 22, 4) + 0;
				$nca++;
			}
		}
	}
}
open(FIT, ">res_fit.list");
for ($i = 0; $i < $nmatch; $i++) {
	# find the epitope residue
	for ($j = 0; $j < $natm1; $j++) {
		if ($rseq1[$j]  ==  $pair[0][$i] and $anam1[$j] eq "CA") {
			$aseq_epi = $aseq1[$j];
			$chid_epi = $chid1[$j];
		}
	}

	# find the corresponding epitope residue
	for ($j = 0; $j < $natm2; $j++) {
		if ($rseq2[$j]  ==  $pair[1][$i] and $anam2[$j] eq "CA") {
			$aseq_scf = $aseq2[$j];
			$chid_scf = $chid2[$j];
		}
	}
	printf FIT ("%1s, %5d, %1s, %5d\n", $chid_epi, $aseq_epi, $chid_scf, $aseq_scf);
}
close(FIT);
system "rmsd_sim $ARGV[0] $ARGV[1] -f res_fit.list > scaffold_fit.pdb";

# read in protein pdb file and store information in global arrays
sub readpdb() {
	my $ipdb, $jpdb, $pdblen, @pdbtxt;
	my $strtmp, @dual, $pdbnam, $find;
	open(PROPDB, "$_[0]");
	@pdbtxt = <PROPDB>;
	$pdblen = @pdbtxt;
	close(PROPDB);
	for ($natm = 0, $ipdb = 0; $ipdb < $pdblen; $ipdb++) {
		if ($pdbtxt[$ipdb] =~ /^ATOM/) {
			chomp $pdbtxt[$ipdb];
			$strtmp = substr($pdbtxt[$ipdb], 11, 10);
			@dual = split(/ +/, $strtmp);

			# correct ile_cd in gromacs-generated pdbs
			if ($dual[2] eq "ILE" and $dual[1] eq "CD") {
				$dual[1] = "CD1";
			}

			# correct c-terminal oxygen in gromacs pdbs
			if ($dual[1] eq "O1") {
				$dual[1] = "O";
			}

			# create a unique pdb name for protein atom
			$pdbnam = $dual[2]."_".$dual[1];
			$find = 0;
			for ($jpdb = 0; $jpdb < $nindex; $jpdb++) {
				if ($pdbnam eq $list[$jpdb]) {
					$find = 1;
					last;
				}
			}
			next if ($find == 0);
			$anam[$natm] = $dual[1];
			$rnam[$natm] = $dual[2];
			$aseq[$natm] = substr($pdbtxt[$ipdb], 6, 5);
			$chid[$natm] = substr($pdbtxt[$ipdb], 21, 1);
			$rseq[$natm] = substr($pdbtxt[$ipdb], 22, 4);
			$xpdb[$natm] = substr($pdbtxt[$ipdb], 30, 8);
			$ypdb[$natm] = substr($pdbtxt[$ipdb], 38, 8);
			$zpdb[$natm] = substr($pdbtxt[$ipdb], 46, 8);
			$natm++;
		}
	}
	# # residues
	$nres = 0;
	$rseq_prev = -1000;
	for ($ipdb = 0; $ipdb < $natm; $ipdb++) {
		if ($rseq[$ipdb] != $rseq_prev) {
			$rseq_prev = $rseq[$ipdb];
			$ratm_star[$nres] = $ipdb;
			$nres++;
		}
	}
	$ratm_star[$nres] = $natm;
	# c-alpha trace
	for ($ipdb = 0; $ipdb < $nres; $ipdb++) {
		$idx_str = $ratm_star[$ipdb];
		$idx_end = $ratm_star[$ipdb + 1];
		for ($jpdb = $idx_str; $jpdb < $idx_end; $jpdb++) {
			if ($anam[$jpdb] eq 'CA') {
				$xca[$ipdb] = $xpdb[$jpdb];
				$yca[$ipdb] = $ypdb[$jpdb];
				$zca[$ipdb] = $zpdb[$jpdb];
				$rca[$ipdb] = $rseq[$jpdb];
			}
		}
	}
}

sub writepdb() {
	my $ipdb;
	open(PROPDB, ">$_[0]");
	for ($ipdb = 0; $ipdb < $natm; $ipdb++) {
		printf PROPDB ("ATOM%7d  %-4s%-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", $ipdb + 1, $anam[$ipdb], $rnam[$ipdb], $chid[$ipdb], $rseq[$ipdb], $xpdb[$ipdb], $ypdb[$ipdb], $zpdb[$ipdb], 1.00, 0.00);
	}
	close(PROPDB);
}
