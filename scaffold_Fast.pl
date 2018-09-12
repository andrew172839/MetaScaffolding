#!/usr/bin/perl

# define reasonable pdb atom name
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

# define sequence in 3-letter and 1-letter formats
@amino3 = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "ASP", "ASN", "GLU", "GLN", "MET", "LYS", "ARG");
@amino1_upper = ("G", "A", "V", "L", "I", "S", "T", "C", "P", "F", "Y", "W", "H", "D", "N", "E", "Q", "M", "K", "R");
@amino1_lower = ("g", "a", "v", "l", "i", "s", "t", "c", "p", "f", "y", "w", "h", "d", "n", "e", "q", "m", "k", "r");
$namino = @amino3;

# surface calculation binary
$surf_bin = "~/programs/jackal/bin/surface";

# surface area of amino acid in the ala-x-ala environment
@area_std = (87.161, 113.564, 156.829, 179.354, 183.745, 127.427, 147.772, 142.608, 147.006, 212.730, 227.144, 247.560, 190.253, 153.640, 154.828, 193.516, 195.611, 209.203, 225.543, 255.639);
$Fast_bin = "~/bin/fast";

# find scaffolds that contain the beta-sheet stem regions
if (@ARGV != 3) {
	print STDERR"scaffold_fast.pl <arg1> <arg2> <arg3>\n";
	print STDERR"<arg1>: segment epitope pdb file, e.g. epitope.pdb\n";
	print STDERR"<arg2>: database dir, e.g. xxx/database/cullpdb-1\n";
	print STDERR"<arg3>: ca-rmsd cutoff of the aligned region, e.g., 2.0\n";
	exit;
}

# check epitope pdb
if (not -f $ARGV[0]) {
	printf "error, scaffold_fast.pl - $ARGV[0] doesn't exist\n";
	exit;
}

# solvent-accessibility for isolated epitope
system "$surf_bin $ARGV[0] >& surf.dat";
open(SRF, "surf.dat");
@srftxt = <SRF>;
$srflen = @srftxt;
close(SRF);
for ($sfscore_ref = 0.0, $j = 0; $j < $srflen; $j++) {
	chomp $srftxt[$j];
	$strlen = length($srftxt[$j]);
	next if ($strlen == 0);
	next if ($srftxt[$j] !~ /!RES          AREA/);
	for ($k = $j + 1; $k < $srflen; $k++) {
		chomp $srftxt[$k];
		$strlen = length($srftxt[$k]);
		last if ($strlen == 0);
		last if ($srftxt[$k] =~ /!Chain/);
		@arr = split / +/, $srftxt[$k];
		$narr = @arr;
		if ($narr == 4) {
			$sfrnam = $arr[$narr - 4];
			$sfrseq = $arr[$narr - 3] + 0;
		}
		else {
			$sfrnam = $arr[$narr - 3];
			$sfrseq = $arr[$narr - 2] + 0;
		}

		# surface area
		$sfarea = $arr[$narr - 1];
		# accessibility
		for ($l = 0; $l < $namino; $l++) {
			if ($sfrnam eq $amino3[$l]) {
				$sfperc = $sfarea / $area_std[$l];
				last;
			}
		}
		$sfscore_ref += $sfperc;
	}
}

# get database dir
opendir(DIR, $ARGV[1]);
@dirtxt = readdir(DIR);
$dirnum = @dirtxt;
closedir(DIR);

# create pdb file list
for ($npdb = 0, $i = 0; $i < $dirnum; $i++) {
	if ($dirtxt[$i] !~ /^\./) {
		next if ($dirtxt[$i] !~ /pdb/);
		$pdblist[$npdb] = $ARGV[1]."/$dirtxt[$i]";
		$pdbname[$npdb] = $dirtxt[$i];
		$npdb++;
	}
}

# align beta-sheet stems onto each protein
for ($i = 0; $i < $npdb; $i++) {
	system "$Fast_bin $ARGV[0] $pdblist[$i] > test.out";

	# 1. get the alignment score - rmsd
	open(OUT, "test.out");
	@outtxt = <OUT>;
	$outlen = @outtxt;
	close(OUT);
	unlink "test.out";

	# get nalign, size 1, size 2 and rmsd
	chomp $outtxt[1];
	@arr1 = split / +/, $outtxt[1];
	@arr2 = split /=/, @arr1[0];
	$nalign = $arr2[1] + 0;
	@arr2 = split /=/, @arr1[3];
	$size1 = $arr2[1] + 0;
	@arr2 = split /=/, @arr1[4];
	$size2 = $arr2[1] + 0;
	@arr2 = split /=/, @arr1[5];
	$rmsd = $arr2[1] + 0.0;

	# skip proteins if the aligned regions is 2-aa shorter
	next if ($nalign < $size1 - 2);
	$alg_seq1 = $alg_corr = $alg_seq2 = "";
	for ($alg_seq2 = "", $j = 0; $j < $outlen; $j++) {
		if ($outtxt[$j] =~ /^ 1: /) {
			chomp $outtxt[$j];
			@arr1 = split /\t + /, $outtxt[$j];
			$narr1 = @arr1;
			$alg_seq1. = $arr1[$narr1 - 1];
		}
		if ($outtxt[$j] =~ /^ 2: /) {
			chomp $outtxt[$j];
			@arr2 = split /\t + /, $outtxt[$j];
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

	# search for first match
	for ($js = 0; $js < $alg_size; $js++) {
		last if (substr($alg_seq1, $js, 1) ne "-" and substr($alg_seq2, $js, 1) ne "-");
	}

	# search for last match (skip *)
	for ($je = $alg_size - 2; $je >= 0; $je--) {
		last if (substr($alg_seq1, $je, 1) ne "-" and substr($alg_seq2, $je, 1) ne "-");
	}
	$strlen = $je-$js + 1;

	# 3 skip proteins if aligned region has high rmsd
	next if  ($rmsd > $ARGV[2]);

	# 2. calculate solvent accessibility for matched region
	# find the correspondence section
	for ($iline = 0, $j = $outlen - 1; $j >= 0; $j--) {
		if ($outtxt[$j] =~ /^ 2:/) {
			$iline = $j + 2;
			last;
		}
	}

	# pair-wise correspondence
	for ($nmatch = 0, $j = $iline; $j < $outlen; $j++) {
		chomp $outtxt[$j];

		# skip empty lines
		next if (length($outtxt[$j]) == 0);

		# get residue correspondence
		@arr1 = split / +/, $outtxt[$j];
		$narr1 = @arr1;
		$pair[0][$nmatch] = $arr1[$narr1 - 6];
		$pair[1][$nmatch] = $arr1[$narr1 - 1];
		$nmatch++;
	}
	# solvent-accessibility
	system "$surf_bin $pdblist[$i] >& surf.dat";
	open(SRF, "surf.dat");
	@srftxt = <SRF>;
	$srflen = @srftxt;
	close(SRF);
	unlink "surf.dat";
	for ($sfscore = 0.0, $j = 0; $j < $srflen; $j++) {
		chomp $srftxt[$j];
		$strlen = length($srftxt[$j]);
		next if ($strlen == 0);
		next if ($srftxt[$j] !~ /!RES          AREA/);
		for ($k = $j + 1; $k < $srflen; $k++) {
			chomp $srftxt[$k];
			$strlen = length($srftxt[$k]);
			last if ($strlen == 0);
			last if ($srftxt[$k] =~ /!Chain/);
			@arr = split / +/, $srftxt[$k];
			$narr = @arr;
			if ($narr == 4) {
				$sfrnam = $arr[$narr - 4];
				$sfrseq = $arr[$narr - 3] + 0;
			}
			else {
				$sfrnam = $arr[$narr - 3];
				$sfrseq = $arr[$narr - 2] + 0;
			}

			# surface area
			$sfarea = $arr[$narr - 1];
			# accessibility
			for ($l = 0; $l < $namino; $l++) {
				if ($sfrnam eq $amino3[$l]) {
					$sfperc = $sfarea / $area_std[$l];
					last;
				}
			}

			# check residue
			for ($l = 0; $l < $nmatch; $l++) {
				if ($sfrseq == $pair[1][$l]) {
					$sfscore += $sfperc;
					last;
				}
			}
		}
	}
	printf ("%5d  %-s   %-5d%-5d%-12.5f%-12.5f\n", $i + 1, $pdbname[$i], $size2, $nalign, $rmsd, $sfscore / $sfscore_ref);
}
