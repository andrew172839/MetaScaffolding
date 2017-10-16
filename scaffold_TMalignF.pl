#!/usr/bin/perl

#define the reasonable pdf atom name
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
#define the sequence in 3-letter and 1-letter formats
@amino3 = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "CYS", "PRO", "PHE", "TYR", "TRP", "HIS", "ASP", "ASN", "GLU", "GLN", "MET", "LYS", "ARG");
@amino1_upper = ("G", "A", "V", "L", "I", "S", "T", "C", "P", "F", "Y", "W", "H", "D", "N", "E", "Q", "M", "K", "R");
@amino1_lower = ("g", "a", "v", "l", "i", "s", "t", "c", "p", "f", "y", "w", "h", "d", "n", "e", "q", "m", "k", "r");
$namino = @amino3;
#surface calculation binary
$surf_bin = "~/programs/jackal/bin/surface";
#surface area of amino acid in the ALA-X-ALA environment
@area_std = (87.161, 113.564, 156.829, 179.354, 183.745, 127.427, 147.772, 142.608, 147.006, 212.730,
	227.144, 247.560, 190.253, 153.640, 154.828, 193.516, 195.611, 209.203, 225.543, 255.639);
$TMalign = "~/bin/TMalignF";
#find scaffolds that contain the beta-sheet stem regions
if (@ARGV != 3) {
	print STDERR "scaffold_tmalignf.pl <arg1> <arg2> <arg3>\n";
	print STDERR "<arg1>: segment epitope pdb file, e.g. epitope.pdb\n";
	print STDERR "<arg2>: database dir, e.g. xxx/database/cullpdb-1\n";
	print STDERR "<arg3>: ca-rmsd cutoff of the aligned region, e.g., 2.0\n";
	exit;
}
#check the epitope pdb
if  (not -f $ARGV[0]) {
	printf "error: scaffold_tmalignf.pl - $ARGV[0] doesn't exist\n";
	exit;
}
#calculate solvent-accessibility for isolated epitope
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
			#residue name
			$sfrnam = $arr[$narr - 4];
			#residue number
			$sfrseq = $arr[$narr - 3] + 0;
		}
		else {
			#residue name
			$sfrnam = $arr[$narr - 3];
			#residue number
			$sfrseq = $arr[$narr - 2] + 0;
		}
		#surface area
		$sfarea = $arr[$narr - 1];
		#calculate accessibility
		for ($l = 0; $l < $namino; $l++) {
			if ($sfrnam eq $amino3[$l]) {
				$sfperc = $sfarea / $area_std[$l];
				last;
			}
		}
		$sfscore_ref += $sfperc;
	}
}
#get database dir
opendir(DIR, $ARGV[1]);
@dirtxt = readdir(DIR);
$dirnum = @dirtxt;
closedir(DIR);
#create pdb file list
for ($npdb = 0, $i = 0; $i < $dirnum; $i++) {
	if ($dirtxt[$i] !~ /^\./) {
		next if ($dirtxt[$i] !~ /pdb/);
		$pdblist[$npdb] = $ARGV[1]."/$dirtxt[$i]";
		$pdbname[$npdb] = $dirtxt[$i];
		$npdb++;
	}
}
#align beta-sheet stems onto each of the proteins
for ($i = 0; $i < $npdb; $i++) {
	system "$TMalign $ARGV[0] $pdblist[$i] > test.out";
	#1. get the specific tm-score
	open(OUT, "test.out");
	@outtxt = <OUT>;
	close(OUT);
	unlink "test.out";
	#get size 1 of template
	chomp $outtxt[8];
	$size1 = substr($outtxt[8], 25, 4) + 0;
	#get size 2 of query
	chomp $outtxt[9];
	$size2 = substr($outtxt[9], 25, 4) + 0;
	#get nalign
	chomp $outtxt[11];
	@arr1 = split /,/, $outtxt[11];
	@arr2 = split /=/, $arr1[0];
	$narr2 = @arr2;
	$nalign = $arr2[$narr2 - 1] + 0;
	#skip proteins if the aligned regions is shorter
	next if ($nalign < $size1 - 2);
	#get tm-score
	@arr1 = split /,/, $outtxt[11];
	@arr2 = split /=/, $arr1[2];
	$narr2 = @arr2;
	$score = $arr2[$narr2 - 1] + 0.0;
	#calculate the average size
	$size_ave = int(($size1 + $size2) / 2);
	#2. manage gaps in the alignment
	$alg_seq1 = $outtxt[14];
	chomp $alg_seq1;
	$alg_corr = $outtxt[15];
	chomp $alg_corr;
	$alg_seq2 = $outtxt[16];
	chomp $alg_seq2;
	$alg_size = length($alg_corr);
	@arr = split / +/, $alg_corr;
	$narr = @arr;
	#search for the first : or .
	for ($js = 0; $js < $alg_size; $js++) {
		last if (substr($alg_corr, $js, 1) eq ":" or substr($alg_corr, $js, 1) eq ".");
	}
	#search for  the last :  or .
	for ($je = $alg_size - 1; $je >= 0; $je--) {
		last if (substr($alg_corr, $je, 1) eq ":" or substr($alg_corr, $je, 1) eq ".");
	}
	$strlen = $je - $js + 1;
	#realign the structures
	system "$TMalign $ARGV[0] $pdblist[$i] -L $size_ave -o TM.sup > test.out";
	open(OUT, "test.out");
	@outtxt = <OUT>;
	close(OUT);
	unlink "test.out";
	#get RMSD
	chomp $outtxt[11];
	@arr1 = split /,/, $outtxt[11];
	@arr2 = split /=/, $arr1[1];
	$narr2 = @arr2;
	$rmsd = $arr2[$narr2 - 1] + 0.0;
	#get tm-score
	@arr1 = split /,/, $outtxt[11];
	@arr2 = split /=/, $arr1[2];
	$narr2 = @arr2;
	$score = $arr2[$narr2 - 1] + 0.0;
	#3.1 skip proteins if the aligned region has high rmsd
	next if ($rmsd > $ARGV[2]);
	#2. calculate solvent accessibility for matched region
	#read in the TM.sup file
	open(SUP, "TM.sup");
	@suptxt = <SUP>;
	$suplen = @suptxt;
	close(SUP);
	unlink "TM.sup";
	$nmatch = 0;
	for ($idx = 0, $j = 0; $j < $suplen; $j++) {
		chomp $suptxt[$j];
		if ($suptxt[$j] =~ /^ATOM/) {
			for ($nca = 0, $k = $j; $k < $suplen; $k++) {
				if  ($suptxt[$k] =~ /^TER/) {
					$nmatch = $nca;
					$nca = 0;
					$j = $k;
					$idx++;
					last;
				}
				$pair[$idx][$nca] = substr($suptxt[$k], 22, 4) + 0;
				$nca++;
			}
		}
	}
	#calculate solvent-accessibility
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
				#residue name
				$sfrnam = $arr[$narr - 4];
				#residue number
				$sfrseq = $arr[$narr - 3] + 0;
			}
			else  {
				#residue name
				$sfrnam = $arr[$narr - 3];
				#residue number
				$sfrseq = $arr[$narr - 2] + 0;
			}
			#surface area
			$sfarea = $arr[$narr - 1];
			#calculate accessibility
			for ($l = 0; $l < $namino; $l++) {
				if ($sfrnam eq $amino3[$l]) {
					$sfperc = $sfarea / $area_std[$l];
					last;
				}
			}
			#check if it is the residue we want
			for ($l = 0; $l < $nmatch; $l++) {
				if ($sfrseq == $pair[1][$l]) {
					$sfscore += $sfperc;
					last;
				}
			}
		}
	}
	#4. print out the output
	printf ("%5d  %-s   %-5d%-5d%-12.5f%-8.3f%-12.5f\n", $i + 1, $pdbname[$i], $size2, $nalign, $rmsd, $score, $sfscore / $sfscore_ref);
}
