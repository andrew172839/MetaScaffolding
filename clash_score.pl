#!/usr/bin/perl

#define the reasonable pdb atom name
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
#set random number seed
srand;
#check the command-line arguments
if (@ARGV != 2) {
	printf "clash_score.pl <arg1> <arg2>\n";
	printf "<arg1>: fitted scaffold pdb directory\n";
	printf "<arg2>: antibdoy structure pdb file\n";
	exit;
}
#check if the pdb structure directory exits
if (not -d $ARGV[0]) {
	printf "clash_score error: cannot find fitted scaffold pdb directory $ARGV[0]\n";
	exit;
}
opendir(DIR, $ARGV[0]);
@dirtxt = readdir(DIR);
$ndirtxt = @dirtxt;
closedir(DIR);
for ($nscff = 0, $i = 0; $i < $ndirtxt; $i++) {
	if ($dirtxt[$i] !~ /^\./) {
		next if ($dirtxt[$i] !~ /pdb/);
		$scfflist[$nscff] = $ARGV[0]."/$dirtxt[$i]";
		$scffname[$nscff] = $dirtxt[$i];
		$nscff++;
	}
}
#check if the antibody pdb exists
if (not -f $ARGV[1]) {
	printf "clash_score error: can not find antibody pdb file $ARGV[1]\n";
	exit;
}
#read in pdb
&readpdb($ARGV[1]);
$natm_nab = $natm;
@aseq_nab = @aseq;
@anam_nab = @anam;
@rnam_nab = @rnam;
@chid_nab = @chid;
@rseq_nab = @rseq;
@xpdb_nab = @xpdb;
@ypdb_nab = @ypdb;
@zpdb_nab = @zpdb;
#calculate clash score
for ($istr = 0; $istr < $nscff; $istr++) {
	#read in a scaffold
	&readpdb($scfflist[$istr]);
	$natm_scf = $natm;
	@aseq_scf = @aseq;
	@anam_scf = @anam;
	@rnam_scf = @rnam;
	@chid_scf = @chid;
	@rseq_scf = @rseq;
	@xpdb_scf = @xpdb;
	@ypdb_scf = @ypdb;
	@zpdb_scf = @zpdb;
	#check for atomic clashes
	for ($ene = 0.0, $iene = 0; $iene < $natm_nab; $iene++) {
		for ($jene = 0; $jene < $natm_scf; $jene++) {
			$xdif  = $xpdb_scf[$iene] - $xpdb_nab[$jene];
			$ydif  = $ypdb_scf[$iene] - $ypdb_nab[$jene];
			$zdif  = $zpdb_scf[$iene] - $zpdb_nab[$jene];
			$dist2 = $xdif * $xdif + $ydif * $ydif + $zdif * $zdif ;
			$dist = sqrt($dist2);
			if ($dist < 1.0) {
				$ene += 1.0 / $dist;
			}
		}
	}
	#print out the clash score
	printf "$scffname[$istr]: $ene\n";
}

#read in protein pdb file and store inforamtion in global arrays
sub readpdb() {
	my $ipdb, $jpdb, $pdblen, @pdbtxt;
	my $strtmp, @dual, $pdbnam, $find;
	open(PROPDB, "$_[0]");
	@pdbtxt = <PROPDB>;
	$pdblen = @pdbtxt;
	close(PROPDB);
	for ($natm = 0, $ipdb = 0; $ipdb < $pdblen; $ipdb++) {
		if ($pdbtxt[$ipdb] = ~/^ATOM/ or $pdbtxt[$ipdb] = ~/^HETATM/) {
			$strtmp = substr($pdbtxt[$ipdb], 11, 10);
			@dual = split(/ +/, $strtmp);
			$pdbnam = $dual[2]."_".$dual[1];
			$aseq[$natm] = substr($pdbtxt[$ipdb], 6, 5);
			$anam[$natm] = $dual[1];
			$rnam[$natm] = $dual[2];
			$chid[$natm] = substr($pdbtxt[$ipdb], 21, 1);
			$rseq[$natm] = substr($pdbtxt[$ipdb], 22, 4);
			$xpdb[$natm] = substr($pdbtxt[$ipdb], 30, 8);
			$ypdb[$natm] = substr($pdbtxt[$ipdb], 38, 8);
			$zpdb[$natm] = substr($pdbtxt[$ipdb], 46, 8);
			$natm++;
		}
	}
	#number of residues
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
}