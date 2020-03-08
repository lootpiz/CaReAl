#!/usr/bin/env python
import os
import sys
import gzip
import time
import string
from distutils.spawn import find_executable
from distutils.version import LooseVersion

#################### Check Python version and modules ####################
if sys.version_info < (2,7,0):
	sys.stderr.write("You need python 2.7 or later to run this script.\n")
	exit(1)

try:
	import tabix
except ImportError:
	print "Please install pytabix module [ pip install pytabix ]"
try:
	import subprocess
	from subprocess import call
	from subprocess import Popen
except ImportError:
	print "Please install subprocess module [ pip install subprocess ]"
try:
	from optparse import OptionParser
except ImportError:
	print "Please install optparse module [ pip install optparse ]"
try:
	import multiprocessing
except ImportError:
	print "Please install multprocessing module [ pip install multiprocessing ]"

samtools_env = find_executable('samtools')
if samtools_env is None:
        sys.stderr.write("SAMTools required. Please visit at http://www.htslib.org/.\n")
        exit(1)
else:
        get_version_cmd = "samtools 2>&1 > /dev/null | grep Version | sed -Ee \'s/Version: ([0-9.]+).*/\\1/;q\'"
        get_version = subprocess.Popen(get_version_cmd,stdout=subprocess.PIPE,stderr=None,shell=True)
        samtools_version = get_version.stdout.read().strip()
	if LooseVersion(samtools_version) < LooseVersion("0.1.19"):
        	sys.stderr.write("SAMTools (version >= 0.1.19) required. Your SAMTools version is %s.\n", samtools_version)
		exit(1)

#################### Arguments #################### 
filteroptions   = []
parser = OptionParser()
parser.add_option("-b", "--bamfile", dest="infile", help="input BAM file or BAM files with directory paths in txt format", default=False)
parser.add_option("-t", "--target", dest="targetfile", help="chromosome and position or chromosome and range, or targets in txt format", default=False)
parser.add_option("-r", "--reference", dest="referencefile", help="reference genome file", default=False)
parser.add_option("-g","--gene", dest="geneinformation", help="gene information file in BED", default=False)
parser.add_option("-v","--with-vcf", dest="vcffile_flag", help="query VCF file and the name of VCF file must be the same as given BAM file (True/False)", default=False)
parser.add_option("-o","--folder", dest="outputfolder", help="out directory name, default is \"OUTPUT\"", default="OUTPUT")
parser.add_option("-n", "--cpu", dest="nr_cpus", help="Number of CPUs to use", default=4)
(options, args) = parser.parse_args()

#################### Codon tables ####################
bases = ["T","C","A","G"]
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons,amino_acids))

#################### R code #################### 
middlecode1="""
if (Collength >= 80) {
	png(file=png_file_name,bg="white",width=(Collength * 18.75), height= 700 + (Rowlength * 24) + vcf_additional, res=200)
} else if (Collength < 80) {
	png(file=png_file_name,bg="white",width=700, height=700 + (Rowlength * 24) + vcf_additional, res=200)
}

layout(matrix(c(1,2,3)), heights = c(lcm(2),2,lcm(ceiling(length(CALL)/10)*3)))
par(mar=c(0,0.2,1.5,0))
barplot(Coverage, xlim=c(1,Collength), ylim=c(0,1), width=1, space=0, axes=F, main=barplot_title)

par(mar=c(0,0,0,0.2))
plot(x=NULL,y=NULL,xlim=c(1, Collength), ylim=c(1,Rowlength2), ann=F, bty='n', type='n', xaxt='n', yaxt='n')

if (length(RasterX) != 0){
	text(RasterX,RasterY,rep("I",length(RasterX)),col="darkorchid4",cex=1)}
if (length(indelX) != 0){
	text(indelX,indelY,rep("D",length(indelX)),col="darkorchid4",cex=1)}
if (length(RperiodX) != 0){
	text(RperiodX,RperiodY,rep(">",length(RperiodX)),col="gray70",cex=0.8)}
if (length(RcommaX) != 0){
	text(RcommaX,RcommaY,rep("<",length(RcommaX)),col="gray70",cex=0.8)}
if (length(RcharAX) != 0){
	text(RcharAX,RcharAY,rep("A",length(RcharAX)),col="springgreen4",cex=1,font=2)}
if (length(RcharTX) != 0){
	text(RcharTX,RcharTY,rep("T",length(RcharTX)),col="red2",cex=1,font=2)}
if (length(RcharGX) != 0){
	text(RcharGX,RcharGY,rep("G",length(RcharGX)),col="orange1",cex=1,font=2)}
if (length(RcharCX) != 0){
	text(RcharCX,RcharCY,rep("C",length(RcharCX)),col="navy",cex=1,font=2)}
"""

middlecode2 = """
if (length(het) != 0){
	text(het,Rowlength+2,rep("1",length(Rowlength+2)),col="black",cex=1)
}
if (length(homo) != 0){
	text(homo,Rowlength+2,rep("2",length(Rowlength+2)),col="black",cex=1)   
}
if (length(vcf_call_ins) != 0){
	text(vcf_call_ins,Rowlength+2,rep("*",length(Rowlength+2)),col="black",cex=1.2, family=\"mono\")	
}
if (length(vcf_call_del) != 0){
	text(vcf_call_del,Rowlength+2,rep("=",length(Rowlength+2)),col="black",cex=1.2, family=\"mono\")	
}
"""

#################### User defined functions #################### 
def printMsg(msg):
	usage= """#####################################################
##      ####   ###   ####   #####   ###   #         #
##     #      #   #  #   #  #      #   #  #         #
##     #      #####  ####   #####  #####  #         #
##     #      #   #  #   #  #      #   #  #         # 
##      ####  #   #  #   #  #####  #   #  #####     # 
#####################################################
Usage: ./careal -b <BAM or BAMs in TXT> -t <terget or targets in TXT> -r <reference genome in FASTA> [-options]

Options:   -g/--gene		<string>	gene information file in BED
	   -v/--with-vcf	<string>	query VCF file [TRUE|FALSE], default is FALSE
	   -o/--folder		<string>	output directory name
	   -n/--cpu		<number>	number of cores to use, default is 4
"""
	
	invalid_range = "Invalid range::end position must be bigger than Start position."
	unknown_error = "Critical error::please contact the author."
	
	if msg == "about":
		print about
	elif msg == "usage":
		print usage
	elif msg == "invalid_range":
		print invalid_range
	elif msg == "unknown_error":
		print unknown_error

def checkArguments(options):
	if not options.infile or not os.path.exists(options.infile):
		print("No such BAM file or '-b' argument required.\n")
		return False
	if not options.targetfile: # or not os.path.exists(options.targetfile):
		print("Invalid target or '-t' argument required.\n")
		return False
	if not options.referencefile or not os.path.exists(options.referencefile):
		print("No such reference genome file or '-r' argument required.\n")
		return False
	global geneinformation_flag
	geneinformation_flag = True
	if not options.geneinformation or not os.path.exists(options.geneinformation):
		geneinformation_flag = False
	try:
		int(options.nr_cpus)
	except Exception, e:
		print("Invalid number for multi-processing: %s\n"%(options.nr_cpus))
		return False
	return True

def istext(filename):
	s=open(filename).read(512)
	text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
	_null_trans = string.maketrans("", "")
	if not s:
	# Empty files are considered text
		return True
	if "\0" in s:
	# Files with null bytes are likely binary
		return False
	# Get the non-text characters (maps a character to itself then
	# use the 'remove' option to get rid of the text characters.)
	t = s.translate(_null_trans, text_characters)
	# If more than 30% non-text characters, then
	# this is considered a binary file
	if float(len(t))/float(len(s)) > 0.30:
		return False
	return True
	
def mean(nucleotide):
	num = 0
	for n in nucleotide:
		num += len(n)
	return float(num / len(nucleotide))
	
def Calculate_plotting_values(input_bam, sample_name, input_range, queue):
	## Input target adjustment
	range_items=[]
	if '-' not in input_range :
		point_items = input_range.strip().split(':')
		chr = point_items[0]
		point_pos = point_items[1]
		start_pos_viewer = int(point_pos) - 40
		end_pos_viewer = int(point_pos) + 40
		gene_pos = int(point_pos)
	elif '-' in input_range :
		range_items.append(input_range.strip().split(':')[0])
		range_items.append(input_range.strip().split(':')[1].split('-')[0])
		range_items.append(input_range.strip().split(':')[1].split('-')[1])
		chr = range_items[0]
		range_start = float(range_items[1])
		range_end = float(range_items[2])
		range = int(range_end) - int(range_start)
		point_pos = int(range_start + int((range_end - range_start)/2))
		gene_pos = int(range_start) + int(round((range_end - range_start)/2))
		if range <= 80 and range >= 0 :
			residue = round((80-range)/2)
			start_pos_viewer = int(range_start) - residue
			end_pos_viewer = int(range_end) + residue
		elif range > 80 :
			start_pos_viewer = int(range_start)
			end_pos_viewer = int(range_end)
		elif range < 0 :
			printMsg("usage")			

	## Extract sequence from the human reference genome and BAM file(s) using SAMTools
	range_in_viewer = int(end_pos_viewer) - int(start_pos_viewer)
	os.environ["COLUMNS"] = str(int(range_in_viewer + 300))

	cmd_out = samtools_env + " tview -d T -p " + chr + ":" + str(int(start_pos_viewer)) + "-" + str(int(end_pos_viewer)) + " " + input_bam + " " + reference
	process_out = subprocess.Popen(cmd_out,stdout=subprocess.PIPE,stderr=None,shell=True)
	output = process_out.stdout.read().strip().split('\n')[1]

	ref_list = []
	for base in output:
		ref_list.append(base)

	idx = 0
	charnum = 0
	for base in ref_list:
		idx = idx + 1
		if base != "*" :
			charnum = charnum + 1
			if charnum > range_in_viewer :
				break

	range_in_viewer_ref = idx
	os.environ["COLUMNS"] = str(range_in_viewer_ref)

	cmd_out = samtools_env + " tview -d T -p " + chr + ":" + str(int(start_pos_viewer)) + "-" + str(int(end_pos_viewer)) + " " + input_bam + " " + reference
	process_out = subprocess.Popen(cmd_out,stdout=subprocess.PIPE,stderr=None,shell=True)
	output = process_out.stdout.read().strip()

	## Get reference sequence and reads in BAM file(s)
	items = output.split('\n')

	refseq = []
	bamseq = []
	linenum = 0
	for seq in items :
		if linenum == 0 :
			linenum = linenum + 1
		elif linenum == 1 :
			refseq.append(seq)
			linenum = linenum + 1
		elif linenum == 2 :
			linenum = linenum + 1
		elif linenum >= 3 :
			bamseq.append(seq)
			linenum = linenum + 1

	num_none = list(refseq[0].upper()).count('N')
	if num_none > 0 :
		ref_start = int(end_pos_viewer) - num_none + 1
		ref_end = int(end_pos_viewer)

		cmd_ref = samtools_env + " faidx " + reference + " " + chr + ":" + str(ref_start) + "-" + str(ref_end)
		ref_out = subprocess.Popen(cmd_ref,stdout=subprocess.PIPE,stderr=None,shell=True)
		outref = ref_out.communicate()[0].strip().split('\n')[1:]

		ref_list_tmp = []
		ref_list_tmp.append(refseq[0].upper().replace('N',''))
		for outref_item in outref :
			ref_list_tmp.append(outref_item)
		ref_final = ''.join(ref_list_tmp)
		refseq[0] = ref_final

	## Declaration of variables
	asterX = []
	asterY = []
	charAX = []
	charAY = []
	charTX = []
	charTY = []
	charGX = []
	charGY = []
	charCX = []
	charCY = []
	commaX = []
	commaY = []
	periodX = []
	periodY = []
	indelX = []
	indelY = []

	rowlength = len(refseq) + len(bamseq)

	numrow = ""
	cat_bin = ""
	if rowlength >0 and rowlength <= 40 :
		numrow = 40
		cat_bin = 1
	elif rowlength > 40 and rowlength <= 80 :
		numrow = 80
		cat_bin = 2
	elif rowlength > 80 :
		numrow = round(rowlength,-1)
		cat_bin = 3

	ref_list = []
	R_idx = 1
	orf_idx = 0
	orf_list1 = []
	orf_list2 = []
	orf_list3 = []
	orf_rest1 = []
	orf_rest1.append(1)
	orf_rest2 = []
	orf_rest3 = []
	ins_list = []
	
	if vcffile_flag:
		vcfline = 2
	else:
		vcfline = 0

	## Reference sequence
	ref_insertion = []
	for base in refseq[0]:
		if base == "*" :
			asterX.append(R_idx)
			asterY.append(numrow + vcfline + 2)
			ins_list.append(R_idx)
			R_idx = R_idx + 1
			ref_list.append(base)
			ref_insertion.append(0)

		elif base == "A" or base == "a":
			charAX.append(R_idx)
			charAY.append(numrow + vcfline + 2)
			R_idx = R_idx + 1
			ref_list.append(base)
			ref_insertion.append(1)
			if orf_idx % 3 == 1 :
				orf_list1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
				if orf_idx != 1 :
					orf_rest3.append(R_idx - 1)
			elif orf_idx % 3 == 2 :
				orf_list2.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest3.append(R_idx - 1)
			elif orf_idx != 0 and orf_idx % 3 == 0 :
				orf_list3.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
			orf_idx = orf_idx + 1

		elif base == "C" or base == "c":
			charCX.append(R_idx)
			charCY.append(numrow + vcfline + 2)
			R_idx = R_idx + 1
			ref_list.append(base)
			ref_insertion.append(1)
			if orf_idx % 3 == 1 :
				orf_list1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
				if orf_idx != 1 :
					orf_rest3.append(R_idx - 1)
			elif orf_idx % 3 == 2 :
				orf_list2.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest3.append(R_idx - 1)
			elif orf_idx != 0 and orf_idx % 3 == 0 :
				orf_list3.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
			orf_idx = orf_idx + 1

		elif base == "T" or base == "t":
			charTX.append(R_idx)
			charTY.append(numrow + vcfline + 2)
			R_idx = R_idx + 1
			ref_list.append(base)
			ref_insertion.append(1)
			if orf_idx % 3 == 1 :
				orf_list1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
				if orf_idx != 1 :
					orf_rest3.append(R_idx - 1)
			elif orf_idx % 3 == 2 :
				orf_list2.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest3.append(R_idx - 1)
			elif orf_idx != 0 and orf_idx % 3 == 0 :
				orf_list3.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
			orf_idx = orf_idx + 1

		elif base == "G" or base == "g":
			charGX.append(R_idx)
			charGY.append(numrow + vcfline + 2)
			R_idx = R_idx + 1
			ref_list.append(base)
			ref_insertion.append(1)
			if orf_idx % 3 == 1 :
				orf_list1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
				if orf_idx != 1 :
					orf_rest3.append(R_idx - 1)
			elif orf_idx % 3 == 2 :
				orf_list2.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest3.append(R_idx - 1)
			elif orf_idx != 0 and orf_idx % 3 == 0 :
				orf_list3.append(R_idx - 1)
				orf_rest1.append(R_idx - 1)
				orf_rest2.append(R_idx - 1)
			orf_idx = orf_idx + 1
		else :
			R_idx = R_idx + 1
			ref_list.append(base)	
			ref_insertion.append(0)
	cov_length = len(ref_list)	
	target_idx = round((int(end_pos_viewer) - int(start_pos_viewer))/2)

	## Six reading frames
	idx = 0
	charnum = 0
	for base in ref_list:
		idx = idx + 1
		if base != "*" :
			charnum = charnum + 1
			if charnum > target_idx :
				break
	viewer_idx  = idx

	seq = refseq[0].replace('*','').upper()
	aa_list1 = []
	aa_list2 = []
	aa_list3 = []
	for i in xrange(0,len(seq),3) :
		codon1 = seq[i:i+3]
		codon2 = seq[i+1:i+4]
		codon3 = seq[i+2:i+5]
		if len(list(codon1)) == 3 :
			amino_acid1 = codon_table[codon1]
			aa_list1.append(amino_acid1)
		if len(list(codon2)) == 3 :
			amino_acid2 = codon_table[codon2]
			aa_list2.append(amino_acid2)
		if len(list(codon3)) == 3 :
			amino_acid3 = codon_table[codon3]
			aa_list3.append(amino_acid3)

	## Calculate NUCLEOTIDE-BASEs position in matrix
	BAM_input = []
	linenum = 0
	var_list = []
	for subbam in bamseq :
		linenum = linenum + 1
		R_idx = 1
		for base in subbam :
			if base == "*" :
				if ref_list[R_idx-1] == "*" :
					asterX.append(R_idx)
					asterY.append(numrow -linenum + 1)
					if R_idx == viewer_idx :
						var_list.append(numrow - linenum + 1)
					R_idx = R_idx + 1
				elif ref_list[R_idx-1] != "*" :
					indelX.append(R_idx)
					indelY.append(numrow -linenum + 1)
					if R_idx == viewer_idx :
						var_list.append(numrow - linenum + 1)
					R_idx = R_idx + 1
			elif base == "A" or base == "a":
				charAX.append(R_idx)
				charAY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == "T" or base == "t":
				charTX.append(R_idx)
				charTY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == "G" or base == "g":
				charGX.append(R_idx)
				charGY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == "C" or base == "c":
				charCX.append(R_idx)
				charCY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == "," :
				commaX.append(R_idx)
				commaY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == "." :
				periodX.append(R_idx)
				periodY.append(numrow - linenum + 1)
				if R_idx == viewer_idx :
					var_list.append(numrow - linenum + 1)
				R_idx = R_idx + 1
			elif base == " " :
				R_idx = R_idx + 1
	numcol = len(ref_list)
	seq_depth=len(var_list)
	if len(var_list) > 0 :
		lower_idx = min(var_list)
	else :
		lower_idx = 1

	## Calculate coverage
	coverage_list = []
	for covnum in xrange(1,cov_length + 1) :
		if covnum not in ins_list :
			pos_cov = charAX.count(covnum) + charTX.count(covnum) + charGX.count(covnum) + charCX.count(covnum) + periodX.count(covnum) + commaX.count(covnum) - 1
			coverage_list.append(pos_cov)
		elif covnum in ins_list :
			pos_cov = charAX.count(covnum) + charTX.count(covnum) + charGX.count(covnum) + charCX.count(covnum) + periodX.count(covnum) + commaX.count(covnum)
			coverage_list.append(pos_cov)

	max_cov = max(coverage_list)
	if max_cov == 0 :
		max_cov = 0.000000001 # pseudo value for barplot ratio

	## Query VCF
	vcf_dir = input_bam.strip().split(sample_name)[0]
	call_list = []
	call_list_plot = {}
	vcffile = vcf_dir + sample_name + ".vcf.gz"
	if vcffile_flag:
		tb = tabix.open(vcffile)
		query = chr + ":" + str(start_pos_viewer) + "-" + str(end_pos_viewer)
		try:
			output_call = tb.querys(query)
			for vcf_call in output_call :
				chr = vcf_call[0]
				pos = vcf_call[1]
				rsid = vcf_call[2]
				ref = vcf_call[3]
				alt = vcf_call[4]
				gt = vcf_call[9].split(':')[0]
				gt2 = gt.replace("|", "/")
				alleles = gt2.split("/")

				# if SNV
				if len(ref) == 1 and mean(alt.split(",")) == 1:
					if alleles[0] == alleles[-1] and int(alleles[0]) != 0 and int(alleles[-1]) != 0:
						gt_info = "homo"
					elif alleles[0] != alleles[-1]:
						gt_info = "hetero"
				# if INDEL
				else:
					if len(ref) < len(alt): # insertion
						gt_info = str(len(ref) - len(alt))
					else:
						gt_info = str(int(len(ref)) - int(len(alt)))

				vcfcall = '{0:17s} {1:12s} {2:6s} {3:6s} {4:15s}'.format(chr+":"+str(pos), rsid, ref, alt, gt)
				call_list.append(vcfcall)
				call_list_plot[pos] = gt_info			
		except tabix.TabixError:
			pass

	## Get gene name
	if geneinformation_flag: 
		gene_ref = geneinformation
		tb = tabix.open(gene_ref)
		query = chr + ":" + str(point_pos) + "-" + str(point_pos)
		try:
			tb_results = tb.querys(query)

			gene_list = []
			for result in tb_results:
				gene_list.append(result[3])

			genename = '_'.join(gene_list)
		except tabix.TabixError:
			genename = ""
	else:
		genename = ""

	## R script, 1. variables
	pos_name=""
	pos_str=""
	pos_outname=chr+":"+str(point_pos)

	if '-' not in input_range :
		pos_name=chr+"_"+point_pos
	elif '-' in input_range :
		pos_name=chr+"_"+str(int(range_start))+'_'+str(int(range_end))
		target_half = round((int(end_pos_viewer) - int(start_pos_viewer))/2)
		target_position = int(start_pos_viewer) + target_half


	lib_str="library(graphics)"
	gene_str="gene <- " + "'" + genename + "'"	
	pos_str="position <- " + "'" + pos_name + "'"
	pos_outstr="outposition <- " + "'" + pos_outname + "'"		
	Ridx_str="Ridx <- " + str(viewer_idx)
	Rowlen_str="Rowlength <- " + str(numrow)
	Collen_str="Collength <- " + str(numcol)
	Lowidx_str="Loweridx <- " + str(lower_idx)
	Bin_str="Binnum <- " + str(cat_bin)
	dep_str="Seqdepth <- " + str(seq_depth)
	filename_str="BAM <- '" + str(sample_name) + "'"
	blank_str="RasterX <- c()\nRasterY <- c() \nRcharAX <- c() \nRcharAY <- c() \nRcharTX <- c() \nRcharTY <- c() \nRcharGX <- c() \nRcharGY <- c() \nRcharCX <- c()\nRcharCY <- c() \nRcommaX <- c() \nRcommaY <- c() \nRperiodX <- c() \nRperiodY <- c() \nindelX <- c() \nindelY <- c()\nCALL <- c()\nCoverage <- c()\nhet <- c()\nhomo <- c()\nvcf_call_ins <- c()\nvcf_call_del <- c()\n"

	f_R=open(output_dir + "/CaReAl_" + sample_name + "_" + pos_name + ".R",'w')
	f_R.write(lib_str + "\n" + gene_str + '\n' + blank_str + '\n' + Ridx_str + '\n' + Rowlen_str + '\n' + Collen_str + '\n' + pos_str + '\n' + pos_outstr + '\n' + Lowidx_str + '\n' + Bin_str + '\n' + dep_str + '\n' + filename_str + '\n\n')

	for subidx in asterX :
		asterX_str="RasterX <- c(RasterX," + str(subidx) + ")"
		f_R.write(asterX_str + '\n')
	for subidx in asterY :
		asterY_str="RasterY <- c(RasterY," + str(subidx) + ")"
		f_R.write(asterY_str + '\n')
	for subidx in indelX :
		indelX_str="indelX <- c(indelX," + str(subidx) + ")"
		f_R.write(indelX_str + '\n')
	for subidx in indelY :
		indelY_str="indelY <- c(indelY," + str(subidx) + ")"
		f_R.write(indelY_str + '\n')
	for subidx in charAX :
		charAX_str="RcharAX <- c(RcharAX," + str(subidx) + ")"
		f_R.write(charAX_str + '\n')
	for subidx in charAY :
		charAY_str="RcharAY <- c(RcharAY," + str(subidx) + ")"
		f_R.write(charAY_str + '\n')
	for subidx in charCX :
		charCX_str="RcharCX <- c(RcharCX," + str(subidx) + ")"
		f_R.write(charCX_str + '\n')
	for subidx in charCY :
		charCY_str="RcharCY <- c(RcharCY," + str(subidx) + ")"
		f_R.write(charCY_str + '\n')
	for subidx in charGX :
		charGX_str="RcharGX <- c(RcharGX," + str(subidx) + ")"
		f_R.write(charGX_str + '\n')
	for subidx in charGY :
		charGY_str="RcharGY <- c(RcharGY," + str(subidx) + ")"
		f_R.write(charGY_str + '\n')
	for subidx in charTX :
		charTX_str="RcharTX <- c(RcharTX," + str(subidx) + ")"
		f_R.write(charTX_str + '\n')
	for subidx in charTY :
		charTY_str="RcharTY <- c(RcharTY," + str(subidx) + ")"
		f_R.write(charTY_str + '\n')
	for subidx in commaX :
		commaX_str="RcommaX <- c(RcommaX," + str(subidx) + ")"
		f_R.write(commaX_str + '\n')
	for subidx in commaY :
		commaY_str="RcommaY <- c(RcommaY," + str(subidx) + ")"
		f_R.write(commaY_str + '\n')
	for subidx in periodX :
		periodX_str="RperiodX <- c(RperiodX," + str(subidx) + ")"
		f_R.write(periodX_str + '\n')
	for subidx in periodY :
		periodY_str="RperiodY <- c(RperiodY," + str(subidx) + ")"
		f_R.write(periodY_str + '\n')

	for cov_item in coverage_list :
		cov_str="Coverage <- c(Coverage," + str(float(cov_item)/max_cov) + ")"
		f_R.write(cov_str + '\n')

	for call_info in call_list :
		call_str="CALL <- c(CALL,'" + str(call_info) + "')"
		f_R.write(call_str + '\n')	

	real_pos = int(start_pos_viewer)-1
	for idx, ref_ins in enumerate(ref_insertion, start=1):
		if int(ref_ins) == 1:
			real_pos += ref_ins
			if str(real_pos) in call_list_plot:
				if call_list_plot[str(real_pos)] == "hetero": # SNV, heterozygous
					hetero="het <- c(het,'" + str(idx) + "')"
					f_R.write(hetero + '\n')	
				elif call_list_plot[str(real_pos)] == "homo": # SNV, homozygous
					homo="homo <- c(homo,'" + str(idx) + "')"
					f_R.write(homo + '\n')
				else: # deletion
					ins_cnt = int(call_list_plot[str(real_pos)])
					if ins_cnt > 0:
						for x in xrange(0,ins_cnt):
							vcf_call_del="vcf_call_del <- c(vcf_call_del,'" + str(idx+x+1) + "')"
							f_R.write(vcf_call_del + '\n')
		elif int(ref_ins) == 0: # insertion
			if str(real_pos) in call_list_plot:
				try:
					ref_ins_idx_value = int(call_list_plot[str(real_pos)])
					ins_cnt = int(call_list_plot[str(real_pos)])
					if ins_cnt < 0:
						for x in xrange(0,abs(ins_cnt)):
							vcf_call_ins="vcf_call_ins <- c(vcf_call_ins,'" + str(idx+x) + "')"
							f_R.write(vcf_call_ins + '\n')
				except ValueError:
					pass				

	## R script, 2. plotting
	if int(numrow) < 80:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline + 6) + "\n"
	elif int(numrow) >= 80 and int(numrow) < 120:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline + 4) + "\n"
	elif int(numrow) >= 120 and int(numrow) < 160:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline + 2) + "\n"
	elif int(numrow) >= 160 and int(numrow) < 200:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline + 1) + "\n"
	elif int(numrow) >= 200 and int(numrow) < 250:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline) + "\n"
	elif int(numrow) >= 250 and int(numrow) < 350:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - 2.5) + "\n"
	elif int(numrow) >= 350 and int(numrow) < 450:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.02)) + "\n"
	elif int(numrow) >= 450 and int(numrow) < 550:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.021)) + "\n"
	elif int(numrow) >= 550 and int(numrow) < 650:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.022)) + "\n"
	elif int(numrow) >= 650 and int(numrow) < 750:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.023)) + "\n"
	elif int(numrow) >= 750 and int(numrow) < 850:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.024)) + "\n"
	elif int(numrow) >= 850 and int(numrow) < 950:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.025)) + "\n"
	elif int(numrow) >= 950 and int(numrow) < 1050:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.027)) + "\n"
	elif int(numrow) >= 1050 and int(numrow) < 1150:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.028)) + "\n"
	elif int(numrow) >= 1150 and int(numrow) < 1250:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.029)) + "\n"
	elif int(numrow) >= 1250 and int(numrow) < 1350:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.03)) + "\n"
	elif int(numrow) >= 1350 and int(numrow) < 1450:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.031)) + "\n"
	elif int(numrow) >= 1450 and int(numrow) < 1550:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.032)) + "\n"
	elif int(numrow) >= 1550 and int(numrow) < 1650:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.033)) + "\n"
	elif int(numrow) >= 1650 and int(numrow) < 1750:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.034)) + "\n"
	elif int(numrow) >= 1750:
		rowlength2 = "Rowlength2 <- " + str(numrow + vcfline - (int(numrow)*0.035)) + "\n"

	if vcffile_flag:
		vcf_additional = "vcf_additional <- (ceiling(length(CALL)/10) * 100)\n"
	else:
		vcf_additional = "vcf_additional <- 0\n"
		
	pre_code = "HOME <- \"" + output_dir + "/\"\n"
	if geneinformation_flag:
		pre_code += "png_file_name <- paste(HOME,position,\"__\",gene,\"__\",BAM,\".png\",sep=\"\")\n"
		pre_code += 'barplot_title <- paste(BAM, ", ", gsub("_", ":", outposition), ", ", gsub("_", "/", gene), ", DP=", Seqdepth, sep="")\n'
	else:
		pre_code += "png_file_name <- paste(HOME,position,\"__\",BAM,\".png\",sep=\"\")\n"
		pre_code += 'barplot_title <- paste(BAM, ", ", gsub("_", ":", outposition), ", ", "DP=", Seqdepth, sep="")\n'
	pre_code += rowlength2
	pre_code += vcf_additional
	pre_code += middlecode1
	f_R.write(pre_code)

	if vcffile_flag:
		f_R.write(middlecode2)
		f_R.write('abline(h=Rowlength+3 ,col = "black", lty = 1)\n')
		f_R.write('abline(h=Rowlength+1 ,col = "black", lty = 1)\n')
		f_R.write('rect(Ridx - 0.5, Loweridx -1, Ridx+0.5, Rowlength+4.7, col=NA, border="black", lty=1, lwd=0.5)\n')
	else:
		f_R.write('abline(h=Rowlength+1 ,col = "black", lty = 1)\n')
		f_R.write('rect(Ridx - 0.5, Loweridx -1, Ridx+0.5, Rowlength+2.7, col=NA, border="black", lty=1, lwd=0.5)\n')

	
	for orf1 in xrange(0,len(aa_list1)) :
		orf1_str="text(" + str(orf_list1[orf1]) + "," + str(numrow+6+vcfline) + ",'" + aa_list1[orf1] + "',cex=1, family=\"mono\", font=2)"
		f_R.write(orf1_str + '\n')
	for orf2 in xrange(0,len(aa_list2)) :
		orf2_str="text(" + str(orf_list2[orf2]) + "," + str(numrow+5+vcfline) + ",'" + aa_list2[orf2] + "',cex=1, family=\"mono\", font=2)"
		f_R.write(orf2_str + '\n')
	for orf3 in xrange(0,len(aa_list3)) :
		orf3_str="text(" + str(orf_list3[orf3]) + "," + str(numrow+4+vcfline) + ",'" + aa_list3[orf3] + "',cex=1, family=\"mono\", font=2)"
		f_R.write(orf3_str + '\n')

	for rest1 in xrange(0,len(orf_rest1)) :
		rest1_str="text(" + str(orf_rest1[rest1]) + "," + str(numrow+6+vcfline) + ",'.',cex=1, family=\"mono\")"
		f_R.write(rest1_str + '\n')
	for rest2 in xrange(0,len(orf_rest2)) :
		rest2_str="text(" + str(orf_rest2[rest2]) + "," + str(numrow+5+vcfline) + ",'.',cex=1, family=\"mono\")"
		f_R.write(rest2_str + '\n')
	for rest3 in xrange(0,len(orf_rest3)) :
		rest3_str="text(" + str(orf_rest3[rest3]) + "," + str(numrow+4+vcfline) + ",'.',cex=1, family=\"mono\")"
		f_R.write(rest3_str + '\n')

	if vcffile_flag:
		post_code = "par(mar=c(0,0,0,0))\n"
		post_code += "plot(x=NULL,y=NULL,xlim=c(0,1), ylim=c(0,ceiling(length(CALL)/10) * 10 + 0.1), ann=F, bty='n', type='n', xaxt='n', yaxt='n')\n"
		post_code += "for (i in length(CALL):1) { text(0, (ceiling(length(CALL)/10) * 10) - length(CALL) + i, CALL[length(CALL)-i+1], pos=4, font=1, family=\"mono\") }\n"
	
		f_R.write(post_code)

	f_R.write("dev.off()")	
	f_R.close()

	## R script, 3. execution Rscript
	os.chdir(output_dir)
	cmd = "R CMD BATCH CaReAl_" + sample_name + "_" + pos_name + ".R"
	os.system(cmd)
	
	return

#################### Main ####################
def main():
	if not checkArguments(options):
		printMsg("usage")
		sys.exit(2)

	global working_dir, reference, vcffile_flag, geneinformation
	working_dir = os.getcwd()
	infile = options.infile
	targetfile = options.targetfile 
	if "/" in options.referencefile:
		reference = options.referencefile
	else:
		reference = working_dir + "/" + options.referencefile
	if geneinformation_flag != False:
		if "/" in options.geneinformation:
			geneinformation = options.geneinformation
		else:
			geneinformation = working_dir + "/" + options.geneinformation
	vcffile_flag = options.vcffile_flag
	
	## BAM file(s)	
	infiles = []
	if istext(infile):
		infile = open(infile)
		lines = infile.readlines()
		for line in lines :
			infiles.append(line.strip())
	else :
		infiles.append(infile)

	## Input target(s)
	targets = []
	if os.path.isfile(targetfile):
		infile = open(targetfile)
		lines = infile.readlines()
		for line in lines :
			targets.append(line.strip())
	else :
		targets.append(targetfile)

	## Set output directory
	global output_dir
	output_dir = working_dir + "/" + options.outputfolder
	
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
		os.chmod(output_dir, 0777)
	os.chdir(working_dir)	
	
	## multi threading
	manager = multiprocessing.Manager()
	queue = manager.Queue()
	pool = multiprocessing.Pool(int(options.nr_cpus))
	
	jobs = []
	for input_bam in infiles:
		sample_name = os.path.basename(input_bam).split('.bam')[0]
		print "Sample: " + str(sample_name)
		for input_range in targets:
			job = pool.apply_async(Calculate_plotting_values, (input_bam,sample_name,input_range,queue))
			jobs.append(job)
			
	for job in jobs:
		job.get()
		
	queue.put('kill')
	pool.close()
	pool.join()

#################### Debug?! ####################	
if __name__ == "__main__":
	start_time = time.time()
        main()
	end_time = time.time()

	print("\nProcess done! : " + str(end_time - start_time) + " second(s)")
	
## CaReAl v1.1 (@'3(o.o);
## EOF 2017.12.29
