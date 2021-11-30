import argparse
import os
import pathlib

import numpy as np

from source.bed_func import get_repetitive_regions, parse_bed_regions, REP_BUFF
from source.corgi_mappability import MappabilityTrack
from source.sv_blacklist import get_common_svs, split_common_svs
from source.vcf_func import filter_vcf_by_readcount, filter_vcf_by_size, filter_vcf_by_type, sort_vcf, split_repetitive_svs
from source.workflow_cmds import EXE_DICT, mv, rm

VALID_REFS = ['hg19', 'hg38']

#
#
#
def main(raw_args=None):
	parser = argparse.ArgumentParser(description='Filt-SV v0.1', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,  metavar='input.vcf',   help="* Input VCF")
	parser.add_argument('-o',  type=str, required=True,  metavar='out_dir/',    help="* Output directory (will be created)")
	parser.add_argument('-r',  type=str, required=False, metavar='',            help="Reference build hg19/hg38", default='hg38')
	#parser.add_argument('-c',  type=str, required=False, metavar='tools.cfg',   help="Config file with exe paths",       default=None)
	#
	parser.add_argument('--readcount',     type=int,   required=False, metavar='', help="Initial readcount filter",          default=2)
	parser.add_argument('--readcount-bnd', type=int,   required=False, metavar='', help="BND readcount filter",              default=5)
	parser.add_argument('--vaf',           type=float, required=False, metavar='', help="Initial VAF filter",                default=0.2)
	parser.add_argument('--vaf-bnd',       type=float, required=False, metavar='', help="BND VAF filter",                    default=0.3)
	parser.add_argument('--control-freq',  type=float, required=False, metavar='', help="Minimum SV frequency in controls",  default=0.5)
	parser.add_argument('--min-size',      type=int,   required=False, metavar='', help="Minimum size (non-repetitive SVs)", default=50)
	parser.add_argument('--min-size-rep',  type=int,   required=False, metavar='', help="Minimum size (repetitive SVs)",     default=300)
	args = parser.parse_args()

	INPUT_VCF = args.i
	OUT_DIR   = args.o
	#
	if OUT_DIR[-1] != '/':
		OUT_DIR += '/'
	#
	OUT_LOG = OUT_DIR + 'log.txt'
	#
	REF_NAME = args.r
	if REF_NAME not in VALID_REFS:
		print('Error: -r must be one of: ' + ', '.join(VALID_REFS))
		exit(1)
	#
	TEMP_DIR = OUT_DIR + 'temp/'
	os.system('mkdir -p ' + OUT_DIR)
	os.system('mkdir -p ' + TEMP_DIR)
	#
	CONFIG = args.c
	if CONFIG != None:
		print('reading config file', CONFIG + '...')
		f = open(CONFIG, 'r')
		for line in f:
			splt = line.strip().split('\t')
			EXE_DICT[splt[0]] = splt[1]
		f.close()
	#
	READCOUNT_INIT = args.readcount
	READCOUNT_BND  = args.readcount_bnd
	VAF_INIT       = args.vaf
	VAF_BND        = args.vaf_bnd
	CONTROL_FREQ   = args.control_freq
	MIN_NON_SV_LEN = args.min_size
	MIN_REP_SV_LEN = args.min_size_rep
	#
	FINAL_VCF      = OUT_DIR + 'filtered_svs.vcf'

	# absolute path to this script
	sim_path = pathlib.Path(__file__).resolve().parent
	res_path = str(sim_path) + '/resources/' + REF_NAME + '/'
	#
	BED_DICT_GENE = parse_bed_regions(res_path + 'gencode_v38_gene.bed.gz', 4)
	BED_DICT_EXON = parse_bed_regions(res_path + 'gencode_v38_exon.bed.gz', 4)
	COMMON_SVS    = get_common_svs(res_path + 'common-svs.bed.gz')
	CONTROL_SVS   = get_common_svs(res_path + 'control-svs.bed.gz', min_control_freq=CONTROL_FREQ)
	REPEAT_BED    = get_repetitive_regions(res_path + 'simple-repeats.bed.gz')
	REPEAT_TRACK  = MappabilityTrack(res_path + 'simple-repeats.bed.gz', bed_buffer=REP_BUFF)

	#
	# [0] INITILIZATION
	#
	TRUVARI_LOG = TEMP_DIR + 'truvari.log'
	TEMP_VCF_0  = TEMP_DIR + 'temp0.vcf'
	TEMP_VCF_1  = TEMP_DIR + 'temp1.vcf'
	TEMP_VCF_2  = TEMP_DIR + 'temp2.vcf'
	#
	NON_REPEAT_VCF = TEMP_DIR + 'unique.vcf'
	REPEAT_VCF     = TEMP_DIR + 'repeats.vcf'
	#
	NR_COMMON = NON_REPEAT_VCF[:-4] + '-common.vcf'
	R_COMMON  = REPEAT_VCF[:-4]     + '-common.vcf'
	NR_NOT_COMMON = NON_REPEAT_VCF[:-4] + '-notcommon.vcf'
	R_NOT_COMMON  = REPEAT_VCF[:-4]     + '-notcommon.vcf'
	#
	NR_NOT_COMMON_BND     = NR_NOT_COMMON[:-4] + '-bnd.vcf'
	NR_NOT_COMMON_NOT_BND = NR_NOT_COMMON[:-4] + '-notbnd.vcf'
	#
	NR_CONTROL = NR_NOT_COMMON[:-4] + '-control.vcf'
	R_CONTROL  = R_NOT_COMMON[:-4]  + '-control.vcf'
	NR_FINAL   = NON_REPEAT_VCF[:-4] + '-final.vcf'
	R_FINAL    = REPEAT_VCF[:-4]     + '-final.vcf'

	#
	# [1] READCOUNT FILTER --> SPLIT REPEATS --> SIZE FILTERS
	#
	print('')
	print('-- initial filtering: --')
	filter_vcf_by_readcount(INPUT_VCF, TEMP_VCF_0, min_readcount=READCOUNT_INIT, min_vaf=VAF_INIT)
	print('')
	print('-- identifying repetitive svs: --')
	split_repetitive_svs(TEMP_VCF_0, TEMP_VCF_1, TEMP_VCF_2, REPEAT_BED, REPEAT_TRACK)
	print('')
	print('-- filtering non-repeats: --')
	filter_vcf_by_size(TEMP_VCF_2, NON_REPEAT_VCF, min_size=MIN_NON_SV_LEN)
	print('')
	print('-- filtering repeats: --')
	filter_vcf_by_size(TEMP_VCF_1, REPEAT_VCF, min_size=MIN_REP_SV_LEN)
	print('')
	rm(TEMP_VCF_0)
	rm(TEMP_VCF_1)
	rm(TEMP_VCF_2)

	#
	# [2] SUBTRACT COMMON SVS
	#
	print('-- subtracting common svs (non-repeats): --')
	split_common_svs(NON_REPEAT_VCF, NR_COMMON, NR_NOT_COMMON, COMMON_SVS)
	print('')
	print('-- subtracting common svs (repeats): --')
	split_common_svs(REPEAT_VCF, R_COMMON, R_NOT_COMMON, COMMON_SVS)
	print('')

	#
	# [3a] NON-REPETITIVE SVS --> SPLIT OFF TRANSLOCATIONS --> READCOUNT FILTER
	# [3b] REPETITIVE SVS --> FILTER BY SIZE
	#
	print('-- identifying translocations: --')
	filter_vcf_by_type(NR_NOT_COMMON, TEMP_VCF_0, NR_NOT_COMMON_NOT_BND, ['BND','TRA'])
	print('')
	print('-- filtering translocations: --')
	filter_vcf_by_readcount(TEMP_VCF_0, NR_NOT_COMMON_BND, min_readcount=READCOUNT_BND, min_vaf=VAF_BND)
	print('')
	rm(TEMP_VCF_0)

	#
	# [4] SUBTRACT CONTROL SVS
	#
	print('-- subtracting control svs (non-repeats): --')
	split_common_svs(NR_NOT_COMMON_NOT_BND, NR_CONTROL, NR_FINAL, CONTROL_SVS)
	print('')
	print('-- subtracting control svs (repeats): --')
	split_common_svs(R_NOT_COMMON, R_CONTROL, R_FINAL, CONTROL_SVS)
	print('')

	#
	# TODO: [5] ANNOTATE SVS THAT INTERSECT GENE REGIONS
	#

	#
	# TODO: [6] CLEANUP
	#
	mv(NR_FINAL, FINAL_VCF)

if __name__ == '__main__':
	main()
