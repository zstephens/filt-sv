import argparse
import os
import pathlib
import re

import numpy as np

from source.bed_func import get_repetitive_regions, parse_bed_regions, REP_BUFF
from source.corgi_mappability import MappabilityTrack
from source.sv_blacklist import get_common_svs, split_common_svs
from source.vcf_func import filter_vcf_by_readcount, sort_vcf, split_repetitive_svs, LEXICO_2_IND
from source.workflow_cmds import EXE_DICT, collapse_pbsv, truvari, bgzip_tabix, rm, exists_and_is_nonzero

#
#
#
def main(raw_args=None):
	parser = argparse.ArgumentParser(description='make-control-vcf v0.1', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,  metavar='input_vcf_dir/', help="* Input directory of VCFs")
	parser.add_argument('-o',  type=str, required=True,  metavar='out_dir/',       help="* Output directory (will be created)")
	parser.add_argument('-r',  type=str, required=True,  metavar='ref.fa',         help="* Reference sequence")
	parser.add_argument('-b',  type=str, required=False, metavar='regions.bed',    help="Only assess SVs in these regions", default=None)
	parser.add_argument('-c',  type=str, required=False, metavar='tools.cfg',      help="Config file with exe paths", default=None)
	args = parser.parse_args()

	INPUT_DIR = args.i
	OUT_DIR   = args.o
	REF       = args.r
	#
	if OUT_DIR[-1] != '/':
		OUT_DIR += '/'
	if args.b != None:
		INCL_BED = ' --includebed ' + args.b
	else:
		INCL_BED = ''
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
	OUT_BED = OUT_DIR + 'control_svs.bed'

	# absolute path to this script
	sim_path = pathlib.Path(__file__).resolve().parent
	res_path = str(sim_path) + '/resources/'

	#
	# [0] INITILIZATION
	#
	TRUVARI_LOG = OUT_DIR + 'truvari.log'
	#
	print('running truvari...')
	INPUT_VCF = sorted([n for n in os.listdir(INPUT_DIR) if n[-4:] == '.vcf'])
	TRUVARI_RDY_SINGLE = []
	TRUVARI_RDY_PAIR   = []
	for i in range(len(INPUT_VCF)):
		current_vcf      = INPUT_DIR + INPUT_VCF[i]
		current_vcf_sort = TEMP_DIR + INPUT_VCF[i][:-4] + '_sort.vcf'
		current_vcf_gz   = TEMP_DIR + INPUT_VCF[i][:-4] + '_sort.vcf.gz'
		if exists_and_is_nonzero(current_vcf_gz) == False:
			#print(i, current_vcf, current_vcf_sort, current_vcf_gz)
			sort_vcf(current_vcf, current_vcf_sort)
			bgzip_tabix(current_vcf_sort, current_vcf_gz, exe_dict=EXE_DICT)
		TRUVARI_RDY_SINGLE.append(current_vcf_gz)
	#
	#
	#
	sv_hits = []	# [sample_i] = {sv_in_sample_i : list_of_samples_we_were_found_in}
	for i in range(len(TRUVARI_RDY_SINGLE)):
		vcf_1 = TRUVARI_RDY_SINGLE[i]
		for j in range(len(TRUVARI_RDY_SINGLE)):
			tv_dir = OUT_DIR + 'truvari_' + str(i) + '_' + str(j)
			vcf_2  = TRUVARI_RDY_SINGLE[j]
			if exists_and_is_nonzero(tv_dir+'/tp-base.vcf') == False or exists_and_is_nonzero(tv_dir+'/tp-call.vcf') == False:
				print('running truvari on samples:',i,j)
				truvari(vcf_1, vcf_2, tv_dir, REF, out_log=TRUVARI_LOG, incl_bed=INCL_BED, exe_dict=EXE_DICT)
		sv_hits.append({})
		for j in range(len(TRUVARI_RDY_SINGLE)):
			tv_dir = OUT_DIR + 'truvari_' + str(i) + '_' + str(j)
			if exists_and_is_nonzero(tv_dir+'/tp-base.vcf'):
				f = open(tv_dir+'/tp-base.vcf', 'r')
				for line in f:
					if line[0] != '#':
						splt = line.strip().split('\t')
						my_info = splt[7]
						re_len = re.findall(r";SVLEN=.*?(?=;)",';'+my_info+';')
						re_end = re.findall(r";END=.*?(?=;)",';'+my_info+';')
						re_typ = re.findall(r";SVTYPE=.*?(?=;)",';'+my_info+';')
						my_typ = re_typ[0][8:]
						if my_typ in ['BND','TRA']:
							my_len = None
							my_end = None
						elif my_typ in ['INV']:
							my_end = int(re_end[0][5:])
							my_len = my_end - int(splt[1])
						else:
							my_len = int(re_len[0][7:])
							if len(re_end):
								my_end = int(re_end[0][5:])
							elif my_typ == 'INS':
								my_end = int(splt[1])
							elif my_typ == 'DEL' and my_len < 0:
								my_end = int(splt[1]) - my_len
							elif my_type == 'DUP' and my_len > 0:
								my_end = int(splt[1]) + my_len
							else:
								# well we tried...
								my_end = None
						my_key = (splt[0], int(splt[1]), my_end, my_typ, my_len)
						if my_key not in sv_hits[i]:
							sv_hits[i][my_key] = {}
						sv_hits[i][my_key][j] = True
				f.close()
	#
	#
	#
	all_svs = {}
	for i in range(len(TRUVARI_RDY_SINGLE)):
		sorted_sv_hits = sorted([(len(sv_hits[i][k]), k) for k in sv_hits[i].keys()])
		for n in sorted_sv_hits:
			if n[1] not in all_svs:
				all_svs[n[1]] = 0
			all_svs[n[1]] = max([all_svs[n[1]], n[0]])
	#
	#
	#
	sorted_keys = [n[2] for n in sorted([(LEXICO_2_IND[k[0]], k[1], k) for k in all_svs.keys()])]
	f = open(OUT_BED, 'w')
	for k in sorted_keys:
		print(k, all_svs[k], float(all_svs[k])/len(INPUT_VCF))
		f.write(k[0] + '\t' + str(k[1]) + '\t' + str(k[2]) + '\t' + k[3] + '\t' + 'CONTROL_FREQ={0:0.3f}'.format(float(all_svs[k])/len(INPUT_VCF)) + '\n')
	f.close()

if __name__ == '__main__':
	main()
