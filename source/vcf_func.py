import pathlib
import re
import sys

from source.bed_func import query_bed_regions

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23, 'chrY' :24, 'chrM' :25}

#
#
#
def get_repeat_hits(my_chr, pos_start, pos_end, rep_dict, rep_track, min_del_frac=0.5, get_repeat_info=False):
	[x1, x2] = sorted([pos_start, pos_end])
	hit_count = rep_track.query_range_faster(my_chr, x1, x2, True)
	# not near repeat annotation, bail!
	if hit_count <= 0:
		return []
	# throw out deletions that mostly affect non-repetitive sequence
	if x2 > x1 and hit_count < min_del_frac*(x2-x1):
		return []
	#
	if get_repeat_info:
		repeat_hit = []
		if my_chr in rep_dict:
			for n in rep_dict[my_chr]:
				[y1,y2] = [n[0],n[1]]
				if x1 <= y2 and y1 <= x2:
					repeat_hit.append(n[2])
			repeat_hit = sorted(list(set(repeat_hit)))
		return repeat_hit
	else:
		return ['who cares']

#
#
#
def parse_info_column_of_sv(my_info, start_pos):
	re_len = re.findall(r";SVLEN=.*?(?=;)",';'+my_info+';')
	re_end = re.findall(r";END=.*?(?=;)",';'+my_info+';')
	re_typ = re.findall(r";SVTYPE=.*?(?=;)",';'+my_info+';')
	my_typ = re_typ[0][8:]
	if my_typ in ['BND','TRA']:
		my_len = None
		my_end = None
	elif my_typ in ['INV']:
		my_end = int(re_end[0][5:])
		my_len = my_end - start_pos
	else:
		my_len = int(re_len[0][7:])
		if len(re_end):
			my_end = int(re_end[0][5:])
		elif my_typ == 'INS':
			my_end = start_pos
		elif my_typ == 'DEL' and my_len < 0:
			my_end = start_pos - my_len
		elif my_typ == 'DUP' and my_len > 0:
			my_end = start_pos + my_len
		else:
			# well we tried...
			my_end = None
	return (my_typ, my_len, my_end)

#
#
#		
def sort_vcf(fn_in, fn_out):
	header = ''
	dat_sort = []
	f = open(fn_in, 'r')
	for line in f:
		if line[0] == '#':
			header += line
		else:
			splt = line.split('\t')
			if splt[0] in LEXICO_2_IND:
				dat_sort.append([LEXICO_2_IND[splt[0]], int(splt[1]), line])
	dat_sort = sorted(dat_sort)
	f.close()
	f = open(fn_out, 'w')
	f.write(header)
	for n in dat_sort:
		f.write(n[2])
	f.close()

#
#
#
def filter_vcf_by_readcount(fn_in, fn_out, min_readcount=2, min_vaf=0.25):
	f  = open(fn_in, 'r')
	f2 = open(fn_out, 'w')
	nskipped = 0	# filtered
	nfailed  = 0	# skipped because we couldn't find the info we wanted, or SV is otherwise weird (e.g. multi-allele)
	nwritten = 0
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
				col_form = splt.index('FORMAT')
				col_samp = len(splt) - 1
			f2.write(line)
		else:
			splt   = line.strip().split('\t')
			splt2  = splt[col_form].split(':')
			splt3  = splt[col_samp].split(':')
			#
			# check for AD first
			if 'AD' in splt2:
				col_ad = splt2.index('AD')
				if '.' in splt3[col_ad]:	# AD values are .
					nfailed += 1
					continue
				readcounts = [int(n) for n in splt3[col_ad].split(',')]
			#
			# if we don't have AD, do we have DR and DV fields instead?
			elif 'DR' in splt2 and 'DV' in splt2:
				col_dr = splt2.index('DR')
				col_dv = splt2.index('DV')
				readcounts = [int(splt3[col_dr]), int(splt3[col_dv])]
			#
			# otherwise we failed to find what we needed.
			else:
				nfailed += 1
				continue
			#
			# why was this SV even called if zero reads support it??
			if sum(readcounts) == 0:
				nskipped += 1
				continue
			my_reads = readcounts[1]
			my_vaf   = float(my_reads)/sum(readcounts)
			#
			# multiple alleles, skipping this for now...
			if len(readcounts) > 2:
				nfailed += 1
				continue
			#
			if len(readcounts) == 2 and my_reads >= min_readcount and my_vaf >= min_vaf:
				f2.write(line)
				nwritten += 1
			else:
				nskipped += 1
	f2.close()
	f.close()
	print('SVS-pass (reads >= '+str(min_readcount)+' & vaf >= '+'{0:0.2f}'.format(min_vaf)+'):', nwritten)
	print('SVS-filt (reads <  '+str(min_readcount)+' | vaf <  '+'{0:0.2f}'.format(min_vaf)+'):', nskipped)
	print('---')
	print('SVs-fail (readcount not found / multi-allelic):', nfailed)

#
#
#
def filter_vcf_by_size(fn_in, fn_out, min_size=300):
	f  = open(fn_in, 'r')
	f2 = open(fn_out, 'w')
	nskipped = 0
	nwritten = 0
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
			f2.write(line)
		else:
			splt = line.strip().split('\t')
			(my_type, my_len, my_end) = parse_info_column_of_sv(splt[col_info], int(splt[1]))
			if my_len == None or abs(my_len) >= min_size:
				f2.write(line)
				nwritten += 1
			else:
				nskipped += 1
	f2.close()
	f.close()
	print('SVS-pass (size >= '+str(min_size)+'):', nwritten)
	print('SVS-filt (size <  '+str(min_size)+'):', nskipped)

#
#
#
def filter_vcf_by_type(fn_in, fn_out_filtered, fn_out_not_filtered, types_to_filter=[]):
	f  = open(fn_in, 'r')
	f2 = open(fn_out_filtered, 'w')
	f3 = open(fn_out_not_filtered, 'w')
	n_filt     = 0
	n_not_filt = 0
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
			f2.write(line)
			f3.write(line)
		else:
			splt = line.strip().split('\t')
			(my_type, my_len, my_end) = parse_info_column_of_sv(splt[col_info], int(splt[1]))
			if my_type in types_to_filter:
				f2.write(line)
				n_filt += 1
			else:
				f3.write(line)
				n_not_filt += 1
	f3.close()
	f2.close()
	f.close()
	out_str = ','.join([str(n) for n in types_to_filter])
	size_diff = len(out_str) - len('other')
	if size_diff >= 0:
		(s1, s2) = ('', ' '*size_diff)
	else:
		(s1, s2) = (' '*size_diff, '')
	print('SVS ('+out_str+'):'+s1, n_filt)
	print('SVS (other):'+s2, n_not_filt)

#
#
#
def split_repetitive_svs(fn_in, fn_out_repeats, fn_out_not_repeats, rep_dict, rep_track):
	f  = open(fn_in, 'r')
	f2 = open(fn_out_repeats, 'w')
	f3 = open(fn_out_not_repeats, 'w')
	n_reps    = 0
	n_notreps = 0
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
				col_form = splt.index('FORMAT')
				col_samp = len(splt) - 1
			f2.write(line)
			f3.write(line)
		else:
			splt = line.strip().split('\t')
			(my_type, my_len, my_end) = parse_info_column_of_sv(splt[col_info], int(splt[1]))
			if my_type == None:
				continue
			if my_type in ['BND', 'INV']:	# assume all svs of these types are not repetitive
				f3.write(line)
				n_notreps += 1
				continue
			elif my_end == None:
				print('what kind of weird SV are you?')
				print(splt)
				exit(1)
			repeat_hit = get_repeat_hits(splt[0], int(splt[1]), my_end, rep_dict, rep_track)
			if len(repeat_hit):
				f2.write(line)
				n_reps += 1
			else:
				f3.write(line)
				n_notreps += 1
	f3.close()
	f2.close()
	f.close()
	print('SVS-repeats:', n_reps)
	print('SVS-other:  ', n_notreps)

#
#
#
def intersect_vcf_with_bed(fn_in, bed_dict, within_dist=10):
	f = open(fn_in, 'r')
	for line in f:
		if line[0] == '#':
			if line[1] != '#':
				splt = line[1:].strip().split('\t')
				col_info = splt.index('INFO')
				col_form = splt.index('FORMAT')
				col_samp = len(splt) - 1
		else:
			splt   = line.strip().split('\t')
			my_chr = splt[0]
			my_pos = int(splt[1])
			(my_type, my_len, my_end) = parse_info_column_of_sv(splt[col_info], my_pos)
			#
			my_hits = query_bed_regions(bed_dict, my_chr, my_pos, q_end=my_end, max_dist=within_dist)
			#
			if len(my_hits):
				print(my_chr, my_pos, my_end, my_len, my_type, my_hits)
	f.close()

