import gzip

SV_TYPE_MAP = {'deletion':                 ['DEL'],
               'herv deletion':            ['DEL'],
               'alu deletion':             ['DEL'],
               'line1 deletion':           ['DEL'],
               'sva deletion':             ['DEL'],
               'duplication':              ['DUP'],
               'copy number loss':         ['DEL'],
               'copy number gain':         ['DUP'],
               'copy number variation':    ['DUP','DEL'],
               'insertion':                ['INS'],
               'line1 insertion':          ['INS'],
               'alu insertion':            ['INS'],
               'sva insertion':            ['INS'],
               'mobile element insertion': ['INS'],
               'DUP': ['DUP'],
               'INS': ['INS'],
               'DEL': ['DEL'],
               'INV': ['INV'],
               'BND': ['BND','TRA'],
               'TRA': ['BND','TRA']}

#
# common SVs
#
# https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd186/
#
def get_common_svs(filename, min_control_freq=-1., collapse_dup_and_ins=True):
	common_svs = {}
	n_svs_read = 0
	f = gzip.open(filename, 'r')
	for line in f:
		splt = line.decode().strip().split('\t')
		if splt[0] not in common_svs:
			common_svs[splt[0]] = []
		my_freq = None
		if len(splt) >= 5:
			if splt[4][:13] == 'CONTROL_FREQ=':
				my_freq = float(splt[4][13:])
		if my_freq == None or my_freq >= min_control_freq:
			n_svs_read += 1
			common_svs[splt[0]].append([int(splt[1]), int(splt[2]), splt[3]])
	f.close()
	#
	if collapse_dup_and_ins:
		for k in SV_TYPE_MAP.keys():
			if 'DUP' in SV_TYPE_MAP[k] and 'INS' not in SV_TYPE_MAP[k]:
				SV_TYPE_MAP[k].append('INS')
			if 'INS' in SV_TYPE_MAP[k] and 'DUP' not in SV_TYPE_MAP[k]:
				SV_TYPE_MAP[k].append('DUP')
	my_fn = filename.split('/')[-1]
	freq_str = ''
	if min_control_freq >= 0.:
		freq_str = ' (freq >= {0:0.3f})'.format(min_control_freq)
	print('read ' + str(n_svs_read) + ' svs from ' + my_fn + freq_str)
	return common_svs

#
#
#
def sv_is_in_blacklist(my_chr, pos_start, pos_end, my_type, blacklist_svs, max_endpoint_dist=50):
	if my_chr in blacklist_svs:
		c2 = sorted([pos_start, pos_end])
		for sv in blacklist_svs[my_chr]:
			if my_type in SV_TYPE_MAP[sv[2]]:	# same type
				c1 = sorted([sv[0], sv[1]])
				d1 = abs(c1[0]-c2[0])
				d2 = abs(c1[1]-c2[1])
				if d1 <= max_endpoint_dist and d2 <= max_endpoint_dist:	# beakpoints are same
					return True
	return False

#
#
#
def split_common_svs(fn_in, fn_out_common, fn_out_not_common, blacklist):
	f  = open(fn_in, 'r')
	f2 = open(fn_out_common, 'w')
	f3 = open(fn_out_not_common, 'w')
	n_common     = 0
	n_not_common = 0
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
			splt  = line.strip().split('\t')
			splt4 = splt[col_info].split(';')
			my_chr   = splt[0]
			my_start = int(splt[1])
			my_type  = None
			my_len   = None
			my_end   = None
			for n in splt4:
				if n[:7] == 'SVTYPE=':
					my_type = n[7:]
				if n[:6] == 'SVLEN=':
					my_len = int(n[6:])
				if n[:4] == 'END=':
					my_end = int(n[4:])
			if my_type == None:
				#print 'SVTYPE not found in info? skipping...'
				continue
			if my_type == 'INS':
				my_end = int(splt[1])
			elif my_type == 'DEL' and my_len < 0:
				my_end = int(splt[1]) - my_len
			elif my_type == 'DUP' and my_len > 0:
				my_end = int(splt[1]) + my_len
			elif my_type in ['BND', 'INV']:	# assume all svs of these types are not in the blacklist
				f3.write(line)
				n_not_common += 1
				continue
			elif my_end == None:
				print('what kind of weird SV are you?')
				print(splt)
				exit(1)
			if sv_is_in_blacklist(my_chr, my_start, my_end, my_type, blacklist):
				f2.write(line)
				n_common += 1
			else:
				f3.write(line)
				n_not_common += 1
	f3.close()
	f2.close()
	f.close()
	print('SVS-pass:', n_not_common)
	print('SVS-filt:', n_common)
