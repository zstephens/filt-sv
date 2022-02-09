import gzip

REP_BUFF = 5

#
# simple repeats
#
def get_repetitive_regions(filename):
	repeat_pos_dict = {}
	f = gzip.open(filename,'r')
	for line in f:
		splt = line.decode().strip().split('\t')
		[p1,p2] = sorted([int(splt[1]),int(splt[2])])
		if splt[0] not in repeat_pos_dict:
			repeat_pos_dict[splt[0]] = []
		repeat_pos_dict[splt[0]].append((p1-REP_BUFF,p2+REP_BUFF,splt[3]))
	f.close()
	return repeat_pos_dict

#
# other bed stuff
#
def parse_bed_regions(fn, anno_col=3):
	dat_out = {}
	f = gzip.open(fn, 'r')
	for line in f:
		splt = line.decode().strip().split('\t')
		if splt[0] not in dat_out:
			dat_out[splt[0]] = []
		if len(splt) <= anno_col:
			dat_out[splt[0]].append((int(splt[1]), int(splt[2])))
		else:
			dat_out[splt[0]].append((int(splt[1]), int(splt[2]), splt[anno_col]))
	f.close()
	return dat_out

#
# bed intersect - work harder, not smarter!
#
def query_bed_regions(bed_dict, q_chr, q_pos, q_end=None, max_dist=10000):
	if q_end == None:
		q_end = q_pos
	[q_pos, q_end] = sorted([q_pos, q_end])
	#
	dat_out = []
	if q_chr in bed_dict:
		for br in bed_dict[q_chr]:
			qr = [q_pos-max_dist, q_end+max_dist]
			if br[0] > qr[1]:						# we've gone too far, stop looking
				break
			if br[0] <= qr[1] and qr[0] <= br[1]:	# ranges overlap
				my_dist = 0
				if q_end < br[0]:
					my_dist = br[0]-q_end
				elif q_pos > br[1]:
					my_dist = q_pos-br[1]
				dat_out.append((br[2], my_dist))
	dat_out = [m[1] for m in sorted(list(set([(n[1], n) for n in dat_out])))]
	return dat_out
