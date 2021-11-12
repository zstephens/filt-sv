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

def query_bed_regions(bed_dict, q_chr, q_pos, max_dist=10000):
	dat_out = []
	if q_chr in bed_dict:
		for br in bed_dict[q_chr]:
			if br[0] > q_pos+max_dist:						# we've gone too far, stop looking
				break
			if q_pos >= br[0]-max_dist and q_pos < br[0]:	# we're close enough, on the right
				dat_out.append((br[2], br[0]-q_pos))
			elif q_pos >= br[0] and q_pos <= br[1]:			# we've got a hit inside a bed region
				dat_out.append((br[2], 0))
			elif q_pos <= br[1]+max_dist and q_pos > br[1]:	# we're close enough, on the left
				dat_out.append((br[2], q_pos-br[1]))
	return dat_out
