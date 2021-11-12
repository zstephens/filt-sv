import os

EXE_DICT = {'bcftools': 'bcftools',
            'bgzip':    'bgzip',
            'tabix':    'tabix',
            'truvari':  'truvari'}

def exists_and_is_nonzero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def rm(fn):
	if os.path.isdir(fn):
		os.system('rm -rf '+fn)
	elif os.path.isfile(fn):
		os.system('rm '+fn)

def mv(src, dest):
	if os.path.isdir(src) or os.path.isfile(src):
		os.system('mv '+src+' '+dest)

# label collapsing one-liners, taken from: https://github.com/PacificBiosciences/sv-benchmark
def collapse_pbsv(fn_in, fn_out, exe_dict=EXE_DICT):
	#
	(BGZIP, TABIX) = (exe_dict['bgzip'], exe_dict['tabix'])
	#
	cmd = BGZIP + ' -c ' + fn_in + ' > ' + fn_out
	os.system(cmd)
	cmd = TABIX + ' ' + fn_out
	os.system(cmd)

def collapse_sniffles(fn_in, fn_out, exe_dict=EXE_DICT):
	#
	(BGZIP, TABIX) = (exe_dict['bgzip'], exe_dict['tabix'])
	#
	cmd  = 'cat ' + fn_in + ' | grep "^#" > ' + fn_out+'.temp_header'
	cmd += '; cat ' + fn_in + ' | grep -vE "^#" | grep \'DUP\|INS\|DEL\' | sed \'s/DUP/INS/g\' | sort -k1,1 -k2,2g > ' + fn_out+'.temp_data'
	os.system(cmd)
	cmd = 'cat ' + fn_out+'.temp_header' + ' ' + fn_out+'.temp_data' + ' | ' + BGZIP + ' -c > ' + fn_out
	os.system(cmd)
	cmd = TABIX + ' ' + fn_out
	os.system(cmd)
	cmd = 'rm ' + fn_out+'.temp_header' + '; rm ' + fn_out+'.temp_data'
	os.system(cmd)

def collapse_svim(fn_in, fn_out, lenFilter=True, exe_dict=EXE_DICT):
	#
	(BGZIP, TABIX) = (exe_dict['bgzip'], exe_dict['tabix'])
	#
	if lenFilter:
		awk_cond = 'if(($5=="<DEL>" || $5=="<INS>") && $6>40)'
	else:
		awk_cond = 'if($5=="<DEL>" || $5=="<INS>")'
	cmd = 'cat ' + fn_in + ' | sed \'s/INS:NOVEL/INS/g\' | sed \'s/DUP:INT/INS/g\' | sed \'s/DUP:TANDEM/INS/g\' | awk \'{ if($1 ~ /^#/) { print $0 } else { ' + awk_cond + ' { print $0 } } }\' > ' + fn_out+'.temp'
	os.system(cmd)
	cmd = BGZIP + ' -c ' + fn_out+'.temp' + ' > ' + fn_out
	os.system(cmd)
	cmd = TABIX + ' ' + fn_out
	os.system(cmd)
	cmd = 'rm ' + fn_out+'.temp'
	os.system(cmd)

#
# oldschool truvari command (v1.0)
#
def truvari_old(vcf_base, vcf_call, out_dir, ref_seq, out_log=None, incl_bed='', exe_dict=EXE_DICT):
	#
	(TRUVARI) = (exe_dict['truvari'])
	#
	if out_dir[-1] == '/':
		out_dir = out_dir[:-1]
	cmd = TRUVARI + ' -f ' + ref_seq + ' -b ' + vcf_base + ' -c ' + vcf_call + ' -o ' + out_dir + ' -r 1000 -p 0.00 --passonly' + incl_bed
	if out_log != None:
		cmd += ' > ' + out_log + ' 2>&1'
	os.system(cmd)

#
# updated truvari command (v3.0)
#
def truvari(vcf_base, vcf_call, out_dir, ref_seq, out_log=None, incl_bed='', exe_dict=EXE_DICT):
	#
	(TRUVARI) = (exe_dict['truvari'])
	#
	if out_dir[-1] == '/':
		out_dir = out_dir[:-1]
	cmd  = TRUVARI + ' bench -f ' + ref_seq + ' -b ' + vcf_base + ' -c ' + vcf_call + ' -o ' + out_dir
	cmd += ' --passonly --multimatch --no-ref a -r 1000 -p 0.0' + incl_bed
	if out_log != None:
		cmd += ' > ' + out_log + ' 2>&1'
	os.system(cmd)

#
#
#
def bgzip_tabix(fn_in, fn_out, exe_dict=EXE_DICT):
	#
	(BGZIP, TABIX) = (exe_dict['bgzip'], exe_dict['tabix'])
	#
	cmd = BGZIP + ' -c ' + fn_in + ' > ' + fn_out
	os.system(cmd)
	cmd = TABIX + ' ' + fn_out
	os.system(cmd)
