import random
import numpy
import math
import pylab as plt
import scipy.stats
import datetime as dtm
import pytz

plt.ion()

def plot_rb_cdfses(rblens=[8, 32, 128, 256, 512, 1024, 2048, 4096, 2*4096], nits=1000, fnum=0):
	do_clf=True
	#
	for i, rblen in enumerate(rblens):
		a=plot_rb_cdfs(rblen=rblen, nits=nits, fnum=fnum, do_clf=do_clf)
		do_clf=False
	#
	return None

def plot_rb_cdfs(cdf_in=None, rblen=100, nits=10000, fnum=0, do_clf=True):
	# ... not properly modularized, but CDFs are easy enough that maybe it's ok.
	if cdf_in==None:
		cdf_in = random_rb_sequence(rblen=rblen, nits=nits)
	#
	log_ratios_gt = [x for x in cdf_in['ratio_lognorm'] if x>0]
	log_ratios_gt.sort()
	log_ratios_lt = [-x for x in cdf_in['ratio_lognorm'] if x<0]
	log_ratios_lt.sort()
	log_ratios_eq = [x for x in cdf_in['ratio_lognorm'] if x==0]
	log_ratios_eq.sort()
	#[x.sort for x in (log_ratios_gt, log_ratios_lt, log_ratios_eq)]
	#Ns=range(1,len(log_ratios)+1)
	#Ps = [float(i+1)/float(nits) for i in xrange(nits)]
	#print "mode: ", scipy.stats.mode(log_ratios)
	for X in (log_ratios_gt, log_ratios_lt, log_ratios_eq):
		print "rblen: %d. max: %f, min: %f, mean: %f, stdev: %f, median: %f/%f, mode: %f" % (rblen, max(X), min(X), numpy.mean(X), numpy.std(X), numpy.median(X), X[int(len(X)/2)], scipy.stats.mode(X)[0])
	#
	plt.figure(fnum)
	plt.ion()
	if do_clf:
		plt.clf()
		plt.legend(loc=0, numpoints=1)
		plt.xlabel('normalized log_ratios')
		plt.ylabel('probability P(l_r)')
		plt.gca().set_yscale('linear')
		plt.gca().set_xscale('log')
	#
	plt.plot(log_ratios_gt, [float(x)/float(nits) for x in xrange(1,len(log_ratios_gt)+1)], '.-', label='rblen_gt=%d' % rblen)
	plt.plot(log_ratios_lt, [float(x)/float(nits) for x in xrange(1,len(log_ratios_lt)+1)], '.-', label='rblen_lt=%d' % rblen)
	plt.plot(log_ratios_eq, [float(x)/float(nits) for x in xrange(1,len(log_ratios_eq)+1)], '.-', label='rblen_eq=%d' % rblen)
	
	#
	return cdf_in
#
def nrb_sequence(sequence_in=None, rb_len=10, seq_len=10000):
	if sequence_in==None:
		R=random.Random()
		sequence_in = [R.random() for j in xrange(seq_len)]
	#
	dtype_names = ['i', 'n_gt', 'n_lt', 'ratio', 'ratio_lognorm']
	rb_vals=[]
	log_N = math.log10(rb_len)
	for i in xrange(rb_len, len(sequence_in)+1):
		rb_gt=[sequence_in[i-rb_len]]
		rb_lt=[sequence_in[i-rb_len]]
		#
		[rb_gt.append(x) for x in sequence_in[i-rb_len:i] if x>rb_gt[-1]]
		[rb_lt.append(x) for x in sequence_in[i-rb_len:i] if x<rb_lt[-1]]
		#
		rb_vals += [[i-1, len(rb_gt), len(rb_lt)]]
		rb_vals[-1]+=[float(rb_vals[-1][-2])/float(rb_vals[-1][-1])]		# record-breaking ratio
		rb_vals[-1]+=[math.log10(rb_vals[-1][-1])/log_N]					# "normalized" record-breaking ratio
	rb_vals = numpy.core.records.fromarrays(zip(*rb_vals), names=dtype_names, formats=[type(x).__name__ for x in rb_vals[0]])
	return rb_vals	
		
#
#def rb_cdfs(rblen=100, nits=10000):
def random_rb_sequence_2(rblen=100, nits=10000):
	# stats on record-breaking ratios of a random sequence.
	# so what this bit does now (after being chopped up a bit) is to produce a sequence of 
	# rb_intervals from nits sequences of length rblen.
	#
	R=random.Random()
	#
	rb_vals = []
	dtype_names = ['i', 'n_gt', 'n_lt', 'ratio', 'ratio_lognorm']
	log_N = math.log10(rblen)
	#
	for i in xrange(nits):
		#vals = [R.random() for j in xrange(nits)]
		vals = [R.random() for j in xrange(rblen)]
		rb_gt = [vals[0]]
		rb_lt = [vals[0]]
		#
		# get greater/lesser record-breaking sub-sequences.
		[rb_gt.append(x) for x in vals[1:] if x>rb_gt[-1]]
		[rb_lt.append(x) for x in vals[1:] if x<rb_lt[-1]] 
		#
		rb_vals += [[i, len(rb_gt), len(rb_lt)]]
		rb_vals[-1]+=[float(rb_vals[-1][-2])/float(rb_vals[-1][-1])]
		rb_vals[-1]+=[math.log10(rb_vals[-1][-1])/log_N]
	#
	rb_vals = numpy.core.records.fromarrays(zip(*rb_vals), names=dtype_names, formats=[type(x).__name__ for x in rb_vals[0]])
	return rb_vals
#
def random_rb_sequence(rblen=100, seq_len=10000):
	# stats on record-breaking ratios of a random sequence.
	# so what this bit does now (after being chopped up a bit) is to produce a sequence of 
	# rb_intervals from nits sequences of length rblen.
	#
	R=random.Random()
	#
	rb_vals = []
	dtype_names = ['i', 'n_gt', 'n_lt', 'ratio', 'ratio_lognorm']
	log_N = math.log10(rblen)
	#
	sequence = [R.random() for j in xrange(seq_len)]
	#
	#for i in xrange(rblen-1, len(sequence)):
	for i, x in enumerate(sequence[rblen-1:]):
		j=i+rblen
		#
		vals = sequence[i:j]
		#
		rb_gt = [vals[0]]
		rb_lt = [vals[0]]
		#
		# get greater/lesser record-breaking sub-sequences.
		[rb_gt.append(x) for x in vals[1:] if x>rb_gt[-1]]
		[rb_lt.append(x) for x in vals[1:] if x<rb_lt[-1]] 
		#
		rb_vals += [[i, len(rb_gt), len(rb_lt)]]
		rb_vals[-1]+=[float(rb_vals[-1][-2])/float(rb_vals[-1][-1])]
		rb_vals[-1]+=[math.log10(rb_vals[-1][-1])/log_N]
	#
	rb_vals = numpy.core.records.fromarrays(zip(*rb_vals), names=dtype_names, formats=[type(x).__name__ for x in rb_vals[0]])
	#return rb_vals, sequence
	return rb_vals
#
def rb_ratio(rb_seq, log_norm=True):
	#print rb_seq
	rb_gt = [rb_seq[0]]
	rb_lt = [rb_seq[0]]
	#
	[rb_gt.append(x) for x in rb_seq if x>rb_gt[-1]]
	[rb_lt.append(x) for x in rb_seq if x<rb_lt[-1]]
	#print "**", rb_gt
	#
	if log_norm==False or len(rb_seq)==1:
		return float(len(rb_gt))/float(len(rb_lt))
	if log_norm:
		# (see now the lognorm() function).
		return math.log10(float(len(rb_gt))/float(len(rb_lt)))/(math.log10(float(len(rb_seq))))
#
def rb_runs(rb_ratios=None, rblen=1, seq_len=100000, log_norm=True):
	# churn up some statistics on rb "runs", aka sequences of 1,0 values (in this case,
	# well use 1 = (ratio>1), -1 = (ratio<1), and 0: ratio=1 ???
	#
	# rb_ratios: a sequence of rb_ratios. note these are effectively in "natural time". if we want to 
	# eximne temporal runs, we'll have to do some integrating.
	#
	rblen = int(rblen)
	seq_len = int(seq_len)
	#
	#if sequence==None:
	if rb_ratios == None:
		R=random.Random()
		sequence = [R.random() for x in xrange(seq_len)]
		#
		rb_ratios = [rb_ratio(sequence[i-rblen:i], log_norm=log_norm) for i in xrange(rblen, seq_len)]
	#runs_gt = []
	#runs_lt = []
	#runs_0  = []
	runs = {-1:[], 1:[], 0:[]}	# parity lists
	#
	parity = tri_parity(rb_ratios[0])
	this_list = [rb_ratios[0]]	
	#while j<seq_len:
	for j, x in enumerate(rb_ratios[1:]):
		#pass
		if tri_parity(x)!=parity:
			# new run.
			if len(this_list)>1: runs[parity]+=[this_list]	# not interested in len(1) runs...
			parity=tri_parity(x)
			this_list=[]
		#
		#print x
		this_list+=[x]
	#
	#for key, rw in runs.iteritems():
	#	print key, len(rw), " :: ", float(len(rw))/float(seq_len)
	#
	return runs
#
def mean_entropy_random(seq_len=10000):
	rb_lens = [4,256]
	#
	mean_H = []
	#
	for rb_len in xrange(rb_lens[0], rb_lens[1]):
		H = shannon_entropy(zip(*random_rb_sequence(rblen=rb_len, seq_len=seq_len))[4])
		mean_H += [[rb_len, numpy.mean(H)]]
		print "<H>: %s" % mean_H[-1]
	#
	plt.figure(0)
	plt.clf()
	plt.plot(zip(*mean_H)[0], zip(*mean_H)[1], '.-')
	#
	return mean_H
		
#
def mean_entropies(seq_len=100000):
	# get a random sequence and random shannon entropy:
	r1, h1 = shannon_entropy(seq_len=1000000)
	rb_random = random_rb_sequence(rblen=220, seq_len=seq_len)
	#
	tohoku = numpy.load('data/tohoku_rb_sequence.pkl')
	tohoku_ratios=zip(*tohoku)[5]
	h_tohoku = shannon_entropy(data_in=tohoku_ratios)
	mean_h_tohoku=numpy.mean(h_tohoku)
	mean_h_tohoku_raw = numpy.mean(shannon_entropy(zip(*tohoku)[4]))
	print numpy.mean(shannon_entropy(zip(*rb_random)[4])), mean_h_tohoku_raw, mean_h_tohoku
	#
	emc = numpy.load('data/emc_rb_sequence.pkl')
	emc_ratios = zip(*emc)[4]
	h_emc = shannon_entropy(data_in=emc_ratios)
	mean_h_emc = numpy.mean(h_emc)
	rb_random_emc = random_rb_sequence(rblen=500, seq_len=seq_len)
	print numpy.mean(shannon_entropy(zip(*rb_random_emc)[4])), mean_h_emc
	#
	pf = numpy.load('data/parkfield-elip-rbsequence.pkl')
	pf_ratios = zip(*pf)[4]
	h_pf = shannon_entropy(data_in=pf_ratios)
	mean_h_pf = numpy.mean(h_pf)
	rb_random_pf = random_rb_sequence(rblen=290, seq_len=seq_len)
	print numpy.mean(shannon_entropy(zip(*rb_random_pf)[4])), mean_h_pf
	#
	chichi=numpy.load('data/chichi-rb_sequence.pkl')
	chichi_ratios = zip(*chichi)[4]
	h_chichi = shannon_entropy(data_in=chichi_ratios)
	mean_h_chichi = numpy.mean(h_chichi)
	rb_random_chichi=random_rb_sequence(rblen=600, seq_len=seq_len)
	print numpy.mean(shannon_entropy(zip(*rb_random_chichi)[4])), mean_h_chichi
#
def permutation_entropies(p_lens=[5,256], seq_len=100000):
	tohoku = numpy.load('data/tohoku_rb_sequence.pkl')
	tohoku_ratios=zip(*tohoku)[4]
	
	emc = numpy.load('data/emc_rb_sequence.pkl')
	emc_ratios = zip(*emc)[4]
	
	pf = numpy.load('data/parkfield-elip-rbsequence.pkl')
	pf_ratios = zip(*pf)[4]
	
	chichi=numpy.load('data/chichi-rb_sequence.pkl')
	chichi_ratios = zip(*chichi)[4]
	
	ratioses = {'random_%d' % x:random_rb_sequence(rblen=x, seq_len=seq_len) for x in [16, 64, 256, 1024]}
	ratioses['tohoku']=tohoku
	ratioses['emc']=emc
	ratioses['pf']=pf
	ratioses['chichi']=chichi
	#
	plt.figure(0)
	plt.clf()
	#
	for lbl, ratios in ratioses.iteritems():
		print "lbl: %s" % lbl
		Y = [numpy.mean(permutation_entropy(data_in=ratios, p_len=p_len)[p_len:]) for p_len in xrange(p_lens[0], p_lens[1])]
		plt.plot(xrange(p_lens[0], p_lens[1]), Y, '.-', lw=2, label='%s'%lbl)
		plt.draw()
		plt.show()
		#
	#
	plt.legend(loc=0, numpoints=1)
		
	
		
#
def shannon_entropy(data_in=None, seq_len = 10000, f_mod=abs, log_base=2.):
	# shannon entropy, H = v-p log(p). note that p requires definition. we'll probably generalize these functions later;
	# for now, let's use a cumulative probability, and for shannon entropy, based on the absolute value of the input.
	# note that in questions of entropy and information, the definition of p is non-trivial.
	#
	if f_mod==None: f_mod = lambda x:x
	#
	if data_in==None: 
		R=random.Random()
		#
		data_in = [R.random() for j in xrange(seq_len)]
	#
	#prob_set = data_in[:]
	prob_set = [f_mod(x) for x in data_in]
	#
	prob_set.sort()
	prob_index = {val:float(i+1.0)/float(len(data_in)) for i, val in enumerate(prob_set)}
	#
	H = [-prob_index[f_mod(x)]*math.log(prob_index[f_mod(x)], log_base) for x in data_in]
	#
	return data_in, H
#
def permutation_entropy(data_in=None, p_len=10, seq_len=10000, f_mod=None, log_base=2.):
	'''
	# permutation entropy:
	# basically, define a sequence S from X_input;
	# S = {s_i} = {x_i<x_(i-1)<x_(i-2)...<x(i-p_len)}
	# then, calc. the probability of each permutation, p(s_i) and calc. p log(p) from that.
	# use f_mod, a function, to modify the data value. example: for some record-breaking ratio sequence analyses,
	# we might want to use abs( log(n_gt/n_lt) ), f_mod() = abs(). for no mod, pass None,
	# and w'll construct a NULL function f(x) = x... of course, for this application, we want the raw values, but see
	# the shannon_entropy, etc. functions. 
	'''
	#
	if f_mod!=None:
		data_in = [f_mod(x) for x in data_in]
	# ... and just in case we use it later:
	if f_mod==None: f_mod = lambda x:x
	#
	if data_in==None: 
		R=random.Random()
		#
		data_in = [R.random() for j in xrange(seq_len)]
	#
	# first, get the permutation set S. we'll also make a permutation dict.
	#permutation_set:
	S=[]
	S_index={}
	n_norm = float(len(data_in))
	#
	for j in xrange(len(data_in)-p_len):
		S+=[tuple([data_in[i+1]>x for i,x in enumerate(data_in[j:j+p_len])])]
		#
		#s_prime = zip(*[data_in[j:j + p_len], range(p_len)])
		#s_prime.sort(key = lambda x:x[0])
		#s_prime = [x for (y,x) in sorted(zip([data_in[j:j + p_len], xrange(p_len)))]
		#S+= [zip(*s_prime)[1]]	# keep the reordered range(), aka, the original positions.
		#
		# in one line:
		#S+= [tuple([x for (y,x) in sorted(zip(*[data_in[j:j + p_len], xrange(p_len)]))])]
		#
		# now, add to index:
		#print "s[-1]: ", S[-1]
		if not S_index.has_key(S[-1]): S_index[S[-1]]=0.0
		#
		S_index[S[-1]]+=1.0/n_norm
	#
	# now we have our permutation/proxy data S and our probability index S_index. lets calculate permutation entropy. assume
	# the sequence is in temporal order with t(x_i)<t(x_{i+1}).
	#
	# return [X, H(x)] ?? do we want to return the dict?
	#return S
	#
	H=[-S_index[s]*math.log(S_index[s], log_base) for s in S]
	H = [None for x in xrange(p_len)] + H
	#
	#return [data_in, H], S_index
	return H

def plot_permutation_entropy(data_file=None, p_len=10, rand_len=10**5):
	# for now, keep it simple and plot from a data file.
	#
	if data_file == None:
		this_N = rand_len
		R=random.Random()
		rb_ratios_raw = [10.**(2.0*R.random()-1.) for i in xrange(this_N)]
		t = range(this_N)
	else:
		#
		rb_data = numpy.load(data_file)
		z_data = zip(*rb_data)
		rb_ratios_raw = z_data[4]
		t= z_data[1]
		#
	H= permutation_entropy(data_in=rb_ratios_raw, p_len=p_len, log_base=2.)
	#print len(H[0]), len(H[1]), len(d), len(t), len(rb_ratios_raw)
	
	#return H
	
	#return H,d
	#
	plt.figure(0)
	plt.ion()
	plt.clf()
	plt.plot(t, H, '.-')
	plt.xlabel('t')
	plt.ylabel('H')
	#
	plt.figure(1)
	plt.clf()
	R_prime, H_prime = zip(*[[r,h] for (r,h) in sorted(zip(rb_ratios_raw,H))])
	plt.plot(numpy.log10(R_prime), H_prime, '.')
	#plt.plot([math.log10(x) for x in rb_ratios_raw], H, '.-')
	#plt.plot(rb_ratios_raw, H, '.')
	plt.xlabel('rb ratio')
	plt.ylabel('H')
	#
	plt.figure(2)
	plt.clf()
	plt.plot(t, numpy.log10(rb_ratios_raw), '.-')
	plt.xlabel('t')
	plt.ylabel('rb_ratios')
	#
	plt.figure(3)
	plt.clf()
	plt.plot(t[p_len:], [math.log10(rb_ratios_raw[i])*H[i] if rb_ratios_raw[i]<1. else 0. for i in xrange(p_len, len(H))], '.-')
	plt.plot([t[0], t[-1]], [0., 0.], 'k-')
	plt.plot(dtm.datetime(2011, 3, 11, 14-9, 46, 24, tzinfo=pytz.timezone('UTC')), 0., 'r*', ms=18)
	return [t, rb_ratios_raw, H]
		

def rb_runs_report(rb_runs_data=None, rblen=128, seq_len=100000, log_norm=True, fignum=0, doplots=True, do_clf=True, cat_name='(test catalog)'):
	# some stats on record-breaking runs.
	if rb_runs_data==None:
		# in default mode, this will produce a record-breaking sequence from a random catalog.
		rb_runs_data = rb_runs(rblen=rblen, seq_len=seq_len, log_norm=log_norm)
	if doplots:
		plt.figure(fignum)
		plt.ion()
		if do_clf: plt.clf()
	#
	# number of runs for gt,lt,0.
	# stats on length of runs.
	run_stats = {key:{} for key in rb_runs_data.iterkeys()}		# makes this dict:  {-1: {}, 0: {}, 1: {}}; we'll populate the keys next.
	total_N_runs = sum([len(x) for x in rb_runs_data.itervalues()])
	#
	for key, datas in rb_runs_data.iteritems():
		# number of runs:
		run_stats[key]['n_runs'] = len(datas)
		#
		run_stats[key]['run_lengths'] = [len(x) for x in datas]
		#run_stats[key]['run_intervals'] = [
		#
		# mean run length:
		run_stats[key]['mean_run_len'] = numpy.mean(run_stats[key]['run_lengths'])
		run_stats[key]['std_run_len']  = numpy.std(run_stats[key]['run_lengths'])
		run_stats[key]['med_run_len']  = numpy.median(run_stats[key]['run_lengths'])
		#
		# mode calc seems to be delicate, like a flower. let's trap and fail gracefully:
		try:
			run_stats[key]['mode_run_len'] = scipy.stats.mode(run_stats[key]['run_lengths'])
		except:
			print "scipy.stats.mode() calc failed. spoof it."
			run_stats[key]['mode_run_len'] = [None, None]
		#
		these_lengths = [x for x in run_stats[key]['run_lengths']]
		norm_run_lengths = [x/float(rblen) for x in these_lengths]
		#
		#Ns = [x/float(total_N_runs) for x in range(1, len(these_lengths)+1)]
		Ns = [x/float(total_N_runs) for x in range(1, len(norm_run_lengths)+1)]
		#these_lengths.sort()
		norm_run_lengths.sort()
		if doplots:
			#plt.gca().set_yscale('log')
			#plt.gca().set_xscale('log')
			plt.gca().set_yscale('linear')
			plt.gca().set_xscale('linear')
			#
			# normalize run lengths by rb_len (basically by GR scaling, aka,. relative to the (expected)
			# seismic sequence length
			#plt.plot(these_lengths, Ns, '.-', label='parity: %d' % key)
			#
			#plt.plot([x/float(rblen) for x in these_lengths], Ns, '.-', label='parity: %d' % key)
			plt.plot(norm_run_lengths, Ns, '.-', label='parity: %d' % key)
		#
	if doplots:
		plt.legend(loc=0, numpoints=1)
		plt.ylabel('Probability $P(\\leq N_{run})$')
		plt.xlabel('Normalized run-length $\Lambda = N_{run}/N_{rb}$')
		plt.title('record-breaking runs report for %s' % cat_name)
	#
	print "summary report:"
	for key, datas in run_stats.iteritems():
		print key, {ky:val for ky,val in datas.iteritems() if ky!='run_lengths'}
	#
	return run_stats
#
# instead of writing a wrapper function, just defind the prams dictionary:
parkfield_rb_report_prams = {'data_file':'data/parkfield-elip-rbsequence.pkl', 'rb_len':310, 'fnum':0, 'ave_len':None, 'cat_name':'Parkfield'}
# alternatively (for parkfield):
#parkfield_rb_report_prams = {'data_file':'data/parkfield-elip-rbsequence.pkl', 'targ_mag':5.96, 'm_c':1.5, 'm_t':7.6, 'fnum':0, 'random_len':100000, 'ave_len':None, 'cat_name':'Parkfield'}

emc_rb_report_prams = {'data_file':'data/emc_rb_sequence.pkl', 'rb_len':500, 'fnum':0, 'ave_len':None, 'cat_name':'EMC'}
tohoku_rb_report_prams = {'data_file':'data/tohoku_rb_sequence.pkl', 'rb_len':220, 'fnum':0, 'ave_len':None, 'cat_name':'Tohoku'}
chi_chi_report_prams = {'data_file':'data/chichi-rb_sequence.pkl', 'rb_len':630, 'fnum':0, 'ave_len':None, 'cat_name':'Chi-chi'}
#
def chi_chi_rb_report(data_file='data/chichi-rb_sequence.pkl', rb_len=630, fnum=0, random_len=10000, ave_len=None):
	return rb_stats_report(data_file=data_file, rb_len=rb_len, fnum=fnum, random_len=random_len)

def tohoku_rb_report(data_file='data/tohoku_rb_sequence.pkl', rb_len=220, fnum=0, random_len=10000, ave_len=None):
	return rb_stats_report(data_file=data_file, rb_len=rb_len, fnum=fnum, random_len=random_len, ave_len=ave_len)
	#return rb_stats_report(*args, **kwargs)
#
def rb_stats_report(data_file='data/tohoku_rb_sequence.pkl', rb_len=220, targ_mag=None, m_t=7.6, m_c=None, fnum=0, random_len=10000, ave_len=None, cat_name='(tohoku)'):
	# we'll need to encode some of these data, but for now, we know that rb_len was 220 (probably).
	# (a good guess for this is the first event number rw_0,item_0=219).
	# ... so this generalized function was orignally written for Tohoku, so "tohoku" probably impies "test data",
	# as opposed to "random control data"
	#
	#
	random_len=int(random_len)
	if targ_mag!=None:
		# calculate the rb_len:
		# note: it may be necessary to specify m_c:
		# this gives a "rb-stats ready" rb_len, including an integer round/truncation so that averaging, etc. will be
		# over integer intervals.
		rb_len = getNsample(m=targ_mag, mc=m_c, mt=m_t, dm=1.0, b1=1.0, b2=1.5, dms0=1.0, doint=True, dmMode=1)
	#
	if ave_len==None:
		ave_len = max(int(rb_len/10), 1)
	if isinstance(ave_len, float):
		# this is a quick way to call fractioal ave_len values, but it could go awry if we make a mistake
		# and submit values like ave_len=1.0 when we mean 1. we could also require odd floats for this sort of work,
		# aka 1.0 --> 1 but 1.1 (or 1.0001) --> 1.0001*max(int(rb_len/10), 1)
		ave_len = int(ave_len*max(int(rb_len/10), 1))
	#
	print "rb_stats_report() for file: %s, rb_len:%d, random_len=%d" % (data_file, rb_len, random_len)
	#	
	rb_data = numpy.load(data_file)
	z_data = zip(*rb_data)
	rb_ratios_raw = z_data[4]
	#try:
	#	rb_ratios_maybe_averaged = z_data[5]		#... and of course, these should be converted to rec_arrays, or manually to dicts?
	#except:
	#	# maybe this hasn't bee created.
	#	rb_ratios_maybe_averaged = [numpy.mean(rb_ratios_raw[max(0, i-ave_len):i]) for i in xrange(1, len(rb_ratios_raw)+1)]
	# let's always do this (so we can vary the ave_len easily).
	rb_ratios_maybe_averaged = [numpy.mean(rb_ratios_raw[max(0, i-ave_len):i]) for i in xrange(1, len(rb_ratios_raw)+1)]
	#
	N_total = float(len(rb_data))		# total number of events...
	log_rb_len = math.log10(rb_len)
	tohoku_rb_gt = [math.log10(x)/log_rb_len for x in rb_ratios_raw if x>1.]
	tohoku_rb_lt = [-math.log10(x)/log_rb_len for x in rb_ratios_raw if x<1.]
	tohoku_rb_eq = [0 for x in rb_ratios_raw if x==1.]
	tohoku_rb_gt.sort()
	tohoku_rb_lt.sort()
	#
	tohoku_rb_gt_mean = [math.log10(x)/log_rb_len for x in rb_ratios_maybe_averaged if x>1.]
	tohoku_rb_lt_mean = [-math.log10(x)/log_rb_len for x in rb_ratios_maybe_averaged if x<1.]
	tohoku_rb_eq_mean = [0 for x in rb_ratios_maybe_averaged if x==1.]
	tohoku_rb_gt_mean.sort()
	tohoku_rb_lt_mean.sort()
	#
	plt.figure(fnum)
	plt.ion()
	plt.clf()
	plt.xlabel('normalized rb_ratio value, $\\frac{\\log (nrb_{gt}/nrb_{lt})}{\\log (N)}$')
	plt.ylabel('Probability $P(abs[r])$')
	#
	#Y = [x/float(len(tohoku_rb_gt)) for x in xrange(1, len(tohoku_rb_gt)+1)]	# prob assuming gt...
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_gt)+1)]						# prob within full sequence.
	plt.plot(tohoku_rb_gt, Y, '.-', label='%s_gt' % cat_name)
	#
	#Y = [x/float(len(tohoku_rb_lt)) for x in xrange(1, len(tohoku_rb_lt)+1)]
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_lt)+1)]
	plt.plot(tohoku_rb_lt, Y, '.-', label='%s_lt' % cat_name)
	#
	#Y = [x/float(len(tohoku_rb_eq)) for x in xrange(1, len(tohoku_rb_eq)+1)]
	#Y = [x/N_total for x in xrange(1, len(tohoku_rb_eq)+1)]
	#plt.plot(tohoku_rb_eq, Y, '.-', label='tohoku_eq')
	#
	# and the means...
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_gt_mean)+1)]						# prob within full sequence.
	plt.plot(tohoku_rb_gt_mean, Y, '.-', label='%s_gt_mean' % cat_name)
	#
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_lt_mean)+1)]
	plt.plot(tohoku_rb_lt_mean, Y, '.-', label='%s_lt_mean' % cat_name)
	#
	#Y = [x/N_total for x in xrange(1, len(tohoku_rb_eq_mean)+1)]
	#plt.plot(tohoku_rb_eq_mean, Y, '.-', label='tohoku_eq')
	#
	############
	# now, plot some random sequences for comparison:
	#rand_plots = plot_rb_cdfs(cdf_in=None, rblen=rb_len, nits=10000, fnum=0, do_clf=True)
	rand_cdf = random_rb_sequence(rblen=rb_len, seq_len=int(random_len))	# returns recarray with: [i, n_gt, n_lt, ratio, ratio_lognorm]
	rand_lognorm = rand_cdf['ratio_lognorm'].copy()
	rand_mean = [numpy.mean(rand_lognorm[max(0, i-ave_len):i]) for i in xrange(1, len(rand_lognorm)+1)]
	#	
	rand_ln_gt = [x for x in rand_lognorm if x>0.]
	rand_ln_lt = [-x for x in rand_lognorm if x<0.]
	rand_ln_eq = [x for x in rand_lognorm if x==0.]
	#
	rand_ln_gt.sort()
	rand_ln_lt.sort()
	#rand_ln_eq.sort()
	#
	rand_ln_gt_mean = [x for x in rand_mean if x>0.]
	rand_ln_lt_mean = [-x for x in rand_mean if x<0.]
	rand_ln_eq_mean = [x for x in rand_mean if x==0.]
	#
	#
	rand_ln_gt_mean.sort()
	rand_ln_lt_mean.sort()
	
	log_N_rand = math.log10(random_len)
	#
	plt.plot(rand_ln_gt, [float(x)/float(random_len) for x in xrange(1,len(rand_ln_gt)+1)], '.-', label='random_gt')
	plt.plot(rand_ln_lt, [x/float(random_len) for x in xrange(1,len(rand_ln_lt)+1)], '.-', label='random_lt')
	#plt.plot(rand_ln_eq, [x/float(random_len) for x in xrange(1,len(rand_ln_eq)+1)], '.-', label='random_eq')
	#
	plt.plot(rand_ln_gt_mean, [float(x)/float(random_len) for x in xrange(1,len(rand_ln_gt_mean)+1)], '.-', label='random_gt_mean')
	plt.plot(rand_ln_lt_mean, [x/float(random_len) for x in xrange(1,len(rand_ln_lt_mean)+1)], '.-', label='random_lt_mean')
	#plt.plot(rand_ln_eq_mean, [x/float(random_len) for x in xrange(1,len(rand_ln_eq_mean)+1)], '.-', label='random_eq_mean')
	#
	plt.legend(loc=0, numpoints=1)
	plt.savefig('output/%s_ratio_distributions.png' % cat_name)
	#
	#
	# "Record-breaking runs" report:
	print "now, generate a run-statistics report..."
	plt.figure(fnum+1)
	plt.clf()
	#
	tohoku_rb_sequence = [math.log10(x)/log_rb_len for x in rb_ratios_raw]
	tohoku_rb_sequence_mean = running_mean(tohoku_rb_sequence, ave_len)
	#tohoku_rb_sequence_mean = [math.log10(x)/log_rb_len for x in rb_ratios_maybe_averaged]
	
	#return [tohoku_rb_sequence, tohoku_rb_sequence_mean]
	
	# now, make a random sequence averaged as if a proper rb sequence:
	#
	random_sequence = nrb_sequence(rb_len=rb_len, seq_len=random_len)['ratio_lognorm']
	random_sequence_mean = [numpy.mean(random_sequence[max(0, i-ave_len):i]) for i in xrange(1,len(random_sequence)+1)]
	#
	#print "rand averaged len: ", len(random_averaged)
	#
	# rb_runs(rb_ratios=tohoku_rb_sequence, rblen=1, seq_len=100000, log_norm=True)
	tohoku_runs = rb_runs(rb_ratios=tohoku_rb_sequence)
	#random_runs = rb_runs(rb_ratios=None, rblen=rb_len, seq_len=random_len, log_norm=True)	# .. but let's use the same random_sequence...
	random_runs = rb_runs(rb_ratios = random_sequence)
	#
	tohoku_runs_mean = rb_runs(rb_ratios=tohoku_rb_sequence_mean)
	random_runs_mean = rb_runs(rb_ratios=random_sequence_mean)
	#
	##
	#random_run_stats = rb_runs_report(rb_runs_data=None, rblen=rb_len, seq_len=random_len, log_norm=True, doplots=False)
	tohoku_run_stats = rb_runs_report(rb_runs_data=tohoku_runs, rblen=rb_len, log_norm=True, doplots=True, do_clf=True, fignum=7, cat_name=cat_name)
	random_run_stats = rb_runs_report(rb_runs_data=random_runs, doplots=True, fignum=7, do_clf=False, cat_name=cat_name )
	plt.title('raw runs %s then random' % cat_name)
	plt.savefig('output/%s_raw_runs.png' % cat_name)
	
	tohoku_means = rb_runs_report(rb_runs_data=tohoku_runs_mean, rblen=rb_len, log_norm=True, doplots=True, fignum=8, do_clf=True, cat_name=cat_name)
	random_run_means = rb_runs_report(rb_runs_data=random_runs_mean, doplots=True, do_clf=False, fignum=8, cat_name=cat_name)
	plt.title('mean runs %s, then random)' % cat_name)
	plt.savefig('output/%s_mean_runs.png' % cat_name)

	#
	#return tohoku2
	return rb_data
#
def running_mean(data_in=None, mean_len=2):
	# return a running mean over mean_len elements. this will return the "best attempt" for the first mean_len
	# elements (aka average them over as many values as possible. if we wan to be more thorough, we can 1) make this optional
	# (though it's easy enough to trim it on the calling side), and 2) account for this by returning the variance.
	#
	if data_in==None: return []
	if len(data_in)<2: return data_in 
	#
	return [numpy.mean(data_in[max(0, i-mean_len):i]) for i in xrange(1,len(data_in)+1)]
#
def lognorm_ratio(num, denom=None, N=1):
	#return math.log10(float(len(rb_gt))/float(len(rb_lt)))/(math.log10(float(len(rb_seq))))
	if denom!=None:
		# calc ratios
		return numpy.log10(float(num)/float(denom))/numpy.log10(N)
	if denom==None:
		# assume "num" are ratios/fractions
		return numpy.log10(num)/numpy.log10(N)
#
def tri_parity(x):
	if x<0: return -1
	if x==0: return 0
	if x>0: return 1
#	
def randtest():
	R1=random.Random()
	R2=random.Random()
	#
	for i in xrange(10):
		print R1.random(), R2.random()
#
def getNsample(m, mc, mt=7.6, dm=1.0, b1=1.0, b2=1.5, dms0=1.0, doint=True, dmMode=1):
	# dmMode: mode for calculating dm for large earthquakes. 0: relative to mainshock, 1: relative to weighted mean 
	# calculated over (m-m_t) and (m_t-m_c). (see code below); default is (should be) to average.
	targmag=m
	if targmag<mt:
		# "small" earthquake
		winlen=10**(targmag-dm-dms0-mc)	# where 2.0 is dmstar + dmprime
	if targmag>=mt:
		# we want to use an average value for both dms and dm, so dm (like dms) is determined from
		# the full sequence, not simply with respect to the mainshock.
		dmfactor = (1.0*(mt-mc) + 1.5*(targmag-mt))/(targmag-mc)
		#
		#dms = dms0*(1.0*(mt-mc) + dms0*1.5*(targmag-mt))/(targmag-mc)
		dms = dmfactor*dms0
		if dmMode==1: thisdm = dmfactor*dm	# relative to full sequence.
		if dmMode==0: thisdm = 1.5*dm		# relative to largest magnitude
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
		#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-dm) - dms)
		#
		winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt) - dms - thisdm)		# as opposed to dm reletive to b=1.5.
	#
	#winlen=int(10*round(winlen/10))
	#print "winlen0: %d" % winlen
	if doint: winlen=int(round(winlen,-1))
	#print "winlen0: %d" % winlen
	if winlen<1: winlen=1
	#
	return winlen
			
		
