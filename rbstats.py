import random
import numpy
import math
import pylab as plt
import scipy.stats

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
		rb_vals[-1]+=[float(rb_vals[-1][-2])/float(rb_vals[-1][-1])]
		rb_vals[-1]+=[math.log10(rb_vals[-1][-1])/log_N]
	rb_vals = numpy.core.records.fromarrays(zip(*rb_vals), names=dtype_names, formats=[type(x).__name__ for x in rb_vals[0]])
	return rb_vals	
		
#
#def rb_cdfs(rblen=100, nits=10000):
def random_rb_sequence(rblen=100, nits=10000):
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
		[rb_gt.append(x) for x in vals if x>rb_gt[-1]]
		[rb_lt.append(x) for x in vals if x<rb_lt[-1]] 
		#
		rb_vals += [[i, len(rb_gt), len(rb_lt)]]
		rb_vals[-1]+=[float(rb_vals[-1][-2])/float(rb_vals[-1][-1])]
		rb_vals[-1]+=[math.log10(rb_vals[-1][-1])/log_N]
	#
	rb_vals = numpy.core.records.fromarrays(zip(*rb_vals), names=dtype_names, formats=[type(x).__name__ for x in rb_vals[0]])
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

def rb_runs_report(rb_runs_data=None, rblen=128, seq_len=100000, log_norm=True, fignum=0, doplots=True, do_clf=True):
	# some stats on record-breaking runs.
	if rb_runs_data==None: rb_runs_data = rb_runs(rblen=rblen, seq_len=seq_len, log_norm=log_norm)
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
		#Ns = [x/float(len(these_lengths)) for x in range(1, len(these_lengths)+1)]
		Ns = [x/float(total_N_runs) for x in range(1, len(these_lengths)+1)]
		these_lengths.sort()
		if doplots:
			#plt.gca().set_yscale('log')
			#plt.gca().set_xscale('log')
			plt.gca().set_yscale('linear')
			plt.gca().set_xscale('linear')
			plt.plot(these_lengths, Ns, '.-', label='parity: %d' % key)
		#
	if doplots:
		plt.legend(loc=0, numpoints=1)
		plt.ylabel('Probability P(L)')
		plt.xlabel('Length L')
	#
	print "summary report:"
	for key, datas in run_stats.iteritems():
		print key, {ky:val for ky,val in datas.iteritems() if ky!='run_lengths'}
	#
	return run_stats
#
# instead of writing a wrapper function, just defind the prams dictionary:
parkfield_rb_report_prams = {'data_file':'data/parkfield-elip-rbsequence.pkl', 'rb_len':310, 'fnum':0, 'random_len':100000, 'ave_len':None}
emc_rb_report_prams = {'data_file':'data/emc_rb_sequence.pkl', 'rb_len':500, 'fnum':0, 'random_len':100000, 'ave_len':None}
tohoku_rb_report_prams = {'data_file':'data/tohoku_rb_sequence.pkl', 'rb_len':220, 'fnum':0, 'random_len':100000, 'ave_len':None}
#
def chi_chi_rb_report(data_file='data/chichi-rb_sequence.pkl', rb_len=630, fnum=0, random_len=10000, ave_len=None):
	return rb_stats_report(data_file=data_file, rb_len=rb_len, fnum=fnum, random_len=random_len)

def tohoku_rb_report(data_file='data/tohoku_rb_sequence.pkl', rb_len=220, fnum=0, random_len=10000, ave_len=None):
	return rb_stats_report(data_file=data_file, rb_len=rb_len, fnum=fnum, random_len=random_len, ave_len=ave_len)
	#return rb_stats_report(*args, **kwargs)
#
def rb_stats_report(data_file='data/tohoku_rb_sequence.pkl', rb_len=220, fnum=0, random_len=10000, ave_len=None):
	# we'll need to encode some of these data, but for now, we know that rb_len was 220 (probably).
	# (a good guess for this is the first event number rw_0,item_0=219).
	# ... so this generalized function was orignally written for Tohoku, so "tohoku" probably impies "test data",
	# as opposed to "random control data"
	#
	if ave_len==None:
		ave_len = max(int(rb_len/10), 1)
	if isinstance(ave_len, float):
		# this is a quick way to call fractioal ave_len values, but it could go awry if we make a mistake
		# and submit values like ave_len=1.0 when we mean 1. we could also require odd floats for this sort of work,
		# aka 1.0 --> 1 but 1.1 (or 1.0001) --> 1.0001*max(int(rb_len/10), 1)
		ave_len = int(ave_len*max(int(rb_len/10), 1))
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
	plt.xlabel('normalized rb_ratio value, $abs(nrb_{gt}/nrb_{lt})$')
	plt.ylabel('Probability $P(abs[r])$')
	#
	#Y = [x/float(len(tohoku_rb_gt)) for x in xrange(1, len(tohoku_rb_gt)+1)]	# prob assuming gt...
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_gt)+1)]						# prob within full sequence.
	plt.plot(tohoku_rb_gt, Y, '.-', label='tohoku_gt')
	#
	#Y = [x/float(len(tohoku_rb_lt)) for x in xrange(1, len(tohoku_rb_lt)+1)]
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_lt)+1)]
	plt.plot(tohoku_rb_lt, Y, '.-', label='tohoku_lt')
	#
	#Y = [x/float(len(tohoku_rb_eq)) for x in xrange(1, len(tohoku_rb_eq)+1)]
	#Y = [x/N_total for x in xrange(1, len(tohoku_rb_eq)+1)]
	#plt.plot(tohoku_rb_eq, Y, '.-', label='tohoku_eq')
	#
	# and the means...
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_gt_mean)+1)]						# prob within full sequence.
	plt.plot(tohoku_rb_gt_mean, Y, '.-', label='tohoku_gt_mean')
	#
	Y = [x/N_total for x in xrange(1, len(tohoku_rb_lt_mean)+1)]
	plt.plot(tohoku_rb_lt_mean, Y, '.-', label='tohoku_lt_mean')
	#
	#Y = [x/N_total for x in xrange(1, len(tohoku_rb_eq_mean)+1)]
	#plt.plot(tohoku_rb_eq_mean, Y, '.-', label='tohoku_eq')
	#
	############
	# now, plot some random sequences for comparison:
	#rand_plots = plot_rb_cdfs(cdf_in=None, rblen=rb_len, nits=10000, fnum=0, do_clf=True)
	rand_cdf = random_rb_sequence(rblen=rb_len, nits=random_len)	# returns recarray with: [i, n_gt, n_lt, ratio, ratio_lognorm]
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
	tohoku_run_stats = rb_runs_report(rb_runs_data=tohoku_runs, rblen=rb_len, log_norm=True, doplots=True, do_clf=True, fignum=7)
	random_run_stats = rb_runs_report(rb_runs_data=random_runs, doplots=True, fignum=7, do_clf=False)
	plt.title('raw runs (tohoku then random')
	
	tohoku_means = rb_runs_report(rb_runs_data=tohoku_runs_mean, rblen=rb_len, log_norm=True, doplots=True, fignum=8, do_clf=True)
	random_run_means = rb_runs_report(rb_runs_data=random_runs_mean, doplots=True, do_clf=False, fignum=8)
	plt.title('mean runs (tohoku, then random)')

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
			
		
