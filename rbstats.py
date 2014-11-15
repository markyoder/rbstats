import random
import numpy
import math
import pylab as plt
import scipy.stats

def plot_rb_cdfses(rblens=[8, 32, 128, 256, 512, 1024, 2048, 4096, 2*4096], nits=1000, fnum=0):
	do_clf=True
	#
	for i, rblen in enumerate(rblens):
		a=plot_rb_cdfs(rblen=rblen, nits=nits, fnum=fnum, do_clf=do_clf)
		do_clf=False
	#
	return None

def plot_rb_cdfs(cdf_in=None, rblen=100, nits=10000, fnum=0, do_clf=True):
	if cdf_in==None:
		cdf_in = rb_cdfs(rblen=rblen, nits=nits)
	#
	log_ratios = [abs(x) for x in cdf_in['ratio_lognorm']]
	log_ratios.sort()
	#Ns=range(1,len(log_ratios)+1)
	Ps = [float(i+1)/float(nits) for i in xrange(nits)]
	#print "mode: ", scipy.stats.mode(log_ratios)
	print "rblen: %d. max: %f, min: %f, mean: %f, stdev: %f, median: %f/%f, mode: %f" % (rblen, max(log_ratios), min(log_ratios), numpy.mean(log_ratios), numpy.std(log_ratios), numpy.median(log_ratios), log_ratios[int(len(log_ratios)/2)], scipy.stats.mode(log_ratios)[0])
	#
	plt.figure(fnum)
	if do_clf: plt.clf()
	plt.gca().set_yscale('linear')
	plt.gca().set_xscale('log')
	plt.plot(log_ratios, Ps, '.-', label='rblen=%d' % rblen)
	plt.legend(loc=0, numpoints=1)
	plt.xlabel('normalized log_ratios')
	plt.ylabel('probability P(l_r)')
	#
	return cdf_in
#
def rb_cdfs(rblen=100, nits=10000):
	R=random.Random()
	#
	rb_vals = []
	dtype_names = ['i', 'n_gt', 'n_lt', 'ratio', 'ratio_lognorm']
	log_N = math.log10(rblen)
	#
	for i in xrange(nits):
		vals = [R.random() for j in xrange(nits)]
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
		return math.log10(float(len(rb_gt))/float(len(rb_lt)))/(math.log10(float(len(rb_seq))))
#
def rb_runs(rblen=1, seq_len=100000, log_norm=True):
	# churn up som statistics on rb "runs", aka sequences of 1,0 values (in this case,
	# well use 1 = (ratio>1), -1 = (ratio<1), and 0: ratio=1 ???
	#
	R=random.Random()
	sequence = [R.random() for x in xrange(seq_len)]
	#
	rb_ratios = [rb_ratio(sequence[i-rblen:i], log_norm=log_norm) for i in xrange(rblen, seq_len)]
	runs_gt = []
	runs_lt = []
	runs_0  = []
	#
	j=0
	while j<seq_len:
		pass
	#
	return rb_ratios
	
def randtest():
	R1=random.Random()
	R2=random.Random()
	#
	for i in xrange(10):
		print R1.random(), R2.random()
			
		
