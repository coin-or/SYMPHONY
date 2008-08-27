#!/usr/bin/python


import array
import os
import commands
import time
import math
import sys
import string
import fileinput
import re

# ---------------------------------------------------------------------------
# IMPORTANT PARAMETERS
# ---------------------------------------------------------------------------
INST_LIST="/home/asm4/instances/deluge.list"
#INST_SET="mittleman_unibo"
# OUTPUT_DIR is overriden -d command line switch
OUTPUT_DIR="/home/asm4/running/d00/deluge"
BEST_BOUND_FILE="/home/asm4/instances/deluge.ub"
#MIP_FILE="jlf.sor"
INFTY=10e10
ABS_GAPTOL=0.001
EPS_TIME = 0.01
EPS_UB=0.01
#TIME_LIMIT = 1800
delim = ""
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------



def find_str(arr0,st0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			return 1
	return 0

def find_float(arr0,st0,fl0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			fl0 = float(re.search('-*[0-9]+\.[0-9]+',line).group())
			return 1,fl0
	return -1,fl0

def find_float_1(arr0,st0,fl0):
	st0 = st0+" "
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			find = re.search('-*[0-9]+\.[0-9]+',line)
			if (find>=0):
				l = re.findall('-*[0-9]+\.[0-9]+',line)
				if (len(l)<1):
					return -1,fl0
				else:
					fl0 = float(l[len(l)-1])
					return 1,fl0
			else:
				return -1,fl0
	return -1,fl0

def find_float_e(arr0,st0,fl0):
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			fl0 = float(re.search('-*[0-9]+\.[0-9]*e*\+*[0-9]*',line).group())
			return 1,fl0
	return -1,fl0

def find_int(arr0,st0,in0):
	in0=0
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			in0 = int(re.search('-*[0-9]+$',line).group())
			return 1,in0
	return -1,in0

def print_usage():
	print "usage: python report_tsp.py -d <path to dir> [-c] [-h] [-n] [-g]",
	print "[-b] [-p] [-m]"
	print "   -c: information about cuts"
	print "   -h: information about heuristics"
	print "   -b: information about branching"
	print "   -n: information about nodes and tree sizes"
	print "   -g: general info about total time, lb, ub, number of solutions"
	print "   -p: general info about presolve (not implemented)"
	print "   -m: memory usage"
	print "opt status and total time taken is always displayed"

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

has_user_dir      = 0
has_user_nodes    = 0
has_user_heurs    = 0
has_user_cuts     = 0
has_user_gen      = 0
has_user_branch   = 0
has_user_presolve = 0
has_user_memory   = 0

if (len(sys.argv)<2):
	print_usage()
	sys.exit(0)

i = 1
while(i<len(sys.argv)):
	if (sys.argv[i]=='-d'):
		has_user_dir = 1
		if (i==len(sys.argv)-1):
			print "Missing argument to '-d'"
			print_usage()
			sys.exit(0)
		else:
			OUTPUT_DIR=sys.argv[i+1]
			if (os.path.exists(OUTPUT_DIR)):
				print "### Reading from directory:", OUTPUT_DIR
				i = i+1
			else:
				print "the specified directory %s is not accessible"%OUTPUT_DIR
				print_usage()
				sys.exit(0)
	elif (sys.argv[i]=='-n'):
		has_user_nodes = 1
	elif (sys.argv[i]=='-h'):
		has_user_heurs = 1
	elif (sys.argv[i]=='-c'):
		has_user_cuts = 1
	elif (sys.argv[i]=='-g'):
		has_user_gen = 1
	elif (sys.argv[i]=='-b'):
		has_user_branch = 1
	elif (sys.argv[i]=='-p'):
		has_user_presolve = 1
	elif (sys.argv[i]=='-m'):
		has_user_memory = 1
	else:
		print "invalid option: %s"%sys.argv[i]
		print_usage()
		sys.exit(0)
	i = i+1
		
if (has_user_dir<1):
	print_usage()
	sys.exit(0)

a = []
fl0 = 0.0
in0 = 0
error = []

flist=open(INST_LIST,'r')
a=flist.read().split()
flist.close()
a.sort()
#print "### Instance set:", INST_SET
print "%18s"%"## Instance","%3s"%"opt","%8s"%"time",
if (has_user_nodes>0):
	print "%8s"%"nodes-c","%8s"%"nodes-a","%8s"%"chains","%8s"%"depth","%8s"%"d-halts",
if (has_user_heurs>0):
	print "%8s"%"fp-time", "%8s"%"fp-sols",
if (has_user_cuts>0):
	print "%8s"%"cuts", "%8s"%"rt-cuts", "%8s"%"bad-coef", "%8s"%"duplicat","%8s"%"time-gom","%8s"%"time-kna","%8s"%"time-odd","%8s"%"time-cli","%8s"%"time-pro","%8s"%"time-flo", "%8s"%"time-che",
	if (has_user_gen<1):
		print "%8s"%"time-cut",
if (has_user_gen>0):
	print "%16s"%"best-ub","%16s"%"ub","%16s"%"lb","%6s"%"gap","%8s"%"lp-time","%8s"%"pre-time","%8s"%"heu-time","%8s"%"cut-time","%8s"%"bra-time","%8s"%"unaccntd",
if (has_user_branch>0):
	print "%8s"%"str-time",
if (has_user_presolve>0):
	print "%8s"%"pre-time",
if (has_user_memory>0):
	print "%8s"%"max-MB",
print ''
for instance in a:
	print "%18s"%instance,

	filename=OUTPUT_DIR+"/"+instance+".condor.out"
	fil=open(filename,'r')
	whole_file=fil.read().split('\n')
	fil.close()

	gave_wrong_result = 0

	# check if claims optimal
	claims_optimal = 0
	find = find_str(whole_file,'Optimal Solution Found')
	if (find>0):
		claims_optimal = 1
	else:
		claims_optimal = 0
	print "%3d"%claims_optimal,

   #==========================================================================
	# total time
	totalTime = INFTY
	find,totalTime=find_float(whole_file,'Total Wallclock Time',totalTime)
	if (find<0 or totalTime >= INFTY):
		print  "%8s"%"-1",
	else:
		print  "%8.2f"%totalTime,

   #==========================================================================
   # find best-ub
	fil=open(BEST_BOUND_FILE,'r')
	whole_file2=fil.read().split('\n')
	fil.close()

	best_ub = INFTY
	find,best_ub=find_float_1(whole_file2,instance,best_ub)
	if (find<0 or best_ub>=INFTY):
		if (has_user_gen>0):
			print  "%16s"%"NF",
	else:
		if (has_user_gen>0):
			print  "%16.2f"%best_ub,

   # find lb, ub, gap
	if (claims_optimal==1):
		lb = INFTY
		ub = INFTY
		gap = INFTY
		find,ub=find_float(whole_file,'Solution Cost:',ub)
		if (find<0 or ub >= INFTY):
			if (has_user_gen>0):
				print  "%16s"%"NF", #ub
				print  "%16s"%"NF", #lb
				print  "%6s"%"-1", #gap
		else:
			gap = ABS_GAPTOL
			lb = ub
			if (has_user_gen>0):
				print  "%16.2f"%ub,
				print  "%16.2f"%lb,
				print  "%6.2f"%gap,
	else:
		ub = INFTY
		find,ub=find_float(whole_file,'Current Upper Bound',ub)
		if (find<0 or ub >= INFTY):
			if (has_user_gen>0):
				print "%16s"%"NF",
		else:
			if (has_user_gen>0):
				print "%16.2f"%ub,

		lb = INFTY
		find,lb=find_float(whole_file,'Current Lower Bound',lb)
		if (find<0 or lb >= INFTY):
			if (has_user_gen>0):
				print  "%16s"%"NF",
		else:
			if (has_user_gen>0):
				print  "%16.2f"%lb,
		
		gap = INFTY
		find,gap=find_float(whole_file,'Gap Percentage',gap)
		if (find<0 or gap >= INFTY):
			if (has_user_gen>0):
				print  "%6s"%"-1",
		else:
			if (gap<=ABS_GAPTOL):
				gap=ABS_GAPTOL
			if (has_user_gen>0):
				print  "%6.2f"%gap,
	
   #==========================================================================
	if (claims_optimal==1 and best_ub < INFTY and abs(best_ub-ub)>EPS_UB):
		error.append(instance)

	if (best_ub < INFTY and ub < INFTY and best_ub-ub>EPS_UB):
		error.append(instance)

	if (best_ub < INFTY and lb-best_ub>EPS_UB):
		error.append(instance)

   #==========================================================================
	if (has_user_gen>0):
		lp_time = INFTY
		find,lp_time=find_float(whole_file,'LP Solution Time',lp_time)
		if (find<0 or lp_time >= INFTY):
			print  "%8s"%"NF",
		else:
			if (lp_time<EPS_TIME):
				lp_time = EPS_TIME
			print  "%8.2f"%lp_time,

		#=======================================================================
		pre_time = INFTY
		find,pre_time=find_float(whole_file,'Presolve Time',pre_time)
		if (find<0 or pre_time >= INFTY):
			print  "%8s"%"NF",
			pre_time = 0
		else:
			print  "%8.2f"%pre_time,

		#=======================================================================
		primal_time = INFTY
		find,primal_time=find_float(whole_file,'Primal Heuristics',primal_time)
		if (find<0 or primal_time >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%primal_time,

		#=======================================================================
		sep_time = INFTY
		find,sep_time=find_float(whole_file,'Separation',sep_time)
		if (find<0 or sep_time >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%sep_time,

		#=======================================================================
		branch_time = INFTY
		find,branch_time=find_float(whole_file,'Strong Branching',branch_time)
		if (find<0 or branch_time >= INFTY):
			print  "%8s"%"NF",
		else:
			if (branch_time<EPS_TIME):
				branch_time = EPS_TIME
			print  "%8.2f"%branch_time,

		#=======================================================================
		if (totalTime>=INFTY):
			print "%8.2s"%"NF",
		else:
			print "%8.2f"%(totalTime-lp_time-pre_time-primal_time-sep_time-
					branch_time),

   #==========================================================================
	if (has_user_cuts>0):
		cuts = INFTY
		find,cuts=find_int(whole_file,'total cuts generated',cuts)
		if (find<0 or cuts >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%cuts,

		cuts = INFTY
		find,cuts=find_int(whole_file,'cuts in root',cuts)
		if (find<0 or cuts >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%cuts,

		cuts = INFTY
		find,cuts=find_int(whole_file,'cuts removed because of bad coeffs:'
				,cuts)
		if (find<0 or cuts >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%cuts,

		cuts = INFTY
		find,cuts=find_int(whole_file,'cuts removed because of duplicacy:'
				,cuts)
		if (find<0 or cuts >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%cuts,

		ctime = INFTY
		find,ctime=find_float(whole_file,'time in gomory cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in knapsack cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in oddhole cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in clique cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in probing cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in flow and cover cuts',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		ctime = INFTY
		find,ctime=find_float(whole_file,'time in checking quality',ctime)
		if (find<0 or ctime >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%ctime,
			
		if (has_user_gen<1):
			#=======================================================================
			sep_time = INFTY
			find,sep_time=find_float(whole_file,'Separation',sep_time)
			if (find<0 or sep_time >= INFTY):
				print  "%8s"%"NF",
			else:
				print  "%8.2f"%sep_time,
   #==========================================================================
	if (has_user_heurs>0):
		fp_time = INFTY
		find,fp_time=find_float(whole_file,'Time spent in feasibility pump',fp_time)
		if (find<0 or fp_time >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%fp_time,

		fp_sols = INFTY
		find,fp_sols=find_int(whole_file,'Number of solutions found by feasibility pump',fp_sols)
		if (find<0 or fp_sols >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%fp_sols,

   #==========================================================================
	if (has_user_branch>0):
		str_time = INFTY
		find,str_time=find_float(whole_file,'Strong Branching',str_time)
		if (find<0 or str_time >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%str_time,

   #==========================================================================
	if (has_user_nodes>0):
		nodes_c = INFTY
		find,nodes_c=find_int(whole_file,'Number of created nodes',nodes_c)
		if (find<0 or nodes_c >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%nodes_c,
		
		nodes_a = INFTY
		find,nodes_a=find_int(whole_file,'Number of analyzed nodes',nodes_a)
		if (find<0 or nodes_a >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%nodes_a,
		
		chains = INFTY
		find,chains=find_int(whole_file,'Number of Chains',chains)
		if (find<0 or chains >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%chains,
		
		depth = INFTY
		find,depth=find_int(whole_file,'Depth of tree',depth)
		if (find<0 or depth >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%depth,
		
		d_halts = INFTY
		find,d_halts=find_int(whole_file,'Number of Diving Halts',d_halts)
		if (find<0 or d_halts >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8d"%d_halts,
		
   #==========================================================================
	if (has_user_memory>0):
		max_mb = INFTY
		find,max_mb=find_float(whole_file,'Virtual memory used',max_mb)
		if (find<0 or max_mb >= INFTY):
			print  "%8s"%"NF",
		else:
			print  "%8.2f"%max_mb,
		
	print ''

print "## errors:",error
