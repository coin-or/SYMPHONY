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
	print "usage: python report_tsp.py -d <path to dir> [-c] [-h] [-n]"

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

has_user_dir = 0
has_user_nodes = 0
has_user_heurs = 0
has_user_cuts = 0

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
	elif (sys.argv[i]=='-c'):
		has_user_cuts = 1
	elif (sys.argv[i]=='-h'):
		has_user_heurs = 1
	elif (sys.argv[i]=='-n'):
		has_user_nodes = 1
	else:
		print "invalid option: %s"%sys.argv[i]
		print_usage()
		sys.exit(0)
	i = i+1
		
if (has_user_dir<1):
	print "usage: python report_tsp.py -d <path to dir> [-c] [-h] [-n]"
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
print "%18s"%"Instance","%16s"%"best ub","%3s"%"opt","%16s"%"ub","%16s"%"lb","%6s"%"gap","%8s"%"time","%8s"%"cut","%8s"%"branch","%8s"%"lp","%8s"%"unaccntd",
if (has_user_nodes>0):
	print "%8s"%"nodes-c","%8s"%"nodes-a",
if (has_user_heurs>0):
	print "%8s"%"fp-time", "%8s"%"fp-sols",
if (has_user_cuts>0):
	print "%8s"%"cuts", "%8s"%"rt-cuts", "%8s"%"thrown", "%8s"%"time-gom","%8s"%"time-kna","%8s"%"time-odd","%8s"%"time-cli","%8s"%"time-pro","%8s"%"time-flo",
print ''
for instance in a:
	print "%18s"%instance,

	fil=open(BEST_BOUND_FILE,'r')
	whole_file=fil.read().split('\n')
	fil.close()

	best_ub = INFTY
	find,best_ub=find_float_1(whole_file,instance,best_ub)
	if (find<0 or best_ub>=INFTY):
		print  "%16s"%"NF",
	else:
		print  "%16.2f"%best_ub,

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

	if (claims_optimal==1):
		lb = INFTY
		ub = INFTY
		gap = INFTY
		find,ub=find_float(whole_file,'Solution Cost:',ub)
		if (find<0 or ub >= INFTY):
			print  "%16s"%"NF", #ub
			print  "%16s"%"NF", #lb
			print  "%6s"%"-1", #gap
		else:
			gap = ABS_GAPTOL
			lb = ub
			print  "%16.2f"%ub,
			print  "%16.2f"%lb,
			print  "%6.2f"%gap,
	else:
		ub = INFTY
		find,ub=find_float(whole_file,'Current Upper Bound',ub)
		if (find<0 or ub >= INFTY):
			print "%16s"%"NF",
		else:
			print "%16.2f"%ub,

		lb = INFTY
		find,lb=find_float(whole_file,'Current Lower Bound',lb)
		if (find<0 or lb >= INFTY):
			print  "%16s"%"NF",
		else:
			print  "%16.2f"%lb,
		
		gap = INFTY
		find,gap=find_float(whole_file,'Gap Percentage',gap)
		if (find<0 or gap >= INFTY):
			print  "%6s"%"-1",
		else:
			if (gap<=ABS_GAPTOL):
				gap=ABS_GAPTOL
			print  "%6.2f"%gap,
	
	if (claims_optimal==1 and best_ub < INFTY and abs(best_ub-ub)>EPS_UB):
		error.append(instance)

	if (best_ub < INFTY and ub < INFTY and best_ub-ub>EPS_UB):
		error.append(instance)

	if (best_ub < INFTY and lb-best_ub>EPS_UB):
		error.append(instance)

	totalTime = INFTY
	find,totalTime=find_float(whole_file,'Total Wallclock Time',totalTime)
	if (find<0 or totalTime >= INFTY):
		print  "%8s"%"-1",
	else:
		print  "%8.2f"%totalTime,

	sep_time = INFTY
	find,sep_time=find_float(whole_file,'Separation',sep_time)
	if (find<0 or sep_time >= INFTY):
		print  "%8s"%"NF",
	else:
		print  "%8.2f"%sep_time,

	branch_time = INFTY
	find,branch_time=find_float(whole_file,'Strong Branching',branch_time)
	if (find<0 or branch_time >= INFTY):
		print  "%8s"%"NF",
	else:
		if (branch_time<EPS_TIME):
			branch_time = EPS_TIME
		print  "%8.2f"%branch_time,

	lp_time = INFTY
	find,lp_time=find_float(whole_file,'LP Solution Time',lp_time)
	if (find<0 or lp_time >= INFTY):
		print  "%8s"%"NF",
	else:
		if (lp_time<EPS_TIME):
			lp_time = EPS_TIME
		print  "%8.2f"%lp_time,

	if (totalTime>=INFTY):
		print "%8.2s"%"NF",
	else:
		print "%8.2f"%(totalTime-sep_time-branch_time-lp_time),

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
		find,cuts=find_int(whole_file,'total cuts discarded',cuts)
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
			
	print ''

print "## errors:",error
