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
ABS_GAPTOL=0.1
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
	for line in arr0:
		find = re.search(st0,line)
		if (find >= 0):
			find = re.search('-*[0-9]+\.[0-9]+',line)
			if (find>=0):
				l = re.findall('-*[0-9]+\.[0-9]+',line)
				if (len(l)<3):
					return -1,fl0
				else:
					fl0 = float(l[2])
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

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

for i in range(0,len(sys.argv)):
	if (sys.argv[i]=='-d'):
		if (i==len(sys.argv)):
			print "Missing argument to '-d'"
			print "usage: python report_tsp.py [-d <path to dir>]"
			sys.exit(0)
		else:
			OUTPUT_DIR=sys.argv[i+1]
			break

a = []
fl0 = 0.0
in0 = 0
error = []

flist=open(INST_LIST,'r')
a=flist.read().split()
flist.close()
a.sort()
#print "### Instance set:", INST_SET
print "### Reading from directory:", OUTPUT_DIR
print "%16s"%"Instance","%14s"%"best ub","%3s"%"opt","%14s"%"ub","%14s"%"lb","%6s"%"gap","%8s"%"time","%8s"%"cut","%8s"%"branch","%8s"%"lp","%7s"%"nodes-c","%7s"%"nodes-a"
for instance in a:
	print "%16s"%instance,

	fil=open(BEST_BOUND_FILE,'r')
	whole_file=fil.read().split('\n')
	fil.close()

	best_ub = INFTY
	find,best_ub=find_float(whole_file,instance,best_ub)
	if (find<0 or best_ub>=INFTY):
		print  "%14s"%"NF",
	else:
		print  "%14.2f"%best_ub,

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
			print  "%14s"%"NF", #ub
			print  "%14s"%"NF", #lb
			print  "%6s"%"NF", #gap
		else:
			gap = 0.0
			lb = ub
			print  "%14.2f"%ub,
			print  "%14.2f"%lb,
			print  "%6.2f"%gap,
	else:
		ub = INFTY
		find,ub=find_float(whole_file,'Current Upper Bound',ub)
		if (find<0 or ub >= INFTY):
			print "%14s"%"NF",
		else:
			print "%14.2f"%ub,

		lb = INFTY
		find,lb=find_float(whole_file,'Current Lower Bound',lb)
		if (find<0 or lb >= INFTY):
			print  "%14s"%"NF",
		else:
			print  "%14.2f"%lb,
		
		gap = INFTY
		find,gap=find_float(whole_file,'Gap Percentage',gap)
		if (find<0 or gap >= INFTY):
			print  "%6s"%"NF",
		else:
			print  "%6.2f"%gap,
	
	if (claims_optimal==1 and best_ub < INFTY and abs(best_ub-ub)>EPS_UB):
		error = error+instance

	if (best_ub < INFTY and ub < INFTY and best_ub-ub>EPS_UB):
		error.append(instance)

	if (best_ub < INFTY and lb-best_ub>EPS_UB):
		error.append(instance)

	totalTime = INFTY
	find,totalTime=find_float(whole_file,'Total Wallclock Time',totalTime)
	if (find<0 or totalTime >= INFTY):
		print  "%8s"%"NF",
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

	nodes_c = INFTY
	find,nodes_c=find_int(whole_file,'Number of created nodes',nodes_c)
	if (find<0 or nodes_c >= INFTY):
		print  "%7s"%"NF",
	else:
		print  "%7d"%nodes_c,
		
	nodes_a = INFTY
	find,nodes_a=find_int(whole_file,'Number of analyzed nodes',nodes_a)
	if (find<0 or nodes_a >= INFTY):
		print  "%7s"%"NF",
	else:
		print  "%7d"%nodes_a,
		
	'''
	noWSTime = INFTY
	find,noWSTime=find_float(whole_file,'Total nonWS solution time:',noWSTime)
	if (find<0 or noWSTime >= INFTY):
		print delim, "%9s"%"NF",
	else:
		if (noWSTime<EPS_TIME):
			noWSTime = EPS_TIME
		print delim, "%9.2f"%noWSTime,

	#sumWSTime=float(noWSTime)+float(WSTime)
	sumWSTime=1
	if (sumWSTime>INFTY):
		print delim, "%9s"%"NF",delim,"%9s"%"NF",
	else:
		print delim,"%9.2f"%(float(WSTime)/sumWSTime*100),
		print delim,"%9.2f"%(float(noWSTime)/sumWSTime*100),

	L1WS = INFTY
	find,L1WS =find_float(whole_file,'L1 WS solution time:',L1WS)
	if (find<0 or L1WS >= INFTY):
		print delim, "%9s"%"NF",
	else:
		print delim, "%9.2f"%L1WS,

	L1NoWS = INFTY
	find,L1NoWS =find_float(whole_file,'L1 nonWS solution time:',L1NoWS)
	if (find<0 or L1NoWS >= INFTY):
		print delim, "%9s"%"NF",
	else:
		print delim, "%9.2f"%L1NoWS,

	L2WS = INFTY
	find,L2WS =find_float(whole_file,'L2 WS solution time:',L2WS)
	if (find<0 or L2WS >= INFTY):
		print delim, "%9s"%"NF",
	else:
		print delim, "%9.2f"%L2WS,

	L2NoWS = INFTY
	find,L2NoWS =find_float(whole_file,'L2 nonWS solution time:',L2NoWS)
	if (find<0 or L2NoWS >= INFTY):
		print delim, "%9s"%"NF",
	else:
		print delim, "%9.2f"%L2NoWS,
	'''
	print ''

print "## errors:",error
