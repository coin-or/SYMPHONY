##############################################################################
##############################################################################
#                                                                            #
# This file is part of a demonstration application for use with the          #
# SYMPHONY Branch, Cut, and Price Library. This application is a solver for  #
# multi-criteria knapsack.                                                   #
#                                                                            #
# (c) Copyright 2004 Lehigh University. All Rights Reserved.                 #
#                                                                            #
# This application was written by Ted Ralphs (tkralphs@lehigh.edu).          #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

##############################################################################
# First, we have to define some variables that will be used
##############################################################################

SYMPHONYROOT = ${HOME}/SYMPHONY
COINROOT = ${HOME}/COIN
MCKPROOT = ${SYMPHONYROOT}/MCKP

##############################################################################
# This section is for listing MATCH source file names
##############################################################################

mckp : mckp_main.c

	g++ -g -L$(SYMPHONYROOT)/lib.LINUX/OSI_CLP -L$(COINROOT)/lib \
	-I$(SYMPHONYROOT)/include -I$(COINROOT)/include \
	mckp_main.c -o mckp -lsym -lCoin -lOsi -lOsiClp -lClp -lCgl

clean : 

	rm mckp mckp_main.o