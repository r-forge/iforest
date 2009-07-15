/*###########################################################################
#
#    Isolation Forest -- Anomaly detection using binary trees
#    Copyright (C) 2009 Fei Tony Liu
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##########################################################################*/

#include <stdbool.h>

#ifndef RT_H
#define RT_H

/* test if the bit at position pos is turned on */
#define isBitOn(x,pos) (((x) & (1 << (pos))) > 0)
/* swap two integers */
#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))


void rTrees(double *x, int *xrow, int *xcol, int *nrnodes, int *ntree, int *hlim, int *rowSamp, int *nRowSamp, int *colSamp, int *nColSamp, int *nmin, double *rFactor, double *colWeight, int *nodeStatus, int *lDaughter, int *rDaughter, int  *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *ntreeSize);

void rTree(double *x, unsigned long int *xref, int xrow, int xcol,  int nsample,int nrnodes, int hlim, int nmin, double rFactor, bool *AttSkip, int *AttPool, int nColSamp, int *nodeStatus, int *lDaughter, int *rDaughter, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *ntreeSize);

void rtDepth(double *x, int *xrow, int *xcol, int *samSize, int *nrnodes, int *ntree, int *hlim, int *nodeStatus, int *lDaughter, int *rDaughter, int  *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int *appRange, double *outF, double *pathLength, int *isoby, double *isoat);

void trvTree(double *x, int xrow, int xcol, int hlim, int *nodeStatus, int *lDaughter, int *rDaughter, int  *splitAtt, double *splitPoint, double *ulim, double *llim, int *nSam, int appRange, double *totDepth, int *isoby, double *isoat);

void SDGainSplit(double *x, unsigned long int *xref, int xrow, int xcol, int ndStart, int ndEnd, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *ndEndl, bool *AttSkip);

void randomSplit(double *x, unsigned long int *xref, int xrow, int xcol, int ndStart, int ndEnd, int *splitAtt, double *splitPoint, double *ulim, double *llim, int *ndEndl, int *AttPool, int nColSamp);

void swapint(int a, int b, int *x);

void swapUnsignedLongInt (int a, int b, unsigned long int *x);

void swapDouble (int a, int b,double *x);

void swapCases (int a, int b, double *x, int n, int m);

void Quicksort (int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol );

int SortAndLocate(int Fp, int Lp, int Att, double *x, unsigned long int *xref, int xrow, int xcol, double splitPoint);

double unif_rand();

//double avgPL(int n);

void zeroInt(int *x, int length);

void zeroDouble(double *x, int length);

void weightedAttPool(double *colWeight, int *AttPool, int nColSamp, int xcol);


#define NODE_TERMINAL -1
#define NODE_TOSPLIT  -2
#define NODE_INTERIOR -3

#define  Max(a,b)               ((a)>(b) ? a : b)
#define	 Min(a,b)               ((a)<(b) ? a : b)
#define	 Round(x)		((int) (x+0.5))
#define  avgPL(n)              (((n-1) <= 0) ? 0.0 : (( 2.0 * (log((double)(n-1)) + 0.5772156649)) - ( 2.0 * (double)(n-1))/( 1.0 + (double)(n-1))))

#endif /* RT_H */
