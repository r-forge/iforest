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

#include <R.h>
#include "rt.h"

void weightedAttPool(double *colWeight, int *AttPool, int nColSamp, int xcol)
{ // This function takes in the column weight and selects nColSamp attributes according to column weight.
  double tempColWeight ;
  int i,j, tempj;

  for (i=0;i<nColSamp;++i){
    tempColWeight = colWeight[i];
    tempj = i;
    for (j = i+1; j<xcol; ++j)
    {
      if(colWeight[j] > tempColWeight) {
        tempColWeight = colWeight[j];
        tempj = j;
      }
    }
    colWeight[tempj] = colWeight[i];
    AttPool[tempj] = AttPool[i];

    colWeight[i] = tempColWeight;
    AttPool[i] = AttPool[tempj];
  }
}


void zeroInt(int *x, int length) {
    memset(x, 0, length * sizeof(int));
}


void zeroDouble(double *x, int length) {
    memset(x, 0, length * sizeof(double));
}

void Quicksort (int Fp,int Lp,int Att,double *x, unsigned long int *xref, int xrow, int xcol ) //adopted from C4.5
{

    register int Lower, Middle;
    register double Thresh;
    register int i;

    if (Fp < Lp)
    {
        Thresh = x[(int)xref[Lp] + Att * xrow];//CVal (Item[Lp], Att)

        /* Isolate all items with values <= threshold */

        Middle = Fp;

        for (i = Fp; i < Lp; i++)
        {
            if (x[(int)xref[i] + Att * xrow] <= Thresh)
            {
                if (i != Middle)
//                    swapCases (Middle, i,x,xrow,xcol);
                    swapUnsignedLongInt(Middle, i, xref);
                Middle++;
            }
        }

        /* Extract all values equal to the threshold */

        Lower = Middle - 1;

        for (i = Lower; i >= Fp; i--)
        {
            if (x[(int)xref[i] + Att* xrow] == Thresh)
            {
                if (i != Lower)
                 //   swapCases (Lower,i,x,xrow,xcol);
                    swapUnsignedLongInt(Lower,i,xref);
                Lower--;
            }
        }

        /* Sort the lower values */

        Quicksort (Fp, Lower, Att, x,xref,xrow,xcol);

        /* Position the middle element */
         swapUnsignedLongInt(Middle, Lp,xref);
        /* Sort the higher values */

        Quicksort (Middle + 1, Lp, Att, x,xref,xrow,xcol);
    }
}

int SortAndLocate(int Fp, int Lp, int Att, double *x, unsigned long int *xref, int xrow, int xcol, double splitPoint)
{
  register int Middle, i;
  Middle = Fp;
  for (i = Fp; i<=Lp; i++)
  {
     if ( x[(int)xref[i]+Att*xrow] <= splitPoint)
     {
       if (i != Middle) swapUnsignedLongInt(Middle, i,xref);
       //swapCases(Middle, i,x, xrow, xcol);
       Middle++;
     }
  }
  return Middle-1 ;
}


void swapUnsignedLongInt (int a, int b, unsigned long int *x)
{
  register unsigned long int hold;
  hold = x[a];
  x[a] = x[b];
  x[b] = hold;
}


void swapint (int a, int b, int *x)
{
  register int hold;
  hold = x[a];
  x[a] = x[b];
  x[b] = hold;
}

void swapDouble (int a, int b,double *x)
{
       register double hold;
       hold = x[a];
       x[a] = x[b];
       x[b] = hold;
}

void swapCases (int a, int b, double *x, int n, int m)
{
  /* n is the number of rows*/
  /* m is the number of cols*/
  register int i;
  for (i=0;i < m;++i)
    swapDouble( a + i * n, b + i * n, x);
}

