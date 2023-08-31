#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "phylip.h"
#include "cont.h"
#include "ml.h"
/* version 3.6. (c) Copyright 1999-2004 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */


node* cont_node_new(node_type type,long index) 
{
  node *n;
  n = Malloc(sizeof(cont_node));
  generic_node_init(n, type, index);
  return n;
}

void cont_node_init(cont_node* n, boolean tip, long index)
{
  generic_node_init((node*)n,tip,index);
  ((node*)n)->copy = cont_node_copy;
}

void cont_node_copy(node* srcn, node* dstn)
{
  cont_node *src, *dst;

  src = (cont_node *)srcn;
  dst = (cont_node *)dstn;
  generic_node_copy(srcn, dstn);
  if ( dst->totalleles != 0 && dst->totalleles != src->totalleles )
    free(dst->view);
  dst->view = NULL;
  dst->totalleles = src->totalleles;
  if ( dst->view == NULL )
    dst->view = Malloc(src->totalleles *sizeof(double));
  memcpy(dst->view, src->view, dst->totalleles*sizeof(double));
}

void alloctree(pointarray *treenode, long nonodes)
{
  /* allocate treenode dynamically */
  /* used in contml & contrast */
  long i, j;
  node *p, *q;

  *treenode = (pointarray)Malloc(nonodes*sizeof(node *));
  for (i = 0; i < spp; i++)
    (*treenode)[i] = functions.new_node(1,i + 1);
  for (i = spp; i < nonodes; i++) {
    q = NULL;
    for (j = 1; j <= 3; j++) {
      p = functions.new_node(0,i+1);
      p->next = q;
      q = p;
    }
    p->next->next->next = p;
    (*treenode)[i] = p;
  }
} /* alloctree */


void freetree(pointarray *treenode, long nonodes)
{
  long i, j;
  node *p, *q;

  for (i = 0; i < spp; i++)
    free((*treenode)[i]);
  for (i = spp; i < nonodes; i++) {
    p = (*treenode)[i];
    for (j = 1; j <= 3; j++) {
      q = p;
      p = p->next;
      if (q != NULL)
        free(q);
    }
  }
  free(*treenode);
} /* freetree */


void setuptree(tree *a, long nonodes)
{
  /* initialize a tree */
  /* used in contml & contrast */
  long i, j;
  node *p;

  for (i = 1; i <= spp; i++) {
    a->nodep[i - 1]->back = NULL;
    a->nodep[i - 1]->iter = true;
  }
  for (i = spp + 1; i <= nonodes; i++) {
    p = a->nodep[i - 1];
    for (j = 1; j <= 3; j++) {
      p->back = NULL;
      p->iter = true;
      p = p->next;
    }
  }
  a->score = -99999.0;
  a->root = a->nodep[0];
}  /* setuptree */


void allocview(tree *a, long nonodes, long totalleles)
{
  /* allocate view */
  /* used in contml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++) {
    ((cont_node*)a->nodep[i])->view = 
      (phenotype3)Malloc(totalleles*sizeof(double));
    ((cont_node*)a->nodep[i])->totalleles = totalleles;
  }
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    /* Assumes bifurcation */
    for (j = 1; j <= 3; j++) {
      ((cont_node*)p)->view = (phenotype3)Malloc(totalleles*sizeof(double));
      ((cont_node*)p)->totalleles = totalleles;
      p = p->next;
    }
  }
}  /* allocview */


void freeview(tree *a, long nonodes)
{
  /* deallocate view */
  /* used in contml */
  long i, j;
  node *p;

  for (i = 0; i < spp; i++)
    free(((cont_node*)a->nodep[i])->view);
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    for (j = 1; j <= 3; j++) {
      free(((cont_node*)p)->view);
      p = p->next;
    }
  }
}  /* freeview */


void standev2(long numtrees, long maxwhich, long a, long b, double maxlogl,
              double *l0gl, double **l0gf, longer seed)
{  /* do paired sites test (KHT or SH) on user-defined trees */
  /* used in contml */
  double **covar, *P, *f, *r;
  long i, j, k;
  double sumw, sum, sum2, sd;
  double temp;

#define SAMPLES 1000
  if (numtrees == 2) {
    fprintf(outfile, "Kishino-Hasegawa-Templeton test\n\n");
    fprintf(outfile, "Tree    logL    Diff logL    Its S.D.");
    fprintf(outfile, "   Significantly worse?\n\n");
    i = 1;
    while (i <= numtrees) {
      fprintf(outfile, "%3ld%10.1f", i, l0gl[i - 1]);
      if (maxwhich == i)
        fprintf(outfile, "  <------ best\n");
      else {
        sumw = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        for (j = a; j <= b; j++) {
          sumw += 1;
          temp = l0gf[i - 1][j] - l0gf[maxwhich - 1][j];
          sum += temp;
          sum2 += temp * temp;
        }
        temp = sum / sumw;
        sd = sqrt(sumw / (sumw - 1.0) * (sum2 - temp * temp));
        fprintf(outfile, "%10.1f%12.4f", (l0gl[i - 1])-maxlogl, sd);
        if (sum > 1.95996 * sd)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
      i++;
    }
    fprintf(outfile, "\n\n");
  } else {           /* Shimodaira-Hasegawa test using normal approximation */
    if(numtrees > MAXSHIMOTREES){
      fprintf(outfile, "Shimodaira-Hasegawa test on first %d of %ld trees\n\n"
              , MAXSHIMOTREES, numtrees);
      numtrees = MAXSHIMOTREES;
    } else {
      fprintf(outfile, "Shimodaira-Hasegawa test\n\n");
    }
    covar = (double **)Malloc(numtrees*sizeof(double *));  
    sumw = b-a+1;
    for (i = 0; i < numtrees; i++)
      covar[i] = (double *)Malloc(numtrees*sizeof(double));  
    for (i = 0; i < numtrees; i++) {        /* compute covariances of trees */
      sum = l0gl[i]/sumw;
      for (j = 0; j <=i; j++) {
        sum2 = l0gl[j]/sumw;
        temp = 0.0;
        for (k = a; k <= b ; k++) {
          temp = temp + (l0gf[i][k]-sum)*(l0gf[j][k]-sum2);
        }
        covar[i][j] = temp;
        if (i != j)
          covar[j][i] = temp;
      }
    }
    for (i = 0; i < numtrees; i++) { /* in-place Cholesky decomposition
                                        of trees x trees covariance matrix */
      sum = 0.0;
      for (j = 0; j <= i-1; j++)
        sum = sum + covar[i][j] * covar[i][j];
      temp = sqrt(covar[i][i] - sum);
      covar[i][i] = temp;
      for (j = i+1; j < numtrees; j++) {
        sum = 0.0;
        for (k = 0; k < i; k++)
          sum = sum + covar[i][k] * covar[j][k];
        if (fabs(temp) < 1.0E-12)
          covar[j][i] = 0.0;
        else
          covar[j][i] = (covar[j][i] - sum)/temp;
      }
    }
    f = (double *)Malloc(numtrees*sizeof(double)); /* resampled likelihoods */
    P = (double *)Malloc(numtrees*sizeof(double)); /* vector of P's of trees */
    r = (double *)Malloc(numtrees*sizeof(double)); /* store Normal variates */
    for (i = 0; i < numtrees; i++)
      P[i] = 0.0;
    for (i = 1; i <= SAMPLES; i++) {           /* loop over resampled trees */
      for (j = 0; j < numtrees; j++)          /* draw Normal variates */
        r[j] = normrand(seed);
      for (j = 0; j < numtrees; j++) {        /* compute vectors */
        sum = 0.0;
        for (k = 0; k <= j; k++)
          sum += covar[j][k]*r[k];
        f[j] = sum;
      }
      sum = f[1];
      for (j = 1; j < numtrees; j++)          /* get max of vector */
        if (f[j] > sum)
          sum = f[j];
      for (j = 0; j < numtrees; j++)          /* accumulate P's */
        if (maxlogl-l0gl[j] < sum-f[j])
          P[j] += 1.0/SAMPLES;
    }
    fprintf(outfile, "Tree    logL    Diff logL    P value");
    fprintf(outfile, "   Significantly worse?\n\n");
    for (i = 0; i < numtrees; i++) {
      fprintf(outfile, "%3ld%10.1f", i+1, l0gl[i]);
      if ((maxwhich-1) == i)
        fprintf(outfile, "  <------ best\n");
      else {
        fprintf(outfile, " %9.1f  %10.3f", l0gl[i]-maxlogl, P[i]);
        if (P[i] < 0.05)
          fprintf(outfile, "           Yes\n");
        else
          fprintf(outfile, "           No\n");
      }
    }
  fprintf(outfile, "\n");
  free(P);             /* free the variables we Malloc'ed */
  free(f);
  free(r);
  for (i = 0; i < numtrees; i++)
    free(covar[i]);
  free(covar);
  }
}  /* standev */
