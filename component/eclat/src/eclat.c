/*----------------------------------------------------------------------
  File    : eclat.c
  Contents: eclat algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2002.06.09 file created from apriori.c
            2002.12.10 option -l (list supporting transactions added)
            2002.08.16 transaction reading improved
            2003.08.18 memory benchmarking functionality added
            2003.08.20 option -t (target type) added
            2003.08.22 based on transaction module from apriori
            2003.08.23 option -q (item sort control) added
            2003.09.12 option -u (sparse representation) added
            2004.08.16 bug concerning option -q0 fixed (min. support)
            2004.11.23 absolute/relative support output changed
            2004.12.09 filter added (binary logarithm of supp. quotient)
            2005.06.20 use of flag for "no item sorting" corrected
            2006.11.26 adapted to new structures ISFMTR and ISEVAL
            2006.11.28 closed and maximal item sets without repository
            2007.02.13 adapted to modified module tabread
            2008.05.02 default limit for maximal number of items removed
            2008.10.26 complete redesign to simplify the structure
            2008.10.28 basic functionality of redesign completed
            2008.10.30 output of transaction ids per item set added
            2008.10.31 closed and maximal item set mining added
            2008.11.13 adapted to changes in transaction management
            2008.12.05 item set reporting order changed (post-order)
            2009.05.28 adapted to modified function tbg_filter()
            2009.10.09 closed/maximal item set check with repository
            2009.10.15 adapted to item set counter in reporter
            2010.03.09 version using transaction ranges added
            2010.03.17 head union tail pruning for maximal sets added
            2010.04.29 bug in memory organization fixed (if n > m)
            2010.06.23 code for transaction id reporting simplified
            2010.07.01 search based on diffsets added (dEclat)
            2010.07.04 bug in tid list setup fixed (after tbg_filter)
            2010.07.09 variant with horizontal processing added
            2010.07.11 filter version of intersection variant added
            2010.07.14 output file made optional (for benchmarking)
            2010.07.15 variant based on an item occurrence table added
            2010.08.05 closedness check based on extensions improved
            2010.08.19 item selection file added as optional input
            2010.08.22 adapted to modified modules tabread and tract
            2010.10.15 adapted to modified interface of module report
            2010.11.24 adapted to modified error reporting (tract)
            2010.11.26 memory handling simplified (list base sizes etc.)
            2010.12.11 adapted to a generic error reporting function
            2010.12.20 adapted to function tbg_icnts() (filter problem)
            2011.03.14 bug in memory allocation in eclat() fixed
            2011.03.19 two range checks for malloc() calls added
            2011.03.20 optional integer transaction weights added
            2011.07.08 adapted to modified function tbg_recode()
            2011.07.27 bug in function eclat_diff() fixed (list length)
            2011.07.29 re-sorting switched off for closed/maximal
            2011.08.15 bit vector version (with masking/reduction) added
            2011.08.16 adapted algorithm variants to finding generators
            2011.08.17 bit vector version modified to use bit map table
            2011.08.28 output of item set counters per size added
            2011.08.31 occurrence deliver version eclat_ocd() added
            2011.09.02 closed/maximal filtering without repo. improved
            2011.09.16 using 16-items machine for trans. ranges added
            2011.09.20 bug in closed/maximal filtering fixed (no repo.)
            2011.09.27 bug in algorithm and mode checking fixed (hut)
            2011.10.01 packing and sorting order for transaction ranges
            2012.04.10 bug in function rec_odfx() fixed (isr_xable())
            2012.05.25 occurrence deliver with item reordering added
            2012.06.13 bug in function rec_odro() fixed (single trans.)
            2012.06.19 function rec_odro() redesigned (delayed 16-items)
            2012.06.20 function fpg_adjust() added (consistency check)
            2012.06.22 use of 16-items machine in rec_odro() improved
            2013.01.24 closed/maximal filtering with vertical database
            2013.02.04 bug in transaction sorting for eclat_trg() fixed
            2013.02.07 check of elim. item support in closed() added
            2013.03.07 direction parameter added to sorting functions
            2013.03.22 adapted to type changes in module tract (SUPP)
            2013.03.26 adapted to type changes in module tract (TID)
            2013.03.28 adapted to type changes in module tract (ITEM)
            2013.10.15 checks of return code of isr_report() added
            2013.10.18 optional pattern spectrum collection added
            2013.10.31 bug in function ecl_adjust fixed (check of mrep)
            2013.11.12 item selection file changed to option -R#
            2013.11.22 bug in function rec_odro() fixed (option -l0)
------------------------------------------------------------------------
  Reference for the Eclat algorithm:
  * M.J. Zaki, S. Parthasarathy, M. Ogihara, and W. Li.
    New Algorithms for Fast Discovery of Association Rules.
    Proc. 3rd Int. Conf. on Knowledge Discovery and Data Mining
    (KDD 1997, Newport Beach, CA), 283-296.
    AAAI Press, Menlo Park, CA, USA 1997
  Reference for the dEclat algorithm (diffsets, option -ad):
  * M.J. Zaki and K. Gouda.
    Fast Vertical Mining Using Diffsets.
    Proc. 9th ACM SIGKDD Int. Conf. on Knowledge Discovery
    and Data Mining (KDD 2003, Washington, DC), 326-335.
    ACM Press, New York, NY, USA 2003
  References for the LCM algorithm (occurrence deliver, option -ao):
  * T. Uno, T. Asai, Y. Uchida, and H. Arimura.
    LCM: An Efficient Algorithm for Enumerating
    Frequent Closed Item Sets.
    Proc. Workshop on Frequent Item Set Mining Implementations
    (FIMI 2003, Melbourne, FL).
    CEUR Workshop Proceedings 90, TU Aachen, Germany 2003
    http://www.ceur-ws.org/Vol-90/
  * T. Uno, M. Kiyomi and H. Arimura.
    LCM ver.2: Efficient Mining Algorithms
    for Frequent/Closed/Maximal Itemsets.
    Proc. Workshop Frequent Item Set Mining Implementations
    (FIMI 2004, Brighton, UK).
    CEUR Workshop Proceedings 126, Aachen, Germany 2004
    http://www.ceur-ws.org/Vol-126/
  * T. Uno, M. Kiyomi, and H. Arimura.
    LCM ver.3: Collaboration of Array, Bitmap and Prefix Tree
    for Frequent Itemset Mining
    Proc. 1st Int. Workshop Open Source Data Mining (OSDM 2005)
    ACM Press, New York, NY, USA 2005
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifndef ISR_PATSPEC
#define ISR_PATSPEC
#endif
#ifdef ECL_MAIN
#ifndef PSP_REPORT
#define PSP_REPORT
#endif
#ifndef TA_READ
#define TA_READ
#endif
#endif
#include "eclat.h"
#include "fim16.h"
#ifdef ECL_MAIN
#include "error.h"
#endif
#ifdef STORAGE
#include "storage.h"
#endif

#define BITMAP_TABLE            /* use a table instead of shifting */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define PRGNAME     "eclat"
#define DESCRIPTION "find frequent item sets with the eclat algorithm"
#define VERSION     "version 4.6 (2014.01.08)         " \
                    "(c) 2002-2014   Christian Borgelt"

/* --- error codes --- */
/* error codes   0 to  -4 defined in tract.h */
#define E_STDIN      (-5)       /* double assignment of stdin */
#define E_OPTION     (-6)       /* unknown option */
#define E_OPTARG     (-7)       /* missing option argument */
#define E_ARGCNT     (-8)       /* too few/many arguments */
#define E_TARGET     (-9)       /* invalid target type */
#define E_SIZE      (-10)       /* invalid item set size */
#define E_SUPPORT   (-11)       /* invalid item set support */
#define E_VARIANT   (-12)       /* invalid algorithm variant */
#define E_MEASURE   (-13)       /* invalid evaluation measure */
/* error codes -15 to -25 defined in tract.h */

#ifndef QUIET                   /* if not quiet version, */
#define MSG         fprintf     /* print messages */
#define XMSG        if (mode & ECL_VERBOSE) fprintf
#else                           /* if quiet version, */
#define MSG(...)                /* suppress messages */
#define XMSG(...)
#endif

#define DIFFSIZE(p,q) ((size_t)((int*)(p)-(int*)(q)) *sizeof(int))
#define SEC_SINCE(t)  ((double)(clock()-(t)) /(double)CLOCKS_PER_SEC)

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- trans. identifier list --- */
  ITEM     item;                /* item identifier (last item in set) */
  SUPP     supp;                /* support of the item (or item set) */
  TID      tids[1];             /* array of transaction identifiers */
} TIDLIST;                      /* (transaction identifier list) */

typedef unsigned int BITBLK;    /* --- bit vector block --- */

typedef struct {                /* --- bit vector --- */
  ITEM     item;                /* item identifier (last item in set) */
  SUPP     supp;                /* support of the item (or item set) */
  BITBLK   bits[1];             /* bit vector over transactions */
} BITVEC;                       /* (bit vector) */

typedef struct {                /* --- transaction id range --- */
  TID      min;                 /* minimum transaction identifier */
  TID      max;                 /* maximum transaction identifier */
  SUPP     wgt;                 /* weight of transactions in range */
} TIDRANGE;                     /* (transaction id range) */

typedef struct {                /* --- transaction range list --- */
  ITEM     item;                /* item identifier (last item in set) */
  SUPP     supp;                /* support of the item (or item set) */
  TIDRANGE trgs[1];             /* array of transaction id ranges */
} TRGLIST;                      /* (transaction id range list) */

typedef struct {                /* --- transaction list --- */
  ITEM     item;                /* item identifier (last item in set) */
  SUPP     supp;                /* support of the item (or item set) */
  TID      cnt;                 /* number of transactions */
  TRACT    *tracts[1];          /* array  of transactions */
} TALIST;                       /* (transaction list) */

typedef struct {                /* --- recursion data --- */
  int      mode;                /* operation mode (e.g. pruning) */
  SUPP     supp;                /* minimum support of an item set */
  ITEM     first;               /* start value for item loops */
  int      dir;                 /* direction   for item loops */
  SUPP     *muls;               /* multiplicity of transactions */
  SUPP     *marks;              /* markers (for item occurrences) */
  ITEM     *cand;               /* to collect candidates (closed()) */
  SUPP     *miss;               /* support still missing (maximal()) */
  BITTA    *btas;               /* array of bit-rep. transactions */
  SUPP     **tab;               /* item occurrence table */
  TRACT    **hash;              /* buffer for hash table */
  TIDLIST  **elim;              /* tra. id lists of eliminated items */
  FIM16    *fim16;              /* 16-items machine */
  TABAG    *tabag;              /* original transaction bag */
  ISREPORT *report;             /* item set reporter */
} RECDATA;                      /* (recursion data) */

typedef TID COMBFN  (TIDLIST *d, TIDLIST *s1, TIDLIST *s2, SUPP *w);
typedef int ECLATFN (TABAG *tabag, int mode, SUPP supp,
                     ISREPORT *report);

/*----------------------------------------------------------------------
  Constants
----------------------------------------------------------------------*/
#if !defined QUIET && defined ECL_MAIN
/* --- error messages --- */
static const char *errmsgs[] = {
  /* E_NONE      0 */  "no error",
  /* E_NOMEM    -1 */  "not enough memory",
  /* E_FOPEN    -2 */  "cannot open file %s",
  /* E_FREAD    -3 */  "read error on file %s",
  /* E_FWRITE   -4 */  "write error on file %s",
  /* E_STDIN    -5 */  "double assignment of standard input",
  /* E_OPTION   -6 */  "unknown option -%c",
  /* E_OPTARG   -7 */  "missing option argument",
  /* E_ARGCNT   -8 */  "wrong number of arguments",
  /* E_TARGET   -9 */  "invalid target type '%c'",
  /* E_SIZE    -10 */  "invalid item set size %"ITEM_FMT,
  /* E_SUPPORT -11 */  "invalid minimum support %g",
  /* E_VARIANT -12 */  "invalid eclat variant '%c'",
  /* E_MEASURE -13 */  "invalid evaluation measure '%c'",
  /*           -14 */  NULL,
  /* E_NOITEMS -15 */  "no (frequent) items found",
  /*           -16 */  "unknown error"
};
#endif

/*----------------------------------------------------------------------
  Global Variables
----------------------------------------------------------------------*/
#ifdef ECL_MAIN
#ifndef QUIET
static CCHAR    *prgname;       /* program name for error messages */
#endif
static TABREAD  *tread  = NULL; /* table/transaction reader */
static ITEMBASE *ibase  = NULL; /* item base */
static TABAG    *tabag  = NULL; /* transaction bag/multiset */
static ISREPORT *report = NULL; /* item set reporter */
static TABWRITE *twrite = NULL; /* table writer for pattern spectrum */
#endif

#ifdef BITMAP_TABLE
static int    bitcnt[256];      /* bit count table */
static BITBLK bitmap[256][256]; /* bit map   table */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions for Debugging
----------------------------------------------------------------------*/
#ifndef NDEBUG

static void indent (int k)
{ while (--k >= 0) printf("   "); }

/*--------------------------------------------------------------------*/

static void show_tid (const char *text, ITEMBASE *base,
                      TIDLIST **lists, ITEM k, int ind)
{                               /* --- show a cond. trans. database */
  ITEM i, j;                    /* item, loop variable */
  TID  *s;                      /* to traverse transaction ids */

  if (text && *text) {          /* print the given text */
    indent(ind); printf("%s\n", text); }
  for (j = 0; j < k; j++) {     /* traverse the items / tid lists */
    indent(ind);                /* indent the output line */
    i = lists[j]->item;         /* print the item name and id */
    if (i < 0) printf("packed  :");
    else       printf("%4s[%2"ITEM_FMT"]:", ib_name(base, i), i);
    for (s = lists[j]->tids; *s >= 0; s++)
      printf(" %"TID_FMT, *s);  /* print the transaction ids */
    printf(" (%"SUPP_FMT")\n", lists[j]->supp);
  }                             /* print the item support */
}  /* show_tid() */

/*--------------------------------------------------------------------*/

static void show_tab (const char *text, ITEMBASE *base,
                      SUPP **tab, TID n, ITEM k)
{                               /* --- show item counter table */
  ITEM i;                       /* loop variable for items */
  TID  r;                       /* loop variable for rows */

  if (text && *text)            /* if it is not empty, */
    printf("%s\n", text);       /* print the given text */
  printf("    ");               /* skip row id/tid column */
  for (r = 0; r < n; r++)       /* print the transaction header */
    printf(" %3"TID_FMT, r);    /* print the row number / tid */
  printf("\n");                 /* terminate the header line */
  for (i = 0; i < k; i++) {     /* traverse the table columns */
    printf("%4s[%2"ITEM_FMT"]:", ib_name(base, i), i);
    for (r = 0; r < n; r++) printf(" %3"SUPP_FMT, tab[i][r]);
    printf("\n");               /* print the item counters */
  }                             /* and terminate the line */
}  /* show_tab() */

/*--------------------------------------------------------------------*/

static void show_trg (const char *text, ITEMBASE *base,
                      TRGLIST **lists, ITEM k, int ind)
{                               /* --- show a cond. trans. database */
  ITEM     i, j;                /* item, loop variable */
  TIDRANGE *r;                  /* to traverse transaction id ranges */

  if (text && *text) {          /* print the given text */
    indent(ind); printf("%s\n", text); }
  for (j = 0; j < k; j++) {     /* traverse the items / range lists */
    indent(ind);                /* indent the output line */
    i = lists[j]->item;         /* get the item identifier */
    r = lists[j]->trgs;         /* and the transaction ranges */
    if (i < 0) {                /* if list for packed items */
      printf("packed:");        /* print special indicator */
      for ( ; r->min >= 0; r++){/* and the transaction ids */
        printf(" %"TID_FMT":%04x", r->min, (unsigned int)r->max);
        printf(":%"SUPP_FMT, r->wgt);
      } }
    else {                      /* if list for a normal item */
      printf("%s[%"ITEM_FMT"]:", ib_name(base, i), i);
      for ( ; r->min >= 0; r++){/* print item name and id */
        printf(" %"TID_FMT"-%"TID_FMT, r->min, r->max);
        printf(":%"SUPP_FMT, r->wgt);
      }                         /* print the transaction ranges */
    }
    printf(" (%"SUPP_FMT")\n", lists[j]->supp);
  }                             /* print the item support */
}  /* show_trg() */

#endif  /* #ifndef NDEBUG */
/*----------------------------------------------------------------------
  Eclat with Transaction Id List Intersection (basic version)
----------------------------------------------------------------------*/

static TID isect (TIDLIST *dst, TIDLIST *src1, TIDLIST *src2,SUPP *muls)
{                               /* --- intersect two trans. id lists */
  TID *s1, *s2, *d;             /* to traverse sources and dest. */

  assert(dst && src1 && src2    /* check the function arguments */
  &&    (src1->tids[0] >= 0) && (src2->tids[0] >= 0) && muls);
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  if (src1->supp > src2->supp) { s2 = src1->tids; s1 = src2->tids; }
  else                         { s1 = src1->tids; s2 = src2->tids; }
  d = dst->tids;                /* get sources and destination */
  while (1) {                   /* trans. id list intersection loop */
    if      (*s1 < *s2) s2++;   /* if one transaction id is larger, */
    else if (*s1 > *s2) s1++;   /* simply skip this transaction id */
    else if (*s1 <   0) break;  /* check for the sentinel */
    else { dst->supp += muls[*d++ = *s1++]; s2++; }
  }                             /* copy equal elements to destination */
  *d++ = (TID)-1;               /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new list */
}  /* isect() */

/*--------------------------------------------------------------------*/

static int rec_base (TIDLIST **lists, ITEM k, size_t x, RECDATA *rd)
{                               /* --- eclat recursion with tid lists */
  int     r;                    /* error status */
  ITEM    i, m, z;              /* loop variables */
  SUPP    pex;                  /* minimum support for perfect exts. */
  TIDLIST *l, *d;               /* to traverse transaction id lists */
  TIDLIST **proj = NULL;        /* trans. id lists of proj. database */
  TID     *p;                   /* to organize the trans. id lists */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for conditional databases) */
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / tid lists */
    r = isr_add(rd->report, l->item, l->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TIDLIST*)(p = (TID*)(proj +k+1));
      for (i = 0; i < k; i++) { /* intersect with preceding lists */
        x = (size_t)isect(d, lists[i], l, rd->muls);
        if      (d->supp >= pex)      /* collect perfect extensions */
          isr_addpex(rd->report, d->item);
        else if (d->supp >= rd->supp) /* collect frequent extensions */
          proj[++m] = d = (TIDLIST*)(p = d->tids +x);
      }                         /* switch to the next output list */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_base(proj, m, DIFFSIZE(p,proj[0]), rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_report(rd->report) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_base() */

/*--------------------------------------------------------------------*/

int eclat_base (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with trans. id lists */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k, m;           /* loop variable, number of items */
  TID        n;                 /* number of transactions */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  SUPP       pex;               /* minimum support for perfect exts. */
  TRACT      *t;                /* to traverse transactions */
  TIDLIST    **lists, *l;       /* to traverse transaction id lists */
  TID        *tids, *p, **next; /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  c = tbg_icnts(tabag, 0);      /* get the number of containing */
  if (!c) return -1;            /* transactions per item */
  lists = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*)
                           +(size_t)k *sizeof(TID*)
                           +(size_t)n *sizeof(SUPP));
  if (!lists) return -1;        /* create initial tid list array */
  next    = (TID**)(lists+k);   /* and split off next position array */
  rd.muls = (SUPP*)(next +k);   /* and transaction multiplicity array */
  x = tbg_extent(tabag);        /* allocate the tid list elements */
  p = tids = (TID*)malloc((size_t)k *sizeof(TIDLIST) +x *sizeof(TID));
  if (!p) { free(lists); return -1; } /* allocate tid list elements */
  for (i = 0; i < k; i++) {     /* traverse the items / tid lists */
    lists[i] = l = (TIDLIST*)p; /* get next transaction id list */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    next[i]  = p = l->tids;     /* note position of next trans. id */
    p += c[i]; *p++ = (TID)-1;  /* skip space for transaction ids */
  }                             /* and store a sentinel at the end */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(tabag, --n);  /* get the next transaction */
    rd.muls[n] = w = ta_wgt(t); /* and store its weight */
    for (s = ta_items(t); *s > TA_END; s++) {
      lists[*s]->supp += w;     /* traverse the transaction's items */
      *next[*s]++      = n;     /* sum the transaction weight and */
    }                           /* collect the transaction ids */
  }
  for (i = m = 0; i < k; i++) { /* traverse the items / tid lists */
    l = lists[i];               /* eliminate all infrequent items and */
    if (l->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (l->supp >= pex) { isr_addpex(report, i); continue; }
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_base(lists, m, DIFFSIZE(p,tids), &rd);
  }                             /* find freq. items sets recursively */
  if (r >= 0)                   /* report the empty item set */
    r = (isr_report(report) < 0) ? -1 : 0;
  free(tids); free(lists);      /* delete the allocated arrays */
  return r;                     /* return the error status */
}  /* eclat_base() */

/*----------------------------------------------------------------------
  Eclat with Transaction Id List Intersection (optimized version)
----------------------------------------------------------------------*/

static int tid_cmp (const void *a, const void *b, void *data)
{                               /* --- compare support of tid lists */
  if (((TIDLIST*)b)->supp > ((TIDLIST*)a)->supp) return  1;
  if (((TIDLIST*)b)->supp < ((TIDLIST*)a)->supp) return -1;
  return 0;                     /* return sign of support difference */
}  /* tid_cmp() */

/*--------------------------------------------------------------------*/

static int tid_cmpx (const void *a, const void *b, void *data)
{                               /* --- compare support of tid lists */
  if (((TIDLIST*)a)->item < 0) return -1;
  if (((TIDLIST*)b)->item < 0) return +1;
  if (((TIDLIST*)b)->supp > ((TIDLIST*)a)->supp) return  1;
  if (((TIDLIST*)b)->supp < ((TIDLIST*)a)->supp) return -1;
  return 0;                     /* return sign of support difference */
}  /* tid_cmpx() */

/*--------------------------------------------------------------------*/

static TID filter (TIDLIST *dst, TIDLIST *src, SUPP *muls)
{                               /* --- filter a transaction id list */
  SUPP m;                       /* multiplicity of transaction */
  TID  *s, *d;                  /* to traverse source and dest. */

  assert(dst && src && muls);   /* check the function arguments */
  dst->item = src->item;        /* copy first item and init. support */
  dst->supp = 0;                /* traverse the source trans. id list */
  for (d = dst->tids, s = src->tids; *s >= 0; s++)
    if ((m = muls[*s]) > 0) {   /* collect the marked trans. ids and */
      dst->supp += m; *d++ = *s; }    /* sum the transaction weights */
  *d++ = (TID)-1;               /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new list */
}  /* filter() */

/*--------------------------------------------------------------------*/

static int closed (TIDLIST *list, RECDATA *rd, ITEM n)
{                               /* --- check for a closed item set */
  TIDLIST    *elim;             /* to traverse eliminated items */
  const ITEM *p;                /* to traverse transaction items */
  TID        *s, *d;            /* to traverse transaction ids */
  ITEM       *t, *r;            /* to traverse items */

  assert(list && rd);           /* check the function arguments */
  if (rd->mode & ECL_VERTICAL){ /* if to use vertical representation */
    while (--n >= 0) {          /* traverse the eliminated items */
      elim = rd->elim[n];       /* skip items with lower support */
      if (elim->supp < list->supp) continue;
      s = list->tids; d = elim->tids;
      while (1) {               /* test for a perfect extension */
        if      (*s < *d) d++;  /* skip missing destination id */
        else if (*s > *d) break;/* if source id is missing, abort */
        else if (*s <  0) return 0;
        else { s++; d++; }      /* check for the sentinel and */
      }                         /* skip matching transaction ids */
    }                           /* (all tids found: perfect ext.) */
    return -1; }                /* return 'item set is closed' */
  else {                        /* if to use horiz. representation */
    p = ta_items(tbg_tract(rd->tabag, list->tids[0]));
    for (r = rd->cand; *p > list->item; p++)
      if (!isr_uses(rd->report, *p))
        *r++ = *p;              /* collect items from a transaction */
    if (r <= rd->cand) return -1;
    *r = TA_END;                /* store a sentinel at the end */
    for (s = list->tids+1; *s >= 0; s++) {
      t = r = rd->cand;         /* traverse the transaction ids */
      p = ta_items(tbg_tract(rd->tabag, *s));
      while (1) {               /* item list intersection loop */
        if      (*t < *p) p++;  /* if one item id is larger, */
        else if (*t > *p) t++;  /* simply skip this item id, */
        else if (*t <  0) break;/* check for the list sentinel */
        else { *r++ = *t++; p++; }
      }                         /* (collect perfect ext. candidates) */
      if (r <= rd->cand) return -1;
      *r = TA_END;              /* if intersection is empty, abort, */
    }                           /* otherwise store a sentinel */
    return 0;                   /* return 'item set is not closed' */
  }
}  /* closed() */

/*--------------------------------------------------------------------*/

static int maximal (TIDLIST *list, RECDATA *rd, ITEM n)
{                               /* --- check for a maximal item set */
  ITEM       i;                 /* loop variable for items */
  SUPP       w;                 /* weight/support buffer */
  const ITEM *p;                /* to traverse transaction items */
  TID        *s, *d;            /* to traverse sources and dest. */

  assert(list && rd);           /* check the function arguments */
  if (rd->mode & ECL_VERTICAL){ /* if to use vertical representation */
    while (--n >= 0) {          /* traverse the eliminated items */
      s = list->tids; d = rd->elim[n]->tids;
      for (w = 0; 1; ) {        /* test for a perfect extension */
        if      (*s < *d) d++;  /* if one transaction id is larger, */
        else if (*s > *d) s++;  /* skip missing destination id */
        else if (*s <  0) break;/* check for the sentinel and */
        else { w += rd->muls[*s++]; d++; }
      }                         /* sum weights of matching trans. ids */
      if (w >= rd->supp) return 0;
    } }                         /* check for a frequent extension */
  else {                        /* if to use horiz. representation */
    for (i = tbg_itemcnt(rd->tabag); --i > list->item; )
      rd->miss[i] = (isr_uses(rd->report, i)) ? list->supp+1 : rd->supp;
    for (s = list->tids; *s >= 0; s++) {
      w = rd->muls[*s];         /* traverse the transactions */
      for (p = ta_items(tbg_tract(rd->tabag, *s)); *p > list->item; p++)
        if ((rd->miss[*p] -= w) <= 0) return 0;
    }                           /* count support of candidate exts.; */
  }                             /* if frequent cand. found, abort */
  return -1;                    /* return 'set is maximal' */
}  /* maximal() */

/*--------------------------------------------------------------------*/

static int rec_tcm (TIDLIST **lists, ITEM k, size_t x, ITEM e,
                    RECDATA *rd)
{                               /* --- eclat recursion with tid lists */
  int     r;                    /* error status */
  ITEM    i, m, z;              /* loop variables */
  SUPP    max;                  /* maximum support of an ext. item */
  SUPP    pex;                  /* minimum support for perfect exts. */
  TIDLIST **proj = NULL;        /* trans. id lists of proj. database */
  TIDLIST *l, *d;               /* to traverse transaction id lists */
  TID     *p;                   /* to traverse transaction ids */
  ITEM    *t;                   /* to collect the tail items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = lists[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for conditional databases) */
  if ((k > 4)                   /* if there are enough items left, */
  &&  (rd->mode & ECL_REORDER)) /* re-sort the items w.r.t. support */
    ptr_qsort(lists, (size_t)k, 1, (rd->fim16) ?tid_cmpx:tid_cmp, NULL);
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / tid lists */
    if (!closed(l, rd, e))      /* if the current set is not closed, */
      continue;                 /* the item need not be processed */
    r = isr_addnc(rd->report, l->item, l->supp);
    if (r < 0) break;           /* add current item to the reporter */
    max = 0;                    /* init. maximal extension support */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TIDLIST*)(proj +k+1);
      if (k < 2) {              /* if there are only few items left */
        /* Benchmark tests showed that this version is faster only */
        /* if there is only one other tid list to intersect with.  */
        if (lists[i = 0]->item < 0) { /* if there are packed items */
          x = (size_t)isect(d, lists[i++], l, rd->muls);
          if (d->supp >= rd->supp) {  /* if they are frequent */
            proj[++m] = d = (TIDLIST*)(d->tids +x); }
        }                       /* add a tid list for packed items */
        for ( ; i < k; i++) {   /* traverse the preceding lists */
          x = (size_t)isect(d, lists[i], l, rd->muls);
          if (d->supp < rd->supp) /* intersect transaction id lists */
            continue;           /* eliminate infrequent items */
          if (d->supp >= pex) { /* collect perfect extensions */
            isr_addpex(rd->report, d->item); continue; }
          if (d->supp > max)    /* find maximal extension support */
            max = d->supp;      /* (for later closed/maximal check) */
          proj[++m] = d = (TIDLIST*)(d->tids +x);
        } }                     /* collect tid lists of freq. items */
      else {                    /* if there are many items left */
        for (p = l->tids; *p >= 0; p++) /* mark transaction ids */
          rd->marks[*p] = rd->muls[*p]; /* in the current list */
        if (lists[i = 0]->item < 0) {   /* if there are packed items */
          x = (size_t)filter(d, lists[i++], rd->marks);
          if (d->supp >= rd->supp) {    /* if they are frequent */
            proj[++m] = d = (TIDLIST*)(d->tids +x); }
        }                       /* add a tid list for packed items */
        for ( ; i < k; i++) {   /* traverse the preceding lists */
          x = (size_t)filter(d, lists[i], rd->marks);
          if (d->supp < rd->supp) /* intersect transaction id lists */
            continue;           /* eliminate infrequent items */
          if (d->supp >= pex) { /* collect perfect extensions */
            isr_addpex(rd->report, d->item); continue; }
          if (d->supp > max)    /* find maximal extension support */
            max = d->supp;      /* (for later closed/maximal check) */
          proj[++m] = d = (TIDLIST*)(d->tids +x);
        }                       /* collect tid lists of freq. items */
        for (p = l->tids; *p >= 0; p++)
          rd->marks[*p] = 0;    /* unmark transaction ids */
      }                         /* in the current list */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_tcm(proj, m, DIFFSIZE(d,proj[0]), e, rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }                           /* (or rather their trans. id lists) */
    if ((rd->mode & ISR_CLOSED) ? (max < l->supp)
    :   ((max < rd->supp) ? maximal(l, rd, e) : 0)) {
      if (isr_reportx(rd->report, l->tids, (diff_t)-l->supp) < 0) {
        r = -1; break; }        /* report the current item set */
    }                           /* and check for an error */
    isr_remove(rd->report, 1);  /* remove the current item and */
    if (rd->mode & ECL_VERTICAL)/* collect the eliminated items */
      rd->elim[e++] = l;        /* (for closed/maximal check) */
  }
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_tcm() */

/*--------------------------------------------------------------------*/

static int rec_tid (TIDLIST **lists, ITEM k, size_t x, RECDATA *rd)
{                               /* --- eclat recursion with tid lists */
  int     r;                    /* error status */
  ITEM    i, m, z;              /* loop variables, error status */
  SUPP    max;                  /* maximum support of an ext. item */
  SUPP    pex;                  /* minimum support for perfect exts. */
  TIDLIST **proj = NULL;        /* trans. id lists of proj. database */
  TIDLIST *l, *d;               /* to traverse transaction id lists */
  TID     *p;                   /* to traverse transaction ids */
  ITEM    *t;                   /* to collect the tail items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = lists[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for conditional databases) */
  if ((k > 4)                   /* if there are enough items left, */
  &&  (rd->mode & ECL_REORDER)) /* re-sort the items w.r.t. support */
    ptr_qsort(lists, (size_t)k, 1, (rd->fim16) ?tid_cmpx:tid_cmp, NULL);
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / tid lists */
    if (l->item < 0) {          /* if this list is for packed items */
      for (p = l->tids; *p >= 0; p++)
        m16_add(rd->fim16, rd->btas[*p], rd->muls[*p]);
      r = m16_mine(rd->fim16);  /* add bit-rep. transaction prefixes */
      if (r >= 0) continue;     /* to the 16-items machine and mine, */
      if (proj) free(proj);     /* then go to the next trans. id list */
      return r;                 /* otherwise free allocated memory */
    }                           /* and abort the function */
    r = isr_add(rd->report, l->item, l->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    max = 0;                    /* init. maximal extension support */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TIDLIST*)(proj +k+1);
      if (k < 2) {              /* if there are only few items left */
        /* Benchmark tests showed that this version is faster only */
        /* if there is only one other tid list to intersect with.  */
        if (lists[i = 0]->item < 0) { /* if there are packed items */
          x = (size_t)isect(d, lists[i++], l, rd->muls);
          if (d->supp >= rd->supp) {  /* if they are frequent */
            proj[++m] = d = (TIDLIST*)(d->tids +x); }
        }                       /* add a tid list for packed items */
        for ( ; i < k; i++) {   /* traverse the preceding lists */
          x = (size_t)isect(d, lists[i], l, rd->muls);
          if (d->supp < rd->supp) /* intersect transaction id lists */
            continue;           /* eliminate infrequent items */
          if (d->supp >= pex) { /* collect perfect extensions */
            isr_addpex(rd->report, d->item); continue; }
          if (d->supp > max)    /* find maximal extension support */
            max = d->supp;      /* (for later closed/maximal check) */
          proj[++m] = d = (TIDLIST*)(d->tids +x);
        } }                     /* collect tid lists of freq. items */
      else {                    /* if there are many items left */
        for (p = l->tids; *p >= 0; p++) /* mark transaction ids */
          rd->marks[*p] = rd->muls[*p]; /* in the current list */
        if (lists[i = 0]->item < 0) {   /* if there are packed items */
          x = (size_t)filter(d, lists[i++], rd->marks);
          if (d->supp >= rd->supp) {    /* if they are frequent */
            proj[++m] = d = (TIDLIST*)(d->tids +x); }
        }                       /* add a tid list for packed items */
        for ( ; i < k; i++) {   /* traverse the preceding lists */
          x = (size_t)filter(d, lists[i], rd->marks);
          if (d->supp < rd->supp) /* intersect transaction id lists */
            continue;           /* eliminate infrequent items */
          if (d->supp >= pex) { /* collect perfect extensions */
            isr_addpex(rd->report, d->item); continue; }
          if (d->supp > max)    /* find maximal extension support */
            max = d->supp;      /* (for later closed/maximal check) */
          proj[++m] = d = (TIDLIST*)(d->tids +x);
        }                       /* collect tid lists of freq. items */
        for (p = l->tids; *p >= 0; p++)
          rd->marks[*p] = 0;    /* unmark transaction ids */
      }                         /* in the current list */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_tid(proj, m, DIFFSIZE(d,proj[0]), rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_reportx(rd->report, l->tids, (diff_t)-l->supp) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_tid() */

/*--------------------------------------------------------------------*/

int eclat_tid (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with trans. id lists */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k, m, e;        /* loop variable, number of items */
  TID        n;                 /* number of transactions */
  size_t     x, z;              /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  SUPP       max;               /* maximum support of an item */
  SUPP       pex;               /* minimum support for perfect exts. */
  TRACT      *t;                /* to traverse transactions */
  TIDLIST    **lists, *l;       /* to traverse transaction id lists */
  TID        *tids, *p, **next; /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (tbg_packcnt(tabag) <= 0) mode &= ~ECL_FIM16;
  if (!(mode & (ISR_CLOSED|ISR_MAXIMAL))) mode &= ~ISR_NOFILTER;
  if (!(mode & ISR_MAXIMAL) || (mode & (ISR_NOFILTER|ECL_FIM16)))
    mode &= ~ECL_TAIL;          /* make search mode consistent */
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  k = tbg_itemcnt(tabag);       /* get and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  c = tbg_icnts(tabag, 0);      /* and the number of containing */
  if (!c) return -1;            /* transactions per item */
  e = ((mode & ISR_NOFILTER) &&  (mode & ECL_VERTICAL)) ? k   : 0;
  m = ((mode & ISR_NOFILTER) && !(mode & ECL_VERTICAL)) ? k+1 : 0;
  x =  (mode & ECL_FIM16) ? (size_t)n *sizeof(BITTA) : 0;
  z = (sizeof(ITEM) > sizeof(SUPP)) ? sizeof(ITEM) : sizeof(SUPP);
  lists = (TIDLIST**)malloc((size_t)(k+e) *sizeof(TIDLIST*)
                           +(size_t) k    *sizeof(TID*)
                           +(size_t)(n+n) *sizeof(SUPP)
                           +(size_t) m    *z +x);
  if (!lists) return -1;        /* create initial tid list array and */
  rd.elim  = lists +k;          /* split off the additional arrays */
  next     = (TID**) (rd.elim +e);
  rd.muls  = (SUPP*) (next    +k);
  rd.miss  = (SUPP*) (rd.muls +n); /* buffer for maximal() */
  rd.cand  = (ITEM*)  rd.miss;     /* buffer for closed() */
  rd.marks = (sizeof(ITEM) > sizeof(SUPP))
           ? (SUPP*)(rd.cand+m) : rd.miss+m;
  rd.btas  = (BITTA*)(rd.marks+n);
  memset(rd.marks, 0, (size_t)n *sizeof(TID));
  for (x = 0, i = 0; i < k; i++)/* traverse the items and sum */
    x += (size_t)c[i];          /* the number of item occurrences */
  /* Do not use tbg_extent(), because it does not take packed items */
  /* properly into account and thus may yield too big a value.      */
  if (x < (size_t)n) x = (size_t)n; /* ensure enough transaction ids */
  p = tids = (TID*)malloc((size_t)k *sizeof(TIDLIST) +x *sizeof(TID));
  if (!p) { free(lists); return -1; } /* allocate tid list elements */
  for (i = 0; i < k; i++) {     /* traverse the items / tid lists */
    lists[i] = l = (TIDLIST*)p; /* get next transaction id list */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    next[i]  = p = l->tids;     /* note position of next trans. id */
    p += c[i]; *p++ = (TID)-1;  /* skip space for transaction ids */
  }                             /* and store a sentinel at the end */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(tabag, --n);  /* get the next transaction */
    rd.muls[n] = w = ta_wgt(t); /* and store its weight */
    for (s = ta_items(t); *s > TA_END; s++) {
      if ((i = *s) < 0) {       /* traverse the transaction's items */
        rd.btas[n] = (BITTA)i; i = 0; }
      lists[i]->supp += w;      /* traverse the transaction's items */
      *next[i]++      = n;      /* sum the transaction weight and */
    }                           /* collect the transaction ids */
  }
  rd.fim16 = NULL;              /* default: no 16-items machine */
  l = lists[i = 0];             /* get the list for packed items */
  if ((mode & ECL_FIM16)        /* if to use a 16-items machine */
  &&  (l->supp >= rd.supp)) {   /* and there are packed items */
    rd.fim16 = m16_create(rd.dir, rd.supp, report);
    if (!rd.fim16) { free(tids); free(lists); return -1; }
    l->item = -1; i = 1;        /* mark list for the packed items */
  }                             /* and add it to the reduced array */
  max = 0;                      /* init. the maximal item support */
  for (m = i; i < k; i++) {     /* traverse the items / tid lists */
    l = lists[i];               /* eliminate all infrequent items and */
    if (l->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (l->supp >= pex) { isr_addpex(report, i); continue; }
    if (l->supp >  max)         /* find the maximal item support */
      max = l->supp;            /* (for later closed/maximal check) */
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = (mode & ISR_NOFILTER)   /* with direct or repo. filtering */
      ? rec_tcm(lists, m, DIFFSIZE(p,tids), 0, &rd)
      : rec_tid(lists, m, DIFFSIZE(p,tids), &rd);
  }                             /* find freq. item sets recursively */
  if (r >= 0) {                 /* if no error occurred */
    i = mode & (ISR_CLOSED|ISR_MAXIMAL);
    w = (i & ISR_MAXIMAL) ? rd.supp : tbg_wgt(tabag);
    if (!i || (max < w)) {      /* if to report the empty set */
      if (!isr_tidfile(report)) /* if not to report transaction ids, */
        r = (isr_report(report) < 0) ? -1 : 0;   /* report empty set */
      else {                    /* if to report transaction ids */
        for (n = tbg_cnt(tabag); n > 0; n--) tids[n] = n;
        r = (isr_reportx(report, tids, (diff_t)n) < 0) ? -1 : 0;
      }                         /* report the empty item set */
    }                           /* with all transaction ids */
  }
  if (rd.fim16)                 /* if a 16-items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16-items machine */
  free(tids); free(lists);      /* delete the allocated arrays */
  return r;                     /* return the error status */
}  /* eclat_tid() */

/*----------------------------------------------------------------------
  Eclat with Bit Vectors
----------------------------------------------------------------------*/

static int bit_cmp (const void *a, const void *b, void *data)
{                               /* --- compare support of tid lists */
  if (((BITVEC*)b)->supp > ((BITVEC*)a)->supp) return  1;
  if (((BITVEC*)b)->supp < ((BITVEC*)a)->supp) return -1;
  return 0;                     /* return sign of support difference */
}  /* bit_cmp() */

/*--------------------------------------------------------------------*/
#ifdef BITMAP_TABLE

static void bit_init (void)
{                               /* --- init. bit count/map tables */
  int i, k, b;                  /* loop variables, bit index */

  if (bitcnt[1] != 0) return;   /* check for an initialized table */
  for (i = 0; ++i < 256; )      /* traverse all byte values */
    for (k = i; k; k >>= 1)     /* traverse the bits in the value */
      bitcnt[i] += k & 1;       /* store their number in the table */
  memset(bitmap[0], 0, sizeof(bitmap[0]));
  for (k = 0; k < 256; ) { bitmap[1][k++] = 0; bitmap[1][k++] = 1; }
  for (i = 1; ++i < 255; ) {    /* traverse the matrix rows (masks) */
    for (b = 8; --b >= 0; ) {   /* traverse set bits of the mask */
      if (((i >> b) & 1) == 0) continue;
      for (k = 0; k < 256; ) {  /* traverse the matrix columns */
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
        bitmap[i][k] = (bitmap[i][k] << 1) | ((k >> b) & 1); k++;
      }                         /* collect the bits of the source */
    }                           /* that are under the mask bits */
  }                             /* for faster bit vector reduction */
  for (k = 0; k < 256; k++) bitmap[255][k] = (BITBLK)k;
}  /* bit_init() */

/*--------------------------------------------------------------------*/

static void bit_isect (BITVEC *dst, BITVEC *src1, BITVEC *src2, TID n)
{                               /* --- intersect two bit vectors */
  BITBLK *s1, *s2, *d;          /* to traverse sources and dest. */
  BITBLK s, m, o, x;            /* source, mask, and output blocks */
  int    b, c;                  /* number of bits in output */

  assert(dst && src1 && src2);  /* check the function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  d = dst->bits; s1 = src1->bits; s2 = src2->bits;
  for (o = 0, b = 0; n > 0; n--) { /* traverse the bit vector blocks */
    s = *s1++; m = *s2++;       /* traverse the bytes of each block */
    for ( ; m != 0; s >>= 8, m >>= 8) {
      dst->supp += (SUPP)bitcnt[x = bitmap[m & 0xff][s & 0xff]];
      o |= x << b;         b += c = bitcnt[m & 0xff];
      if (b < 32) continue;     /* add output bits for current byte */
      b -= 32; *d++ = o;        /* if a bit block is full, store it */
      o = x >> (c-b-1) >> 1;    /* store remaining bits in buffer, */
    }                           /* but note that x >> 32 == x >> 0, */
  }                             /* so simply o = x >> (c-b) fails */
  if (b > 0) *d = o;            /* store the last bit vector block */
}  /* bit_isect() */

/*--------------------------------------------------------------------*/
#else

#define bit_init()              /* no initialization needed */

static void bit_isect (BITVEC *dst, BITVEC *src1, BITVEC *src2, TID n)
{                               /* --- intersect two bit vectors */
  BITBLK *s1, *s2, *d;          /* to traverse sources and dest. */
  BITBLK s, m, o;               /* source, mask, and output block */
  int    b;                     /* number of bits in output */

  assert(dst && src1 && src2);  /* check the function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  d = dst->bits; s1 = src1->bits; s2 = src2->bits;
  for (o = 0, b = 0; n > 0; n--) { /* traverse the bit vector blocks */
    for (s = *s1++, m = *s2++; m != 0; ) {
      if (m & 1) {              /* if first of four mask bits is set */
        dst->supp += (SUPP)(s & 1); o |= (s & 1) << b;
        if (++b >= 32) { *d++ = o; o = 0; b = 0; }
      }                         /* copy the source bit to the output */
      s >>= 1; m >>= 1;         /* get the next source and mask bit */
      if (m & 1) {              /* if first of four mask bits is set */
        dst->supp += (SUPP)(s & 1); o |= (s & 1) << b;
        if (++b >= 32) { *d++ = o; o = 0; b = 0; }
      }                         /* copy the source bit to the output */
      s >>= 1; m >>= 1;         /* get the next source and mask bit */
      if (m & 1) {              /* if first of four mask bits is set */
        dst->supp += (SUPP)(s & 1); o |= (s & 1) << b;
        if (++b >= 32) { *d++ = o; o = 0; b = 0; }
      }                         /* copy the source bit to the output */
      s >>= 1; m >>= 1;         /* get the next source and mask bit */
      if (m & 1) {              /* if first of four mask bits is set */
        dst->supp += (SUPP)(s & 1); o |= (s & 1) << b;
        if (++b >= 32) { *d++ = o; o = 0; b = 0; }
      }                         /* copy the source bit to the output */
      s >>= 1; m >>= 1;         /* get the next source and mask bit */
    }                           /* collect the source bits */
  }                             /* for which a mask bit is set */
  if (b > 0) *d = o;            /* store the last bit block */
}  /* bit_isect() */

#endif
/*--------------------------------------------------------------------*/

static int rec_bit (BITVEC **vecs, ITEM k, TID n, RECDATA *rd)
{                               /* --- eclat recursion with bit vecs. */
  int    r;                     /* error status */
  ITEM   i, m, z;               /* loop variables */
  SUPP   pex;                   /* minimum support for perf. exts. */
  TID    len;                   /* length of (reduced) bit vectors */
  BITVEC **proj = NULL;         /* bit vectors of projected database */
  BITVEC *v, *d;                /* to traverse bit vectors */
  BITBLK *p;                    /* to traverse bit vector blocks */
  ITEM   *t;                    /* to collect the tail items */

  assert(vecs && (k > 0) && rd);/* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = vecs[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (BITVEC**)malloc((size_t)k                *sizeof(BITVEC*)
                          + (size_t)k                *sizeof(BITVEC)
                          +((size_t)k*(size_t)(n-1)) *sizeof(BITBLK));
    if (!proj) return -1;       /* allocate bit vectors and array */
  }                             /* (memory for conditional databases) */
  if ((k > 4)                   /* if there are enough items left, */
  &&  (rd->mode & ECL_REORDER)) /* re-sort the items w.r.t. support */
    ptr_qsort(vecs, (size_t)k, +1, bit_cmp, NULL);
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    v = vecs[k];                /* traverse the remaining items */
    r = isr_add(rd->report, v->item, v->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (proj && (k > 0)) {      /* if another item can be added */
      len = (TID)(v->supp+31) >> 5;    /* get new vector length */
      pex = (rd->mode & ECL_PERFECT) ? v->supp : SUPP_MAX;
      proj[m = 0] = d = (BITVEC*)(p = (BITBLK*)(proj +k+1));
      for (i = 0; i < k; i++) { /* traverse preceding vectors */
        bit_isect(d, vecs[i], v, n);
        if (d->supp < rd->supp) /* intersect transaction bit vectors */
          continue;             /* eliminate infrequent items */
        if (d->supp >= pex) {   /* collect perfect extensions */
          isr_addpex(rd->report, d->item); continue; }
        proj[++m] = d = (BITVEC*)(p = d->bits +len);
      }                         /* collect the remaining bit vectors */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_bit(proj, m, len, rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_report(rd->report) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete bit vectors and array */
  return r;                     /* return the error status */
}  /* rec_bit() */

/*--------------------------------------------------------------------*/

int eclat_bit (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with bit vectors */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k, m;           /* loop variable, number of items */
  TID        n;                 /* number of transactions */
  TID        x;                 /* number of item instances */
  SUPP       pex;               /* minimum support for perfect exts. */
  TRACT      *t;                /* to traverse transactions */
  BITVEC     **vecs, *v;        /* to traverse bit vectors */
  BITBLK     *p;                /* to traverse bit vector blocks */
  const ITEM *s;                /* to traverse transaction items */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  bit_init();                   /* initialize the bit count table */
  x = (n + 31) >> 5;            /* and compute the bit vector size */
  vecs = (BITVEC**)malloc((size_t)k                *sizeof(BITVEC*)
                        + (size_t)k                *sizeof(BITVEC)
                        +((size_t)k*(size_t)(x-1)) *sizeof(BITBLK));
  if (!vecs) return -1;         /* create initial bit vector array */
  p = (BITBLK*)(vecs+k);        /* and get the bit vector memory */
  for (i = 0; i < k; i++) {     /* traverse the items / bit vectors */
    vecs[i] = v = (BITVEC*)p;   /* get/create the next bit vector */
    v->item = i;                /* initialize the bit vector item */
    v->supp = 0;                /* and the support counter */
    memset(v->bits, 0, (size_t)x *sizeof(BITBLK));
    p = v->bits +x;             /* clear all transaction bits and */
  }                             /* skip them to get the next vector */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(tabag, --n);  /* retrieve the next transaction */
    assert(ta_wgt(t) == 1);     /* transaction weight must be 1 */
    for (s = ta_items(t); *s > TA_END; s++) {
      v = vecs[*s];             /* traverse the transaction's items */
      v->supp += 1;             /* sum/count the transaction weight */
      v->bits[n >> 5] |= (BITBLK)(1 << (n & 0x1f));
    }                           /* set the bit for the current trans. */
  }                             /* to indicate that item is contained */
  for (i = m = 0; i < k; i++) { /* traverse the items / bit vectors */
    v = vecs[i];                /* eliminate all infrequent items and */
    if (v->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (v->supp >= pex) { isr_addpex(report, i); continue; }
    vecs[m++] = v;              /* collect vectors for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_bit(vecs, m, x, &rd);
  }                             /* find freq. items sets recursively */
  if (r >= 0)                   /* report the empty item set */
    r = (isr_report(report) < 0) ? -1 : 0;
  free(vecs);                   /* delete the allocated bit vectors */
  return r;                     /* return the error status */
}  /* eclat_bit() */

/*----------------------------------------------------------------------
  Eclat with an Occurrence Indicator Table
----------------------------------------------------------------------*/

static int rec_tab (TIDLIST **lists, ITEM k, size_t x, RECDATA *rd)
{                               /* --- eclat recursion with table */
  int     r;                    /* error status */
  ITEM    i, m, z;              /* loop variables */
  SUPP    pex;                  /* minimum support for perfect exts. */
  TIDLIST *l, *d;               /* to traverse transaction id lists */
  TIDLIST **proj = NULL;        /* trans. id lists of proj. database */
  ITEM    *t;                   /* to collect the tail items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = lists[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for projected database) */
  if ((k > 4)                   /* if there are enough items left, */
  &&  (rd->mode & ECL_REORDER)) /* re-sort the items w.r.t. support */
    ptr_qsort(lists, (size_t)k, +1, tid_cmp, NULL);
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / tid lists */
    r = isr_add(rd->report, l->item, l->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TIDLIST*)(proj +k+1);
      for (i = 0; i < k; i++) { /* traverse the preceding lists */
        x = (size_t)filter(d, lists[i], rd->tab[l->item]);
        if (d->supp < rd->supp) /* filter transaction id list */
          continue;             /* eliminate infrequent items */
        if (d->supp >= pex) {   /* collect perfect extensions */
          isr_addpex(rd->report, d->item); continue; }
        proj[++m] = d = (TIDLIST*)(d->tids +x);
      }                         /* collect tid lists of freq. items */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_tab(proj, m, DIFFSIZE(d,proj[0]), rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_reportx(rd->report, l->tids, (diff_t)-l->supp) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_tab() */

/*--------------------------------------------------------------------*/

int eclat_tab (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with occurrence table */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k, m;           /* loop variable, number of items */
  TID        n;                 /* number of transactions */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  SUPP       max;               /* maximum support of an item */
  SUPP       pex;               /* minimum support for perfect exts. */
  SUPP       *d;                /* to traverse occurrence table rows */
  TRACT      *t;                /* to traverse transactions */
  TIDLIST    **lists, *l;       /* to traverse transaction id lists */
  TID        *tids, *p, **next; /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  x = tbg_extent(tabag);        /* get the number of item instances */
  c = tbg_icnts(tabag, 0);      /* and the number of containing */
  if (!c) return -1;            /* transactions per item */
  if ((SIZE_MAX/sizeof(SUPP) -x) / (size_t)(n+4) < (size_t)k)
    return -1;                  /* check the table size */
  lists = (TIDLIST**)malloc((size_t) k              *sizeof(TIDLIST*)
                           +(size_t) k              *sizeof(TID*)
                           +(size_t) k              *sizeof(SUPP*)
                           +(size_t)(k+1)*(size_t)n *sizeof(SUPP));
  if (!lists) return -1;        /* create initial tid list array */
  next     = (TID**) (lists +k);/* and split off arrays */
  rd.tab   = (SUPP**)(next  +k);/* get item occ. table header */
  rd.muls  = (SUPP*) (rd.tab+k);/* split off trans. weight array */
  d = (SUPP*)memset(rd.muls +n, 0, (size_t)k*(size_t)n *sizeof(SUPP));
  if (x < (size_t)n) x = (size_t)n; /* ensure enough transaction ids */
  p = tids = (TID*)malloc((size_t)k *sizeof(TIDLIST) +x *sizeof(TID));
  if (!p) { free(lists); return -1; }
  for (i = 0; i < k; i++) {     /* traverse the items / tid lists */
    rd.tab[i] = d; d += n;      /* organize the table rows */
    lists[i] = l = (TIDLIST*)p; /* get/create the next trans. id list */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    next[i]  = p = l->tids;     /* note position of next trans. id */
    p += c[i]; *p++ = (TID)-1;  /* skip space for transaction ids */
  }                             /* and store a sentinel at the end */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(tabag, --n);  /* get the next transaction */
    rd.muls[n] = w = ta_wgt(t); /* and store its weight */
    for (s = ta_items(t); *s > TA_END; s++) {
      rd.tab[*s][n]    = w;     /* traverse the transaction's items */
      lists[*s]->supp += w;     /* and set the item occurrence flags */
      *next[*s]++      = n;     /* sum the transaction weight and */
    }                           /* collect the transaction ids */
  }
  max = 0;                      /* init. the maximal item support */
  for (i = m = 0; i < k; i++) { /* traverse the items / tid lists */
    l = lists[i];               /* eliminate all infrequent items and */
    if (l->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (l->supp >= pex) { isr_addpex(report, i); continue; }
    if (l->supp >  max)         /* find the maximal item support */
      max = l->supp;            /* (for later closed/maximal check) */
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_tab(lists, m, DIFFSIZE(p,tids), &rd);
  }                             /* find freq. item sets recursively */
  if (r >= 0) {                 /* if no error occurred */
    i = mode & (ISR_CLOSED|ISR_MAXIMAL);
    w = (i & ISR_MAXIMAL) ? rd.supp : tbg_wgt(tabag);
    if (!i || (max < w)) {      /* if to report the empty set */
      if (!isr_tidfile(report)) /* if not to report transaction ids, */
        r = (isr_report(report) < 0) ? -1 : 0;   /* report empty set */
      else {                    /* if to report transaction ids */
        for (n = tbg_cnt(tabag); n > 0; n--) tids[n] = n;
        r = (isr_reportx(report, tids, (diff_t)n) < 0) ? -1 : 0;
      }                         /* report the empty item set */
    }                           /* with all transaction ids */
  }
  free(tids); free(lists);      /* delete the allocated arrays */
  return r;                     /* return the error status */
}  /* eclat_tab() */

/*----------------------------------------------------------------------
  Eclat with an Occurrence Indicator Table (Simplified)
----------------------------------------------------------------------*/

static int rec_simp (TID *tids, SUPP n, ITEM k, RECDATA *rd)
{                               /* --- eclat recursion (table based) */
  int  r;                       /* error status */
  ITEM z;                       /* loop variable */
  SUPP s, w;                    /* item set support, weight buffer */
  SUPP pex;                     /* trans. count for perfect exts. */
  TID  *dst, *d, *p;            /* to traverse transaction ids */
  SUPP *row;                    /* to traverse occurrence table rows */

  assert(tids                   /* check the function arguments */
  &&    (n > 0) && (k > 0) && rd);
  pex = (rd->mode & ECL_PERFECT) ? n : SUPP_MAX;
  dst = tids +n+1;              /* get destination for intersections */
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir){ /* traverse the remaining items */
    row = rd->tab[k]; s = 0;    /* filter tids with item's table row */
    for (d = dst, p = tids; *p >= 0; p++)
      if ((w = row[*p]) > 0) {  /* compute the item set support */
        s += w; *d++ = *p; }    /* and the reduced trans. id list */
    if (s < rd->supp) continue; /* skip infrequent items and */
    if ((w = (SUPP)(d-dst)) >= pex) { /* collect perfect extensions */
      isr_addpex(rd->report, k); continue; }
    *d = -1;                    /* store a sentinel at the list end */
    r  = isr_add(rd->report, k, s);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if ((k > 0)                 /* if another item can be added */
    &&  isr_xable(rd->report,1) /* and upper size limit not reached */
    &&  ((r = rec_simp(dst, w, k, rd)) < 0))
      break;                    /* recursively find freq. item sets */
    if (isr_reportx(rd->report, tids, (diff_t)-s) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  return r;                     /* return the error status */
}  /* rec_simp() */

/*----------------------------------------------------------------------
Note that no memory is allocated in the above function; all processing
is done in the single memory block that is allocated in the function
below. The size of this memory block is O(n*k), where n is the number
of items and k the number of transactions. Additional memory is only
allocated in the item set reporter if closed or maximal item sets are
to be found, since this requires setting up an item set repository.
----------------------------------------------------------------------*/

int eclat_simp (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with occurrence table */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k;              /* loop variable, number of items */
  TID        n, m;              /* number of transactions */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  TID        *tids;             /* transaction identifier array */
  SUPP       *p;                /* to traverse occurrence table rows */
  const ITEM *s;                /* to traverse transaction items */
  TRACT      *t;                /* to traverse the transactions */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = -1;                 /* (supports only downward currently) */
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  if (tbg_wgt(tabag) < rd.supp) /* check the total transaction weight */
    return 0;                   /* against the minimum support */
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  x = tbg_extent(tabag);        /* get the number of item instances */
  if ((SIZE_MAX/sizeof(TID) -x-(size_t)n-1) / (size_t)(n+2) < (size_t)k)
    return -1;                  /* check the database/table size */
  x += (size_t)n+1+(size_t)k;   /* compute the database/table size */
  rd.tab = (SUPP**)malloc((size_t)k           *sizeof(SUPP*)
                        + (size_t)k*(size_t)n *sizeof(SUPP)
                        +         x           *sizeof(TID));
  if (!rd.tab) return -1;       /* allocate working memory */
  p = (SUPP*)memset(rd.tab +k, 0, (size_t)k*(size_t)n *sizeof(SUPP));
  for (i = 0; i < k; i++) {     /* init and organize the table rows */
    rd.tab[i] = p; p += n; }    /* (one table row per item) */
  tids = (TID*)p;               /* get the transaction id array */
  for (m = 0; m < n; m++) {     /* traverse the transactions */
    tids[m] = m;                /* set the initial (full) tid set */
    t = tbg_tract(tabag, m);    /* get the next transaction */
    w = ta_wgt(t);              /* and note its weight */
    for (s = ta_items(t); *s > TA_END; s++)
      rd.tab[*s][m] = w;        /* set the item occurrence flags */
  }                             /* (item *s occurs in transaction i) */
  tids[n] = (TID)-1;            /* store a sentinel at the end */
  if (isr_xable(report, 1)) {   /* if item set may be extended */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_simp(tids, (SUPP)n, k, &rd);
  }                             /* recursively find freq. item sets */
  if (r >= 0)                   /* report the empty item set */
    r = (isr_report(report) < 0) ? -1 : 0;
  free(rd.tab);                 /* delete the allocated table/arrays */
  return r;                     /* return the error status */
}  /* eclat_simp() */

/*----------------------------------------------------------------------
  Eclat with Transaction Ranges
----------------------------------------------------------------------*/

static void build_trg (TRGLIST **lists, TIDRANGE **next,
                       TABAG *tabag, TID min, TID max, ITEM off)
{                               /* --- build the trans. range lists */
  ITEM     i;                   /* loop variable */
  TID      k;                   /* loop variable */
  SUPP     w;                   /* weight buffer */
  ITEM     item;                /* to traverse items at offset */
  TRGLIST  *l;                  /* to access the trans. range lists */
  TIDRANGE *r;                  /* to access the transaction ranges */
  TRACT    *t;                  /* to traverse the transactions */

  assert(lists && tabag         /* check the function arguments */
  &&    (min >= 0) && (max < (TID)tbg_cnt(tabag)) && (off >= 0));

  /* --- skip short transactions --- */
  while ((min <= max)           /* traverse the transactions */
  &&     (ta_items(tbg_tract(tabag, min))[off] <= TA_END))
    ++min;                      /* skip trans. that are too short */
  if (min > max) return;        /* check for an empty trans. range */

  /* --- handle packed items --- */
  if (off <= 0) {               /* if at first item in transactions */
    l = lists[0];               /* get the list for packed items */
    for (k = min; min <= max; min++) {
      t = tbg_tract(tabag,min); /* traverse the transactions */
      i = ta_items(t)[off];     /* get the first item from them */
      if (i >= 0) break;        /* if it is not packed, abort loop */
      r = next[0]++;            /* get the current range in list */
      r->min = min;             /* store the transaction id and */
      r->max = (TID)(BITTA)i;   /* the bit repr. of the items */
      l->supp += r->wgt = ta_wgt(t);
    }                           /* store and sum transaction weight */
    if (min > k) {              /* if the trans. range is not empty */
      build_trg(lists, next, tabag, k, min-1, off+1);
      if (min > max) return;    /* recursively build trans. ranges, */
    }                           /* check whether an empty range */
  }                             /* is left to be processed */

  /* --- handle normal items --- */
  t = tbg_tract(tabag, min);    /* get the first transaction */
  i = item = ta_items(t)[off];  /* and from it the first item */
  do {                          /* traverse the longer transactions */
    w = ta_wgt(t);              /* init. the transaction weight */
    for (k = min; ++min <= max; ) {  /* while not at end of section */
      t = tbg_tract(tabag,min); /* get the next transaction and */
      i = ta_items(t)[off];     /* from it the item at the offset */
      if (i != item) break;     /* if the item differs, abort loop */
      w += ta_wgt(t);           /* otherwise sum the trans. weight */
    }                           /* (collect trans. with same item) */
    l = lists[item];            /* get list for the current item */
    r = next[item]++; item = i; /* and create a new trans. id range */
    l->supp += r->wgt = w;      /* store the transaction weights */
    build_trg(lists, next, tabag, r->min = k, r->max = min-1, off+1);
  } while (min <= max);         /* create the children recursively */
}  /* build_trg() */            /* while the range is not empty */

/*--------------------------------------------------------------------*/

static TID isect_trg (TRGLIST *dst, TRGLIST *src1, TRGLIST *src2)
{                               /* --- intersect two range lists */
  TIDRANGE *s1, *s2, *d, *p;    /* to traverse sources and dest. */

  assert(dst && src1 && src2);  /* check the function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  s1 = src1->trgs; s2 = src2->trgs; d = dst->trgs-1; p = NULL;
  while (1) {                   /* range list intersection loop */
    if      (s1->max < s2->min) {    /* skip transaction ranges */
      if ((++s1)->min < 0) break; }  /* that do not overlap */
    else if (s2->max < s1->min) {    /* check for the sentinel */
      if ((++s2)->min < 0) break; }  /* after advancing the range */
    else {                      /* if the transaction ranges overlap */
      if (s1 == p)              /* if there was a previous overlap, */
        d->wgt += s2->wgt;      /* only add the additional support */
      else {                    /* if this is a new overlap, */
        p = s1; ++d;            /* get/create a new trans. range */
        d->min = s1->min;       /* note the minimum and the */
        d->max = s1->max;       /* maximum transaction identifier */
        d->wgt = s2->wgt;       /* and the corresponding support */
      }
      dst->supp += s2->wgt;     /* sum the support for the item */
      if ((++s2)->min < 0) break;
    }                           /* skip the processed trans. range */
  }                             /* and check for the sentinel */
  (++d)->min = -1;              /* store a sentinel at the list end */
  return (TID)(++d -dst->trgs); /* return the size of the new list */
}  /* isect_trg() */

/*--------------------------------------------------------------------*/

static TID filter_trg (TRGLIST *dst, TRGLIST *src1, TRGLIST *src2)
{                               /* --- filter tids with a range list */
  TIDRANGE *s1, *s2, *d;        /* to traverse sources and dest. */

  assert(dst && src1 && src2);  /* check the function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = 0;                /* initialize the support */
  s1 = src1->trgs; s2 = src2->trgs; d = dst->trgs-1;
  while (1) {                   /* transaction id filtering loop */
    if      (s1->min < s2->min) {    /* skip transaction ids */
      if ((++s1)->min < 0) break; }  /* that are not in range */
    else if (s1->min > s2->max) {    /* check for the sentinel */
      if ((++s2)->min < 0) break; }  /* after advancing the range */
    else {                      /* if transaction id is in range */
      *++d = *s1;               /* copy the entry to the dest. */
      dst->supp += s1->wgt;     /* sum the support for the item */
      if ((++s1)->min < 0) break;
    }                           /* check for the sentinel */
  }
  (++d)->min = -1;              /* store a sentinel at the list end */
  return (TID)(++d -dst->trgs); /* return the size of the new list */
}  /* filter_trg() */

/*--------------------------------------------------------------------*/

static int rec_trg (TRGLIST **lists, ITEM k, size_t x, RECDATA *rd)
{                               /* --- eclat recursion with ranges */
  int      r;                   /* error status */
  ITEM     i, m, z;             /* loop variables */
  SUPP     pex;                 /* minimum support for perfect exts. */
  TRGLIST  **proj = NULL;       /* range lists of projected database */
  TRGLIST  *l, *d;              /* to traverse trans. range lists */
  TIDRANGE *p;                  /* to traverse transaction ranges */
  ITEM     *t;                  /* to collect the tail items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = lists[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TRGLIST**)malloc((size_t)k *sizeof(TRGLIST*) +x);
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for projected database) */
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / range lists */
    if (l->item < 0) {          /* if this list is for packed items */
      for (p = l->trgs; p->min >= 0; p++)
        m16_add(rd->fim16, (BITTA)p->max, p->wgt);
      r = m16_mine(rd->fim16);  /* add bit-rep. transaction prefixes */
      if (r >= 0) continue;     /* to the 16-items machine and mine, */
      if (proj) free(proj);     /* then go to the next tid range list */
      return r;                 /* otherwise free allocated memory */
    }                           /* and abort the function */
    r = isr_add(rd->report, l->item, l->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TRGLIST*)(proj +k+1);
      if (lists[i = 0]->item < 0) { /* if there are packed items */
        x = (size_t)filter_trg(d, lists[i++], l);
        if (d->supp >= rd->supp)    /* if they are frequent */
          proj[++m] = d = (TRGLIST*)(d->trgs +x);
      }                         /* add a range list for packed items */
      for ( ; i < k; i++) {     /* traverse the preceding lists */
        x = (size_t)isect_trg(d, lists[i], l);
        if (d->supp < rd->supp) /* intersect tid range lists */
          continue;             /* eliminate infrequent items */
        if (d->supp >= pex) {   /* collect perfect extensions */
          isr_addpex(rd->report, d->item); continue; }
        proj[++m] = d = (TRGLIST*)(d->trgs +x);
      }                         /* collect the trans. range lists */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_trg(proj, m, DIFFSIZE(d,proj[0]), rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_report(rd->report) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_trg() */

/*--------------------------------------------------------------------*/

int eclat_trg (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with transaction ranges */
  int       r = 0;              /* result of recursion/error status */
  ITEM      i, k, m;            /* loop variable, number of items */
  TID       n;                  /* number of transactions */
  size_t    x;                  /* number of item instances */
  SUPP      pex;                /* minimum support for perfect exts. */
  TRGLIST   **lists, *l;        /* to traverse trans. range lists */
  TIDRANGE  *trgs, *p, **next;  /* to traverse transaction ranges */
  const TID *c;                 /* item occurrence counters */
  RECDATA   rd;                 /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (tbg_packcnt(tabag) <= 0) mode &= ~ECL_FIM16;
  if (!(mode & ISR_MAXIMAL) || (mode & ECL_FIM16)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  c = tbg_icnts(tabag, 0);      /* get the number of containing */
  if (!c) return -1;            /* transactions per item */
  lists = (TRGLIST**)malloc((size_t)k *sizeof(TRGLIST*)
                           +(size_t)k *sizeof(TIDRANGE*));
  if (!lists) return -1;        /* create initial range lists array */
  next = (TIDRANGE**)(lists+k); /* and split off next range array */
  for (x = 0, i = 0; i < k; i++)/* traverse the items and sum */
    x += (size_t)c[i];          /* the number of item occurrences */
  /* Do not use tbg_extent(), because it does not take packed items */
  /* properly into account and thus may yield too big a value.      */
  p = trgs = (TIDRANGE*)malloc((size_t)k *sizeof(TRGLIST)
                                      +x *sizeof(TIDRANGE));
  if (!p) { free(lists); return -1; }
  for (i = 0; i < k; i++) {     /* allocate range list elements */
    lists[i] = l = (TRGLIST*)p; /* and organize the range lists */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    next[i]  = p = l->trgs;     /* note position of next trans. id */
    p += c[i];                  /* skip space for transaction ids */
    (p++)->min = (TID)-1;       /* and store a sentinel at the end */
  }
  build_trg(lists, next, tabag, 0, (TID)(n-1), 0);
  rd.fim16 = NULL;              /* build the transaction ranges */
  l = lists[i = 0];             /* get the list for packed items */
  if ((l->supp >= rd.supp)      /* if there are packed items */
  &&  (mode & ECL_FIM16)) {     /* and to use a 16-items machine */
    rd.fim16 = m16_create(rd.dir, rd.supp, report);
    if (!rd.fim16) { free(trgs); free(lists); return -1; }
    next[i++]->min = (TID)-1;   /* store a sentinel at the list end */
    l->item        = -1;        /* mark list for the packed items */
  }                             /* and store it in the reduced array */
  for (m = i; i < k; i++) {     /* traverse the trans. range lists */
    l = lists[i];               /* eliminate all infrequent items and */
    if (l->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (l->supp >= pex) { isr_addpex(report, i); continue; }
    next[i]->min = (TID)-1;     /* store a sentinel at the list end */
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_trg(lists, m, DIFFSIZE(p,trgs), &rd);
  }                             /* find freq. items sets recursively */
  if (r >= 0)                   /* if no error occurred, report */
    r = (isr_report(report) < 0) ? -1 : 0;/* the empty item set */
  if (rd.fim16)                 /* if a 16-items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16-items machine */
  free(trgs); free(lists);      /* delete the allocated arrays */
  return r;                     /* return the error status */
}  /* eclat_trg() */

/*----------------------------------------------------------------------
  Eclat with Occurrence Deliver (LCM-style)
----------------------------------------------------------------------*/

static int rec_odfx (TALIST **lists, ITEM k, RECDATA *rd)
{                               /* --- occ. deliver w/o reordering */
  int        r;                 /* error status */
  ITEM       i, m;              /* loop variables */
  TID        n;                 /* loop variable for transactions */
  SUPP       w;                 /* weight/support buffer */
  SUPP       pex;               /* minimum support for perfect exts. */
  TALIST     *l, *p;            /* to traverse transaction lists */
  TRACT      *t;                /* to traverse transactions */
  const ITEM *s;                /* to traverse items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  l = lists[k];                 /* collate equal transactions */
  taa_collate(l->tracts, l->cnt, k);
  for (n = 0; n < l->cnt; n++){ /* traverse the transactions, */
    t = l->tracts[n];           /* but skip collated transactions */
    if ((w = ta_wgt(t)) <= 0) continue;
    s = ta_items(t);            /* if there are packed items, */
    if (ispacked(*s))           /* add them to the 16-items machine */
      m16_add(rd->fim16, (BITTA)*s++, w);
    for ( ; (UITEM)*s < (UITEM)k; s++) {
      p = lists[*s]; p->supp += w; p->tracts[p->cnt++] = t; }
  }                             /* deliver the item occurrences */
  pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
  for (m = 0, i = rd->first; i < k; i++) {
    p = lists[i];               /* traverse the items / trans. lists */
    if (p->supp <  rd->supp) {  /* eliminate infrequent items */
      p->supp = 0; p->cnt = 0; continue; }
    if (p->supp >= pex) {       /* collect perfect extension items */
      p->supp = 0; p->cnt = 0; isr_addpex(rd->report, i); continue; }
    m++;                        /* count the frequent items */
  }                             /* (to see whether loop is needed) */
  r = (rd->fim16)               /* if there is a 16-items machine, */
    ? m16_mine(rd->fim16) : 0;  /* execute the 16-items machine */
  if (m <= 0) {                 /* if no frequent items found, abort */
    taa_uncoll(l->tracts, l->cnt); return r; }
  m = isr_xable(rd->report, 2) ? 0 : ITEM_MAX;
  for (i = rd->first; i < k; i++) {
    p = lists[i];               /* traverse the items / trans. lists, */
    if (p->supp <= 0) continue; /* but skip all eliminated items */
    r = isr_add(rd->report, i, p->supp);
    if (r < 0) break;           /* add current item to the reporter */
    if (r > 0) {                /* if the item needs processing */
      if (i > m) {              /* if to compute a projection, */
        r = rec_odfx(lists, i, rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* and check for a recursion error */
      if (isr_report(rd->report) < 0) {
        r = -1; break; }        /* report the current item set */
      isr_remove(rd->report,1); /* remove the current item */
    }                           /* from the item set reporter */
    p->supp = 0; p->cnt = 0;    /* reinitialize the transaction list */
  }
  taa_uncoll(l->tracts,l->cnt); /* uncollate the transactions */
  return r;                     /* and return the error status */
}  /* rec_odfx() */

/*--------------------------------------------------------------------*/

static int rec_odro (TALIST **lists, ITEM k, RECDATA *rd)
{                               /* --- occ. deliver with reordering */
  int        r;                 /* error status */
  ITEM       i, m, b;           /* loop variables */
  TID        n;                 /* number of transactions */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  SUPP       pex;               /* minimum support for perfect exts. */
  TALIST     *l, *p;            /* to traverse transaction lists */
  TRACT      *t;                /* to traverse transactions */
  const ITEM *s;                /* to traverse items */
  SUPP       *supp;             /* item support array */
  ITEM       *map, *inv;        /* item identifier maps */
  TALIST     **dst;             /* destination for reduction */
  void       *mem = NULL;       /* memory for reduction */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  supp = (SUPP*)memset(rd->muls, 0, (size_t)k *sizeof(SUPP));
  l    = lists[k];              /* initialize the item support array */
  for (x = 0, n = 0; n < l->cnt; n++) {
    t = l->tracts[n];           /* traverse the transactions */
    w = ta_wgt(t);              /* get the transaction weight */
    for (s = ta_items(t); (UITEM)*s < (UITEM)k; s++)
      supp[*s] += w;            /* determine the support of the items */
    x += (size_t)(s -ta_items(t));
  }                             /* compute the size of the database */
  pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
  for (inv = rd->cand, i = m = 0; i < k; i++) {
    if (supp[i] <  rd->supp) {  /* traverse items and their support */
      supp[i] = -1; continue; } /* eliminate infrequent items an */
    if (supp[i] >= pex) {       /* collect perfect extension items */
      supp[i] = -1; isr_addpex(rd->report, lists[i]->item); continue; }
    inv[m++] = i;               /* collect the remaining items */
  }                             /* (map from new to old identifiers) */
  if (m <= 0) return 0;         /* if no frequent items found, abort */
  i = inv[m-1];                 /* get the highest frequent item */
  if (++i < k) k = i;           /* and compute limit for item loop */
  i2s_sort(inv, (size_t)m, -1, supp);
  map = inv+m;                  /* sort items by their support */
  if ((m <= 16) && rd->fim16) { /* if at most 16-items left */
    for (i = 0; i < k; i++) map[i] = -1;
    for (i = 0; i < m; i++) {   /* build a map from the old */
      map[inv[i]] = i;          /* to the new item identifiers */
      m16_setmap(rd->fim16, i, lists[inv[i]]->item);
    }                           /* set the item identifier map */
    for (n = 0; n < l->cnt; n++) {
      t = l->tracts[n]; b = 0;  /* traverse the transactions */
      for (s = ta_items(t); (UITEM)*s < (UITEM)k; s++)
        if ((i = map[*s]) >= 0) b |= 1 << i;
      m16_add(rd->fim16, (BITTA)b, ta_wgt(t));
    }                           /* add bit-represented transactions */
    return m16_mine(rd->fim16); /* to the 16-items machine and */
  }                             /* mine frequent item sets */
  for (i = 0; i < k; i++) {     /* copy support to trans. lists */
    if (supp[i] > 0) lists[i]->supp = supp[i];
    else           { lists[i]->supp = 0; map[i] = -1; }
  }                             /* set map for eliminated items */
  dst = lists;                  /* get the trans. lists to process */
  if ((l->cnt >= 6)             /* if there are enough transactions */
  &&  (m >= ((rd->fim16) ? 20 : 6))) {     /* and enough items left */
    dst = (TALIST**)malloc((size_t)m *sizeof(TALIST*));
    if (!dst) return -1;        /* allocate memory for projection */
    for (i = 0; i < m; i++) {   /* copy the transaction lists that */
      dst[i] = lists[inv[i]];   /* are needed for the projection */
      map[inv[i]] = i;          /* and build a map from the old */
    }                           /* to the new item identifiers */
    mem = malloc(taa_dstsize(l->cnt, x)); /* allocate memory for */
    if (!mem) { free(dst); return -1; }   /* destination trans. */
    l->cnt = taa_reduce(l->tracts, l->cnt, k, map, rd->hash, &mem);
    k = m;                      /* reduce the transactions */
  }                             /* (remove items, collate trans.) */
  if (rd->fim16                 /* if to use a 16-items machine */
  && (rd->dir > 0)) {           /* and forward processing direction */
    for (n = 0; n < l->cnt; n++) {
      t = l->tracts[n]; b = 0;  /* traverse the transactions */
      for (s = ta_items(t); (UITEM)*s < (UITEM)16; s++)
        b |= 1 << *s;           /* add trans. to 16-items machine */
      m16_add(rd->fim16, (BITTA)b, ta_wgt(t));
      for ( ; (UITEM)*s < (UITEM)k; s++) {
        p = dst[*s]; p->tracts[p->cnt++] = t; }
    }                           /* deliver the item occurrences */
    for (i = 0; i < 16; i++) {  /* traverse the first 16 items */
      l = dst[i]; l->supp = 0; l->cnt = 0; /* and clear support */
      m16_setmap(rd->fim16, i, l->item);
    }                           /* set the item identifier map */
    r = m16_mine(rd->fim16);    /* mine with 16-items machine */
    if (r < 0) return r; }      /* and check for an error */
  else {                        /* if not to use a 16-items machine */
    for (n = 0; n < l->cnt; n++) {
      t = l->tracts[n];         /* traverse the transactions */
      for (s = ta_items(t); (UITEM)*s < (UITEM)k; s++) {
        p = dst[*s]; p->tracts[p->cnt++] = t; }
    }                           /* deliver the item occurrences */
    i = 0;                      /* to the transaction lists and */
  }                             /* get first item index to process */
  m = isr_xable(rd->report, 2) ? 0 : ITEM_MAX;
  for (r = 0; i < k; i++) {     /* traverse the items/trans. lists, */
    l = dst[i];                 /* but skip all eliminated items */
    if (l->supp <= 0) { l->cnt = 0; continue; }
    r = isr_add(rd->report, l->item, l->supp);
    if (r < 0) break;           /* add current item to the reporter */
    if (r > 0) {                /* if the item needs processing */
      if (i > m) {              /* if to compute a projection */
        r = rec_odro(dst, i, rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* and check for a recursion error */
      if (isr_report(rd->report) < 0) {
        r = -1; break; }        /* report the current item set */
      isr_remove(rd->report,1); /* remove the current item */
    }                           /* from the item set reporter */
    l->supp = 0; l->cnt = 0;    /* reinitialize the transaction list */
  }
  if (mem) {                    /* delete projection lists */
    free(mem); free(dst); }     /* and destination memory */
  return r;                     /* return the error status */
}  /* rec_odro() */

/*--------------------------------------------------------------------*/

int eclat_ocd (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with occurrence deliver */
  int       r = 0;              /* result of recursion/error status */
  ITEM      i, k;               /* loop variable, number of items */
  TID       n, m;               /* number of transactions */
  size_t    x, h;               /* extent, hash table size */
  SUPP      pex;                /* minimum support for perfect exts. */
  TALIST    **lists, *l;        /* to traverse transaction lists */
  TRACT     **tras, **p;        /* to traverse transactions */
  const TID *c;                 /* item occurrence counters */
  RECDATA   rd;                 /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  c = tbg_icnts(tabag, 0);      /* get the number of containing */
  if (!c) return -1;            /* transactions per item */
  lists = (TALIST**)malloc((size_t)(k+1) *sizeof(TALIST*));
  if (!lists) return -1;        /* create the trans. id list array */
  for (x = 0, i = 0; i < k; i++)/* traverse the items and sum */
    x += (size_t)c[i];          /* the numbers of item occurrences */
  /* Do not use tbg_extent(), because it does not take packed items */
  /* properly into account and thus may yield too big a value.      */
  h = (size_t)taa_tabsize(n);   /* get the hash table size and */
  x = (size_t)n +x -(size_t)(k+1);        /* the database size */
  p = tras = (TRACT**)malloc((size_t)(k+1) *sizeof(TALIST)
                           +         (x+h) *sizeof(TRACT*)
                           + (size_t) k    *sizeof(SUPP)
                           + (size_t)(k+k) *sizeof(ITEM));
  if (!p) { free(lists); return -1; }
  for (i = 0; i < k; i++) {     /* allocate the list elements and */
    lists[i] = l = (TALIST*)p;  /* traverse the items / trans. lists */
    l->item  = i;               /* set the item identifier */
    l->supp  = 0;               /* clear the item support */
    l->cnt   = 0;               /* and the transaction counter */
    p = l->tracts +c[i];        /* skip space for transactions */
  }                             /* and a sentinel at the end */
  lists[k] = l = (TALIST*)p;    /* set last list (all transactions) */
  l->item  = k; l->cnt = n;     /* for a dummy item (> all items) */
  l->supp  = tbg_wgt(tabag);    /* with the full database support */
  for (m = 0; m < n; m++)       /* copy the transactions */
    l->tracts[m] = tbg_tract(tabag, m);
  rd.hash  = (TRACT**)memset(l->tracts+n, 0, h *sizeof(TRACT*));
  rd.muls  = (SUPP*)  memset(rd.hash  +h, 0, (size_t)k *sizeof(SUPP));
  rd.cand  = (ITEM*)(rd.muls+k);/* get the auxiliary arrays */
  rd.fim16 = NULL; rd.first = 0;/* default: no 16-items machine */
  if (mode & ECL_FIM16) {       /* if to use a 16-items machine */
    rd.fim16 = m16_create(rd.dir, rd.supp, report);
    if (!rd.fim16) { free(tras); free(lists); return -1; }
    rd.first = tbg_packcnt(tabag);
  }                             /* get the number of packed items */
  rd.report = report;           /* initialize the recursion data */
  rd.tabag  = tabag;            /* (store reporter and transactions) */
  r = (mode & ECL_REORDER)      /* execute the eclat recursion */
    ? rec_odro(lists, k, &rd)   /* with    item reordering */
    : rec_odfx(lists, k, &rd);  /* without item reordering */
  if (r >= 0)                   /* report the empty item set */
    r = (isr_report(report) < 0) ? -1 : 0;
  if (rd.fim16)                 /* if a 16-items machine was used, */
    m16_delete(rd.fim16);       /* delete the 16-items machine */
  free(tras); free(lists);      /* deallocate the transaction array */
  return r;                     /* return the error status */
}  /* eclat_ocd() */

/*----------------------------------------------------------------------
  Eclat with Diffsets
----------------------------------------------------------------------*/

static TID cmpl (TIDLIST *dst, TIDLIST *src1, TIDLIST *src2, SUPP *muls)
{                               /* --- complement two trans. id lists */
  TID *s1, *s2, *d;             /* to traverse sources and dest. */

  assert(dst && src1 && src2 && muls); /* check function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = src1->supp;       /* initialize the support */
  d = dst->tids; s1 = src1->tids; s2 = src2->tids;
  while (1) {                   /* trans. id list difference loop */
    if      (*s1 > *s2) dst->supp -= muls[*s1++];
    else if (*s1 < *s2) *d++ = *s2++;
    else if (*s1 < 0) break;    /* collect elements of second source */
    else { s1++; s2++; }        /* that are not in the first source */
  }                             /* (form complement of first source) */
  *d++ = -1;                    /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new lists */
}  /* cmpl() */

/*--------------------------------------------------------------------*/

static TID diff (TIDLIST *dst, TIDLIST *src1, TIDLIST *src2, SUPP *muls)
{                               /* --- subtract two trans. id lists */
  TID *s1, *s2, *d;             /* to traverse sources and dest. */

  assert(dst && src1 && src2 && muls); /* check function arguments */
  dst->item = src1->item;       /* copy the first item and */
  dst->supp = src1->supp;       /* initialize the support */
  d = dst->tids; s1 = src1->tids; s2 = src2->tids;
  while (1) {                   /* trans. id list difference loop */
    if      (*s1 > *s2) *d++ = *s1++;
    else if (*s1 < *s2) dst->supp -= muls[*s2++];
    else if (*s1 < 0) break;    /* remove all elements of the second */
    else { s1++; s2++; }        /* source from the first source */
  }                             /* (form difference of tid lists) */
  *d++ = -1;                    /* store a sentinel at the list end */
  return (TID)(d -dst->tids);   /* return the size of the new lists */
}  /* diff() */

/*--------------------------------------------------------------------*/

static int rec_diff (TIDLIST **lists, ITEM k, TID x,
                     COMBFN comb, RECDATA *rd)
{                               /* --- eclat recursion with diffsets */
  int     r;                    /* error status */
  ITEM    i, m, z;              /* loop variables */
  TID     c;                    /* size of combined lists */
  SUPP    pex;                  /* minimum support for perfect exts. */
  TIDLIST *l, *d;               /* to traverse transaction id lists */
  TIDLIST **proj = NULL;        /* trans. id lists of proj. database */
  ITEM    *t;                   /* to collect the tail items */

  assert(lists && (k > 0) && rd);  /* check the function arguments */
  if (rd->mode & ECL_TAIL) {    /* if to use tail to prune w/ repo. */
    t = isr_buf(rd->report);    /* collect the tail items in buffer */
    for (m = 0, i = k; --i >= 0; ) t[m++] = lists[i]->item;
    r = isr_tail(rd->report, t, m);
    if (r) return r;            /* if tail need not be processed, */
  }                             /* abort the recursion */
  if ((k > 1)                   /* if there is more than one item */
  &&  isr_xable(rd->report,2)){ /* and another item can be added */
    proj = (TIDLIST**)malloc((size_t)k           *sizeof(TIDLIST*)
                            +(size_t)k           *sizeof(TIDLIST)
                            +(size_t)k*(size_t)x *sizeof(TID));
    if (!proj) return -1;       /* allocate list and element arrays */
  }                             /* (memory for projected database) */
  if ((k > 4)                   /* if there are enough items left, */
  &&  (rd->mode & ECL_REORDER)) /* re-sort the items w.r.t. support */
    ptr_qsort(lists, (size_t)k, +1, tid_cmp, NULL);
  if (rd->dir > 0) { z =  k; k  = 0; }
  else             { z = -1; k -= 1; }
  for (r = 0; k != z; k += rd->dir) {
    l = lists[k];               /* traverse the items / tid lists */
    r = isr_add(rd->report, l->item, l->supp);
    if (r <  0) break;          /* add current item to the reporter */
    if (r <= 0) continue;       /* check if item needs processing */
    if (proj && (k > 0)) {      /* if another item can be added */
      pex = (rd->mode & ECL_PERFECT) ? l->supp : SUPP_MAX;
      proj[m = 0] = d = (TIDLIST*)(proj +k+1); x = 0;
      for (i = 0; i < k; i++) { /* traverse the preceding lists */
        c = comb(d, lists[i], l, rd->muls);
        if (d->supp < rd->supp) /* combine transaction id lists */
          continue;             /* eliminate infrequent items */
        if (d->supp >= pex) {   /* collect perfect extensions */
          isr_addpex(rd->report, d->item); continue; }
        proj[++m] = d = (TIDLIST*)(d->tids +c);
        if (c > x) x = c;       /* collect the trans. id lists and */
      }                         /* determine their maximum length */
      if (m > 0) {              /* if the projection is not empty */
        r = rec_diff(proj, m, x, diff, rd);
        if (r < 0) break;       /* recursively find freq. item sets */
      }                         /* in the created projection */
    }
    if (isr_report(rd->report) < 0) {
      r = -1; break; }          /* report the current item set */
    isr_remove(rd->report, 1);  /* remove the current item */
  }                             /* from the item set reporter */
  if (proj) free(proj);         /* delete the list and element arrays */
  return r;                     /* return the error status */
}  /* rec_diff() */

/*--------------------------------------------------------------------*/

int eclat_diff (TABAG *tabag, int mode, SUPP supp, ISREPORT *report)
{                               /* --- eclat with difference sets */
  int        r = 0;             /* result of recursion/error status */
  ITEM       i, k, m;           /* loop variable, number of items */
  TID        n, z;              /* (maximum) number of transactions */
  size_t     x;                 /* number of item instances */
  SUPP       w;                 /* weight/support buffer */
  SUPP       pex;               /* minimum support for perfect exts. */
  TRACT      *t;                /* to traverse transactions */
  TIDLIST    **lists, *l;       /* to traverse transaction id lists */
  TID        *tids, *p, **next; /* to traverse transaction ids */
  const ITEM *s;                /* to traverse transaction items */
  const TID  *c;                /* item occurrence counters */
  RECDATA    rd;                /* recursion data */

  assert(tabag && report);      /* check the function arguments */
  if (!(mode & ISR_MAXIMAL)) mode &= ~ECL_TAIL;
  rd.mode = mode;               /* store search mode and item dir. */
  rd.dir  = (mode & (ISR_CLOSED|ISR_MAXIMAL)) ? -1 : +1;
  rd.supp = (supp > 0) ? supp : 1;  /* check and adapt the support */
  pex     = tbg_wgt(tabag);     /* check the total transaction weight */
  if (rd.supp > pex) return 0;  /* and get support for perfect exts. */
  if (!(mode & ECL_PERFECT)) pex = SUPP_MAX;
  n = tbg_cnt(tabag);           /* get the number of transactions */
  k = tbg_itemcnt(tabag);       /* and check the number of items */
  if (k <= 0) return (isr_report(report) < 0) ? -1 : 0;
  c = tbg_icnts(tabag, 0);      /* get the number of containing */
  if (!c) return -1;            /* transactions per item */
  lists = (TIDLIST**)malloc((size_t)k *sizeof(TIDLIST*)
                           +(size_t)k *sizeof(TID*)
                           +(size_t)n *sizeof(SUPP));
  if (!lists) return -1;        /* create initial tid list array */
  next     = (TID**)(lists+k);  /* and split off next position array */
  rd.muls  = (SUPP*)(next +k);  /* and transaction multiplicity array */
  x = tbg_extent(tabag);        /* get the number of item occurrences */
  p = tids = (TID*)malloc((size_t)k *sizeof(TIDLIST) +x *sizeof(TID));
  if (!p) { free(lists); return -1; }  /* allocate tid list elements */
  for (i = 0; i < k; i++) {     /* traverse the items / tid lists */
    lists[i] = l = (TIDLIST*)p; /* get/create the next trans. id list */
    l->item  = i;               /* initialize the list item */
    l->supp  = 0;               /* and the support counter */
    next[i]  = p = l->tids;     /* note position of next trans. id */
    p += c[i]; *p++ = (TID)-1;  /* skip space for transaction ids */
  }                             /* and store a sentinel at the end */
  while (n > 0) {               /* traverse the transactions */
    t = tbg_tract(tabag, --n);  /* get the next transaction */
    rd.muls[n] = w = ta_wgt(t); /* and store its weight */
    for (s = ta_items(t); *s > TA_END; s++) {
      lists[*s]->supp += w;     /* traverse the transaction's items */
      *next[*s]++      = n;     /* sum the transaction weight and */
    }                           /* collect the transaction ids */
  }
  z = 0;                        /* init. the maximal list length */
  for (i = m = 0; i < k; i++) { /* traverse the items / tid lists */
    l = lists[i];               /* eliminate all infrequent items and */
    if (l->supp <  rd.supp) continue;   /* collect perfect extensions */
    if (l->supp >= pex) { isr_addpex(report, i); continue; }
    n = (TID)(next[i] -l->tids);
    if (n > z) z = n;           /* find maximum trans. id list length */
    lists[m++] = l;             /* collect lists for frequent items */
  }                             /* (eliminate infrequent items) */
  if (m > 0) {                  /* if there are frequent items */
    rd.report = report;         /* initialize the recursion data */
    rd.tabag  = tabag;          /* (store reporter and transactions) */
    r = rec_diff(lists, m, z, cmpl, &rd);
  }                             /* find freq. items sets recursively */
  if (r >= 0)                   /* report the empty item set */
    r = (isr_report(report) < 0) ? -1 : 0;
  free(tids); free(lists);      /* delete the allocated arrays */
  return r;                     /* return the error status */
}  /* eclat_diff() */

/*----------------------------------------------------------------------
  Eclat (generic)
----------------------------------------------------------------------*/

void ecl_adjust (int target, int eval,
                 int *algo, int *mode, int *pack, int *mrep)
{                               /* --- adjust algorithm and modes */
  assert(algo && mode);         /* check the function arguments */
  if ((*mode & ECL_TIDOUT)      /* if trans. identifiers requested */
  &&  (*algo != ECL_LISTS) && (*algo != ECL_TABLE))
    *algo  =  ECL_LISTS;        /* trans. identifier lists are needed */
  if (*algo != ECL_LISTS)       /* only tid lists allow for filtering */
    *mode &= ~ISR_NOFILTER;     /* for closed/maximal without a repo. */
  if (target & (ISR_CLOSED|ISR_MAXIMAL)) {
    *mode &= ~ECL_REORDER;      /* cannot reorder for closed/maximal */
    if (*algo == ECL_OCCDLV) *algo = ECL_LISTS; }
  else if (target & ISR_GENERA){/* if to filter for generators, */
    *mode |= ECL_PERFECT;       /* need perfect extension pruning */
    if (*algo == ECL_SIMPLE) *algo = ECL_TABLE;
  }                             /* cannot use simple table variant */
  if (((*algo != ECL_LISTS) && (*algo != ECL_RANGES)
  &&   (*algo != ECL_OCCDLV)) || (*mode & ISR_NOFILTER))
    *mode &= ~ECL_FIM16;        /* restrict use of 16-items machine */
  if ((*algo == ECL_RANGES) || (*algo == ECL_SIMPLE))
    *mode &= ~ECL_REORDER;      /* no all variants allow reordering */
  if (mrep)                     /* copy filter flag to reporting */
    *mrep = (*mrep & ~ISR_NOFILTER) | (*mode & ISR_NOFILTER);
  if ((target == ISR_GENERA) && (*mode & ECL_REORDER)) {
    if (mrep) *mrep |= ISR_SORT; } /* reordering requires set sorting */
  if (pack) {                   /* if a packing parameter is given */
    if (*pack <  0) *pack =  0; /* clamp the number of items */
    if (*pack > 16) *pack = 16; /* for the k-items machine */
    if (mode && (*pack == 0)) *mode &= ~ECL_FIM16;
  }                             /* k-items machine requires packing */
  if ((*algo == ECL_OCCDLV) && (*mode & ECL_REORDER)) {
    if (pack) *pack = 0; }      /* delayed packing if reordering */
  if (mrep && (eval == 'b'))    /* if to evaluate found item sets, */
    *mrep |= ISR_LOGS;          /* logarithms need to be computed */
  if (mrep) *mrep |= target;    /* store target in report mode */
}  /* ecl_adjust() */

/*--------------------------------------------------------------------*/

static ECLATFN* eclatvars[] = { /* --- table of eclat variants */
  eclat_base,                   /* trans. id lists (basic) */
  eclat_tid,                    /* trans. id lists (improved) */
  eclat_bit,                    /* bit vector over transactions */
  eclat_tab,                    /* item occurrence table */
  eclat_simp,                   /* simplified version with table */
  eclat_trg,                    /* transaction identifier ranges */
  eclat_ocd,                    /* occurrence deliver (LCM-style) */
  eclat_diff,                   /* difference sets (diffsets) */
};

/*--------------------------------------------------------------------*/

int eclat (TABAG *tabag, int target, int algo, int mode,
           SUPP supp, int eval, double minval, ISREPORT *report)
{                               /* --- eclat algorithm */
  int     r;                    /* result of eclat algorithm */
  clock_t t;                    /* timer for measurements */

  assert(tabag && report);      /* check the function arguments */
  t = clock();                  /* start the timer */
  if (eval == 'b')              /* if to compute add. evaluation */
    isr_seteval(report, isr_logrto, NULL, +1, minval);
  XMSG(stderr, "writing %s ... ", isr_name(report));
  r = eclatvars[algo](tabag, target|mode, (SUPP)supp, report);
  if (r < 0) return -1;         /* search for frequent item sets */
  XMSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(report));
  XMSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  return 0;                     /* return 'ok' */
}  /* eclat() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/
#ifdef ECL_MAIN

static void help (void)
{                               /* --- print add. option information */
  #ifndef QUIET
  fprintf(stderr, "\n");        /* terminate startup message */
  printf("eclat algorithm variants (option -a#)\n");
  printf("  e   transaction id lists intersection (basic)\n");
  printf("  i   transaction id lists intersection (improved)\n");
  printf("  b   transaction id lists represented as bit vectors\n");
  printf("  t   filtering with item occurrence table (standard)\n");
  printf("  s   filtering with item occurrence table (simplified)\n");
  printf("  r   transaction id range lists intersection\n");
  printf("  o   occurrence deliver from transaction lists (default)\n");
  printf("  d   transaction id difference sets (diffsets/dEclat)\n");
  printf("Currently the algorithm variant 'o' cannot be used to find ");
  printf("closed\nor maximal item sets. Therefore it is automatically ");
  printf("changed to 'i'\nif one of the options -tc or -tm is given. ");
  printf("These restrictions may\nbe removed in future versions of ");
  printf("this program.\n");
  printf("\n");
  printf("additional evaluation measures (option -e#)\n");
  printf("  x   no measure (default)\n");
  printf("  b   binary logarithm of support quotient\n");
  printf("\n");
  printf("information output format characters (option -v#)\n");
  printf("  %%%%  a percent sign\n");
  printf("  %%i  number of items (item set size)\n");
  printf("  %%a  absolute item set support\n");
  printf("  %%s  relative item set support as a fraction\n");
  printf("  %%S  relative item set support as a percentage\n");
  printf("  %%e  additional evaluation measure\n");
  printf("  %%E  additional evaluation measure as a percentage\n");
  printf("All format characters can be preceded by the number\n");
  printf("of significant digits to be printed (at most 32 digits),\n");
  printf("even though this value is ignored for integer numbers.\n");
  #endif                        /* print help information */
  exit(0);                      /* abort the program */
}  /* help() */

/*--------------------------------------------------------------------*/

#ifndef NDEBUG                  /* if debug version */
  #undef  CLEANUP               /* clean up memory and close files */
  #define  CLEANUP \
  if (twrite) twr_delete(twrite, 1); \
  if (report) isr_delete(report, 0); \
  if (tabag)  tbg_delete(tabag,  0); \
  if (tread)  trd_delete(tread,  1); \
  if (ibase)  ib_delete (ibase);
#endif

GENERROR(error, exit)           /* generic error reporting function */

/*--------------------------------------------------------------------*/

int main (int argc, char *argv[])
{                               /* --- main function */
  int     i, k = 0;             /* loop variables, counters */
  char    *s;                   /* to traverse the options */
  CCHAR   **optarg = NULL;      /* option argument */
  CCHAR   *fn_inp  = NULL;      /* name of input  file */
  CCHAR   *fn_out  = NULL;      /* name of output file */
  CCHAR   *fn_sel  = NULL;      /* name of item selection file */
  CCHAR   *fn_tid  = NULL;      /* name of transaction ids file */
  CCHAR   *fn_psp  = NULL;      /* name of pattern spectrum file */
  CCHAR   *recseps = NULL;      /* record  separators */
  CCHAR   *fldseps = NULL;      /* field   separators */
  CCHAR   *blanks  = NULL;      /* blank   characters */
  CCHAR   *comment = NULL;      /* comment characters */
  CCHAR   *hdr     = "";        /* record header  for output */
  CCHAR   *sep     = " ";       /* item separator for output */
  CCHAR   *dflt    = " (%S)";   /* default format for check */
  CCHAR   *format  = dflt;      /* format for information output */
  int     target   = 's';       /* target type (closed/maximal) */
  ITEM    min      = 1;         /* minimum size of an item set */
  ITEM    max      = ITEM_MAX;  /* maximum size of an item set */
  double  supp     = 10;        /* minimum support (in percent) */
  int     eval     = 'x';       /* additional evaluation measure */
  double  minval   = 10;        /* minimum evaluation measure value */
  int     sort     =  2;        /* flag for item sorting and recoding */
  int     algo     = 'o';       /* variant of eclat algorithm */
  int     cmfilt   = -1;        /* mode for closed/maximal filtering */
  int     pack     = 16;        /* number of bit-packed items */
  int     mode     = ECL_DEFAULT;  /* search mode (e.g. pruning) */
  int     mtar     =  0;        /* mode for transaction reading */
  int     mrep     =  0;        /* mode for item set reporting */
  int     stats    =  0;        /* flag for item set statistics */
  PATSPEC *psp;                 /* collected pattern spectrum */
  ITEM    m;                    /* number of items */
  TID     n;                    /* number of transactions */
  SUPP    w;                    /* total transaction weight */
  clock_t t;                    /* timer for measurements */

  #ifndef QUIET                 /* if not quiet version */
  prgname = argv[0];            /* get program name for error msgs. */

  /* --- print usage message --- */
  if (argc > 1) {               /* if arguments are given */
    fprintf(stderr, "%s - %s\n", argv[0], DESCRIPTION);
    fprintf(stderr, VERSION); } /* print a startup message */
  else {                        /* if no argument is given */
    printf("usage: %s [options] infile [outfile]\n", argv[0]);
    printf("%s\n", DESCRIPTION);
    printf("%s\n", VERSION);
    printf("-t#      target type                              "
                    "(default: %c)\n", target);
    printf("         (s: frequent, c: closed, m: maximal item sets, "
                     "g: generators)\n");
    printf("-m#      minimum number of items per item set     "
                    "(default: %"ITEM_FMT")\n", min);
    printf("-n#      maximum number of items per item set     "
                    "(default: no limit)\n");
    printf("-s#      minimum support of an item set           "
                    "(default: %g%%)\n", supp);
    printf("         (positive: percentage, "
                     "negative: absolute number)\n");
    printf("-e#      additional evaluation measure            "
                    "(default: none)\n");
    printf("-d#      minimum value of add. evaluation measure "
                    "(default: %g%%)\n", minval);
    printf("-q#      sort items w.r.t. their frequency        "
                    "(default: %d)\n", sort);
    printf("         (1: ascending, -1: descending, 0: do not sort,\n"
           "          2: ascending, -2: descending w.r.t. "
                    "transaction size sum)\n");
    printf("-a#      variant of the eclat algorithm to use    "
                    "(default: %c)\n", algo);
    printf("-x       do not prune with perfect extensions     "
                    "(default: prune)\n");
    printf("-l#      number of items for k-items machine      "
                    "(default: %d)\n", pack);
    printf("         (only for algorithm variants i,r,o,   "
                    "options -ai/-ar/-ao)\n");
    printf("-p       do not sort items w.r.t. cond. support   "
                    "(default: sort)\n");
    printf("         (only for algorithm variants i,b,t,d, "
                    "options -ai/-ab/-at/-ad)\n");
    printf("-y#      check extensions for closed/maximal sets "
                    "(default: repository)\n");
    printf("         (0: horizontal, > 0: vertical representation)\n");
    printf("         (only with improved tid lists variant, "
                    "option -ai)\n");
    printf("-z       do not use head union tail (hut) pruning "
                    "(default: use hut)\n");
    printf("         (only for maximal item sets, option -tm, "
                    "not with option -ab)\n");
    printf("-R#      read an item selection from a file\n");
    printf("-P#      write a pattern spectrum to a file\n");
    printf("-Z       print item set statistics "
                    "(number of item sets per size)\n");
    printf("-g       write output in scanable form "
                    "(quote certain characters)\n");
    printf("-h#      record header  for output                "
                    "(default: \"%s\")\n", hdr);
    printf("-k#      item separator for output                "
                    "(default: \"%s\")\n", sep);
    printf("-v#      output format for item set information   "
                    "(default: \"%s\")\n", format);
    printf("-w       transaction weight in last field         "
                    "(default: only items)\n");
    printf("-r#      record/transaction separators            "
                    "(default: \"\\n\")\n");
    printf("-f#      field /item        separators            "
                    "(default: \" \\t,\")\n");
    printf("-b#      blank   characters                       "
                    "(default: \" \\t\\r\")\n");
    printf("-C#      comment characters                       "
                    "(default: \"#\")\n");
    printf("-T#      file to write transaction identifiers to "
                    "(default: none)\n");
    printf("-!       print additional option information\n");
    printf("infile   file to read transactions from           "
                    "[required]\n");
    printf("outfile  file to write frequent item sets to      "
                    "[optional]\n");
    return 0;                   /* print a usage message */
  }                             /* and abort the program */
  #endif  /* #ifndef QUIET */
  /* free option characters: cijlou [A-Z]\[CPRTZ] */

  /* --- evaluate arguments --- */
  for (i = 1; i < argc; i++) {  /* traverse arguments */
    s = argv[i];                /* get option argument */
    if (optarg) { *optarg = s; optarg = NULL; continue; }
    if ((*s == '-') && *++s) {  /* -- if argument is an option */
      while (*s) {              /* traverse options */
        switch (*s++) {         /* evaluate switches */
          case '!': help();                          break;
          case 't': target = (*s) ? *s++ : 's';      break;
          case 'm': min    = (ITEM)strtol(s, &s, 0); break;
          case 'n': max    = (ITEM)strtol(s, &s, 0); break;
          case 's': supp   =       strtod(s, &s);    break;
          case 'e': eval   = (*s) ? *s++ : 0;        break;
          case 'd': minval =       strtod(s, &s);    break;
          case 'q': sort   = (int) strtol(s, &s, 0); break;
          case 'a': algo   = (*s) ? *s++ : 0;        break;
          case 'x': mode  &= ~ECL_PERFECT;           break;
          case 'l': pack   = (int) strtol(s, &s, 0); break;
          case 'p': mode  &= ~ECL_REORDER;           break;
          case 'y': cmfilt = (int) strtol(s, &s, 0); break;
          case 'z': mode  &= ~ECL_TAIL;              break;
          case 'R': optarg = &fn_sel;                break;
          case 'P': optarg = &fn_psp;                break;
          case 'Z': stats  = 1;                      break;
          case 'g': mrep   = ISR_SCAN;               break;
          case 'h': optarg = &hdr;                   break;
          case 'k': optarg = &sep;                   break;
          case 'v': optarg = &format;                break;
          case 'w': mtar  |= TA_WEIGHT;              break;
          case 'r': optarg = &recseps;               break;
          case 'f': optarg = &fldseps;               break;
          case 'b': optarg = &blanks;                break;
          case 'C': optarg = &comment;               break;
          case 'T': optarg = &fn_tid;                break;
          default : error(E_OPTION, *--s);           break;
        }                       /* set option variables */
        if (optarg && *s) { *optarg = s; optarg = NULL; break; }
      } }                       /* get option argument */
    else {                      /* -- if argument is no option */
      switch (k++) {            /* evaluate non-options */
        case  0: fn_inp = s;      break;
        case  1: fn_out = s;      break;
        default: error(E_ARGCNT); break;
      }                         /* note filenames */
    }
  }
  if (optarg)     error(E_OPTARG);    /* check (option) arguments */
  if (k    < 1)   error(E_ARGCNT);    /* and number of arguments */
  if (min  < 0)   error(E_SIZE, min); /* check the size limits */
  if (max  < 0)   error(E_SIZE, max); /* and the minimum support */
  if (supp > 100) error(E_SUPPORT, supp);
  if ((!fn_inp || !*fn_inp) && (fn_sel && !*fn_sel))
    error(E_STDIN);             /* stdin must not be used twice */
  switch (target) {             /* check and translate target type */
    case 's': target = ISR_ALL;              break;
    case 'c': target = ISR_CLOSED;           break;
    case 'm': target = ISR_MAXIMAL;          break;
    case 'g': target = ISR_GENERA;           break;
    default : error(E_TARGET, (char)target); break;
  }                             /* (get target type code) */
  switch (algo) {               /* check and translate alg. variant */
    case 'e': algo = ECL_BASIC;              break;
    case 'i': algo = ECL_LISTS;              break;
    case 't': algo = ECL_TABLE;              break;
    case 'r': algo = ECL_RANGES;             break;
    case 'b': algo = ECL_BITS;               break;
    case 's': algo = ECL_SIMPLE;             break;
    case 'o': algo = ECL_OCCDLV;             break;
    case 'd': algo = ECL_DIFFS;              break;
    default : error(E_VARIANT, (char)algo);  break;
  }                             /* (get eclat algorithm code) */
  if ((eval != 'x') && (eval != 'b'))
    error(E_MEASURE,(char)eval);/* check evaluation measure */
  if (fn_tid) {                 /* turn "-" into "" for consistency */
    if (strcmp(fn_tid, "-") == 0) fn_tid = ""; }
  if ((format == dflt) && (supp < 0))
    format = " (%a)";           /* adapt the default info. format */
  if ((cmfilt >= 0) && (target & (ISR_CLOSED|ISR_MAXIMAL)))
    mode |= ISR_NOFILTER | ((cmfilt > 0) ? ECL_VERTICAL : 0);
  MSG(stderr, "\n");            /* terminate the startup message */

  /* --- make algorithm and modes consistent --- */
  if (fn_tid)                   /* if file for transaction ids given, */
    mode |= ECL_TIDOUT;         /* set transaction identifier flag */
  ecl_adjust(target, eval, &algo, &mode, &pack, &mrep);
  if (mode & ECL_REORDER)       /* simplified sorting if reordering */
    sort = (sort < 0) ? -1 : (sort > 0) ? +1 : 0;

  /* --- read item selection --- */
  ibase = ib_create(0, 0);      /* create an item base */
  if (!ibase) error(E_NOMEM);   /* to manage the items */
  tread = trd_create();         /* create a transaction reader */
  if (!tread) error(E_NOMEM);   /* and configure the characters */
  trd_allchs(tread, recseps, fldseps, blanks, "", comment);
  if (fn_sel) {                 /* if item appearances are given */
    t = clock();                /* start timer, open input file */
    if (trd_open(tread, NULL, fn_sel) != 0)
      error(E_FOPEN, trd_name(tread));
    MSG(stderr, "reading %s ... ", trd_name(tread));
    m = ib_readsel(ibase,tread);/* read the given item selection */
    if (m < 0) error((int)-m, ib_errmsg(ibase, NULL, 0));
    trd_close(tread);           /* close the input file */
    MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
    MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  }                             /* print a log message */

  /* --- read transaction database --- */
  tabag = tbg_create(ibase);    /* create a transaction bag */
  if (!tabag) error(E_NOMEM);   /* to store the transactions */
  t = clock();                  /* start timer, open input file */
  if (trd_open(tread, NULL, fn_inp) != 0)
    error(E_FOPEN, trd_name(tread));
  MSG(stderr, "reading %s ... ", trd_name(tread));
  k = tbg_read(tabag, tread, mtar);
  if (k < 0)                    /* read the transaction database */
    error(-k, tbg_errmsg(tabag, NULL, 0));
  trd_delete(tread, 1);         /* close the input file and */
  tread = NULL;                 /* delete the table reader */
  m = ib_cnt(ibase);            /* get the number of items, */
  n = tbg_cnt(tabag);           /* the number of transactions, */
  w = tbg_wgt(tabag);           /* the total transaction weight */
  MSG(stderr, "[%"ITEM_FMT" item(s), %"TID_FMT, m, n);
  if (w != (SUPP)n) MSG(stderr, "/%"SUPP_FMT, w);
  MSG(stderr, " transaction(s)] done [%.2fs].", SEC_SINCE(t));
  if ((m <= 0) || (n <= 0))     /* check for at least one item */
    error(E_NOITEMS);           /* and at least one transaction */
  MSG(stderr, "\n");            /* compute absolute support value */
  supp = ceilsupp((supp >= 0) ? 0.01 *supp *(double)w : -supp);

  /* --- sort and recode items --- */
  t = clock();                  /* start timer, print log message */
  MSG(stderr, "filtering, sorting and recoding items ... ");
  m = tbg_recode(tabag, (SUPP)supp, -1, -1, -sort);
  if (m <  0) error(E_NOMEM);   /* recode items and transactions */
  if (m <= 0) error(E_NOITEMS); /* and check the number of items */
  MSG(stderr, "[%"ITEM_FMT" item(s)]", m);
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));

  /* --- sort and reduce transactions --- */
  t = clock();                  /* start timer, print log message */
  s = (fn_tid || (algo == ECL_BITS)) ? "filter" : "sorting and reduc";
  MSG(stderr, "%sing transactions ... ", s);
  tbg_filter(tabag,min,NULL,0); /* remove items of short transactions */
  i = ((algo == ECL_RANGES) || (algo == ECL_OCCDLV)) ? +1 : -1;
  tbg_itsort(tabag, i, 0);      /* sort items in transactions */
  if ((algo == ECL_RANGES)      /* if to use transaction ranges */
  &&  (mode &  ECL_FIM16)) {    /* together with a 16-items machine */
    tbg_pack(tabag, pack);      /* pack the most frequent items */
    tbg_sort(tabag,1,TA_EQPACK);/* sort trans. lexicographically and */
    n = tbg_reduce(tabag, 0); } /* reduce transactions to unique ones */
  else if (!fn_tid              /* if not to report transaction ids */
  &&      (algo != ECL_BITS)) { /* and not to use bit vectors */
    tbg_sort(tabag, i, 0);      /* sort trans. lexicographically and */
    n = tbg_reduce(tabag, 0);   /* reduce transactions to unique ones */
    if (mode & ECL_FIM16)       /* if to use a 16-items machine, */
      tbg_pack(tabag, pack);    /* pack the most frequent items */
  }                             /* (bit-represented transactions) */
  MSG(stderr, "[%"TID_FMT, n);  /* print number of transactions */
  if (w != (SUPP)n) MSG(stderr, "/%"SUPP_FMT, w);
  MSG(stderr, " transaction(s)] done [%.2fs].\n", SEC_SINCE(t));

  /* --- find frequent item sets --- */
  t = clock();                  /* start the timer */
  report = isr_create(ibase, target|mrep, -1, hdr, sep, NULL);
  if (!report) error(E_NOMEM);  /* create an item set reporter */
  isr_setfmt (report, format);  /* and configure it: set flags, */
  isr_setsize(report, min, max);/* info. format and size range, */
  if (fn_psp && (isr_addpsp(report, NULL) < 0))
    error(E_NOMEM);             /* add a pattern spectrum if req. */
  if (isr_open   (report, NULL, fn_out) != 0)
    error(E_FOPEN, isr_name(report));    /* open the item set file */
  if (isr_tidopen(report, NULL, fn_tid) != 0)
    error(E_FOPEN, isr_tidname(report)); /* open the trans. id file */
  MSG(stderr, "writing %s ... ", isr_name(report));
  k = eclat(tabag, target, algo, mode, (SUPP)supp,
            eval, 0.01*minval, report);
  if (k < 0) error(E_NOMEM);    /* search for frequent item sets */
  if (isr_tidclose(report) != 0)/* close the trans. id output file */
    error(E_FWRITE, isr_tidname(report));
  if (isr_close   (report) != 0)/* close the item set output file */
    error(E_FWRITE, isr_name(report));
  MSG(stderr, "[%"SIZE_FMT" set(s)]", isr_repcnt(report));
  MSG(stderr, " done [%.2fs].\n", SEC_SINCE(t));
  if (stats) isr_prstats(report, stdout, 0);

  /* --- write pattern spectrum --- */
  if (fn_psp) {                 /* if to write a pattern spectrum */
    psp    = isr_getpsp(report);/* get the pattern spectrum */
    twrite = twr_create();      /* create a table writer and */
    if (!twrite) error(E_NOMEM);/* open the output file */
    if (twr_open(twrite, NULL, fn_psp) != 0)
      error(E_FOPEN,  twr_name(twrite));
    if (psp_report(psp, twrite) != 0)
      error(E_FWRITE, twr_name(twrite));
    twr_delete(twrite, 1);      /* write the pattern spectrum, */
    twrite = NULL;              /* delete the table writer, and */
  }                             /* clear the writer variable */

  /* --- clean up --- */
  CLEANUP;                      /* clean up memory and close files */
  SHOWMEM;                      /* show (final) memory usage */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif
