/*----------------------------------------------------------------------
  File    : eclat.h
  Contents: eclat algorithm for finding frequent item sets
  Author  : Christian Borgelt
  History : 2011.08.22 file created
            2011.08.31 occurrence deliver variant ECL_OCCDLV added
            2012.04.17 interface for external call added (e.g. python)
            2012.06.20 function ecl_adjust() added (consistency check)
----------------------------------------------------------------------*/
#ifndef __ECLAT__
#define __ECLAT__
#include "tract.h"
#ifndef ISR_CLOMAX
#define ISR_CLOMAX
#endif
#include "report.h"

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
/* --- eclat variants --- */
#define ECL_BASIC     0         /* tid lists intersection (basic) */
#define ECL_LISTS     1         /* tid lists intersection (improved) */
#define ECL_BITS      2         /* bit vectors over transactions */
#define ECL_TABLE     3         /* item occurrence table (standard) */
#define ECL_SIMPLE    4         /* item occurrence table (simplified) */
#define ECL_RANGES    5         /* tid range lists intersection */
#define ECL_OCCDLV    6         /* occurrence deliver (LCM-style) */
#define ECL_DIFFS     7         /* tid difference sets (diffsets) */

/* --- operation modes --- */
#define ECL_FIM16     0x0100    /* use 16 items machine (bit rep.) */
#define ECL_PERFECT   0x0200    /* perfect extension pruning */
#define ECL_REORDER   0x0400    /* reorder items in cond. databases */
#define ECL_TAIL      0x0800    /* use head union tail pruning */
#define ECL_VERTICAL  0x1000    /* vertical extensions tests */
#define ECL_TIDOUT    0x2000    /* flag for trans. identifier output */
#define ECL_DEFAULT   (ECL_FIM16|ECL_PERFECT|ECL_REORDER|ECL_TAIL)
#define ECL_VERBOSE   INT_MIN   /* verbose message output */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern void ecl_adjust (int target, int eval,
                        int *algo, int *mode, int *pack, int *mrep);
extern int  eclat      (TABAG *tabag, int target, int algo, int mode,
                        SUPP supp, int eval, double minval,
                        ISREPORT *report);
#endif
