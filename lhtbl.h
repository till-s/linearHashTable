#ifndef LHTBL_LINPROBE_H
#define LHTBL_LINPROBE_H
/* $Id$ */

/* Linear Hashing with open addressing and linear probing                                              */

/* Author: Till Straumann <strauman@slac.stanford.edu>, 2009                                           */

/* C.f. 'LINEAR HASHING; A NEW TOOL FOR FILE AND TABLE ADDRESSING',
 *      Witold Litwin (INRIA), 1980.
 */
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle for a hash table                                                                      */
typedef struct    LHTblRec_ *LHTbl;
/* A key; always use this type since definition may change in the future...                            */
typedef uint32_t  LHTblKey;
/* Opaque handle for table entry; layout is defined by the user but all entries must have the 'key'
 * member at the same offset.
 */
typedef void      *LHTblEntry;

/* Error return codes                                                                                  */
#define LHTBL_FULL			(-1) /* no more memory could be allocated or all indirect blocks are full  */
#define LHTBL_KEY_EXISTS	(-2) /* entry with 'key' already exists (lhtblAdd)                         */
#define LHTBL_KEY_NOTFND	(-3) /* entry with 'key' not found                                         */
#define LHTBL_ARG_TOOSMALL  (-4) /* supplied argument too small                                        */
#define LHTBL_ARG_TOOBIG    (-5) /* supplied argument too big                                          */

/* Configuration parameters for lhtblCreate().
 *
 * auto_manage_enable:           If automatic management is enabled then the insert and replace
 *                               operations are allowed to 'split' the table if needed and the delete
 *                               operation is allowed to 'group' the table if needed.
 *
 *                               The criteria for implicit splitting and grouping are as follows:
 *
 *                               IF ('insert' OR ( 'replace' AND 'new_entry_is_added' ) ) THEN
 *                                   IF ('auto_manage_enable') THEN
 *                                      IF (    ( 'current_load' > 'load_percentage_hi' / 100 ) ) 
 *                                           OR ( 'new_overflow_bucket_needed' ) ) THEN
 *                                        lhtblSplit();
 *                                      END
 *                                   END
 *                               END
 *
 *
 *                               IF ( 'delete' ) THEN
 *                                   IF ('auto_manage_enable') THEN
 *                                      IF (    ( 'current_load' < 'load_percentage_lo' / 100 ) ) THEN
 *                                        lhtblGroup();
 *                                      END
 *                                   END
 *                               END
 *
 *                               'current_load' is defined as the ratio of 'used/occupied' buckets
 *                               to 'primary buckets'.
 *
 *                               NOTE: if 'auto_manage_enable' is true but the thresholds
 *                                     'load_percentage_thres_hi == 100' and 
 *                                     'load_percentage_thres_lo == 0' then automatic splitting
 *                                     and grouping for load-balancing is disabled. However,
 *                                     automatic expansion of the overflow area is still possible.
 *
 * load_percentage_thres_hi:     Threshold for automatic load management in percent (0..100).
 *                               Note that '0' is equivalent to '100', i.e., implicit splitting
 *                               is DISABLED.
 *                               See 'auto_manage_enable' for how the algorithm works.
 *
 * load_percentage_thres_lo:     Threshold for automatic load management in percent (0..100).
 *                               See 'auto_manage_enable' for how the algorithm works.
 *
 * ovfl_buckets_preserve:        At initialization, make sure at least 'ovfl_buckets_preserve'
 *                               buckets are available in the overflow area w/o the need for
 *                               table expansion.
 *
 * NOTE: Automatic load management distributes the computational effort of 'split/group' 
 *       operations, i.e., each 'insert/replace' operation only does O(1) split if necessary
 *       and each 'delete' operation does O(1) group if necessary.
 *       Since the table size is managed in blocks this does not optimize the use of
 *       available space: once a table expansion is necessary a whole block of new buckets
 *       is added but they only gradually become available as 'splits' are executed.
 *       It would be better to split the whole block (while reserving a certain amount
 *       of buckets for the overflow area) as soon as it becomes available but that
 *       would be computationally more expensive.
 *
 *       A solution would be a separate (low-priority) thread which takes care of load-
 *       management. Such a thread could monitor the load and expand the table as
 *       necessary. After expanding the table the thread would execute a number of
 *       split operations to make use of the additional space (while reserving some
 *       buckets for the overflow area).
 *       In a similar fashion the thread would group/collapse the table.
 */                                                                   

typedef struct LHTblConfigRec_ {
	int  auto_manage_enable;		/* enable automatic splitting (ins/rpl) and grouping (del)         */
	int  load_percentage_thres_hi;  /* load-threshold (in percent 0..100) for splitting (0 == disable) */
	int  load_percentage_thres_lo;  /* load-threshold (in percent 0..100) for grouping  (0 == disable) */
	int  ovfl_buckets_preserve;     /* how much buckets to preserve as a overflow area                 */
} LHTblConfigRec, *LHTblConfig;

/*
 * Create a new hash-table with (initially) 'n_buckets'
 * primary buckets.
 *
 * 'key_off' is the byte-offset where the LHTblKey can be
 * located in the (user-defined) table-entry struct.
 *
 * 'cfg' points to a table of configuration parameters
 * (see above). 'cfg' may be NULL in which case a set of
 * defaults is used:
 *
 *      auto_manage_enable       = yes;
 *      load_percentage_thres_hi =  80;
 *      load_percentage_thres_lo =   0; // never shrinks
 *      ovfl_buckets_preserve    =  20;
 *
 * RETURNS: new hash table on success, NULL on failure. 
 */
LHTbl
lhtblCreate(unsigned n_buckets, unsigned long key_off, LHTblConfig cfg);

/*
 * Destroy a hash table releasing all associated resources.
 *
 * If the 'cleanup' function-pointer is non-NULL then
 * All buckets are scanned and the cleanup is executed on
 * all non-empty buckets.
 *
 * NOTE: table destruction is not thread-safe.
 */
void
lhtblDestroy(LHTbl t, void (*cleanup)(LHTblEntry, void *closure), void *closure);

/*
 * Locate the entry associated with a given key
 * in a hash table.
 *
 * RETURNS: entry or NULL if the key was not found.
 */
LHTblEntry
lhtblFind(LHTbl t, LHTblKey key);

/*
 * Add a new entry to a table.
 *
 *
 * RETURNS: 0 on success, nonzero on failure.
 *          Reasons for failure can e.g., be
 *
 *            LHTBL_KEY_EXISTS: entry with key already present.
 *            LHTBL_FULL:       no more memory available or
 *                              no buckets available (if automatic
 *                              management is disabled).
 *
 *          The current table entry is returned in *pe,
 *          i.e., if the routine returns LHTBL_KEY_EXISTS then
 *          the old, existing entry is returned in *pe.
 * 
 * NOTES:   If 'automatic' table management is enabled
 *          and the 'high' threshold is set > 0 then 
 *          this routine splits the table if it would
 *          raise the load above the 'high' threshold.
 *
 *          The table is also split if 'automatic' management
 *          is enabled and new overflow-buckets past
 *          the current end of the table are needed to
 *          add the new entry.
 */
int
lhtblFindAdd(LHTbl t, LHTblEntry *pe);

/*
 * Helper entry point for lhtblFindAdd()
 *
 *  lhtblAdd(LHTbl t, LHTblEntry e)
 *  {
 *     return lhtblFindAdd(t, &e);
 *  }
 */

int
lhtblAdd(LHTbl t, LHTblEntry e);

#define LHTBL_ADD_FAIL  1

/*
 * Replace a table entry.
 * 
 * The key of entry (*pe) is looked up and if found
 * then the current entry is replaced by (*pe) and
 * the old entry passed back to the user in (*pe).
 *
 * If the key is not found and 'add_fail' is nonzero
 * then the routine RETURNS SHTBL_KEY_NOTFND.
 * If 'add_fail' is false then the routine behaves
 * like 'lhtblAdd()' (including automatic table
 * splitting/expansion, load control etc.)
 *
 * RETURNS: zero on success, nonzero on error.
 *          Possible errors are SHTBL_KEY_NOTFND
 *          (if 'add_fail' is true) of the errors
 *          returned by lhtblAdd().
 */
int
lhtblRpl(LHTbl t, LHTblEntry *pe, int add_fail);

/*
 * Remove an entry from the table.
 *
 * RETURNS: 0 on success, LHTBL_KEY_NOTFND if
 *          the entry could not be located.
 *
 * NOTES:   If automatic table management is enabled
 *          and the 'low' threshold is set > 0 then
 *          this routine 'groups' the table thus
 *          reducing its size so that the load
 *          increases above the 'low' threshold.
 *          
 *          If the size is reduced sufficiently then
 *          memory can be free()d.
 */
int
lhtblDel(LHTbl t, LHTblEntry e);

/*
 * Compute and RETURN an upper-bound for the number
 * of buckets that would have to be added to the
 * table if it were split.
 */
int
lhtblCntUpperBnd(LHTbl t);

/*
 * Split operation. Increase the number of primary buckets
 * by one; as a side-effect, more memory might have to be
 * allocated and added to the table.
 *
 *
 * RETURNS: zero on success, LHTBL_FULL if the size of the
 *          table could not be increased.
 *
 * NOTES:   The table is unlocked while memory is allocated.
 */
int
lhtblSplit(LHTbl t);

/*
 * Group table, i.e., decrement the number of primary slots
 * by one; as a side-effect, memory might be removed from
 * the table and released.
 *
 * RETURNS: zero on success. Current implementation cannot fail.
 *
 * NOTES:   The table is unlocked while memory is released.
 *
 */
int
lhtblGroup(LHTbl t);

/*
 * Expand the table by allocating memory so that it can
 * hold at least 'x' buckets.
 *
 * RETURNS: 0 on success, LHTBL_FULL if no memory could
 *          be allocated or the block array is full.
 *
 * NOTES:   If allocation is necessary then this routine
 *          RELEASES THE LOCK during allocation thus
 *          allowing other tasks to access the hash
 *          table.
 *          Only after re-acquiring the lock is the new
 *          memory actually attached to the table.
 *
 *          Under pathological conditions two threads
 *          t1 and t2 could probably 'starve' each other
 *          if t1 tries to expand while t2 is collapsing
 *          the table because the table shrinks while
 *          the lock is released forcing more while-loop
 *          iterations (t2 would loop in a similar way).
 */
int
lhtblExpand(LHTbl t, int n_buckets);

/*
 * Collapse the table freeing up memory to 'x'
 * buckets (or more -- as much memory as possible
 * is released to reduce the number of available
 * slots to 'x').

 * RETURNS: The new max. # of buckets.
 *
 * NOTES:   The algorithm collects unused blocks
 *          while holding the lock but 'free()' is
 *          only called after the lock has been
 *          released.
 */
int
lhtblCollapse(LHTbl t, int n_buckets);

/* Statistics */
typedef struct LHTblStatsRec_ {
	int         tble_level;   /* current 'level' of table                                        */
	int         prim_buckets; /* number of primary buckets in table                              */
	int         ovfl_buckets; /* number of 'overflow' buckets past table holding collision chain */
	int         avlb_buckets; /* total number of available buckets; table can grow to this number
                               * w/o having to allocate more memory.
                               */ 
	int         used_buckets; /* used/occupied buckets                                           */
} LHTblStatsRec, *LHTblStats;

/* Obtain statistics; user must supply a pointer to
 * a LHTblStatsRec area.
 */
void
lhtblGetStats(LHTbl t, LHTblStats stats);

/* Obtain statistics and dump to a file; 
 * 'f' may be NULL in which case 'stdout'
 * is used.
 */
void
lhtblDumpStats(LHTbl t, FILE *f);
#ifdef __cplusplus
}
#endif

#endif
