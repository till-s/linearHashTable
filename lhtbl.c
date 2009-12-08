/* $Id$ */

/* Author: Till Straumann <strauman@slac.stanford.edu>, 2009    */

/* C.f. 'LINEAR HASHING; A NEW TOOL FOR FILE AND TABLE ADDRESSING',
 *      Witold Litwin (INRIA), 1980.
 */
#include <inttypes.h>

#ifdef DEBUG
#undef NDEBUG
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lhtbl.h>

/* Printf format for key               */
#define KEYu	PRIu32

/* Locking -- implementation dependent */
#ifdef USE_EPICS

#include <epicsMutex.h>
#include <epicsInterrupt.h>

#define  LOCK_T	              epicsMutexId
#define  INT_LEVEL            int
#define  NULL_MUTEX           (0)

#define  __INT_LOCK(lvl)      do { lvl = epicsInterruptLock();        } while (0)
#define  __INT_UNLOCK(lvl)    do { epicsInterruptUnlock( lvl );       } while (0)

#define __LHT_LOCK_CHECK(m)   ( epicsMutexLockOK != epicsMutexLock(m) )
#define __LHT_UNLOCK_CHECK(m) ( epicsMutexUnlock(m), 0                )

#define __LHT_LOCK_CREATE()   (epicsMutexCreate()                     )
#define __LHT_LOCK_DELETE(m)  do { epicsMutexDestroy((m));            } while (0)

#elif defined(__rtems__)

#include <rtems.h>

#define LOCK_T               rtems_id
#define INT_LEVEL            rtems_interrupt_level
#define NULL_MUTEX           (0)

#define  __INT_LOCK(lvl)     do { rtems_interrupt_disable(lvl);       } while (0)
#define  __INT_UNLOCK(lvl)   do { rtems_interrupt_enable(lvl);        } while (0)

static __inline__ int
__LHT_LOCK_CHECK(rtems_id m)
{
	return (   RTEMS_SUCCESSFUL
            != rtems_semaphore_obtain(m, RTEMS_WAIT, RTEMS_NO_TIMEOUT)
           );
}

static __inline__ int
__LHT_UNLOCK_CHECK(rtems_id m)
{
	return (   RTEMS_SUCCESSFUL
            != rtems_semaphore_release(m)
           );
}

static __inline__ rtems_id
__LHT_LOCK_CREATE()
{
rtems_id           m;
rtems_status_code st;
	st = rtems_semaphore_create(
            rtems_build_name('L','H','T','M'),
            1,
            RTEMS_PRIORITY            |
            RTEMS_BINARY_SEMAPHORE    |
            RTEMS_INHERIT_PRIORITY    |
            RTEMS_NO_PRIORITY_CEILING |
            RTEMS_LOCAL,
            0,
            &m);
	return RTEMS_SUCCESSFUL == st ? m : NULL_MUTEX;
}

static inline void
__LHT_LOCK_DELETE(rtems_id m)
{
	rtems_semaphore_delete(m);
}

#else
 
#define LOCK_T                int
#define NULL_MUTEX            (0)

#define __LHT_LOCK_CHECK(m)   (0)
#define __LHT_UNLOCK_CHECK(m) (0)

#define __LHT_LOCK_CREATE()   (1)
#define __LHT_LOCK_DELETE(m)  do { } while (0)

#endif

typedef int      Hash;

#if ! defined(HFUNC)
/* use Knuth 'golden-ratio' and try to get upper bits to lower position */
#define HFUNC(k) __builtin_bswap32( (k) * 2654435769U & 0xffffffffU )
#endif

#if ! defined(HFUNC) || defined(TESTING)
#warning "Using trivial / testing hash function"
#undef  HFUNC
#define HFUNC(k)  (k)
#endif

/* Memory management for the hash table.
 *
 * We use a simple approach, a fixed-size
 * array of pointers to blocks of memory
 * which are all of the same size which
 * must be a power of two.
 *
 * A more sophisticated approach could be
 * implemented -- this is transparent to
 * most of the code and can be changed by
 * modifying the lhtblExpand() and lhtblCollapse()
 * routines as well as the inlines for
 * indexing the table (GENTRY() and GENTRY_P())
 */

#ifdef TESTING
#define LDBLOCKSZ   1
#define COLL_LOWAT	0
#else
#define LDBLOCKSZ   10   /* log2 of block size */
#define COLL_LOWAT	16   /* 'low-water' mark for the 'collapse' routine:
                          * memory is only released if the available
                          * number of buckets - COLL_LOWAT is bigger
                          * than the block size.
                          */
#endif

#define TEBLOCKSZ	(1<<(LDBLOCKSZ))
#define MAXBLOCKS   1024 /* Max. number of blocks; this limits the max.
                          * number of buckets to MAXBLOCKS * TEBLOCKSZ - 1
                          */

/* Index of the last block computed by up-aligning the
 * current max. table size, which may be less than a multiple
 * of blocks because we need to store a sentinel NULL pointer
 * (see below).
 */
#define LSTBLKIDX(t) (((t)->bckts_max + TEBLOCKSZ - 1) >> LDBLOCKSZ)

/* We always want space for one NULL pointer
 * 'behind' the table which marks the end of
 * a possible chain off the last primary
 * slot:
 *
 *  bckts_prm         bckts_tot
 *
 *  last              last   
 *  primary           bucket
 *  bucket                  sentinel
 *    |   
 * |  v  |     |     |  v  |     |
 * |  w  |  x  |  y  |  z  | NUL |
 *
 *    ^.................^
 *     chain in overflow area
 */
#define NULLSENTINEL 1 /* reserve one slot for NULL sentinel                       */

/* A block of buckets */
typedef LHTblEntry TEBlock[TEBLOCKSZ];

/* Data structure describing the hash table */
typedef struct LHTblRec_ {
	int          table_lvl;  /* table 'level'                                       */
	Hash         split_idx;  /* split index for dynamic hash function               */
	Hash         mod2l_msk;  /* bitmask to compute modulo 2^level                   */
	unsigned     offst_key;  /* byte-offset of key in entry records                 */
	int          bckts_tot;  /* total # of buckets in table (incl. overflow at end) */
	int          bckts_prm;  /* # of primary slots/buckets in table                 */
	int          bckts_max;  /* available memory for table                          */
	int          bckts_usd;  /* # of used/occupied buckets                          */
	int          thres_hig;  /* Hi load threshold for splitting table (percent)     */
	int          thres_low;  /* Lo load threshold for grouping table (percent)      */
	int          mngmt_aut;  /* Enable automatic split/group/expand/collapse        */
	volatile int write_prg;  /* A table-modifying operation is in progress          */
	LOCK_T       mutex_hdl;  /* Mutex to protect table from concurrent modification */
	TEBlock     *block_arr[MAXBLOCKS];
} LHTblRec;

/* Macro to locate the key in the (user-defined) entry struct                       */
#define KEY(t, e) (*(LHTblKey*)((void*)(e) + (t)->offst_key))

static __inline__ void __LHT_LOCK_R(LHTbl t)
{
	if ( __LHT_LOCK_CHECK(t->mutex_hdl) ) {
		fprintf(stderr,"lhtbl: FATAL ERROR; unable to lock mutex\n");
		abort();
	}
}

static __inline__ void __LHT_UNLOCK_R(LHTbl t)
{
	if ( __LHT_UNLOCK_CHECK(t->mutex_hdl) ) {
		fprintf(stderr,"lhtbl: FATAL ERROR; unable to unlock mutex\n");
		abort();
	}
}

static __inline__ void __LHT_LOCK(LHTbl t)
{
	/* lock and mark */
	__LHT_LOCK_R(t);
	t->write_prg = 1;
	__asm__ __volatile__("":::"memory");
}

static __inline__ void __LHT_UNLOCK(LHTbl t)
{
	/* unmark and unlock */
	__asm__ __volatile__("":::"memory");
	t->write_prg = 0;
	__LHT_UNLOCK_R(t);
}


#ifndef NDEBUG
static  int      debug  = 
#ifdef TESTING
	1
#else
	0
#endif
;
#endif

/* Inlines to locate bucket number 'idx' in memory (the
 * layout of which is opaque to most of the code).
 */

/* Locate entry 'idx'                                                              */
static __inline__ LHTblEntry GENTRY(LHTbl t, Hash idx)
{
	return (*t->block_arr[idx>>LDBLOCKSZ])[idx & (TEBLOCKSZ-1)];
}

/* Locate a pointer to entry 'idx'                                                 */
static __inline__ LHTblEntry *GENTRY_P(LHTbl t, Hash idx)
{
	return & (*t->block_arr[idx>>LDBLOCKSZ])[idx & (TEBLOCKSZ-1)];
}

/* binary search for MSB */
static int ld2up(unsigned x)
{
int rval  = 0;
int bitno = sizeof(x) * 8 / 2;
unsigned m = -1;

	if ( 0 == x )
		return -1;

	m = m ^ (m>>bitno);	/* 0xffffffff ^ 0x0000ffff = 0xffff0000 */

	while ( bitno ) {
		if ( m & x ) {
			rval += bitno;
			x >>= bitno;
		}
		bitno>>=1;
		m>>=bitno;
	}
	return rval;
}

/*
 * Expand the table by allocating memory so that it can
 * hold at least 'x' buckets.
 *
 * RETURNS: 0 on success, LHTBL_FULL if no memory could
 *          be allocated or the block array (MAXBLOCKS)
 *          is full.
 *
 * NOTES:   The caller is supposed to hold the lock when
 *          executing this routine.
 *          If allocation is necessary then this routine
 *          RELEASES THE LOCK during allocation thus
 *          allowing other tasks to access the hash
 *          table.
 *          Only after re-acquiring the lock is the new
 *          memory actually attached to the table.
 *
 *          This routine holds the lock on return.
 *          
 *          Under pathological conditions two threads
 *          t1 and t2 could probably 'starve' each other
 *          if t1 tries to expand while t2 is collapsing
 *          the table because the table shrinks while
 *          the lock is released forcing more while-loop
 *          iterations (t2 would loop in a similar way).
 */
static int
expand_unlocked(LHTbl t, int x)
{
int      blkidx, nblks;
void    *blk_list = 0;
TEBlock *blk;
int      blks     = 0;

#ifndef NDEBUG
	if ( debug )
		printf("New table size; expanding to %u\n",x);
#endif

	while ( x > t->bckts_max ) {
		/* Need to obtain more memory; calculate the total # of
         * blocks needed (including the ones already existing).
         * The calculation takes into account the space needed
		 * for the 'sentinel' NULL pointer.
         */
		nblks  = ((x + NULLSENTINEL) + TEBLOCKSZ - 1) >> LDBLOCKSZ;

		/* index of last (currently existing) block           */
		blkidx = LSTBLKIDX(t);

		__LHT_UNLOCK(t);

		/* pre-allocate w/o holding the lock                  */
		for ( ; blkidx < nblks; blkidx++ ) {
			if (    blkidx  >= MAXBLOCKS 
					|| ( ! (blk = calloc(1, sizeof(*blk))) ) ) {
				fprintf(stderr,"lhtbl: unable to allocate more memory to expand table\n");
				/* force the outmost 'while' loop to terminate */
				x = -2;
				break;
			}
			/* link into a list */
			(*blk)[0] = blk_list;
			blk_list  = blk;
			blks++;
		}

		__LHT_LOCK(t);

		/* Things might have changed during the time we relinquished
         * the lock; must recompute ...
		 */
		blkidx = LSTBLKIDX(t);

		/* Append what we have */
		for ( ; blkidx < MAXBLOCKS && blk_list; blkidx++ ) {
			blk       = t->block_arr[blkidx] = blk_list;
			blk_list  = (*blk)[0];
			(*blk)[0] = 0;
			t->bckts_max += TEBLOCKSZ;
			blks--;
		}

		/* It is possible that someone else has expanded so that
		 * not all the blocks we allocated are really needed.
		 */
		if ( blk_list ) {
			__LHT_UNLOCK(t);
			while ( blk_list ) {
				blk = blk_list;
				blk_list = (*blk)[0];
				free(blk);
			}
			__LHT_LOCK(t);
		}

		/* might have to try over again -- this could
		 * be necessary if someone else collapsed while
		 * we were allocating blocks (with the lock released).
		 * NOTE: There could be pathological scenarios where
		 *       two tasks 'starve' each other with one collapsing
		 *       and the other one expanding at the same time...
	     */
	}

	return -2 == x ? LHTBL_FULL : 0;
}

/* Locked wrapper around expand_unlocked()
 *
 * See comments for expand_unlocked() for
 * RETURN value and NOTES.
 */
int
lhtblExpand(LHTbl t, int x)
{
int rval;
	__LHT_LOCK(t);
		rval = expand_unlocked(t, x);
	__LHT_UNLOCK(t);
	return rval;
}

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
lhtblCollapse(LHTbl t, int x)
{
int  blkidx, rval;
void *to_free = 0;

#ifndef NDEBUG
	if ( debug )
		printf("New table size; collapsing to %u\n",x);
#endif

	__LHT_LOCK(t);

	if ( x < t->bckts_tot ) {
		/* Cannot collapse the used area */
		x = t->bckts_tot;
	}

	while ( x <= (t->bckts_max - TEBLOCKSZ - COLL_LOWAT ) ) {
		/* can release a block */
		t->bckts_max -= TEBLOCKSZ;
		blkidx = LSTBLKIDX(t);
		/* build a linked list we can then free once the
		 * lock is released.
		 */
		(*t->block_arr[blkidx])[0] = to_free;
		to_free = t->block_arr[blkidx];
		t->block_arr[blkidx] = 0;
	}

	rval = t->bckts_max;

	__LHT_UNLOCK(t);

	while ( to_free ) {
		TEBlock *bp = to_free;
		to_free = (*bp)[0];
		free(bp);
	}

	return rval;
}

/* Initialize essential (and some redundant) parameters
 * of the table that can be derived from the initial
 * bucket count.
 */
static int
tblrst(LHTbl t, int nelm)
{
	/* make sure we have enough buckets by allocating them
     * if necessary.
     */
	if ( lhtblExpand(t, nelm) )
		return -1;

	/* The initial table level is the smallest power of two
     * that is bigger than the bucket count.
     * Note that if 'nelm' is itself a power of two then
     * the initial 2^level is *twice* the bucket count (nelm).
     */
	t->table_lvl = ld2up(nelm) + 1;
	/* Bitmask for fast computation of MOD 2^level
     * used e.g. by the 'dynhash' function.
     */
	t->mod2l_msk = (1<<t->table_lvl) - 1;
	/* Number of total buckets equals the number
     * of primary buckets plus any overflow buckets
	 * past the table which are appended due to
	 * collisions.
	 * Initially the table is empty and hence
	 * bckts_tot == bckts_prm
	 */
	t->bckts_tot = 
	/* Number of primary buckets = 'nelm' */
    t->bckts_prm = nelm;
	/* No buckets are used yet            */
	t->bckts_usd = 0;
	/* The initial split index is equal to
	 * the number of primary buckets past
	 * the first 'half' of the table 0 .. (2^lvl / 2 - 1)
     * minus one (0-based index).
	 *
	 * Since nelm < 2^lvl and nelm >= 2^(lvl-1) [otherwise
     * 'lvl-1' would already be the initial level] the
	 * split index 
	 *
	 *     -1 <= split_idx < 2^(lvl-1) - 1
	 *
	 * If nelm is a power of two then the split index
	 * is '-1' (no split).
	 */
	t->split_idx = t->bckts_prm - (1 << (t->table_lvl - 1)) - 1;
	return 0;
}

/*
 * Destroy a hash table releasing all associated resources.
 *
 * If the 'cleanup' function-pointer is non-NULL then
 * All buckets are scanned and the cleanup is executed on
 * all non-empty buckets.
 */
void
lhtblDestroy(LHTbl t, void (*cleanup)(LHTblEntry, void *), void *closure)
{
int        i;
LHTblEntry e;

	if ( cleanup ) {
		for ( i=0; i < t->bckts_tot; i++ ) {
			if ( (e = GENTRY(t, i)) ) {
				cleanup(e, closure);
			}
		}
	}
	for ( i=0; i < LSTBLKIDX(t); i++ ) {
		free( t->block_arr[i] );
	}

	if ( NULL_MUTEX != t->mutex_hdl ) {
		__LHT_LOCK_DELETE(t->mutex_hdl);
		t->mutex_hdl = NULL_MUTEX;
	}
	
	free(t);
}

static LHTblConfigRec dflt_cfg = {
	auto_manage_enable          :  1,
	load_percentage_thres_hi    : 80,
	load_percentage_thres_lo    :  0, /* never shrink */
	ovfl_buckets_preserve       : 20,
};

/*
 * Create a new hash-table with (initially) 'nelms' primary
 * buckets.
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
lhtblCreate(unsigned nelms, unsigned long key_off, LHTblConfig cfg)
{
LHTbl rval = calloc(1, sizeof(*rval));
int   res;

	if ( rval ) {

		if ( !cfg ) {
			cfg = &dflt_cfg;
		}

		rval->mngmt_aut = cfg->auto_manage_enable;

		rval->thres_hig = cfg->load_percentage_thres_hi;
		if ( rval->thres_hig < 0 || rval->thres_hig >= 100 )
			rval->thres_hig = 0;

		rval->thres_low = cfg->load_percentage_thres_lo;
		if ( rval->thres_low < 0 || rval->thres_low >= 100 )
			rval->thres_low = 0;

		if ( !rval->mngmt_aut ) {
			/* avoid computation of load */
			rval->thres_low = rval->thres_hig = 0;
		}

		res = cfg->ovfl_buckets_preserve;

		if ( res < 0 )
			res = 0;
		if ( res > TEBLOCKSZ - NULLSENTINEL - 1 )
			res = TEBLOCKSZ - NULLSENTINEL - 1;

		if ( TEBLOCKSZ - ((nelms + NULLSENTINEL) % TEBLOCKSZ) < res ) {
			nelms = ((nelms + (TEBLOCKSZ - 1)) & ~(TEBLOCKSZ - 1)) + 1;
		}

		rval->bckts_max = - NULLSENTINEL;

		/* Must have the lock before we can allocate memory */
		if ( NULL_MUTEX == (rval->mutex_hdl = __LHT_LOCK_CREATE()) ) {
			fprintf(stderr,"lhtbl: ERROR; unable to create mutex\n");
			lhtblDestroy(rval, 0, 0);
			return 0;
		}

		if ( tblrst(rval, nelms) ) {
			lhtblDestroy(rval, 0, 0);
			return 0;
		}

		rval->offst_key = key_off;

	}
	return rval;
}

/* The dynamic hash function:
 *  - compute h1 = hash modulo 2^(level)
 *  - compute h0 = hash modulo 2^(level-1) (strip MSB)
 *  - if h0 is in the area already split (h0 <= split_idx)
 *    then use h1 (0 .. 2^(level-1) + split_idx)
 */
static __inline__ Hash
dynhash(LHTblKey k, Hash mod2l_msk, Hash split_idx)
{
register Hash h1 = HFUNC(k) & mod2l_msk;
register Hash h0 = h1 & (mod2l_msk>>1);
	return h0 <= split_idx ? h1 : h0;
}

/*
 * Locate the entry associated with a given key
 * in a hash table.
 *
 * RETURNS: entry or NULL if the key was not found.
 */
LHTblEntry
lhtblFind(LHTbl t, LHTblKey key)
{
Hash          h;
LHTblEntry    e;
#ifdef __INT_LOCK
INT_LEVEL     lvl;
int           have_lock;
#endif

	/* We want to make this routine really fast; just disabling interrupts may
	 * be cheaper than taking and releasing a mutex.
	 * However, if we execute this while some other thread holds the mutex
	 * then we must be more careful: re-enable interrupts and do a 'normal'
	 * mutex lock/unlock.
	 */
#ifdef __INT_LOCK
	__INT_LOCK(lvl);

	if ( (have_lock = t->write_prg) ) {
		/* write is in progress -- we must wait for the 'normal' lock */
		__INT_UNLOCK(lvl);

		__LHT_LOCK_R(t);
	}
#else
	__LHT_LOCK_R(t);
#endif

	for ( h = dynhash(key,t->mod2l_msk,t->split_idx); (e = GENTRY(t, h)); h++ )
		if ( key == KEY(t, e) )
			break;

#ifdef __INT_LOCK
	if ( ! have_lock )
		__INT_UNLOCK(lvl);
	else
#endif
		__LHT_UNLOCK_R(t);

	return e;
}

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
lhtblFindAdd(LHTbl t, LHTblEntry *pe)
{
Hash       h;
LHTblEntry e_h, e = *pe;
int        split_necessary;
int        st;


	do {

		/* Load - controlled split;
		 * ASSUME: we can atomically read these variables and obtaining
		 *         a slightly inaccurate 'usd' count due to race-conditions
		 *         doesn't matter. Cheaper than a LOCK - UNLOCK sequence.
		 */

		split_necessary = t->mngmt_aut && t->thres_hig > 0 && (t->bckts_usd + 1) * 100 > t->thres_hig * t->bckts_prm;

		if ( split_necessary ) {
			if ( (st = lhtblSplit(t)) )
				return st;
			split_necessary = 0;
		}

		__LHT_LOCK(t);
		for ( h = dynhash(KEY(t, e),t->mod2l_msk,t->split_idx); (e_h = GENTRY(t, h)); h++ ) {
			if ( KEY(t, e_h) == KEY(t, e) ) {
				__LHT_UNLOCK(t);
				*pe = e;
				return LHTBL_KEY_EXISTS;
			}
		}

		if ( h >= t->bckts_tot ) {
			assert( h == t->bckts_tot );

			/* would need to add the new entry past the table;
			 * see if we have space.
			 */

			/* if new bucket requires more memory than is available
			 * then split and try again
			 */
			if (  h >= t->bckts_max ) {

				split_necessary = t->mngmt_aut;

				__LHT_UNLOCK(t);

				if ( split_necessary ) {
					if ( (st = lhtblSplit(t)) )
						return st;
				} else {
					return LHTBL_FULL;
				}

			} else {
				/* we get here if h < t->bckts_max
				 * which means that we have memory available.
				 * We don't need to split and just write into
				 * the overflow area.
				 */

				t->bckts_tot = h+1;
			}
		}
	} while ( split_necessary );

	*GENTRY_P(t,h) = e;
	t->bckts_usd++;
	
	__LHT_UNLOCK(t);

	return 0;
}

int
lhtblAdd(LHTbl t, LHTblEntry e)
{
	return lhtblFindAdd(t, &e);
}

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
lhtblRpl(LHTbl t, LHTblEntry *pe, int add_fail)
{
LHTblEntry  e = *pe;
LHTblEntry *eh_p;
Hash        h;
int         split_necessary;
int         st;

	do {

		__LHT_LOCK(t);

		for ( h = dynhash(KEY(t, e),t->mod2l_msk,t->split_idx); *(eh_p = GENTRY_P(t, h)); h++ ) {
			if ( KEY(t, *eh_p) == KEY(t, e) ) {
				/* found; replace */
				*pe   = *eh_p;
				*eh_p = e;
				__LHT_UNLOCK(t);
				return 0;
			}
		}

		if ( add_fail ) {
			__LHT_UNLOCK(t);
			return LHTBL_KEY_NOTFND;
		}

		/* Load - controlled split */
		split_necessary = ( t->mngmt_aut && t->thres_hig > 0 && (t->bckts_usd + 1) * 100 > t->thres_hig * t->bckts_prm );

		if ( split_necessary ) {
			__LHT_UNLOCK(t);
			if ( (st = lhtblSplit(t)) )
				return st;
		} else {
			if ( h >= t->bckts_tot ) {
				assert( h == t->bckts_tot );

				/* would need to add the new entry past the table;
				 * see if we have space.
				 */

				/* if new bucket requires more memory than is available
				 * then split and try again
				 */
				if (  h >= t->bckts_max ) {
					split_necessary = t->mngmt_aut;

					__LHT_UNLOCK(t);

					if ( split_necessary )
						lhtblSplit(t);
					else
						return LHTBL_FULL;
				} else {
					/* we get here if h < t->bckts_max
					 * which means that we have memory available.
					 * We don't need to split and just write into
					 * the overflow area.
					 */

					t->bckts_tot = h+1;
				}
			}
		}
	} while ( split_necessary );

	*eh_p = e;
	*pe   = 0;

	t->bckts_usd++;

	__LHT_UNLOCK(t);

	return 0;
}

/* Update a linear chain of colliding entries
 * by moving eligible entries upstream first
 * filling bucket 'i' and then trying to fill
 * the new vacancy etc.
 */
static Hash
upd_chain(LHTbl t, Hash i)
{
Hash        j;
LHTblEntry  e_j;

	/* move something to entry i */
	for ( j = i + 1; (e_j = GENTRY(t, j)); j++ ) {

		Hash h     = dynhash(KEY(t, e_j),t->mod2l_msk,t->split_idx);

		/* when updating a chain during a split then

		Hash h_old = h & (mod2l_msk>>1);

		 * we have to be careful since we cannot blindly
		 * use the dynhash function for elements that are
		 * subject to the split (old hash == split_idx).
		 * We can have the following cases:
         *
         *   h == h_old
         *     ==> OK to move element
         *
         *   h != h_old; element already in slot >= h
         *     ==> MUST leave element in place
         *   
         *   h != h_old; element (still) in slot < h
         *     ==> OK to move element
         *
         * Hence we could adjust 'h:=h_old' in this
		 * special case:

		if ( h_old == split_idx && j < h )
			h = h_old;

         *
         *
         * However: we note that 
         *   a) for any key not currently being split
	     *      the dynhash function is correct.
         *
         *   b) position j >= dynhash(k) for valid table entries
         *
         *   c) for any key under a split which still
         *      hashes to the 'old' key the current
         *      position, 'j' is >= h_old.
         *
         * Hence 'j<h' implies 'key being split' AND 'h = hold + (1<<(table_lvl-1))'
         */

		/* we can only move elements if their hash is <= i
		 * otherwise, we break the chain into two pieces.
		 * If the chain is updated during a split operation
		 * then we have to be more careful. The possible
		 * scenarios are explained above...
		 */

		if ( h <= i || j < h ) {
			*GENTRY_P(t, i) = e_j;
			i = j;
		}
	}
	*GENTRY_P(t, i) = e_j;

	/* if we just moved the last bucket in the overflow area then
	 * the number 'total buckets' needs to be decremented.
	 */
	if ( t->bckts_tot > t->bckts_prm && 0 == GENTRY(t, t->bckts_tot-1) ) {
		t->bckts_tot--;
	}

	return i;
}

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
lhtblDel(LHTbl t, LHTblEntry e)
{
Hash       i,rval;
LHTblKey   key = KEY(t, e);
LHTblEntry e_i;

	__LHT_LOCK(t);

	for ( i = dynhash(key,t->mod2l_msk,t->split_idx); (e_i = GENTRY(t, i)); i++ ) {
		if ( KEY(t, e_i) == key ) {
			rval = upd_chain(t, i);

			/* Load - controlled grouping */
			t->bckts_usd--;

			i = t->mngmt_aut && t->thres_low > 0 && t->bckts_usd * 100 < t->thres_low * t->bckts_prm;

			__LHT_UNLOCK(t);

			if ( i )
				lhtblGroup(t);

			return 0;
		}
	}

	__LHT_UNLOCK(t);

	return LHTBL_KEY_NOTFND;
}

/* It is non-trivial to find out by how much
 * the table size needs to be increased to
 * accommodate all buckets that overflow the
 * end of the table as a result of the next
 * split.
 * However, we can easily compute an upper-
 * bound...
 */
static int
cntUpperBnd_unlocked(LHTbl t)
{
Hash       n = t->split_idx + 1;
Hash       l = t->table_lvl;
Hash       m = t->mod2l_msk;
Hash       i,h;
int        nup, sz;
LHTblEntry e_i;

	/* compute temporary table_lvl + mod 2^l bitmask */
	if ( n == (1<<(l-1))) {
		/* split would lead to a level-increase */
		l++;
		m = (m<<1) + 1;
		n = 0;
	}

	/* 'upper' position/hash of split bucket */
	h  = n + (1<<(l-1));

	/* Find number of elements that need to
	 * be moved up (using the 'new' table
     * parameters).
	 */
	for ( nup = 0, i = n; (e_i = GENTRY(t, i)) && i < h; i++ ) {
		if ( h == dynhash(KEY(t, e_i), m, n) )
			nup++;
	}

	/* new table size */
	sz = t->bckts_tot > t->bckts_prm ? t->bckts_tot : t->bckts_prm + 1;

	/* reduce by number of currently available primary buckets
	 * in table at & beyond 'upper/new' hash.
	 */
	for ( i = h; nup > 0 && i < sz; i++ ) {
		if ( ! GENTRY(t, i) ) {
			nup--;
		}
	}

	return sz + nup - t->bckts_tot;
}

/*
 * Compute an upper-bound for the number of
 * buckets that would have to be added to
 * the table if it were split.
 */
int
lhtblCntUpperBnd(LHTbl t)
{
int rval;

	__LHT_LOCK(t);
		rval = cntUpperBnd_unlocked(t);
	__LHT_UNLOCK(t);
	return rval;
}

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
lhtblSplit(LHTbl t)
{
Hash       i,j;
Hash       nsplit_idx;
int        holes    = 0;
int        append;
LHTblEntry e_i, *pe_j;

	__LHT_LOCK(t);

	/* Find out if we need to obtain more memory */
	while ( (i = cntUpperBnd_unlocked(t) + t->bckts_tot) > t->bckts_max ) {

#ifndef NDEBUG
		if ( debug ) {
			j = (t->split_idx + 1) & (t->mod2l_msk >> 1 );
			printf("Will need to increment table size by <= %u elements (splitting %u)\n",i - t->bckts_tot,j);
		}
#endif

		if ( expand_unlocked(t, i) ) {
			__LHT_UNLOCK(t);
			fprintf(stderr,"lhtbl: table expansion failed; cannot SPLIT\n");
			return LHTBL_FULL;
		}
		/* 'expand_unlock()' released the lock during allocation;
		 * therefore we must start over again...
		 */
	}

	/* new split index */
	nsplit_idx = t->split_idx + 1;

	/* Increment the number of primary buckets.
	 * If this results in the overflow area being
     * all covered by primary buckets then adjust
	 * the overflow-area index 'bckts_tot' as well.
	 */
	if ( ++t->bckts_prm > t->bckts_tot ) {
		assert( t->bckts_prm == t->bckts_tot + 1 );
		t->bckts_tot = t->bckts_prm;
	}

	/* If split index wraps around (would reach
	 * the 'upper-half' of the table) then the
	 * level must be incremented.
	 */
	t->split_idx = nsplit_idx;
	if ( t->split_idx == (1<<(t->table_lvl-1))) {
		/* increment level, recompute bitmask */
		t->table_lvl++;
		t->mod2l_msk = (t->mod2l_msk<<1) + 1;
		/* we must now split the first bucket */
		t->split_idx = nsplit_idx = 0;
	}

	/* Scan the collision-chain starting at the bucket to-be split */
	for ( i = nsplit_idx; (e_i = GENTRY(t, i)); i++) {
		LHTblKey  key   = KEY(t, e_i);
		Hash      h     = dynhash(key,t->mod2l_msk,nsplit_idx);
		Hash      h_old = h & (t->mod2l_msk>>1);
		Hash      hole;

		if ( h_old == nsplit_idx ) {
			/* This is a candidate. */ 

			if ( h != h_old ) {
				/* Yep. Using the 'new' hash function this element
				 * needs to be moved into the 'upper-half' of the
				 * table.
				 * However, we must do that only if it is not already
				 * there (because the collision-chain already overlaps
				 * the 'upper-half' of the table.
				 */
				assert( h == h_old + (1<<(t->table_lvl-1)) );
				/* move up ? */
				if ( i < h ) {
					/* yes; ('normal' case) element is in the 'lower-half'
					 * of the table but now hashes to the 'upper-half'.
					 */
					hole = upd_chain(t, i);

					/* must be moved; find a new slot in the upper-half */
					for ( append=h; GENTRY(t, append); append++)
						;

#ifndef NDEBUG
					if ( debug )
						printf("Moving %"PRIu32" from %u up to %u\n", key, i, append);
#endif

					/* we should have pre-allocated enough memory but do
					 * a paranoia check here...
					 */
					if ( append >= t->bckts_tot ) {
						assert( t->bckts_tot == append );
						if ( (t->bckts_tot = append + 1) > t->bckts_max ) {
							fprintf(stderr,"lhtbl: FATAL INTERNAL ERROR, incorrect pre-allocation for SPLIT\n");
							abort();
						}
					}

					*GENTRY_P(t, append) = e_i;

					/* count the number of buckets we free-ed up in the
					 * lower-half
					 */
					if ( hole < h )
						holes++;

					/* If we just moved an element downstream the
					 * collision-chain into bucket 'i' then we
					 * must check that bucket again.
					 * If 'upd_chain()' created a hole in bucket 'i'
					 * then we must continue checking 'i+1' (normal
					 * for-loop operation).
					 */
					if ( (e_i = GENTRY(t, i)) ) {
						/* check this slot again; decrement 'i' so that
						 * for-loop increment results finds the same
						 * bucket again.
						 */
						i--;
					}
				}
			} else {
				/* This element still hashes to the 'old' bucket. However,
				 * it could be that it currently is in the 'upper-half'
				 * of the table because it is member of a collision-chain
				 * that overlaps the 'upper-half'.
				 * Since there could now be free buckets in the 'lower-half'
				 * we must check if we can move this element down.
				 *
				 * We must do that if 
				 *  - element is at or beyond 'h'
				 *  - there is at least one new hole (as a result of moving
				 *    a predecessor of 'i' into the 'upper-half')
				 */
				if ( (i >= /*h_old + (1<<(table_lvl-1)) ==> */ h ) && holes ) {
					/* yes */
					upd_chain(t, i);

					/* find a new bucket for the element in the 'lower-half' */
					for ( j = h_old; *(pe_j = GENTRY_P(t,j)); j++ )
						;

					*pe_j = e_i;

#ifndef NDEBUG
					if ( debug )
						printf("Moving %"PRIu32" from %u down to %u\n", key, i, j);
#endif

					/* We ate up a hole */
					holes--;

					/* See comments above for why we check this slot again   */
					if ( (e_i = GENTRY(t, i)) ) {
						/* check this slot again */
						i--;
					}
				}
			}
		}
	}

	__LHT_UNLOCK(t);

	return 0;
}

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
lhtblGroup(LHTbl t)
{
Hash       i,j;
Hash       tgt;
int        rval;
LHTblEntry e_i, *pe_j;

	__LHT_LOCK(t);

	if( t->table_lvl < 1 ) {
		/* cannot reduce further */
		__LHT_UNLOCK(t);
		return 0;
	}

	if ( t->split_idx < 0 ) {
		t->table_lvl--;
		t->mod2l_msk >>= 1;
		t->split_idx = (1 << (t->table_lvl - 1)) - 1;
	}

	tgt = t->split_idx + (1<<(t->table_lvl-1));

	for ( i = tgt; (e_i = GENTRY(t, i)); i++ ) {
		LHTblKey  key   = KEY(t, e_i);
		Hash      h     = dynhash(key,t->mod2l_msk,t->split_idx);
		if ( h == tgt ) {
			*GENTRY_P(t, i) = 0;
			/* look for empty slot in lower section of table */
			for ( j = t->split_idx; *(pe_j = GENTRY_P(t, j)); j++ )
				;
#ifndef NDEBUG
			if ( debug )
				printf("Moving %"KEYu" from %u down to %u\n", key, i, j);
#endif
			*pe_j = e_i;
			if ( j == i ) {
				/* lower section overlaps upper section and
				 * the first available entry is the original
				 * position of 'e'.
				 * No further work to do...
				 */
				rval = 0;
				break;
			}
			upd_chain(t, i);
			if ( (e_i = GENTRY(t,i)) ) {
				i--;
			}
		}
	}

	if ( 0 == t->split_idx ) {
		t->table_lvl--;
		t->mod2l_msk >>= 1;
		t->split_idx = 1 << (t->table_lvl - 1);
	}
	t->split_idx--;
	t->bckts_prm--;

	/* count number of entries freed up in overflow area */
	for ( rval=0, i=t->bckts_tot-1; i >= t->bckts_prm && 0 == GENTRY(t,i); i-- )
		rval++;

	/* upd_chain should have found all these buckets */
	assert( rval <= 1 && t->bckts_tot >= t->bckts_prm );

	t->bckts_tot -= rval;
	i = t->bckts_tot;

	__LHT_UNLOCK(t);

	if ( rval ) {
		lhtblCollapse( t, i );
	}

	return 0;
}

void
lhtblGetStats(LHTbl t, LHTblStats stats)
{
	__LHT_LOCK(t);
		stats->tble_level    = t->table_lvl;
		stats->prim_buckets  = t->bckts_prm;
		stats->ovfl_buckets  = t->bckts_tot - t->bckts_prm;
		stats->avlb_buckets  = t->bckts_max;
		stats->used_buckets  = t->bckts_usd;
	__LHT_UNLOCK(t);
}

void
lhtblDumpStats(LHTbl t, FILE *f)
{
LHTblStatsRec stats;
	if ( ! f )
		f = stdout;
	lhtblGetStats(t, &stats);
	fprintf(f, "LH-Table Statistics:\n");
	fprintf(f, " Table Level:              %5i\n",   stats.tble_level);
	fprintf(f, " Primary Buckets:          %5i\n",   stats.prim_buckets);
	fprintf(f, " Overflow Buckets:         %5i\n",   stats.ovfl_buckets);
	fprintf(f, " Available Buckets (incl.\n");
    fprintf(f, "   buckets not yet attached\n");
    fprintf(f, "   to table):              %5i\n",   stats.avlb_buckets);
	fprintf(f, " Used/Occupied Buckets:    %5i\n",   stats.used_buckets);
	fprintf(f, " Table Load                  %5.1f%%\n", (float)stats.used_buckets/(float)stats.prim_buckets*100.0);
}

#ifdef TEST_MAIN

#ifdef TESTING
#define MAX 4
#else
#define MAX 32
#endif

#include <stdarg.h>

static LHTblConfigRec ex_cfg = {
	auto_manage_enable          :  1,
	load_percentage_thres_hi    :  0,
	load_percentage_thres_lo    :  0, /* never shrink */
	ovfl_buckets_preserve       :  0,
};

typedef struct TEntryRec_ {
	char     *str;
	LHTblKey k;	
	char     buf[];
} TEntryRec, *TEntry;

void prtble(LHTbl t, TEntry e)
{
	if ( !e )
		printf("<EMPTY>\n");
	else
		printf("0x%04x -> %s\n",KEY(t, e), e->str ? e->str : "<null>");
}

void prtbl(LHTbl t)
{
int i;
	for (i=0; i<t->bckts_tot; i++ ) {
		if ( 1<<(t->table_lvl-1) == i )
			printf("- - - - - -\n");
		printf("%c", i == t->split_idx ? '*' : ' ');
		prtble(t,GENTRY(t,i));
	}
	printf("==============\n");
}

static void
prnum(TEntry ents, int n)
{
int i;
	if ( n > 0 )
		printf("%u",ents[0].k);
	for ( i=1; i<n; i++ )
		printf(".%u",ents[i].k);
}

void clnte(LHTblEntry te)
{
	free(te);
}

void
purge(LHTbl t)
{
int i;
LHTblEntry *p_e;
	for ( i=0; i<t->bckts_max; i++ ) {
		p_e = GENTRY_P(t, i);
		free(*p_e);
		*p_e = 0;
	}
	tblrst(t, MAX);
}

int newperm(TEntry entries, int n)
{
int i=n;
	
	for ( i=0; i<n; i++ ) {
		if ( ((++entries[i].k) & 0xff) < n ) {
			return 0;
		}
		entries[i].k &= ~0xff;
	}
	return -1;
}

int
htck(LHTbl t, char *comment, TEntry ents, int n)
{
int  i;
Hash h;

	if ( t->bckts_tot > t->bckts_max ) {
		printf("**** HTAB MAX too small\n");
		abort();
	}
	for ( i=0; i<t->bckts_tot; i++ ) {
		LHTblEntry e;
		if ( (e = GENTRY(t, i)) ) {
			if ( (h = dynhash(KEY(t, e),t->mod2l_msk,t->split_idx)) > i ) {
				printf("***** HTAB INCONSISTENCY FOUND (%s) ******", comment);
				if ( n )
					prnum(ents,n);
				printf("\n");
				prtbl(t);
				printf("Slot %u (key 0x%04x) hashes to %u!\n", i, KEY(t, e), h);
				abort();
			}
		}
	}
	for ( i=t->bckts_prm; i<t->bckts_tot; i++ ) {
		if ( ! GENTRY(t, i) ) {
			printf("***** HTAB INCONSISTENCY FOUND (%s) *****", comment);
			if ( n )
				prnum(ents,n);
			printf("\n");
			prtbl(t);
			printf("Not all overflow elements used\n");
			abort();
		}
	}
	for ( i=t->bckts_tot; i<t->bckts_max; i++ ) {
		LHTblEntry e;
		if ( (e = GENTRY(t, i)) ) {
			printf("***** HTAB INCONSISTENCY FOUND (%s) *****", comment);
			if ( n )
				prnum(ents,n);
			printf("\n");
			prtbl(t);
			printf("Non-NULL element in overflow area found: (i=%u)\n", i);
			prtble(t,e);
			abort();
		}
	}
	return 0;
}

static void
doex(LHTbl t, void (*ex)(LHTbl t, TEntry, int), TEntry ents, int n)
{
int i;
	ex(t, ents, n);
	for ( i=0; i<t->bckts_max; i++ ) {
		assert( GENTRY(t, i) == 0);
	}
	tblrst(t, 1);
}

static void
chksplit(LHTbl t)
{
int	nbckts_tot = lhtblCntUpperBnd(t) + t->bckts_tot;
	assert( 0 == lhtblSplit(t));
	assert( nbckts_tot >= t->bckts_tot );
}

static void
catcher(TEntry ents, int n, ...) __attribute__((unused));

static void
catcher(TEntry ents, int n, ...)
{
va_list ap;
int     i;

	va_start(ap, n);

	for ( i=0; i<n; i++ ) {
		if ( ents[i].k != va_arg(ap, int) ) {
			va_end(ap);
			return;
		}
	}
	printf("FOUND\n");
}

/* forallentries(split-insert)
 * forallentries(del)
 */
static void
ex1(LHTbl t, TEntry ents, int n)
{
int i;


	for ( i=0; i<n; i++ ) {
		chksplit(t);
		htck(t, "ex1-split", ents, n);
		lhtblAdd(t, &ents[i]);
		htck(t, "ex1-add", ents, n);
	}

	for ( i=0; i<n; i++ ) {
		if ( lhtblDel(t, &ents[i]) ) {
			prtbl(t);
			printf("****EXERCISE 1 FAILURE*****");
			prnum(ents,n);
			printf("\n");
			abort();
		}
		htck(t, "ex1-del", ents, n);
	}
}

/* forallentries(insert),
 * forallentries(split)
 * forallentries(del)
 * forallentries(group)
 */
static void
ex2(LHTbl t, TEntry ents, int n)
{
int i;
int to = t->bckts_tot;

	for ( i=0; i<n; i++ ) {
		lhtblAdd(t, &ents[i]);
		htck(t, "ex2-add", ents, n);
		assert(t->bckts_tot >= t->bckts_prm);
	}

	for ( i=0; i<n; i++ ) {
		chksplit(t);
		htck(t, "ex2-split", ents, n);
		assert(t->bckts_tot >= t->bckts_prm);
	}

	for ( i=0; i<n; i++ ) {
		if ( lhtblDel(t, &ents[i]) ) {
			prtbl(t);
			printf("****EXERCISE 2 FAILURE*****");
			prnum(ents,n);
			printf("\n");
			abort();
		}
		htck(t, "ex2-del", ents, n);
	}
	for ( i=0; i<n; i++ ) {
		lhtblGroup(t);
		htck(t, "ex2-group", ents, n);
	}

	assert( to == t->bckts_tot );
}

/* forallentries(split-insert)
 * forallentries(del - group)
 */
static void
ex3(LHTbl t, TEntry ents, int n)
{
int i;
int to = t->bckts_tot;

	for ( i=0; i<n; i++ ) {
		chksplit(t);
		htck(t, "ex3-split", ents, n);
		lhtblAdd(t, &ents[i]);
		htck(t, "ex3-add", ents, n);
	}

	for ( i=0; i<n; i++ ) {
		if ( lhtblDel(t, &ents[i]) ) {
			prtbl(t);
			printf("****EXERCISE 3 FAILURE*****");
			prnum(ents,n);
			printf("\n");
			abort();
		}
		htck(t, "ex3-del", ents, n);
		lhtblGroup(t);
		htck(t, "ex3-group", ents, n);
	}
	assert( to == t->bckts_tot );
}

int
exercise(int n)
{
struct TEntryRec_ ents[n];
int    i,p;
LHTbl  t;

	t = lhtblCreate( n, (void*)&(((TEntry)0)->k) - (void*)0, &ex_cfg);
	assert( t );
	tblrst(t, 1);

	for ( i=0; i<n; i++ ) {
		ents[i].k   = i<<8;
		ents[i].str = 0;
	}
	
	p = 0;
	do { 
		doex(t, ex1,ents,n);
		doex(t, ex2,ents,n);
		/* Repeat ex1 but pre-allocate a few more empty buckets */
		for (i=t->bckts_tot; i<n/4; i++)
			chksplit(t);
		doex(t, ex1,ents,n);
		doex(t, ex3,ents,n);
		p++;
	} while ( 0 == newperm(ents,n) );

	lhtblDestroy(t, 0, 0);
	printf("Exercise Done (%u permutations)\n", p);
	return 0;
}

static LHTblConfigRec test_cfg = {
	auto_manage_enable          :  0,
	load_percentage_thres_hi    :  0,
	load_percentage_thres_lo    :  0, /* never shrink */
	ovfl_buckets_preserve       :  0,
};

int
main(int argc, char **argv)
{
char     buf[1024];
char     cmd;
int      st;
LHTblKey key;
TEntry   te;
int      i;
LHTbl    t = 0;

FILE *f = stdin;

t = lhtblCreate(MAX, (void*)&(((TEntry)0)->k) - (void*)0, &test_cfg);

printf("Debug mode: %s\n", debug ? "ON" : "OFF");
	do {
		switch ( (cmd = getc(f)) ) {
			case -1:
			break;

			case 'a': case 'A':
			case 'r': case 'R':
				if  ( fscanf(f, "%i%s",&key, buf) > 0 ) {
					te = calloc(1, sizeof(*te) + strlen(buf)+1);
					te->str = te->buf;
					strcpy(te->str,buf);
					te->k   = key;
					if ( 'r' != cmd && 'R' != cmd ) {
						if ( (st = lhtblAdd(t, (LHTblEntry)te)) ) {
							fprintf(stderr,"Adding failed: %i\n",st);
							free(te);
						}
					} 
				else {
						if ( (st = lhtblRpl(t, (LHTblEntry*)&te, 0)) ) {
							fprintf(stderr,"Replacing failed: %i\n",st);
						} else {
							printf("Replaced:\n");
							prtble(t, te);
						}
						free(te);
					}
				}
			break;

			case 'c': case 'C':
				htck(t, "usr", 0, 0);
			break;

			case 'd': case 'D':
				if  ( fscanf(f, "%i",&key) > 0 ) {
					if ( 0 == (te  = lhtblFind(t, key)) ) {
						fprintf(stderr,"Key %i not found\n", key);
					} else {
						lhtblDel(t, te);
						free(te);
					}
				}
			break;

			case 'e': case 'E':
				if ( fscanf(f, "%i", &i) > 0 ) {
					exercise(i);
				}
			break;

			case 'f': case 'F':
				if  ( fscanf(f, "%i",&key) > 0 ) {
					te = lhtblFind(t, key);
					if ( 0 == te )
						printf("<Not Found>\n");
					else
						prtble(t, te);
				}
			break;

			case 'g': case 'G':
				lhtblGroup(t);
			break;

			case 'p': case 'P': prtbl(t);
			break;

			case 'q': case 'Q':
			goto bail;

			/* 'r', 'R': see 'a'/'A' */

			case 's': case 'S':
				if ( (st=lhtblSplit(t)) )
					fprintf(stderr,"Split failed: %i\n",st);
			break;

#ifndef NDEBUG
			case 't': case 'T':
				debug = !debug;
			break;
#endif

			case 'x': case 'X':
				purge(t);
			break;

			case '#':
				printf("Table size would have to be incremented by <= %u elements (splitting %u)\n", lhtblCntUpperBnd(t), (t->split_idx + 1) & (t->mod2l_msk >> 1));
			break;

			case '%':
				lhtblDumpStats(t,stdout);
			break;
		}

		do { 
			cmd = getchar();
		} while ( -1 != cmd && '\n' != cmd );
	} while ( -1 != cmd );
bail:

	return 0;
}
#endif
