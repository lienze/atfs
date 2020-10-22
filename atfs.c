#include <linux/init.h>
#include <linux/module.h>
#include <linux/fs.h>
#include <linux/mount.h>
#include <linux/slab.h>
#include <linux/pagemap.h>
#include <linux/iversion.h>
#include <linux/buffer_head.h>
#include <linux/mpage.h>
#include <linux/quotaops.h>
#include <linux/cred.h>
#include <linux/blkdev.h>

#include "atfs.h"

static struct vfsmount *atfs_mnt;

const struct inode_operations atfs_dir_inode_operations;
const struct file_operations atfs_dir_operations;
const struct file_operations atfs_file_operations;
const struct address_space_operations atfs_aops;


struct atfs_dir_entry {
	__le32	inode;			/* Inode number */
	__le16	rec_len;		/* Directory entry length */
	__u8	name_len;		/* Name length */
	__u8	file_type;
	char	name[];			/* File name */
};

typedef struct {
	__le32	*p;
	__le32	key;
	struct buffer_head *bh;
} Indirect;

static inline void add_chain(Indirect *p, struct buffer_head *bh, __le32 *v)
{
	p->key = *(p->p = v);
	p->bh = bh;
}

static inline int verify_chain(Indirect *from, Indirect *to)
{
	while (from <= to && from->key == *from->p)
		from++;
	return (from > to);
}

typedef struct atfs_dir_entry atfs_dirent;

static inline unsigned atfs_chunk_size(struct inode *inode)
{
	return inode->i_sb->s_blocksize;
}

static inline void atfs_put_page(struct page *page)
{
	kunmap(page);
	put_page(page);
}

static struct page * atfs_get_page(struct inode *dir, unsigned long n,
				   int quiet)
{
	struct address_space *mapping = dir->i_mapping;
	struct page *page = read_mapping_page(mapping, n, NULL);
	if (!IS_ERR(page)) {
		kmap(page);
		if (unlikely(!PageChecked(page))) {
			if (PageError(page))
				goto fail;
		}
	}
	return page;

fail:
	atfs_put_page(page);
	return ERR_PTR(-EIO);
}

/*
 * ATFS_DIR_PAD defines the directory entries boundaries
 *
 * NOTE: It must be a multiple of 4
 */
#define ATFS_DIR_PAD		 	4
#define ATFS_DIR_ROUND 			(ATFS_DIR_PAD - 1)
#define ATFS_DIR_REC_LEN(name_len)	(((name_len) + 8 + ATFS_DIR_ROUND) & \
					 ~ATFS_DIR_ROUND)
#define ATFS_MAX_REC_LEN		((1<<16)-1)
static inline unsigned atfs_rec_len_from_disk(__le16 dlen)
{
	unsigned len = le16_to_cpu(dlen);

#if (PAGE_SIZE >= 65536)
	if (len == ATFS_MAX_REC_LEN)
		return 1 << 16;
#endif
	return len;
}

static int atfs_set_super(struct super_block *sb, void *data)
{
	sb->s_fs_info = data;
	return set_anon_super(sb, data);
}

static struct dentry *atfs_lookup(struct inode * dir,
		struct dentry *dentry, unsigned int flags)
{
	printk(KERN_INFO "==atfs== atfs_lookup invoked");
	return NULL;
}

ssize_t atfs_file_read(struct file *filp, char __user *buf,
		size_t count, loff_t *ppos)
{
	printk(KERN_INFO "==atfs== atfs_file_read invoked");
	return 0;
}

ssize_t atfs_file_write(struct file *filp, const char __user *buf,
		size_t count, loff_t *ppos)
{
	printk(KERN_INFO "==atfs== atfs_file_write invoked");
	return 0;
}

static int atfs_file_open(struct inode *inode, struct file *filp)
{
	//BUG();
	printk(KERN_INFO "==atfs== atfs_file_open invoked...");
	return 0;
}

static ssize_t atfs_file_read_iter(struct kiocb *iocb, struct iov_iter *to)
{
	return generic_file_read_iter(iocb, to);
}

static int atfs_mkdir(struct inode * dir, struct dentry * dentry, umode_t mode)
{
	struct inode * inode;

	inode = new_inode(dir->i_sb);
	inode->i_mode |= S_IFDIR;
	inode->i_op = &atfs_dir_inode_operations;
	inode->i_fop = &atfs_file_operations;
	inode->i_state = I_NEW;

	d_instantiate_new(dentry, inode);
	return 0;
}

/**
 *	atfs_block_to_path - parse the block number into array of offsets
 *	@inode: inode in question (we are only interested in its superblock)
 *	@i_block: block number to be parsed
 *	@offsets: array to store the offsets in
 *      @boundary: set this non-zero if the referred-to block is likely to be
 *             followed (on disk) by an indirect block.
 *	To store the locations of file's data atfs uses a data structure common
 *	for UNIX filesystems - tree of pointers anchored in the inode, with
 *	data blocks at leaves and indirect blocks in intermediate nodes.
 *	This function translates the block number into path in that tree -
 *	return value is the path length and @offsets[n] is the offset of
 *	pointer to (n+1)th node in the nth one. If @block is out of range
 *	(negative or too large) warning is printed and zero returned.
 *
 *	Note: function doesn't find node addresses, so no IO is needed. All
 *	we need to know is the capacity of indirect blocks (taken from the
 *	inode->i_sb).
 */

/*
 * Portability note: the last comparison (check that we fit into triple
 * indirect block) is spelled differently, because otherwise on an
 * architecture with 32-bit longs and 8Kb pages we might get into trouble
 * if our filesystem had 8Kb blocks. We might use long long, but that would
 * kill us on x86. Oh, well, at least the sign propagation does not matter -
 * i_block would have to be negative in the very beginning, so we would not
 * get there at all.
 */

static int atfs_block_to_path(struct inode *inode,
			long i_block, int offsets[4], int *boundary)
{
	int ptrs = ATFS_ADDR_PER_BLOCK(inode->i_sb);
	int ptrs_bits = ATFS_ADDR_PER_BLOCK_BITS(inode->i_sb);

	const long direct_blocks = ATFS_NDIR_BLOCKS,
		indirect_blocks = ptrs,
		double_blocks = (1 << (ptrs_bits * 2));
	int n = 0;
	int final = 0;

	printk(KERN_INFO "==atfs== atfs_block_to_path ptrs_bits:%d",ptrs_bits);

	if (i_block < 0) {
		/*
		 * ext2_msg(inode->i_sb, KERN_WARNING,
		 *     "warning: %s: block < 0", __func__);
		 */
	} else if (i_block < direct_blocks) {
		offsets[n++] = i_block;
		final = direct_blocks;
	} else if ( (i_block -= direct_blocks) < indirect_blocks) {
		offsets[n++] = ATFS_IND_BLOCK;
		offsets[n++] = i_block;
		final = ptrs;
	} else if ((i_block -= indirect_blocks) < double_blocks) {
		offsets[n++] = ATFS_DIND_BLOCK;
		offsets[n++] = i_block >> ptrs_bits;
		offsets[n++] = i_block & (ptrs - 1);
		final = ptrs;
	} else if (((i_block -= double_blocks) >> (ptrs_bits * 2)) < ptrs) {
		offsets[n++] = ATFS_TIND_BLOCK;
		offsets[n++] = i_block >> (ptrs_bits * 2);
		offsets[n++] = (i_block >> ptrs_bits) & (ptrs - 1);
		offsets[n++] = i_block & (ptrs - 1);
		final = ptrs;
	} else {
		/*
		 * ext2_msg(inode->i_sb, KERN_WARNING,
		 *     "warning: %s: block is too big", __func__);
		 */
	}
	if (boundary)
		*boundary = final - 1 - (i_block & (ptrs - 1));

	return n;
}

/**
 *	atfs_get_branch - read the chain of indirect blocks leading to data
 *	@inode: inode in question
 *	@depth: depth of the chain (1 - direct pointer, etc.)
 *	@offsets: offsets of pointers in inode/indirect blocks
 *	@chain: place to store the result
 *	@err: here we store the error value
 *
 *	Function fills the array of triples <key, p, bh> and returns %NULL
 *	if everything went OK or the pointer to the last filled triple
 *	(incomplete one) otherwise. Upon the return chain[i].key contains
 *	the number of (i+1)-th block in the chain (as it is stored in memory,
 *	i.e. little-endian 32-bit), chain[i].p contains the address of that
 *	number (it points into struct inode for i==0 and into the bh->b_data
 *	for i>0) and chain[i].bh points to the buffer_head of i-th indirect
 *	block for i>0 and NULL for i==0. In other words, it holds the block
 *	numbers of the chain, addresses they were taken from (and where we can
 *	verify that chain did not change) and buffer_heads hosting these
 *	numbers.
 *
 *	Function stops when it stumbles upon zero pointer (absent block)
 *		(pointer to last triple returned, *@err == 0)
 *	or when it gets an IO error reading an indirect block
 *		(ditto, *@err == -EIO)
 *	or when it notices that chain had been changed while it was reading
 *		(ditto, *@err == -EAGAIN)
 *	or when it reads all @depth-1 indirect blocks successfully and finds
 *	the whole chain, all way to the data (returns %NULL, *err == 0).
 */
static Indirect *atfs_get_branch(struct inode *inode,
				 int depth,
				 int *offsets,
				 Indirect chain[4],
				 int *err)
{
	struct super_block *sb = inode->i_sb;
	Indirect *p = chain;
	struct buffer_head *bh;

	*err = 0;
	/* i_data is not going away, no lock needed */
	add_chain (chain, NULL, ATFS_I(inode)->i_data + *offsets);
	if (!p->key)
		goto no_block;
	while (--depth) {
		bh = sb_bread(sb, le32_to_cpu(p->key));
		if (!bh)
			goto failure;
		read_lock(&ATFS_I(inode)->i_meta_lock);
		if (!verify_chain(chain, p))
			goto changed;
		add_chain(++p, bh, (__le32*)bh->b_data + *++offsets);
		read_unlock(&ATFS_I(inode)->i_meta_lock);
		if (!p->key)
			goto no_block;
	}
	return NULL;

changed:
	read_unlock(&ATFS_I(inode)->i_meta_lock);
	brelse(bh);
	*err = -EAGAIN;
	goto no_block;
failure:
	*err = -EIO;
no_block:
	return p;
}

/**
 * atfs_init_block_alloc_info()
 * @inode:		file inode structure
 *
 * Allocate and initialize the  reservation window structure, and
 * link the window to the atfs inode structure at last
 *
 * The reservation window structure is only dynamically allocated
 * and linked to atfs inode the first time the open file
 * needs a new block. So, before every atfs_new_block(s) call, for
 * regular files, we should check whether the reservation window
 * structure exists or not. In the latter case, this function is called.
 * Fail to do so will result in block reservation being turned off for that
 * open file.
 *
 * This function is called from atfs_get_blocks_handle(), also called
 * when setting the reservation window size through ioctl before the file
 * is open for write (needs block allocation).
 *
 * Needs truncate_mutex protection prior to calling this function.
 */
void atfs_init_block_alloc_info(struct inode *inode)
{
	struct atfs_inode_info *ei = ATFS_I(inode);
	struct atfs_block_alloc_info *block_i;
	struct super_block *sb = inode->i_sb;

	block_i = kmalloc(sizeof(*block_i), GFP_NOFS);
	if (block_i) {
		struct atfs_reserve_window_node *rsv = &block_i->rsv_window_node;

		rsv->rsv_start = ATFS_RESERVE_WINDOW_NOT_ALLOCATED;
		rsv->rsv_end = ATFS_RESERVE_WINDOW_NOT_ALLOCATED;

	 	/*
		 * if filesystem is mounted with NORESERVATION, the goal
		 * reservation window size is set to zero to indicate
		 * block reservation is off
		 */
		if (!test_opt(sb, RESERVATION))
			rsv->rsv_goal_size = 0;
		else
			rsv->rsv_goal_size = ATFS_DEFAULT_RESERVE_BLOCKS;
		rsv->rsv_alloc_hit = 0;
		block_i->last_alloc_logical_block = 0;
		block_i->last_alloc_physical_block = 0;
	}
	ei->i_block_alloc_info = block_i;
}

/**
 *	atfs_find_near - find a place for allocation with sufficient locality
 *	@inode: owner
 *	@ind: descriptor of indirect block.
 *
 *	This function returns the preferred place for block allocation.
 *	It is used when heuristic for sequential allocation fails.
 *	Rules are:
 *	  + if there is a block to the left of our position - allocate near it.
 *	  + if pointer will live in indirect block - allocate near that block.
 *	  + if pointer will live in inode - allocate in the same cylinder group.
 *
 * In the latter case we colour the starting block by the callers PID to
 * prevent it from clashing with concurrent allocations for a different inode
 * in the same block group.   The PID is used here so that functionally related
 * files will be close-by on-disk.
 *
 *	Caller must make sure that @ind is valid and will stay that way.
 */

static atfs_fsblk_t atfs_find_near(struct inode *inode, Indirect *ind)
{
	struct atfs_inode_info *ei = ATFS_I(inode);
	__le32 *start = ind->bh ? (__le32 *) ind->bh->b_data : ei->i_data;
	__le32 *p;
	atfs_fsblk_t bg_start;
	atfs_fsblk_t colour;

	/* Try to find previous block */
	for (p = ind->p - 1; p >= start; p--)
		if (*p)
			return le32_to_cpu(*p);

	/* No such thing, so let's try location of indirect block */
	if (ind->bh)
		return ind->bh->b_blocknr;

	/*
	 * It is going to be referred from inode itself? OK, just put it into
	 * the same cylinder group then.
	 */
	bg_start = atfs_group_first_block_no(inode->i_sb, ei->i_block_group);
	colour = (current->pid % 16) *
			(ATFS_BLOCKS_PER_GROUP(inode->i_sb) / 16);
	return bg_start + colour;
}


/**
 *	atfs_find_goal - find a preferred place for allocation.
 *	@inode: owner
 *	@block:  block we want
 *	@partial: pointer to the last triple within a chain
 *
 *	Returns preferred place for a block (the goal).
 */

static inline atfs_fsblk_t atfs_find_goal(struct inode *inode, long block,
					  Indirect *partial)
{
	struct atfs_block_alloc_info *block_i;

	block_i = ATFS_I(inode)->i_block_alloc_info;

	/*
	 * try the heuristic for sequential allocation,
	 * failing that at least try to get decent locality.
	 */
	if (block_i && (block == block_i->last_alloc_logical_block + 1)
		&& (block_i->last_alloc_physical_block != 0)) {
		return block_i->last_alloc_physical_block + 1;
	}

	return atfs_find_near(inode, partial);
}

/**
 *	atfs_blks_to_allocate: Look up the block map and count the number
 *	of direct blocks need to be allocated for the given branch.
 *
 * 	@branch: chain of indirect blocks
 *	@k: number of blocks need for indirect blocks
 *	@blks: number of data blocks to be mapped.
 *	@blocks_to_boundary:  the offset in the indirect block
 *
 *	return the total number of blocks to be allocate, including the
 *	direct and indirect blocks.
 */
static int
atfs_blks_to_allocate(Indirect * branch, int k, unsigned long blks,
		int blocks_to_boundary)
{
	unsigned long count = 0;

	/*
	 * Simple case, [t,d]Indirect block(s) has not allocated yet
	 * then it's clear blocks on that path have not allocated
	 */
	if (k > 0) {
		/* right now don't hanel cross boundary allocation */
		if (blks < blocks_to_boundary + 1)
			count += blks;
		else
			count += blocks_to_boundary + 1;
		return count;
	}

	count++;
	while (count < blks && count <= blocks_to_boundary
		&& le32_to_cpu(*(branch[0].p + count)) == 0) {
		count++;
	}
	return count;
}

/**
 * atfs_splice_branch - splice the allocated branch onto inode.
 * @inode: owner
 * @block: (logical) number of block we are adding
 * @where: location of missing link
 * @num:   number of indirect blocks we are adding
 * @blks:  number of direct blocks we are adding
 *
 * This function fills the missing link and does all housekeeping needed in
 * inode (->i_blocks, etc.). In case of success we end up with the full
 * chain to new block and return 0.
 */
static void atfs_splice_branch(struct inode *inode,
			long block, Indirect *where, int num, int blks)
{
	int i;
	struct atfs_block_alloc_info *block_i;
	atfs_fsblk_t current_block;

	block_i = ATFS_I(inode)->i_block_alloc_info;

	/* XXX LOCKING probably should have i_meta_lock ?*/
	/* That's it */

	*where->p = where->key;

	/*
	 * Update the host buffer_head or inode to point to more just allocated
	 * direct blocks blocks
	 */
	if (num == 0 && blks > 1) {
		current_block = le32_to_cpu(where->key) + 1;
		for (i = 1; i < blks; i++)
			*(where->p + i ) = cpu_to_le32(current_block++);
	}

	/*
	 * update the most recently allocated logical & physical block
	 * in i_block_alloc_info, to assist find the proper goal block for next
	 * allocation
	 */
	if (block_i) {
		block_i->last_alloc_logical_block = block + blks - 1;
		block_i->last_alloc_physical_block =
				le32_to_cpu(where[num].key) + blks - 1;
	}

	/* We are done with atomic stuff, now do the rest of housekeeping */

	/* had we spliced it onto indirect block? */
	if (where->bh)
		mark_buffer_dirty_inode(where->bh, inode);

	inode->i_ctime = current_time(inode);
	mark_inode_dirty(inode);
}

/*
 * Returns 1 if the passed-in block region is valid; 0 if some part overlaps
 * with filesystem metadata blocks.
 */
int atfs_data_block_valid(struct atfs_sb_info *sbi, atfs_fsblk_t start_blk,
			  unsigned int count)
{
	if ((start_blk <= le32_to_cpu(sbi->s_es->s_first_data_block)) ||
	    (start_blk + count - 1 < start_blk) ||
	    (start_blk + count - 1 >= le32_to_cpu(sbi->s_es->s_blocks_count)))
		return 0;

	/* Ensure we do not step over superblock */
	if ((start_blk <= sbi->s_sb_block) &&
	    (start_blk + count - 1 >= sbi->s_sb_block))
		return 0;

	return 1;
}

struct atfs_group_desc * atfs_get_group_desc(struct super_block * sb,
					     unsigned int block_group,
					     struct buffer_head ** bh)
{
	unsigned long group_desc;
	unsigned long offset;
	struct atfs_group_desc * desc;
	struct atfs_sb_info *sbi = ATFS_SB(sb);

	if (block_group >= sbi->s_groups_count) {
/*
 *         ext2_error (sb, "ext2_get_group_desc",
 *                 "block_group >= groups_count - "
 *                 "block_group = %d, groups_count = %lu",
 *                 block_group, sbi->s_groups_count);
 * 
 */
		return NULL;
	}

	group_desc = block_group >> ATFS_DESC_PER_BLOCK_BITS(sb);
	offset = block_group & (ATFS_DESC_PER_BLOCK(sb) - 1);
	if (!sbi->s_group_desc[group_desc]) {
		/*
		 * ext2_error (sb, "ext2_get_group_desc",
		 *         "Group descriptor not loaded - "
		 *         "block_group = %d, group_desc = %lu, desc = %lu",
		 *          block_group, group_desc, offset);
		 */
		return NULL;
	}

	desc = (struct atfs_group_desc *) sbi->s_group_desc[group_desc]->b_data;
	if (bh)
		*bh = sbi->s_group_desc[group_desc];
	return desc + offset;
}

static int atfs_valid_block_bitmap(struct super_block *sb,
					struct atfs_group_desc *desc,
					unsigned int block_group,
					struct buffer_head *bh)
{
	atfs_grpblk_t offset;
	atfs_grpblk_t next_zero_bit;
	atfs_fsblk_t bitmap_blk;
	atfs_fsblk_t group_first_block;

	group_first_block = atfs_group_first_block_no(sb, block_group);

	/* check whether block bitmap block number is set */
	bitmap_blk = le32_to_cpu(desc->bg_block_bitmap);
	offset = bitmap_blk - group_first_block;
	if (!atfs_test_bit(offset, bh->b_data))
		/* bad block bitmap */
		goto err_out;

	/* check whether the inode bitmap block number is set */
	bitmap_blk = le32_to_cpu(desc->bg_inode_bitmap);
	offset = bitmap_blk - group_first_block;
	if (!atfs_test_bit(offset, bh->b_data))
		/* bad block bitmap */
		goto err_out;

	/* check whether the inode table block number is set */
	bitmap_blk = le32_to_cpu(desc->bg_inode_table);
	offset = bitmap_blk - group_first_block;
	next_zero_bit = atfs_find_next_zero_bit(bh->b_data,
				offset + ATFS_SB(sb)->s_itb_per_group,
				offset);
	if (next_zero_bit >= offset + ATFS_SB(sb)->s_itb_per_group)
		/* good bitmap for inode tables */
		return 1;

err_out:
	/*
	 * ext2_error(sb, __func__,
	 *         "Invalid block bitmap - "
	 *         "block_group = %d, block = %lu",
	 *         block_group, bitmap_blk);
	 */
	return 0;
}

/*
 * Read the bitmap for a given block_group,and validate the
 * bits for block/inode/inode tables are set in the bitmaps
 *
 * Return buffer_head on success or NULL in case of failure.
 */
static struct buffer_head *
read_block_bitmap(struct super_block *sb, unsigned int block_group)
{
	struct atfs_group_desc * desc;
	struct buffer_head * bh = NULL;
	atfs_fsblk_t bitmap_blk;

	desc = atfs_get_group_desc(sb, block_group, NULL);
	if (!desc)
		return NULL;
	bitmap_blk = le32_to_cpu(desc->bg_block_bitmap);
	bh = sb_getblk(sb, bitmap_blk);
	if (unlikely(!bh)) {
		/*
		 * ext2_error(sb, __func__,
		 *         "Cannot read block bitmap - "
		 *         "block_group = %d, block_bitmap = %u",
		 *         block_group, le32_to_cpu(desc->bg_block_bitmap));
		 */
		return NULL;
	}
	if (likely(bh_uptodate_or_lock(bh)))
		return bh;

	if (bh_submit_read(bh) < 0) {
		brelse(bh);
		/*
		 * ext2_error(sb, __func__,
		 *         "Cannot read block bitmap - "
		 *         "block_group = %d, block_bitmap = %u",
		 *         block_group, le32_to_cpu(desc->bg_block_bitmap));
		 */
		return NULL;
	}

	atfs_valid_block_bitmap(sb, desc, block_group, bh);
	/*
	 * file system mounted not to panic on error, continue with corrupt
	 * bitmap
	 */
	return bh;
}

static void group_adjust_blocks(struct super_block *sb, int group_no,
	struct atfs_group_desc *desc, struct buffer_head *bh, int count)
{
	if (count) {
		struct atfs_sb_info *sbi = ATFS_SB(sb);
		unsigned free_blocks;

		spin_lock(sb_bgl_lock(sbi, group_no));
		free_blocks = le16_to_cpu(desc->bg_free_blocks_count);
		desc->bg_free_blocks_count = cpu_to_le16(free_blocks + count);
		spin_unlock(sb_bgl_lock(sbi, group_no));
		mark_buffer_dirty(bh);
	}
}


#define in_range(b, first, len)	((b) >= (first) && (b) <= (first) + (len) - 1)

/**
 * atfs_free_blocks() -- Free given blocks and update quota and i_blocks
 * @inode:		inode
 * @block:		start physical block to free
 * @count:		number of blocks to free
 */
void atfs_free_blocks (struct inode * inode, unsigned long block,
		       unsigned long count)
{
	struct buffer_head *bitmap_bh = NULL;
	struct buffer_head * bh2;
	unsigned long block_group;
	unsigned long bit;
	unsigned long i;
	unsigned long overflow;
	struct super_block * sb = inode->i_sb;
	struct atfs_sb_info * sbi = ATFS_SB(sb);
	struct atfs_group_desc * desc;
	struct atfs_super_block * es = sbi->s_es;
	unsigned freed = 0, group_freed;

	if (!atfs_data_block_valid(sbi, block, count)) {
		/*
		 * ext2_error (sb, "ext2_free_blocks",
		 *         "Freeing blocks not in datazone - "
		 *         "block = %lu, count = %lu", block, count);
		 */
		goto error_return;
	}

	/*
	 * ext2_debug ("freeing block(s) %lu-%lu\n", block, block + count - 1);
	 */

do_more:
	overflow = 0;
	block_group = (block - le32_to_cpu(es->s_first_data_block)) /
		      ATFS_BLOCKS_PER_GROUP(sb);
	bit = (block - le32_to_cpu(es->s_first_data_block)) %
		      ATFS_BLOCKS_PER_GROUP(sb);
	/*
	 * Check to see if we are freeing blocks across a group
	 * boundary.
	 */
	if (bit + count > ATFS_BLOCKS_PER_GROUP(sb)) {
		overflow = bit + count - ATFS_BLOCKS_PER_GROUP(sb);
		count -= overflow;
	}
	brelse(bitmap_bh);
	bitmap_bh = read_block_bitmap(sb, block_group);
	if (!bitmap_bh)
		goto error_return;

	desc = atfs_get_group_desc(sb, block_group, &bh2);
	if (!desc)
		goto error_return;

	if (in_range (le32_to_cpu(desc->bg_block_bitmap), block, count) ||
	    in_range (le32_to_cpu(desc->bg_inode_bitmap), block, count) ||
	    in_range (block, le32_to_cpu(desc->bg_inode_table),
		      sbi->s_itb_per_group) ||
	    in_range (block + count - 1, le32_to_cpu(desc->bg_inode_table),
		      sbi->s_itb_per_group)) {
		/*
		 * ext2_error (sb, "ext2_free_blocks",
		 *         "Freeing blocks in system zones - "
		 *         "Block = %lu, count = %lu",
		 *         block, count);
		 */
		goto error_return;
	}

	for (i = 0, group_freed = 0; i < count; i++) {
		if (!atfs_clear_bit_atomic(sb_bgl_lock(sbi, block_group),
						bit + i, bitmap_bh->b_data)) {
			/*
			 * ext2_error(sb, __func__,
			 *     "bit already cleared for block %lu", block + i);
			 */
		} else {
			group_freed++;
		}
	}

	mark_buffer_dirty(bitmap_bh);
	if (sb->s_flags & SB_SYNCHRONOUS)
		sync_dirty_buffer(bitmap_bh);

	group_adjust_blocks(sb, block_group, desc, bh2, group_freed);
	freed += group_freed;

	if (overflow) {
		block += count;
		count = overflow;
		goto do_more;
	}
error_return:
	brelse(bitmap_bh);
	if (freed) {
		percpu_counter_add(&sbi->s_freeblocks_counter, freed);
		dquot_free_block_nodirty(inode, freed);
		mark_inode_dirty(inode);
	}
}

/**
 * atfs_has_free_blocks()
 * @sbi:		in-core super block structure.
 *
 * Check if filesystem has at least 1 free block available for allocation.
 */
static int atfs_has_free_blocks(struct atfs_sb_info *sbi)
{
	atfs_fsblk_t free_blocks, root_blocks;

	free_blocks = percpu_counter_read_positive(&sbi->s_freeblocks_counter);
	root_blocks = le32_to_cpu(sbi->s_es->s_r_blocks_count);
	if (free_blocks < root_blocks + 1 && !capable(CAP_SYS_RESOURCE) &&
		!uid_eq(sbi->s_resuid, current_fsuid()) &&
		(gid_eq(sbi->s_resgid, GLOBAL_ROOT_GID) ||
		 !in_group_p (sbi->s_resgid))) {
		return 0;
	}
	return 1;
}

/*
 * rsv_is_empty() -- Check if the reservation window is allocated.
 * @rsv:		given reservation window to check
 *
 * returns 1 if the end block is ATFS_RESERVE_WINDOW_NOT_ALLOCATED.
 */
static inline int rsv_is_empty(struct atfs_reserve_window *rsv)
{
	/* a valid reservation end block could not be 0 */
	return (rsv->_rsv_end == ATFS_RESERVE_WINDOW_NOT_ALLOCATED);
}

/**
 * bitmap_search_next_usable_block()
 * @start:		the starting block (group relative) of the search
 * @bh:			bufferhead contains the block group bitmap
 * @maxblocks:		the ending block (group relative) of the reservation
 *
 * The bitmap search --- search forward through the actual bitmap on disk until
 * we find a bit free.
 */
static atfs_grpblk_t
bitmap_search_next_usable_block(atfs_grpblk_t start, struct buffer_head *bh,
					atfs_grpblk_t maxblocks)
{
	atfs_grpblk_t next;

	next = atfs_find_next_zero_bit(bh->b_data, maxblocks, start);
	if (next >= maxblocks)
		return -1;
	return next;
}


/**
 * find_next_usable_block()
 * @start:		the starting block (group relative) to find next
 * 			allocatable block in bitmap.
 * @bh:			bufferhead contains the block group bitmap
 * @maxblocks:		the ending block (group relative) for the search
 *
 * Find an allocatable block in a bitmap.  We perform the "most
 * appropriate allocation" algorithm of looking for a free block near
 * the initial goal; then for a free byte somewhere in the bitmap;
 * then for any free bit in the bitmap.
 */
static atfs_grpblk_t
find_next_usable_block(int start, struct buffer_head *bh, int maxblocks)
{
	atfs_grpblk_t here, next;
	char *p, *r;

	if (start > 0) {
		/*
		 * The goal was occupied; search forward for a free 
		 * block within the next XX blocks.
		 *
		 * end_goal is more or less random, but it has to be
		 * less than ATFS_BLOCKS_PER_GROUP. Aligning up to the
		 * next 64-bit boundary is simple..
		 */
		atfs_grpblk_t end_goal = (start + 63) & ~63;
		if (end_goal > maxblocks)
			end_goal = maxblocks;
		here = atfs_find_next_zero_bit(bh->b_data, end_goal, start);
		if (here < end_goal)
			return here;
		/*
		 * ext2_debug("Bit not found near goal\n");
		 */
	}

	here = start;
	if (here < 0)
		here = 0;

	p = ((char *)bh->b_data) + (here >> 3);
	r = memscan(p, 0, ((maxblocks + 7) >> 3) - (here >> 3));
	next = (r - ((char *)bh->b_data)) << 3;

	if (next < maxblocks && next >= here)
		return next;

	here = bitmap_search_next_usable_block(here, bh, maxblocks);
	return here;
}


/**
 * atfs_try_to_allocate()
 * @sb:			superblock
 * @group:		given allocation block group
 * @bitmap_bh:		bufferhead holds the block bitmap
 * @grp_goal:		given target block within the group
 * @count:		target number of blocks to allocate
 * @my_rsv:		reservation window
 *
 * Attempt to allocate blocks within a give range. Set the range of allocation
 * first, then find the first free bit(s) from the bitmap (within the range),
 * and at last, allocate the blocks by claiming the found free bit as allocated.
 *
 * To set the range of this allocation:
 * 	if there is a reservation window, only try to allocate block(s)
 * 	from the file's own reservation window;
 * 	Otherwise, the allocation range starts from the give goal block,
 * 	ends at the block group's last block.
 *
 * If we failed to allocate the desired block then we may end up crossing to a
 * new bitmap.
 */
static int
atfs_try_to_allocate(struct super_block *sb, int group,
			struct buffer_head *bitmap_bh, atfs_grpblk_t grp_goal,
			unsigned long *count,
			struct atfs_reserve_window *my_rsv)
{
	atfs_fsblk_t group_first_block;
	atfs_grpblk_t start, end;
	unsigned long num = 0;

	/* we do allocation within the reservation window if we have a window */
	if (my_rsv) {
		group_first_block = atfs_group_first_block_no(sb, group);
		if (my_rsv->_rsv_start >= group_first_block)
			start = my_rsv->_rsv_start - group_first_block;
		else
			/* reservation window cross group boundary */
			start = 0;
		end = my_rsv->_rsv_end - group_first_block + 1;
		if (end > ATFS_BLOCKS_PER_GROUP(sb))
			/* reservation window crosses group boundary */
			end = ATFS_BLOCKS_PER_GROUP(sb);
		if ((start <= grp_goal) && (grp_goal < end))
			start = grp_goal;
		else
			grp_goal = -1;
	} else {
		if (grp_goal > 0)
			start = grp_goal;
		else
			start = 0;
		end = ATFS_BLOCKS_PER_GROUP(sb);
	}

	BUG_ON(start > ATFS_BLOCKS_PER_GROUP(sb));

repeat:
	if (grp_goal < 0) {
		grp_goal = find_next_usable_block(start, bitmap_bh, end);
		if (grp_goal < 0)
			goto fail_access;
		if (!my_rsv) {
			int i;

			for (i = 0; i < 7 && grp_goal > start &&
					!atfs_test_bit(grp_goal - 1,
								bitmap_bh->b_data);
			     		i++, grp_goal--)
				;
		}
	}
	start = grp_goal;

	if (atfs_set_bit_atomic(sb_bgl_lock(ATFS_SB(sb), group), grp_goal,
			       				bitmap_bh->b_data)) {
		/*
		 * The block was allocated by another thread, or it was
		 * allocated and then freed by another thread
		 */
		start++;
		grp_goal++;
		if (start >= end)
			goto fail_access;
		goto repeat;
	}
	num++;
	grp_goal++;
	while (num < *count && grp_goal < end
		&& !atfs_set_bit_atomic(sb_bgl_lock(ATFS_SB(sb), group),
					grp_goal, bitmap_bh->b_data)) {
		num++;
		grp_goal++;
	}
	*count = num;
	return grp_goal - num;
fail_access:
	*count = num;
	return -1;
}

/**
 * __rsv_window_dump() -- Dump the filesystem block allocation reservation map
 * @rb_root:		root of per-filesystem reservation rb tree
 * @verbose:		verbose mode
 * @fn:			function which wishes to dump the reservation map
 *
 * If verbose is turned on, it will print the whole block reservation
 * windows(start, end). Otherwise, it will only print out the "bad" windows,
 * those windows that overlap with their immediate neighbors.
 */
#if 1
static void __rsv_window_dump(struct rb_root *root, int verbose,
			      const char *fn)
{
	struct rb_node *n;
	struct atfs_reserve_window_node *rsv, *prev;
	int bad;

restart:
	n = rb_first(root);
	bad = 0;
	prev = NULL;

	printk("Block Allocation Reservation Windows Map (%s):\n", fn);
	while (n) {
		rsv = rb_entry(n, struct atfs_reserve_window_node, rsv_node);
		if (verbose)
			printk("reservation window 0x%p "
				"start: %lu, end: %lu\n",
				rsv, rsv->rsv_start, rsv->rsv_end);
		if (rsv->rsv_start && rsv->rsv_start >= rsv->rsv_end) {
			printk("Bad reservation %p (start >= end)\n",
			       rsv);
			bad = 1;
		}
		if (prev && prev->rsv_end >= rsv->rsv_start) {
			printk("Bad reservation %p (prev->end >= start)\n",
			       rsv);
			bad = 1;
		}
		if (bad) {
			if (!verbose) {
				printk("Restarting reservation walk in verbose mode\n");
				verbose = 1;
				goto restart;
			}
		}
		n = rb_next(n);
		prev = rsv;
	}
	printk("Window map complete.\n");
	BUG_ON(bad);
}
#define rsv_window_dump(root, verbose) \
	__rsv_window_dump((root), (verbose), __func__)
#else
#define rsv_window_dump(root, verbose) do {} while (0)
#endif

/**
 * goal_in_my_reservation()
 * @rsv:		inode's reservation window
 * @grp_goal:		given goal block relative to the allocation block group
 * @group:		the current allocation block group
 * @sb:			filesystem super block
 *
 * Test if the given goal block (group relative) is within the file's
 * own block reservation window range.
 *
 * If the reservation window is outside the goal allocation group, return 0;
 * grp_goal (given goal block) could be -1, which means no specific
 * goal block. In this case, always return 1.
 * If the goal block is within the reservation window, return 1;
 * otherwise, return 0;
 */
static int
goal_in_my_reservation(struct atfs_reserve_window *rsv, atfs_grpblk_t grp_goal,
			unsigned int group, struct super_block * sb)
{
	atfs_fsblk_t group_first_block, group_last_block;

	group_first_block = atfs_group_first_block_no(sb, group);
	group_last_block = group_first_block + ATFS_BLOCKS_PER_GROUP(sb) - 1;

	if ((rsv->_rsv_start > group_last_block) ||
	    (rsv->_rsv_end < group_first_block))
		return 0;
	if ((grp_goal >= 0) && ((grp_goal + group_first_block < rsv->_rsv_start)
		|| (grp_goal + group_first_block > rsv->_rsv_end)))
		return 0;
	return 1;
}

/**
 * rsv_window_remove() -- unlink a window from the reservation rb tree
 * @sb:			super block
 * @rsv:		reservation window to remove
 *
 * Mark the block reservation window as not allocated, and unlink it
 * from the filesystem reservation window rb tree. Must be called with
 * rsv_lock held.
 */
static void rsv_window_remove(struct super_block *sb,
			      struct atfs_reserve_window_node *rsv)
{
	rsv->rsv_start = ATFS_RESERVE_WINDOW_NOT_ALLOCATED;
	rsv->rsv_end = ATFS_RESERVE_WINDOW_NOT_ALLOCATED;
	rsv->rsv_alloc_hit = 0;
	rb_erase(&rsv->rsv_node, &ATFS_SB(sb)->s_rsv_window_root);
}

/*
 * atfs_rsv_window_add() -- Insert a window to the block reservation rb tree.
 * @sb:			super block
 * @rsv:		reservation window to add
 *
 * Must be called with rsv_lock held.
 */
void atfs_rsv_window_add(struct super_block *sb,
		    struct atfs_reserve_window_node *rsv)
{
	struct rb_root *root = &ATFS_SB(sb)->s_rsv_window_root;
	struct rb_node *node = &rsv->rsv_node;
	atfs_fsblk_t start = rsv->rsv_start;

	struct rb_node ** p = &root->rb_node;
	struct rb_node * parent = NULL;
	struct atfs_reserve_window_node *this;

	while (*p)
	{
		parent = *p;
		this = rb_entry(parent, struct atfs_reserve_window_node, rsv_node);

		if (start < this->rsv_start)
			p = &(*p)->rb_left;
		else if (start > this->rsv_end)
			p = &(*p)->rb_right;
		else {
			rsv_window_dump(root, 1);
			BUG();
		}
	}

	rb_link_node(node, parent, p);
	rb_insert_color(node, root);
}


/**
 * 	find_next_reservable_window():
 *		find a reservable space within the given range.
 *		It does not allocate the reservation window for now:
 *		alloc_new_reservation() will do the work later.
 *
 * 	@search_head: the head of the searching list;
 *		This is not necessarily the list head of the whole filesystem
 *
 *		We have both head and start_block to assist the search
 *		for the reservable space. The list starts from head,
 *		but we will shift to the place where start_block is,
 *		then start from there, when looking for a reservable space.
 *
 * 	@size: the target new reservation window size
 *
 * 	@group_first_block: the first block we consider to start
 *			the real search from
 *
 * 	@last_block:
 *		the maximum block number that our goal reservable space
 *		could start from. This is normally the last block in this
 *		group. The search will end when we found the start of next
 *		possible reservable space is out of this boundary.
 *		This could handle the cross boundary reservation window
 *		request.
 *
 * 	basically we search from the given range, rather than the whole
 * 	reservation double linked list, (start_block, last_block)
 * 	to find a free region that is of my size and has not
 * 	been reserved.
 *
 */
static int find_next_reservable_window(
				struct atfs_reserve_window_node *search_head,
				struct atfs_reserve_window_node *my_rsv,
				struct super_block * sb,
				atfs_fsblk_t start_block,
				atfs_fsblk_t last_block)
{
	struct rb_node *next;
	struct atfs_reserve_window_node *rsv, *prev;
	atfs_fsblk_t cur;
	int size = my_rsv->rsv_goal_size;

	/* TODO: make the start of the reservation window byte-aligned */
	/* cur = *start_block & ~7;*/
	cur = start_block;
	rsv = search_head;
	if (!rsv)
		return -1;

	while (1) {
		if (cur <= rsv->rsv_end)
			cur = rsv->rsv_end + 1;

		/* TODO?
		 * in the case we could not find a reservable space
		 * that is what is expected, during the re-search, we could
		 * remember what's the largest reservable space we could have
		 * and return that one.
		 *
		 * For now it will fail if we could not find the reservable
		 * space with expected-size (or more)...
		 */
		if (cur > last_block)
			return -1;		/* fail */

		prev = rsv;
		next = rb_next(&rsv->rsv_node);
		rsv = rb_entry(next,struct atfs_reserve_window_node,rsv_node);

		/*
		 * Reached the last reservation, we can just append to the
		 * previous one.
		 */
		if (!next)
			break;

		if (cur + size <= rsv->rsv_start) {
			/*
			 * Found a reserveable space big enough.  We could
			 * have a reservation across the group boundary here
		 	 */
			break;
		}
	}
	/*
	 * we come here either :
	 * when we reach the end of the whole list,
	 * and there is empty reservable space after last entry in the list.
	 * append it to the end of the list.
	 *
	 * or we found one reservable space in the middle of the list,
	 * return the reservation window that we could append to.
	 * succeed.
	 */

	if ((prev != my_rsv) && (!rsv_is_empty(&my_rsv->rsv_window)))
		rsv_window_remove(sb, my_rsv);

	/*
	 * Let's book the whole available window for now.  We will check the
	 * disk bitmap later and then, if there are free blocks then we adjust
	 * the window size if it's larger than requested.
	 * Otherwise, we will remove this node from the tree next time
	 * call find_next_reservable_window.
	 */
	my_rsv->rsv_start = cur;
	my_rsv->rsv_end = cur + size - 1;
	my_rsv->rsv_alloc_hit = 0;

	if (prev != my_rsv)
		atfs_rsv_window_add(sb, my_rsv);
	return 0;
}

/**
 * search_reserve_window()
 * @rb_root:		root of reservation tree
 * @goal:		target allocation block
 *
 * Find the reserved window which includes the goal, or the previous one
 * if the goal is not in any window.
 * Returns NULL if there are no windows or if all windows start after the goal.
 */
static struct atfs_reserve_window_node *
search_reserve_window(struct rb_root *root, atfs_fsblk_t goal)
{
	struct rb_node *n = root->rb_node;
	struct atfs_reserve_window_node *rsv;

	if (!n)
		return NULL;

	do {
		rsv = rb_entry(n, struct atfs_reserve_window_node, rsv_node);

		if (goal < rsv->rsv_start)
			n = n->rb_left;
		else if (goal > rsv->rsv_end)
			n = n->rb_right;
		else
			return rsv;
	} while (n);
	/*
	 * We've fallen off the end of the tree: the goal wasn't inside
	 * any particular node.  OK, the previous node must be to one
	 * side of the interval containing the goal.  If it's the RHS,
	 * we need to back up one.
	 */
	if (rsv->rsv_start > goal) {
		n = rb_prev(&rsv->rsv_node);
		rsv = rb_entry(n, struct atfs_reserve_window_node, rsv_node);
	}
	return rsv;
}


/**
 * 	alloc_new_reservation()--allocate a new reservation window
 *
 *		To make a new reservation, we search part of the filesystem
 *		reservation list (the list that inside the group). We try to
 *		allocate a new reservation window near the allocation goal,
 *		or the beginning of the group, if there is no goal.
 *
 *		We first find a reservable space after the goal, then from
 *		there, we check the bitmap for the first free block after
 *		it. If there is no free block until the end of group, then the
 *		whole group is full, we failed. Otherwise, check if the free
 *		block is inside the expected reservable space, if so, we
 *		succeed.
 *		If the first free block is outside the reservable space, then
 *		start from the first free block, we search for next available
 *		space, and go on.
 *
 *	on succeed, a new reservation will be found and inserted into the list
 *	It contains at least one free block, and it does not overlap with other
 *	reservation windows.
 *
 *	failed: we failed to find a reservation window in this group
 *
 *	@rsv: the reservation
 *
 *	@grp_goal: The goal (group-relative).  It is where the search for a
 *		free reservable space should start from.
 *		if we have a goal(goal >0 ), then start from there,
 *		no goal(goal = -1), we start from the first block
 *		of the group.
 *
 *	@sb: the super block
 *	@group: the group we are trying to allocate in
 *	@bitmap_bh: the block group block bitmap
 *
 */
static int alloc_new_reservation(struct atfs_reserve_window_node *my_rsv,
		atfs_grpblk_t grp_goal, struct super_block *sb,
		unsigned int group, struct buffer_head *bitmap_bh)
{
	struct atfs_reserve_window_node *search_head;
	atfs_fsblk_t group_first_block, group_end_block, start_block;
	atfs_grpblk_t first_free_block;
	struct rb_root *fs_rsv_root = &ATFS_SB(sb)->s_rsv_window_root;
	unsigned long size;
	int ret;
	spinlock_t *rsv_lock = &ATFS_SB(sb)->s_rsv_window_lock;

	group_first_block = atfs_group_first_block_no(sb, group);
	group_end_block = group_first_block + (ATFS_BLOCKS_PER_GROUP(sb) - 1);

	if (grp_goal < 0)
		start_block = group_first_block;
	else
		start_block = grp_goal + group_first_block;

	size = my_rsv->rsv_goal_size;

	if (!rsv_is_empty(&my_rsv->rsv_window)) {
		/*
		 * if the old reservation is cross group boundary
		 * and if the goal is inside the old reservation window,
		 * we will come here when we just failed to allocate from
		 * the first part of the window. We still have another part
		 * that belongs to the next group. In this case, there is no
		 * point to discard our window and try to allocate a new one
		 * in this group(which will fail). we should
		 * keep the reservation window, just simply move on.
		 *
		 * Maybe we could shift the start block of the reservation
		 * window to the first block of next group.
		 */

		if ((my_rsv->rsv_start <= group_end_block) &&
				(my_rsv->rsv_end > group_end_block) &&
				(start_block >= my_rsv->rsv_start))
			return -1;

		if ((my_rsv->rsv_alloc_hit >
		     (my_rsv->rsv_end - my_rsv->rsv_start + 1) / 2)) {
			/*
			 * if the previously allocation hit ratio is
			 * greater than 1/2, then we double the size of
			 * the reservation window the next time,
			 * otherwise we keep the same size window
			 */
			size = size * 2;
			if (size > ATFS_MAX_RESERVE_BLOCKS)
				size = ATFS_MAX_RESERVE_BLOCKS;
			my_rsv->rsv_goal_size= size;
		}
	}

	spin_lock(rsv_lock);
	/*
	 * shift the search start to the window near the goal block
	 */
	search_head = search_reserve_window(fs_rsv_root, start_block);

	/*
	 * find_next_reservable_window() simply finds a reservable window
	 * inside the given range(start_block, group_end_block).
	 *
	 * To make sure the reservation window has a free bit inside it, we
	 * need to check the bitmap after we found a reservable window.
	 */
retry:
	ret = find_next_reservable_window(search_head, my_rsv, sb,
						start_block, group_end_block);

	if (ret == -1) {
		if (!rsv_is_empty(&my_rsv->rsv_window))
			rsv_window_remove(sb, my_rsv);
		spin_unlock(rsv_lock);
		return -1;
	}

	/*
	 * On success, find_next_reservable_window() returns the
	 * reservation window where there is a reservable space after it.
	 * Before we reserve this reservable space, we need
	 * to make sure there is at least a free block inside this region.
	 *
	 * Search the first free bit on the block bitmap.  Search starts from
	 * the start block of the reservable space we just found.
	 */
	spin_unlock(rsv_lock);
	first_free_block = bitmap_search_next_usable_block(
			my_rsv->rsv_start - group_first_block,
			bitmap_bh, group_end_block - group_first_block + 1);

	if (first_free_block < 0) {
		/*
		 * no free block left on the bitmap, no point
		 * to reserve the space. return failed.
		 */
		spin_lock(rsv_lock);
		if (!rsv_is_empty(&my_rsv->rsv_window))
			rsv_window_remove(sb, my_rsv);
		spin_unlock(rsv_lock);
		return -1;		/* failed */
	}

	start_block = first_free_block + group_first_block;
	/*
	 * check if the first free block is within the
	 * free space we just reserved
	 */
	if (start_block >= my_rsv->rsv_start && start_block <= my_rsv->rsv_end)
		return 0;		/* success */
	/*
	 * if the first free bit we found is out of the reservable space
	 * continue search for next reservable space,
	 * start from where the free block is,
	 * we also shift the list head to where we stopped last time
	 */
	search_head = my_rsv;
	spin_lock(rsv_lock);
	goto retry;
}

/**
 * try_to_extend_reservation()
 * @my_rsv:		given reservation window
 * @sb:			super block
 * @size:		the delta to extend
 *
 * Attempt to expand the reservation window large enough to have
 * required number of free blocks
 *
 * Since atfs_try_to_allocate() will always allocate blocks within
 * the reservation window range, if the window size is too small,
 * multiple blocks allocation has to stop at the end of the reservation
 * window. To make this more efficient, given the total number of
 * blocks needed and the current size of the window, we try to
 * expand the reservation window size if necessary on a best-effort
 * basis before atfs_new_blocks() tries to allocate blocks.
 */
static void try_to_extend_reservation(struct atfs_reserve_window_node *my_rsv,
			struct super_block *sb, int size)
{
	struct atfs_reserve_window_node *next_rsv;
	struct rb_node *next;
	spinlock_t *rsv_lock = &ATFS_SB(sb)->s_rsv_window_lock;

	if (!spin_trylock(rsv_lock))
		return;

	next = rb_next(&my_rsv->rsv_node);

	if (!next)
		my_rsv->rsv_end += size;
	else {
		next_rsv = rb_entry(next, struct atfs_reserve_window_node, rsv_node);

		if ((next_rsv->rsv_start - my_rsv->rsv_end - 1) >= size)
			my_rsv->rsv_end += size;
		else
			my_rsv->rsv_end = next_rsv->rsv_start - 1;
	}
	spin_unlock(rsv_lock);
}

/**
 * atfs_try_to_allocate_with_rsv()
 * @sb:			superblock
 * @group:		given allocation block group
 * @bitmap_bh:		bufferhead holds the block bitmap
 * @grp_goal:		given target block within the group
 * @count:		target number of blocks to allocate
 * @my_rsv:		reservation window
 *
 * This is the main function used to allocate a new block and its reservation
 * window.
 *
 * Each time when a new block allocation is need, first try to allocate from
 * its own reservation.  If it does not have a reservation window, instead of
 * looking for a free bit on bitmap first, then look up the reservation list to
 * see if it is inside somebody else's reservation window, we try to allocate a
 * reservation window for it starting from the goal first. Then do the block
 * allocation within the reservation window.
 *
 * This will avoid keeping on searching the reservation list again and
 * again when somebody is looking for a free block (without
 * reservation), and there are lots of free blocks, but they are all
 * being reserved.
 *
 * We use a red-black tree for the per-filesystem reservation list.
 */
static atfs_grpblk_t
atfs_try_to_allocate_with_rsv(struct super_block *sb, unsigned int group,
			struct buffer_head *bitmap_bh, atfs_grpblk_t grp_goal,
			struct atfs_reserve_window_node * my_rsv,
			unsigned long *count)
{
	atfs_fsblk_t group_first_block, group_last_block;
	atfs_grpblk_t ret = 0;
	unsigned long num = *count;

	/*
	 * we don't deal with reservation when
	 * filesystem is mounted without reservation
	 * or the file is not a regular file
	 * or last attempt to allocate a block with reservation turned on failed
	 */
	if (my_rsv == NULL) {
		return atfs_try_to_allocate(sb, group, bitmap_bh,
						grp_goal, count, NULL);
	}
	/*
	 * grp_goal is a group relative block number (if there is a goal)
	 * 0 <= grp_goal < ATFS_BLOCKS_PER_GROUP(sb)
	 * first block is a filesystem wide block number
	 * first block is the block number of the first block in this group
	 */
	group_first_block = atfs_group_first_block_no(sb, group);
	group_last_block = group_first_block + (ATFS_BLOCKS_PER_GROUP(sb) - 1);

	/*
	 * Basically we will allocate a new block from inode's reservation
	 * window.
	 *
	 * We need to allocate a new reservation window, if:
	 * a) inode does not have a reservation window; or
	 * b) last attempt to allocate a block from existing reservation
	 *    failed; or
	 * c) we come here with a goal and with a reservation window
	 *
	 * We do not need to allocate a new reservation window if we come here
	 * at the beginning with a goal and the goal is inside the window, or
	 * we don't have a goal but already have a reservation window.
	 * then we could go to allocate from the reservation window directly.
	 */
	while (1) {
		if (rsv_is_empty(&my_rsv->rsv_window) || (ret < 0) ||
			!goal_in_my_reservation(&my_rsv->rsv_window,
						grp_goal, group, sb)) {
			if (my_rsv->rsv_goal_size < *count)
				my_rsv->rsv_goal_size = *count;
			ret = alloc_new_reservation(my_rsv, grp_goal, sb,
							group, bitmap_bh);
			if (ret < 0)
				break;			/* failed */

			if (!goal_in_my_reservation(&my_rsv->rsv_window,
							grp_goal, group, sb))
				grp_goal = -1;
		} else if (grp_goal >= 0) {
			int curr = my_rsv->rsv_end -
					(grp_goal + group_first_block) + 1;

			if (curr < *count)
				try_to_extend_reservation(my_rsv, sb,
							*count - curr);
		}

		if ((my_rsv->rsv_start > group_last_block) ||
				(my_rsv->rsv_end < group_first_block)) {
			rsv_window_dump(&ATFS_SB(sb)->s_rsv_window_root, 1);
			BUG();
		}
		ret = atfs_try_to_allocate(sb, group, bitmap_bh, grp_goal,
					   &num, &my_rsv->rsv_window);
		if (ret >= 0) {
			my_rsv->rsv_alloc_hit += num;
			*count = num;
			break;				/* succeed */
		}
		num = *count;
	}
	return ret;
}


/*
 * atfs_new_blocks() -- core block(s) allocation function
 * @inode:		file inode
 * @goal:		given target block(filesystem wide)
 * @count:		target number of blocks to allocate
 * @errp:		error code
 *
 * atfs_new_blocks uses a goal block to assist allocation.  If the goal is
 * free, or there is a free block within 32 blocks of the goal, that block
 * is allocated.  Otherwise a forward search is made for a free block; within 
 * each block group the search first looks for an entire free byte in the block
 * bitmap, and then for any free bit if that fails.
 * This function also updates quota and i_blocks field.
 */
atfs_fsblk_t atfs_new_blocks(struct inode *inode, atfs_fsblk_t goal,
		    unsigned long *count, int *errp)
{
	struct buffer_head *bitmap_bh = NULL;
	struct buffer_head *gdp_bh;
	int group_no;
	int goal_group;
	atfs_grpblk_t grp_target_blk;	/* blockgroup relative goal block */
	atfs_grpblk_t grp_alloc_blk;	/* blockgroup-relative allocated block*/
	atfs_fsblk_t ret_block;		/* filesyetem-wide allocated block */
	int bgi;			/* blockgroup iteration index */
	int performed_allocation = 0;
	atfs_grpblk_t free_blocks;	/* number of free blocks in a group */
	struct super_block *sb;
	struct atfs_group_desc *gdp;
	struct atfs_super_block *es;
	struct atfs_sb_info *sbi;
	struct atfs_reserve_window_node *my_rsv = NULL;
	struct atfs_block_alloc_info *block_i;
	unsigned short windowsz = 0;
	unsigned long ngroups;
	unsigned long num = *count;
	int ret;

	*errp = -ENOSPC;
	sb = inode->i_sb;

	/*
	 * Check quota for allocation of this block.
	 */
	ret = dquot_alloc_block(inode, num);
	if (ret) {
		*errp = ret;
		return 0;
	}

	sbi = ATFS_SB(sb);
	es = ATFS_SB(sb)->s_es;
	/*
	 * ext2_debug("goal=%lu.\n", goal);
	 */
	/*
	 * Allocate a block from reservation only when
	 * filesystem is mounted with reservation(default,-o reservation), and
	 * it's a regular file, and
	 * the desired window size is greater than 0 (One could use ioctl
	 * command ATFS_IOC_SETRSVSZ to set the window size to 0 to turn off
	 * reservation on that particular file)
	 */
	block_i = ATFS_I(inode)->i_block_alloc_info;
	if (block_i) {
		windowsz = block_i->rsv_window_node.rsv_goal_size;
		if (windowsz > 0)
			my_rsv = &block_i->rsv_window_node;
	}

	if (!atfs_has_free_blocks(sbi)) {
		*errp = -ENOSPC;
		goto out;
	}

	/*
	 * First, test whether the goal block is free.
	 */
	if (goal < le32_to_cpu(es->s_first_data_block) ||
	    goal >= le32_to_cpu(es->s_blocks_count))
		goal = le32_to_cpu(es->s_first_data_block);
	group_no = (goal - le32_to_cpu(es->s_first_data_block)) /
			ATFS_BLOCKS_PER_GROUP(sb);
	goal_group = group_no;
retry_alloc:
	gdp = atfs_get_group_desc(sb, group_no, &gdp_bh);
	if (!gdp)
		goto io_error;

	free_blocks = le16_to_cpu(gdp->bg_free_blocks_count);
	/*
	 * if there is not enough free blocks to make a new resevation
	 * turn off reservation for this allocation
	 */
	if (my_rsv && (free_blocks < windowsz)
		&& (free_blocks > 0)
		&& (rsv_is_empty(&my_rsv->rsv_window)))
		my_rsv = NULL;

	if (free_blocks > 0) {
		grp_target_blk = ((goal - le32_to_cpu(es->s_first_data_block)) %
				ATFS_BLOCKS_PER_GROUP(sb));
		bitmap_bh = read_block_bitmap(sb, group_no);
		if (!bitmap_bh)
			goto io_error;
		grp_alloc_blk = atfs_try_to_allocate_with_rsv(sb, group_no,
					bitmap_bh, grp_target_blk,
					my_rsv, &num);
		if (grp_alloc_blk >= 0)
			goto allocated;
	}

	ngroups = ATFS_SB(sb)->s_groups_count;
	smp_rmb();

	/*
	 * Now search the rest of the groups.  We assume that
	 * group_no and gdp correctly point to the last group visited.
	 */
	for (bgi = 0; bgi < ngroups; bgi++) {
		group_no++;
		if (group_no >= ngroups)
			group_no = 0;
		gdp = atfs_get_group_desc(sb, group_no, &gdp_bh);
		if (!gdp)
			goto io_error;

		free_blocks = le16_to_cpu(gdp->bg_free_blocks_count);
		/*
		 * skip this group (and avoid loading bitmap) if there
		 * are no free blocks
		 */
		if (!free_blocks)
			continue;
		/*
		 * skip this group if the number of
		 * free blocks is less than half of the reservation
		 * window size.
		 */
		if (my_rsv && (free_blocks <= (windowsz/2)))
			continue;

		brelse(bitmap_bh);
		bitmap_bh = read_block_bitmap(sb, group_no);
		if (!bitmap_bh)
			goto io_error;
		/*
		 * try to allocate block(s) from this group, without a goal(-1).
		 */
		grp_alloc_blk = atfs_try_to_allocate_with_rsv(sb, group_no,
					bitmap_bh, -1, my_rsv, &num);
		if (grp_alloc_blk >= 0)
			goto allocated;
	}
	/*
	 * We may end up a bogus earlier ENOSPC error due to
	 * filesystem is "full" of reservations, but
	 * there maybe indeed free blocks available on disk
	 * In this case, we just forget about the reservations
	 * just do block allocation as without reservations.
	 */
	if (my_rsv) {
		my_rsv = NULL;
		windowsz = 0;
		group_no = goal_group;
		goto retry_alloc;
	}
	/* No space left on the device */
	*errp = -ENOSPC;
	goto out;

allocated:

	/*
	 * ext2_debug("using block group %d(%d)\n",
	 *         group_no, gdp->bg_free_blocks_count);
	 */

	ret_block = grp_alloc_blk + atfs_group_first_block_no(sb, group_no);

	if (in_range(le32_to_cpu(gdp->bg_block_bitmap), ret_block, num) ||
	    in_range(le32_to_cpu(gdp->bg_inode_bitmap), ret_block, num) ||
	    in_range(ret_block, le32_to_cpu(gdp->bg_inode_table),
		      ATFS_SB(sb)->s_itb_per_group) ||
	    in_range(ret_block + num - 1, le32_to_cpu(gdp->bg_inode_table),
		      ATFS_SB(sb)->s_itb_per_group)) {
		/*
		 * ext2_error(sb, "ext2_new_blocks",
		 *         "Allocating block in system zone - "
		 *         "blocks from "E2FSBLK", length %lu",
		 *         ret_block, num);
		 */

		/*
		 * ext2_try_to_allocate marked the blocks we allocated as in
		 * use.  So we may want to selectively mark some of the blocks
		 * as free
		 */
		goto retry_alloc;
	}

	performed_allocation = 1;

	if (ret_block + num - 1 >= le32_to_cpu(es->s_blocks_count)) {
		/*
		 * ext2_error(sb, "ext2_new_blocks",
		 *         "block("E2FSBLK") >= blocks count(%d) - "
		 *         "block_group = %d, es == %p ", ret_block,
		 *     le32_to_cpu(es->s_blocks_count), group_no, es);
		 */
		goto out;
	}

	group_adjust_blocks(sb, group_no, gdp, gdp_bh, -num);
	percpu_counter_sub(&sbi->s_freeblocks_counter, num);

	mark_buffer_dirty(bitmap_bh);
	if (sb->s_flags & SB_SYNCHRONOUS)
		sync_dirty_buffer(bitmap_bh);

	*errp = 0;
	brelse(bitmap_bh);
	if (num < *count) {
		dquot_free_block_nodirty(inode, *count-num);
		mark_inode_dirty(inode);
		*count = num;
	}
	return ret_block;

io_error:
	*errp = -EIO;
out:
	/*
	 * Undo the block allocation
	 */
	if (!performed_allocation) {
		dquot_free_block_nodirty(inode, *count);
		mark_inode_dirty(inode);
	}
	brelse(bitmap_bh);
	return 0;
}

atfs_fsblk_t atfs_new_block(struct inode *inode, unsigned long goal, int *errp)
{
	unsigned long count = 1;

	return atfs_new_blocks(inode, goal, &count, errp);
}


/**
 *	atfs_alloc_blocks: multiple allocate blocks needed for a branch
 *	@indirect_blks: the number of blocks need to allocate for indirect
 *			blocks
 *
 *	@new_blocks: on return it will store the new block numbers for
 *	the indirect blocks(if needed) and the first direct block,
 *	@blks:	on return it will store the total number of allocated
 *		direct blocks
 */
static int atfs_alloc_blocks(struct inode *inode,
			atfs_fsblk_t goal, int indirect_blks, int blks,
			atfs_fsblk_t new_blocks[4], int *err)
{
	int target, i;
	unsigned long count = 0;
	int index = 0;
	atfs_fsblk_t current_block = 0;
	int ret = 0;

	/*
	 * Here we try to allocate the requested multiple blocks at once,
	 * on a best-effort basis.
	 * To build a branch, we should allocate blocks for
	 * the indirect blocks(if not allocated yet), and at least
	 * the first direct block of this branch.  That's the
	 * minimum number of blocks need to allocate(required)
	 */
	target = blks + indirect_blks;

	while (1) {
		count = target;
		/* allocating blocks for indirect blocks and direct blocks */
		current_block = atfs_new_blocks(inode,goal,&count,err);
		if (*err)
			goto failed_out;

		target -= count;
		/* allocate blocks for indirect blocks */
		while (index < indirect_blks && count) {
			new_blocks[index++] = current_block++;
			count--;
		}

		if (count > 0)
			break;
	}

	/* save the new block number for the first direct block */
	new_blocks[index] = current_block;

	/* total number of blocks allocated for direct blocks */
	ret = count;
	*err = 0;
	return ret;
failed_out:
	for (i = 0; i <index; i++)
		atfs_free_blocks(inode, new_blocks[i], 1);
	if (index)
		mark_inode_dirty(inode);
	return ret;
}

/**
 *	atfs_alloc_branch - allocate and set up a chain of blocks.
 *	@inode: owner
 *	@indirect_blks: depth of the chain (number of blocks to allocate)
 *	@blks: number of allocated direct blocks
 *	@goal: preferred place for allocation
 *	@offsets: offsets (in the blocks) to store the pointers to next.
 *	@branch: place to store the chain in.
 *
 *	This function allocates @num blocks, zeroes out all but the last one,
 *	links them into chain and (if we are synchronous) writes them to disk.
 *	In other words, it prepares a branch that can be spliced onto the
 *	inode. It stores the information about that chain in the branch[], in
 *	the same format as atfs_get_branch() would do. We are calling it after
 *	we had read the existing part of chain and partial points to the last
 *	triple of that (one with zero ->key). Upon the exit we have the same
 *	picture as after the successful atfs_get_block(), except that in one
 *	place chain is disconnected - *branch->p is still zero (we did not
 *	set the last link), but branch->key contains the number that should
 *	be placed into *branch->p to fill that gap.
 *
 *	If allocation fails we free all blocks we've allocated (and forget
 *	their buffer_heads) and return the error value the from failed
 *	atfs_alloc_block() (normally -ENOSPC). Otherwise we set the chain
 *	as described above and return 0.
 */

static int atfs_alloc_branch(struct inode *inode,
			int indirect_blks, int *blks, atfs_fsblk_t goal,
			int *offsets, Indirect *branch)
{
	int blocksize = inode->i_sb->s_blocksize;
	int i, n = 0;
	int err = 0;
	struct buffer_head *bh;
	int num;
	atfs_fsblk_t new_blocks[4];
	atfs_fsblk_t current_block;

	num = atfs_alloc_blocks(inode, goal, indirect_blks,
				*blks, new_blocks, &err);
	if (err)
		return err;

	branch[0].key = cpu_to_le32(new_blocks[0]);
	/*
	 * metadata blocks and data blocks are allocated.
	 */
	for (n = 1; n <= indirect_blks;  n++) {
		/*
		 * Get buffer_head for parent block, zero it out
		 * and set the pointer to new one, then send
		 * parent to disk.
		 */
		bh = sb_getblk(inode->i_sb, new_blocks[n-1]);
		if (unlikely(!bh)) {
			err = -ENOMEM;
			goto failed;
		}
		branch[n].bh = bh;
		lock_buffer(bh);
		memset(bh->b_data, 0, blocksize);
		branch[n].p = (__le32 *) bh->b_data + offsets[n];
		branch[n].key = cpu_to_le32(new_blocks[n]);
		*branch[n].p = branch[n].key;
		if ( n == indirect_blks) {
			current_block = new_blocks[n];
			/*
			 * End of chain, update the last new metablock of
			 * the chain to point to the new allocated
			 * data blocks numbers
			 */
			for (i=1; i < num; i++)
				*(branch[n].p + i) = cpu_to_le32(++current_block);
		}
		set_buffer_uptodate(bh);
		unlock_buffer(bh);
		mark_buffer_dirty_inode(bh, inode);
		/* We used to sync bh here if IS_SYNC(inode).
		 * But we now rely upon generic_write_sync()
		 * and b_inode_buffers.  But not for directories.
		 */
		if (S_ISDIR(inode->i_mode) && IS_DIRSYNC(inode))
			sync_dirty_buffer(bh);
	}
	*blks = num;
	return err;

failed:
	for (i = 1; i < n; i++)
		bforget(branch[i].bh);
	for (i = 0; i < indirect_blks; i++)
		atfs_free_blocks(inode, new_blocks[i], 1);
	atfs_free_blocks(inode, new_blocks[i], num);
	return err;
}


/*
 * Allocation strategy is simple: if we have to allocate something, we will
 * have to go the whole way to leaf. So let's do it before attaching anything
 * to tree, set linkage between the newborn blocks, write them if sync is
 * required, recheck the path, free and repeat if check fails, otherwise
 * set the last missing link (that will protect us from any truncate-generated
 * removals - all blocks on the path are immune now) and possibly force the
 * write on the parent block.
 * That has a nice additional property: no special recovery from the failed
 * allocations is needed - we simply release blocks and do not touch anything
 * reachable from inode.
 *
 * `handle' can be NULL if create == 0.
 *
 * return > 0, # of blocks mapped or allocated.
 * return = 0, if plain lookup failed.
 * return < 0, error case.
 */
static int atfs_get_blocks(struct inode *inode,
			   sector_t iblock, unsigned long maxblocks,
			   u32 *bno, bool *new, bool *boundary,
			   int create)
{
	int err;
	int offsets[4];
	Indirect chain[4];
	Indirect *partial;
	atfs_fsblk_t goal;
	int indirect_blks;
	int blocks_to_boundary = 0;
	int depth;
	struct atfs_inode_info *ei = ATFS_I(inode);
	int count = 0;
	atfs_fsblk_t first_block = 0;

	BUG_ON(maxblocks == 0);

	depth = atfs_block_to_path(inode,iblock,offsets,&blocks_to_boundary);

	if (depth == 0)
		return -EIO;

	partial = atfs_get_branch(inode, depth, offsets, chain, &err);
	/* Simplest case - block found, no allocation needed */
	if (!partial) {
		first_block = le32_to_cpu(chain[depth - 1].key);
		count++;
		/*map more blocks*/
		while (count < maxblocks && count <= blocks_to_boundary) {
			atfs_fsblk_t blk;

			if (!verify_chain(chain, chain + depth - 1)) {
				/*
				 * Indirect block might be removed by
				 * truncate while we were reading it.
				 * Handling of that case: forget what we've
				 * got now, go to reread.
				 */
				err = -EAGAIN;
				count = 0;
				partial = chain + depth - 1;
				break;
			}
			blk = le32_to_cpu(*(chain[depth-1].p + count));
			if (blk == first_block + count)
				count++;
			else
				break;
		}
		if (err != -EAGAIN)
			goto got_it;
	}

	/* Next simple case - plain lookup or failed read of indirect block */
	if (!create || err == -EIO)
		goto cleanup;

	mutex_lock(&ei->truncate_mutex);
	/*
	 * If the indirect block is missing while we are reading
	 * the chain(atfs_get_branch() returns -EAGAIN err), or
	 * if the chain has been changed after we grab the semaphore,
	 * (either because another process truncated this branch, or
	 * another get_block allocated this branch) re-grab the chain to see if
	 * the request block has been allocated or not.
	 *
	 * Since we already block the truncate/other get_block
	 * at this point, we will have the current copy of the chain when we
	 * splice the branch into the tree.
	 */
	if (err == -EAGAIN || !verify_chain(chain, partial)) {
		while (partial > chain) {
			brelse(partial->bh);
			partial--;
		}
		partial = atfs_get_branch(inode, depth, offsets, chain, &err);
		if (!partial) {
			count++;
			mutex_unlock(&ei->truncate_mutex);
			goto got_it;
		}

		if (err) {
			mutex_unlock(&ei->truncate_mutex);
			goto cleanup;
		}
	}

	/*
	 * Okay, we need to do block allocation.  Lazily initialize the block
	 * allocation info here if necessary
	*/
	if (S_ISREG(inode->i_mode) && (!ei->i_block_alloc_info))
		atfs_init_block_alloc_info(inode);

	goal = atfs_find_goal(inode, iblock, partial);

	/* the number of blocks need to allocate for [d,t]indirect blocks */
	indirect_blks = (chain + depth) - partial - 1;
	/*
	 * Next look up the indirect map to count the total number of
	 * direct blocks to allocate for this branch.
	 */
	count = atfs_blks_to_allocate(partial, indirect_blks,
					maxblocks, blocks_to_boundary);
	/*
	 * XXX ???? Block out atfs_truncate while we alter the tree
	 */
	err = atfs_alloc_branch(inode, indirect_blks, &count, goal,
				offsets + (partial - chain), partial);

	if (err) {
		mutex_unlock(&ei->truncate_mutex);
		goto cleanup;
	}

	if (IS_DAX(inode)) {
		/*
		 * We must unmap blocks before zeroing so that writeback cannot
		 * overwrite zeros with stale data from block device page cache.
		 */
		clean_bdev_aliases(inode->i_sb->s_bdev,
				   le32_to_cpu(chain[depth-1].key),
				   count);
		/*
		 * block must be initialised before we put it in the tree
		 * so that it's not found by another thread before it's
		 * initialised
		 */
		err = sb_issue_zeroout(inode->i_sb,
				le32_to_cpu(chain[depth-1].key), count,
				GFP_NOFS);
		if (err) {
			mutex_unlock(&ei->truncate_mutex);
			goto cleanup;
		}
	}
	*new = true;

	atfs_splice_branch(inode, iblock, partial, indirect_blks, count);
	mutex_unlock(&ei->truncate_mutex);
got_it:
	if (count > blocks_to_boundary)
		*boundary = true;
	err = count;
	/* Clean up and exit */
	partial = chain + depth - 1;	/* the whole chain */
cleanup:
	while (partial > chain) {
		brelse(partial->bh);
		partial--;
	}
	if (err > 0)
		*bno = le32_to_cpu(chain[depth-1].key);
	return err;
}

int atfs_get_block(struct inode *inode, sector_t iblock,
		struct buffer_head *bh_result, int create)
{
	unsigned max_blocks = bh_result->b_size >> inode->i_blkbits;
	bool new = false, boundary = false;
	u32 bno;
	int ret;
	ret = atfs_get_blocks(inode, iblock, max_blocks, &bno, &new, &boundary,create);
	if (ret <= 0)
		return ret;

	map_bh(bh_result, inode->i_sb, bno);
	bh_result->b_size = (ret << inode->i_blkbits);
	if (new)
		set_buffer_new(bh_result);
	if (boundary)
		set_buffer_boundary(bh_result);
	return 0;
}

static int atfs_readpage(struct file *file, struct page *page)
{
	return mpage_readpage(page, atfs_get_block);
}

const struct address_space_operations atfs_aops = {
	.readpage		= atfs_readpage,
};

/*
 * p is at least 6 bytes before the end of page
 */
static inline atfs_dirent *atfs_next_entry(atfs_dirent *p)
{
	return (atfs_dirent *)((char *)p +
			atfs_rec_len_from_disk(p->rec_len));
}

static inline unsigned 
atfs_validate_entry(char *base, unsigned offset, unsigned mask)
{
	atfs_dirent *de = (atfs_dirent*)(base + offset);
	atfs_dirent *p = (atfs_dirent*)(base + (offset&mask));
	while ((char*)p < (char*)de) {
		if (p->rec_len == 0)
			break;
		p = atfs_next_entry(p);
	}
	return (char *)p - base;
}

/*
 * Return the offset into page `page_nr' of the last valid
 * byte in that page, plus one.
 */
static unsigned
atfs_last_byte(struct inode *inode, unsigned long page_nr)
{
	unsigned last_byte = inode->i_size;

	last_byte -= page_nr << PAGE_SHIFT;
	if (last_byte > PAGE_SIZE)
		last_byte = PAGE_SIZE;
	return last_byte;
}

static int atfs_readdir(struct file *file, struct dir_context *ctx)
{
	loff_t pos = ctx->pos;
	struct inode *inode = file_inode(file);
	//struct super_block *sb = inode->i_sb;
	unsigned int offset = pos & ~PAGE_MASK;
	unsigned long n = pos >> PAGE_SHIFT;
	unsigned long npages = dir_pages(inode);
	unsigned chunk_mask = ~(atfs_chunk_size(inode)-1);
	bool need_revalidate = !inode_eq_iversion(inode, file->f_version);
	bool has_filetype;

	/*
	 * printk(KERN_INFO "==atfs== atfs_readdir pos:%lld npages:%d inode_size:%lld",
	 *         pos, npages, inode->i_size);
	 */

	//if (pos > inode->i_size - ATFS_DIR_REC_LEN(1))
	//	return 0;

	//has_filetype =
	//	ATFS_HAS_INCOMPAT_FEATURE(sb, ATFS_FEATURE_INCOMPAT_FILETYPE);

	for ( ; n < npages; n++, offset = 0) {
		char *kaddr, *limit;
		atfs_dirent *de;
		struct page *page = atfs_get_page(inode, n, 0);

		if (IS_ERR(page)) {
			/*
			 * ext2_error(sb, __func__,
			 *        "bad page in #%lu",
			 *        inode->i_ino);
			 */
			ctx->pos += PAGE_SIZE - offset;
			return PTR_ERR(page);
		}
		kaddr = page_address(page);
		if (unlikely(need_revalidate)) {
			if (offset) {
				offset = atfs_validate_entry(kaddr, offset, chunk_mask);
				ctx->pos = (n<<PAGE_SHIFT) + offset;
			}
			file->f_version = inode_query_iversion(inode);
			need_revalidate = false;
		}
		de = (atfs_dirent *)(kaddr+offset);
		limit = kaddr + atfs_last_byte(inode, n) - ATFS_DIR_REC_LEN(1);
		for ( ;(char*)de <= limit; de = atfs_next_entry(de)) {
			if (de->rec_len == 0) {
				/*
				 * ext2_error(sb, __func__,
				 *     "zero-length directory entry");
				 */
				atfs_put_page(page);
				return -EIO;
			}
			if (de->inode) {
				unsigned char d_type = DT_UNKNOWN;

				if (has_filetype)
					d_type = fs_ftype_to_dtype(de->file_type);

				if (!dir_emit(ctx, de->name, de->name_len,
						le32_to_cpu(de->inode),
						d_type)) {
					atfs_put_page(page);
					return 0;
				}
			}
			ctx->pos += atfs_rec_len_from_disk(de->rec_len);
		}
		atfs_put_page(page);
	}
	return 0;
}

const struct inode_operations atfs_dir_inode_operations = {
	.lookup	= atfs_lookup,
	.mkdir	= atfs_mkdir,
};

const struct file_operations atfs_dir_operations = {
	.iterate_shared	= atfs_readdir,
};

const struct file_operations atfs_file_operations = {
	.open		= atfs_file_open,
	.write		= atfs_file_write,
	.read		= atfs_file_read,
	.read_iter	= atfs_file_read_iter,
};

static struct dentry *atfs_mount(struct file_system_type *fs_type, int flags,
		   const char *dev_name, void *data)
{
	struct block_device *bdev;
	struct super_block *sb = NULL;
	struct inode *inode;
	struct dentry *root;
	struct atfs_sb_info * sb_info;
	fmode_t mode = FMODE_READ | FMODE_EXCL;

	if (atfs_mnt)
		return ERR_PTR(-EINVAL);
	printk(KERN_INFO "==atfs== atfs_mount invoked...");

	bdev = blkdev_get_by_path(dev_name, mode, fs_type);
	if (IS_ERR(bdev))
		return ERR_CAST(bdev);

	sb_info = kzalloc(sizeof(*sb_info), GFP_KERNEL);
	if (!sb_info)
		goto failed;

	sb = sget(fs_type, NULL, atfs_set_super, flags, sb_info);
	if (!sb) {
		printk(KERN_INFO "==atfs== super_block is NULL");
		return ERR_PTR(-EINVAL);
	}
	sb->s_bdev = bdev;
	sb_set_blocksize(sb, block_size(bdev));

	sb_info->s_addr_per_block_bits = 
		ilog2(ATFS_ADDR_PER_BLOCK(sb));

	inode = new_inode(sb);
	inode->i_mode |= S_IFDIR;
	inode->i_op = &atfs_dir_inode_operations;
	inode->i_fop = &atfs_dir_operations;
	inode->i_size = 1024;
	inode->i_mapping->a_ops = &atfs_aops;
	root = d_make_root(inode);
	sb->s_root = root;
	return dget(sb->s_root);
failed:
	return NULL;
}

static struct file_system_type atfs_fs_type = {
	.name	= "atfs",
	.owner	= THIS_MODULE,
	.mount	= atfs_mount,
};

struct dentry *atfs_create_file(const char* name)
{
	struct dentry *root_dentry = NULL;
	struct dentry *new_dentry = NULL;
	struct inode *inode = NULL;
	struct qstr dentry_name;
	char *tmp_name;
	tmp_name = kmalloc(NAME_MAX + 1, GFP_KERNEL);
	if (!tmp_name)
		return ERR_PTR(-ENOMEM);
	memcpy(tmp_name, name, 6);
	dentry_name.name = tmp_name;
	dentry_name.len = 6;

	printk(KERN_INFO "==atfs== atfs_create_file start...");
	root_dentry = atfs_mnt->mnt_sb->s_root;
	new_dentry = d_alloc(root_dentry, &dentry_name);
	inode = new_inode(root_dentry->d_inode->i_sb);
	inode->i_mode = S_IFDIR;
	if (inode) {
		d_instantiate(new_dentry, inode);
		printk(KERN_INFO "==atfs== root_dentry half...");
		dget(new_dentry);
		printk(KERN_INFO "==atfs== root_dentry dget...");
	}
	printk(KERN_INFO "==atfs== atfs_create_file end.");
	return ERR_PTR(-EINVAL);
}

struct dentry *atfs_create_dir(const char *name)
{
	return atfs_create_file(name);
}

/*
 * register atfs
 */
int atfs_register(void)
{
	int ret;
	printk(KERN_INFO "==atfs== atfs_register start...");
	ret = register_filesystem(&atfs_fs_type);
	if (ret < 0) {
		printk(KERN_ERR "==atfs== register file system fail %d", ret);
		return ret;
	} else {
		printk(KERN_INFO "==atfs== register file system succeed %d", ret);
	}
	printk(KERN_INFO "==atfs== atfs_register end...");
	return 0;
}

int atfs_unregister(void)
{
	int ret;
	printk(KERN_INFO "==atfs== atfs_unregister start...");
	ret = unregister_filesystem(&atfs_fs_type);
	if (ret < 0) {
		printk(KERN_ERR "==atfs== unregister file system fail %d", ret);
		return ret;
	}
	return 0;
}

static int __init atfs_init(void)
{
	printk(KERN_INFO "==atfs== atfs_init start...\n");
	atfs_register();
	return 0;
}
module_init(atfs_init);

static void __exit atfs_exit(void)
{
	printk(KERN_INFO "==atfs== atfs_exit start...\n");
	atfs_unregister();
	return;
}
module_exit(atfs_exit);

MODULE_AUTHOR("Enze Li");
MODULE_LICENSE("GPL v2");
MODULE_DESCRIPTION("A Tiny File System");
MODULE_ALIAS("atfs");
