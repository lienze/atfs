
#include <linux/blockgroup_lock.h>

/* data type for block offset of block group */
typedef int atfs_grpblk_t;

/* data type for filesystem-wide blocks number */
typedef unsigned long atfs_fsblk_t;


struct atfs_reserve_window {
	atfs_fsblk_t		_rsv_start;	/* First byte reserved */
	atfs_fsblk_t		_rsv_end;	/* Last byte reserved or 0 */
};

struct atfs_reserve_window_node {
	struct rb_node	 	rsv_node;
	__u32			rsv_goal_size;
	__u32			rsv_alloc_hit;
	struct atfs_reserve_window	rsv_window;
};

#define rsv_start rsv_window._rsv_start
#define rsv_end rsv_window._rsv_end

struct atfs_block_alloc_info {
	/* information about reservation window */
	struct atfs_reserve_window_node	rsv_window_node;
	/*
	 * was i_next_alloc_block in atfs_inode_info
	 * is the logical (file-relative) number of the
	 * most-recently-allocated block in this file.
	 * We use this for detecting linearly ascending allocation requests.
	 */
	__u32			last_alloc_logical_block;
	/*
	 * Was i_next_alloc_goal in atfs_inode_info
	 * is the *physical* companion to i_next_alloc_block.
	 * it the the physical block number of the block which was most-recentl
	 * allocated to this file.  This give us the goal (target) for the next
	 * allocation when we detect linearly ascending requests.
	 */
	atfs_fsblk_t		last_alloc_physical_block;
};

/*
 * Structure of the super block
 */
struct atfs_super_block {
	__le32	s_inodes_count;		/* Inodes count */
	__le32	s_blocks_count;		/* Blocks count */
	__le32	s_r_blocks_count;	/* Reserved blocks count */
	__le32	s_free_blocks_count;	/* Free blocks count */
	__le32	s_free_inodes_count;	/* Free inodes count */
	__le32	s_first_data_block;	/* First Data Block */
	__le32	s_log_block_size;	/* Block size */
	__le32	s_log_frag_size;	/* Fragment size */
	__le32	s_blocks_per_group;	/* # Blocks per group */
	__le32	s_frags_per_group;	/* # Fragments per group */
	__le32	s_inodes_per_group;	/* # Inodes per group */
	__le32	s_mtime;		/* Mount time */
	__le32	s_wtime;		/* Write time */
	__le16	s_mnt_count;		/* Mount count */
	__le16	s_max_mnt_count;	/* Maximal mount count */
	__le16	s_magic;		/* Magic signature */
	__le16	s_state;		/* File system state */
	__le16	s_errors;		/* Behaviour when detecting errors */
	__le16	s_minor_rev_level; 	/* minor revision level */
	__le32	s_lastcheck;		/* time of last check */
	__le32	s_checkinterval;	/* max. time between checks */
	__le32	s_creator_os;		/* OS */
	__le32	s_rev_level;		/* Revision level */
	__le16	s_def_resuid;		/* Default uid for reserved blocks */
	__le16	s_def_resgid;		/* Default gid for reserved blocks */
	/*
	 * These fields are for ATFS_DYNAMIC_REV superblocks only.
	 *
	 * Note: the difference between the compatible feature set and
	 * the incompatible feature set is that if there is a bit set
	 * in the incompatible feature set that the kernel doesn't
	 * know about, it should refuse to mount the filesystem.
	 * 
	 * e2fsck's requirements are more strict; if it doesn't know
	 * about a feature in either the compatible or incompatible
	 * feature set, it must abort and not try to meddle with
	 * things it doesn't understand...
	 */
	__le32	s_first_ino; 		/* First non-reserved inode */
	__le16	s_inode_size; 		/* size of inode structure */
	__le16	s_block_group_nr; 	/* block group # of this superblock */
	__le32	s_feature_compat; 	/* compatible feature set */
	__le32	s_feature_incompat; 	/* incompatible feature set */
	__le32	s_feature_ro_compat; 	/* readonly-compatible feature set */
	__u8	s_uuid[16];		/* 128-bit uuid for volume */
	char	s_volume_name[16]; 	/* volume name */
	char	s_last_mounted[64]; 	/* directory where last mounted */
	__le32	s_algorithm_usage_bitmap; /* For compression */
	/*
	 * Performance hints.  Directory preallocation should only
	 * happen if the ATFS_COMPAT_PREALLOC flag is on.
	 */
	__u8	s_prealloc_blocks;	/* Nr of blocks to try to preallocate*/
	__u8	s_prealloc_dir_blocks;	/* Nr to preallocate for dirs */
	__u16	s_padding1;
	/*
	 * Journaling support valid if EXT3_FEATURE_COMPAT_HAS_JOURNAL set.
	 */
	__u8	s_journal_uuid[16];	/* uuid of journal superblock */
	__u32	s_journal_inum;		/* inode number of journal file */
	__u32	s_journal_dev;		/* device number of journal file */
	__u32	s_last_orphan;		/* start of list of inodes to delete */
	__u32	s_hash_seed[4];		/* HTREE hash seed */
	__u8	s_def_hash_version;	/* Default hash version to use */
	__u8	s_reserved_char_pad;
	__u16	s_reserved_word_pad;
	__le32	s_default_mount_opts;
 	__le32	s_first_meta_bg; 	/* First metablock block group */
	__u32	s_reserved[190];	/* Padding to the end of the block */
};


/*
 * atfs file system inode data in memory
 */
struct atfs_inode_info {
	__le32	i_data[15];
	__u32	i_flags;
	__u32	i_faddr;
	__u8	i_frag_no;
	__u8	i_frag_size;
	__u16	i_state;
	__u32	i_file_acl;
	__u32	i_dir_acl;
	__u32	i_dtime;

	/*
	 * i_block_group is the number of the block group which contains
	 * this file's inode.  Constant across the lifetime of the inode,
	 * it is used for making block allocation decisions - we try to
	 * place a file's data blocks near its inode block, and new inodes
	 * near to their parent directory's inode.
	 */
	__u32	i_block_group;

	/* block reservation info */
	struct atfs_block_alloc_info *i_block_alloc_info;

	__u32	i_dir_start_lookup;
	rwlock_t i_meta_lock;
	/*
	 * truncate_mutex is for serialising atfs_truncate() against
	 * atfs_getblock().  It also protects the internals of the inode's
	 * reservation data structures: atfs_reserve_window and
	 * atfs_reserve_window_node.
	 */
	struct mutex truncate_mutex;
	struct inode	vfs_inode;
	struct list_head i_orphan;	/* unlinked but open inodes */
};

/*
 * atfs file system super-block data in memory
 */
struct atfs_sb_info {
	unsigned long s_frag_size;	/* Size of a fragment in bytes */
	unsigned long s_frags_per_block;/* Number of fragments per block */
	unsigned long s_inodes_per_block;/* Number of inodes per block */
	unsigned long s_frags_per_group;/* Number of fragments in a group */
	unsigned long s_blocks_per_group;/* Number of blocks in a group */
	unsigned long s_inodes_per_group;/* Number of inodes in a group */
	unsigned long s_itb_per_group;	/* Number of inode table blocks per group */
	unsigned long s_gdb_count;	/* Number of group descriptor blocks */
	unsigned long s_desc_per_block;	/* Number of group descriptors per block */
	unsigned long s_groups_count;	/* Number of groups in the fs */
	unsigned long s_overhead_last;  /* Last calculated overhead */
	unsigned long s_blocks_last;    /* Last seen block count */
	struct buffer_head * s_sbh;	/* Buffer containing the super block */
	struct atfs_super_block * s_es;	/* Pointer to the super block in the buffer */
	struct buffer_head ** s_group_desc;
	unsigned long  s_mount_opt;
	unsigned long s_sb_block;
	kuid_t s_resuid;
	kgid_t s_resgid;
	unsigned short s_mount_state;
	unsigned short s_pad;
	int s_addr_per_block_bits;
	int s_desc_per_block_bits;
	int s_inode_size;
	int s_first_ino;
	spinlock_t s_next_gen_lock;
	u32 s_next_generation;
	unsigned long s_dir_count;
	u8 *s_debts;
	struct percpu_counter s_freeblocks_counter;
	struct percpu_counter s_freeinodes_counter;
	struct percpu_counter s_dirs_counter;
	struct blockgroup_lock *s_blockgroup_lock;
	/* root of the per fs reservation window tree */
	spinlock_t s_rsv_window_lock;
	struct rb_root s_rsv_window_root;
	struct atfs_reserve_window_node s_rsv_window_head;
	/*
	 * s_lock protects against concurrent modifications of s_mount_state,
	 * s_blocks_last, s_overhead_last and the content of superblock's
	 * buffer pointed to by sbi->s_es.
	 *
	 * Note: It is used in atfs_show_options() to provide a consistent view
	 * of the mount options.
	 */
	spinlock_t s_lock;
	struct mb_cache *s_ea_block_cache;
	struct dax_device *s_daxdev;
};

static inline struct atfs_sb_info *ATFS_SB(struct super_block *sb)
{
	return sb->s_fs_info;
}

/*
 * Structure of a blocks group descriptor
 */
struct atfs_group_desc
{
	__le32	bg_block_bitmap;		/* Blocks bitmap block */
	__le32	bg_inode_bitmap;		/* Inodes bitmap block */
	__le32	bg_inode_table;		/* Inodes table block */
	__le16	bg_free_blocks_count;	/* Free blocks count */
	__le16	bg_free_inodes_count;	/* Free inodes count */
	__le16	bg_used_dirs_count;	/* Directories count */
	__le16	bg_pad;
	__le32	bg_reserved[3];
};


/*
 * Macro-instructions used to manage group descriptors
 */
#define ATFS_BLOCKS_PER_GROUP(s)	(ATFS_SB(s)->s_blocks_per_group)
#define ATFS_DESC_PER_BLOCK(s)		(ATFS_SB(s)->s_desc_per_block)
#define ATFS_INODES_PER_GROUP(s)	(ATFS_SB(s)->s_inodes_per_group)
#define ATFS_DESC_PER_BLOCK_BITS(s)	(ATFS_SB(s)->s_desc_per_block_bits)

/*
 * Macro-instructions used to manage several block sizes
 */
#define ATFS_MIN_BLOCK_SIZE		1024
#define	ATFS_MAX_BLOCK_SIZE		4096
#define ATFS_MIN_BLOCK_LOG_SIZE		  10
#define ATFS_BLOCK_SIZE(s)		((s)->s_blocksize)
#define	ATFS_ADDR_PER_BLOCK(s)		(ATFS_BLOCK_SIZE(s) / sizeof (__u32))
#define ATFS_BLOCK_SIZE_BITS(s)		((s)->s_blocksize_bits)
#define	ATFS_ADDR_PER_BLOCK_BITS(s)	(ATFS_SB(s)->s_addr_per_block_bits)
#define ATFS_INODE_SIZE(s)		(ATFS_SB(s)->s_inode_size)
#define ATFS_FIRST_INO(s)		(ATFS_SB(s)->s_first_ino)

/*
 * Macro-instructions used to manage fragments
 */
#define ATFS_MIN_FRAG_SIZE		1024
#define	ATFS_MAX_FRAG_SIZE		4096
#define ATFS_MIN_FRAG_LOG_SIZE		  10
#define ATFS_FRAG_SIZE(s)		(ATFS_SB(s)->s_frag_size)
#define ATFS_FRAGS_PER_BLOCK(s)		(ATFS_SB(s)->s_frags_per_block)

/*
 * Constants relative to the data blocks
 */
#define	ATFS_NDIR_BLOCKS		12
#define	ATFS_IND_BLOCK			ATFS_NDIR_BLOCKS
#define	ATFS_DIND_BLOCK			(ATFS_IND_BLOCK + 1)
#define	ATFS_TIND_BLOCK			(ATFS_DIND_BLOCK + 1)
#define	ATFS_N_BLOCKS			(ATFS_TIND_BLOCK + 1)

/*
 * Define ATFS_RESERVATION to reserve data blocks for expanding files
 */
#define ATFS_DEFAULT_RESERVE_BLOCKS     8
/*max window size: 1024(direct blocks) + 3([t,d]indirect blocks) */
#define ATFS_MAX_RESERVE_BLOCKS         1027
#define ATFS_RESERVE_WINDOW_NOT_ALLOCATED 0
/*
 * ATFS file system version
 */
#define ATFSFS_DATE		"2020/10/20"
#define ATFSFS_VERSION		"0.1b"

/*
 * Mount flags
 */
#define ATFS_MOUNT_CHECK		0x000001  /* Do mount-time checks */
#define ATFS_MOUNT_OLDALLOC		0x000002  /* Don't use the new Orlov allocator */
#define ATFS_MOUNT_GRPID		0x000004  /* Create files with directory's group */
#define ATFS_MOUNT_DEBUG		0x000008  /* Some debugging messages */
#define ATFS_MOUNT_ERRORS_CONT		0x000010  /* Continue on errors */
#define ATFS_MOUNT_ERRORS_RO		0x000020  /* Remount fs ro on errors */
#define ATFS_MOUNT_ERRORS_PANIC		0x000040  /* Panic on errors */
#define ATFS_MOUNT_MINIX_DF		0x000080  /* Mimics the Minix statfs */
#define ATFS_MOUNT_NOBH			0x000100  /* No buffer_heads */
#define ATFS_MOUNT_NO_UID32		0x000200  /* Disable 32-bit UIDs */
#define ATFS_MOUNT_XATTR_USER		0x004000  /* Extended user attributes */
#define ATFS_MOUNT_POSIX_ACL		0x008000  /* POSIX Access Control Lists */
#define ATFS_MOUNT_XIP			0x010000  /* Obsolete, use DAX */
#define ATFS_MOUNT_USRQUOTA		0x020000  /* user quota */
#define ATFS_MOUNT_GRPQUOTA		0x040000  /* group quota */
#define ATFS_MOUNT_RESERVATION		0x080000  /* Preallocation */
#define ATFS_MOUNT_DAX			0x100000  /* Direct Access */


#define clear_opt(o, opt)		o &= ~ATFS_MOUNT_##opt
#define set_opt(o, opt)			o |= ATFS_MOUNT_##opt
#define test_opt(sb, opt)		(ATFS_SB(sb)->s_mount_opt & \
					 ATFS_MOUNT_##opt)

static inline spinlock_t *
sb_bgl_lock(struct atfs_sb_info *sbi, unsigned int block_group)
{
	return bgl_lock_ptr(sbi->s_blockgroup_lock, block_group);
}

/*
 * Ok, these declarations are also in <linux/kernel.h> but none of the
 * atfs source programs needs to include it so they are duplicated here.
 */

static inline struct atfs_inode_info *ATFS_I(struct inode *inode)
{
	return container_of(inode, struct atfs_inode_info, vfs_inode);
}

static inline atfs_fsblk_t
atfs_group_first_block_no(struct super_block *sb, unsigned long group_no)
{
	return group_no * (atfs_fsblk_t)ATFS_BLOCKS_PER_GROUP(sb) +
		le32_to_cpu(ATFS_SB(sb)->s_es->s_first_data_block);
}

#define atfs_set_bit	__test_and_set_bit_le
#define atfs_clear_bit	__test_and_clear_bit_le
#define atfs_test_bit	test_bit_le
#define atfs_find_first_zero_bit	find_first_zero_bit_le
#define atfs_find_next_zero_bit		find_next_zero_bit_le

#define atfs_set_bit_atomic(l, nr, addr)	test_and_set_bit_le(nr, addr)
#define atfs_clear_bit_atomic(l, nr, addr)	test_and_clear_bit_le(nr, addr)


