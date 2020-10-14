#include <linux/init.h>
#include <linux/module.h>
#include <linux/fs.h>
#include <linux/mount.h>
#include <linux/slab.h>

static struct vfsmount *atfs_mnt;

const struct inode_operations atfs_dir_inode_operations;
const struct file_operations atfs_dir_operations;

struct atfs_dir_entry {
	__le32	inode;			/* Inode number */
	__le16	rec_len;		/* Directory entry length */
	__u8	name_len;		/* Name length */
	__u8	file_type;
	char	name[];			/* File name */
};

typedef struct atfs_dir_entry atfs_dirent;

static inline void atfs_put_page(struct page *page)
{
	//kunmap(page);
	//put_page(page);
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
	BUG();
	printk(KERN_INFO "==atfs== atfs_file_open invoked...");
	return 0;
}

static int atfs_mkdir(struct inode * dir, struct dentry * dentry, umode_t mode)
{
	struct inode * inode;

	inode = new_inode(dir->i_sb);
	inode->i_mode |= S_IFDIR;
	inode->i_op = &atfs_dir_inode_operations;
	inode->i_fop = &atfs_dir_operations;
	inode->i_state = I_NEW;

	d_instantiate_new(dentry, inode);
	return 0;
}

const struct inode_operations atfs_dir_inode_operations = {
	.lookup	= atfs_lookup,
	.mkdir	= atfs_mkdir,
};

const struct file_operations atfs_dir_operations = {
};

const struct file_operations atfs_file_operations = {
	.open	= atfs_file_open,
	.write	= atfs_file_write,
	.read	= atfs_file_read,
};

static struct dentry *atfs_mount(struct file_system_type *fs_type, int flags,
		   const char *dev_name, void *data)
{
	struct super_block *sb = NULL;
	struct inode *inode;
	struct dentry *root;
	if (atfs_mnt)
		return ERR_PTR(-EINVAL);
	printk(KERN_INFO "==atfs== atfs_mount invoked...");
	sb = sget(fs_type, NULL, atfs_set_super, flags, NULL);
	printk(KERN_INFO "==atfs== super_block:%p", sb);
	if (!sb) {
		printk(KERN_INFO "==atfs== super_block is NULL");
		return ERR_PTR(-EINVAL);
	}
	inode = new_inode(sb);
	inode->i_mode |= S_IFDIR;
	inode->i_op = &atfs_dir_inode_operations;
	inode->i_fop = &atfs_file_operations;
	root = d_make_root(inode);
	sb->s_root = root;
	return dget(sb->s_root);
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
