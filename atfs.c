#include <linux/init.h>
#include <linux/module.h>
#include <linux/fs.h>
#include <linux/mount.h>

static struct vfsmount *atfs_mnt;

static int atfs_set_super(struct super_block *sb, void *data)
{
	sb->s_fs_info = data;
	return set_anon_super(sb, data);
}

static struct dentry *atfs_mount(struct file_system_type *fs_type, int flags,
		   const char *dev_name, void *data)
{
	struct super_block *sb = NULL;
	struct inode *inode;
	struct dentry *root;
	printk(KERN_INFO "atfs_mount...");
	sb = sget(fs_type, NULL, atfs_set_super, flags, NULL);
	if (!sb) {
		printk(KERN_INFO "atfs super_block is NULL");
		return ERR_PTR(-EINVAL);
	}
	inode = new_inode(sb);
	root = d_make_root(inode);		
	sb->s_root = root;
	return dget(sb->s_root);
}

static struct file_system_type atfs_fs_type = {
	.name	= "atfs",
	.owner	= THIS_MODULE,
	.mount	= atfs_mount,
};

ssize_t atfs_file_read(struct file *filp, char __user *buf, 
		size_t count, loff_t *ppos)
{
	return 0;
}

ssize_t atfs_file_write(struct file *filp, const char __user *buf, 
		size_t count, loff_t *ppos)
{
	return 0;
}

int atfs_open_file(struct inode *inode, struct file *filp)
{
	return 0;
}

struct dentry *atfs_create_file(const char* name)
{
	struct dentry *root_dentry = NULL;
	struct inode *inode = NULL;
	BUG_ON(!atfs_mnt);
	BUG_ON(!atfs_mnt->mnt_sb);
	BUG_ON(!atfs_mnt->mnt_sb->s_root);
	root_dentry = atfs_mnt->mnt_sb->s_root;
	inode = new_inode(root_dentry->d_inode->i_sb);
	if (root_dentry) {
		d_instantiate(root_dentry, inode);
		dget(root_dentry);
	}
	return ERR_PTR(-EINVAL);
}

struct dentry *atfs_create_dir(const char *name)
{
	return atfs_create_file(name);
}
struct file_operations atfs_file_operations = {
	read:		atfs_file_read,
	write:		atfs_file_write,
	open:		atfs_open_file,
};

/*
 * register atfs
 */
int atfs_register(void)
{
	int ret;
	printk(KERN_INFO "atfs_register...");
	ret = register_filesystem(&atfs_fs_type);
	if (ret < 0) {
		printk(KERN_ERR "register file system fail %d", ret);
		return ret;
	} else {
		printk(KERN_INFO "register file system succeed %d", ret);
		atfs_mnt = kern_mount(&atfs_fs_type);
		if (IS_ERR(atfs_mnt)) {
			printk(KERN_ERR "mount filesystem error %ld", PTR_ERR(atfs_mnt));
			return PTR_ERR(atfs_mnt);
		}
	}
	atfs_create_dir("testDir");
	return 0;
}

int atfs_unregister(void)
{
	int ret;
	printk(KERN_INFO "atfs_unregister...");
	ret = unregister_filesystem(&atfs_fs_type);
	if (ret < 0) {
		printk(KERN_ERR "unregister file system fail %d", ret);
		return ret;
	}
	return 0;
}

static int __init atfs_init(void)
{
	printk(KERN_INFO "atfs init...\n");
	atfs_register();
	return 0;
}
module_init(atfs_init);

static void __exit atfs_exit(void)
{
	printk(KERN_INFO "atfs exit...\n");
	atfs_unregister();
	return;
}
module_exit(atfs_exit);

MODULE_AUTHOR("Enze Li");
MODULE_LICENSE("GPL v2");
MODULE_DESCRIPTION("A Tiny File System");
MODULE_ALIAS("atfs");
