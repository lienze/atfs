#include <linux/init.h>
#include <linux/module.h>
#include <linux/fs.h>

static struct file_system_type atfs_fs_type = {
	.name	= "atfs",
	.owner	= THIS_MODULE,
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
	}
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
