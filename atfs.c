#include <linux/init.h>
#include <linux/module.h>
#include <linux/fs.h>

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
	printk(KERN_INFO "atfs_register...");
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
}
module_exit(atfs_exit);

MODULE_AUTHOR("Hacker");
MODULE_LICENSE("GPL v2");
MODULE_DESCRIPTION("A Test File System");
MODULE_ALIAS("atfs module");
