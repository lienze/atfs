#include <linux/init.h>
#include <linux/module.h>

static int __init atfs_init(void)
{
	printk(KERN_INFO "atfs init...\n");
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
