EXTRA_CFLAGS := -I$(PWD) -I$(PWD)/../../include -I$(PWD)/../../models/common  -I$(PWD)/../../models/int -I/usr/realtime/include -I/usr/include  -I/usr/local/comedi/include  -ffast-math -mhard-float -lm
KDIR         := /lib/modules/$(shell uname -r)/build
PWD          := $(shell pwd)
KERNEL_24    := $(if $(wildcard $(KDIR)/Rules.make),1,0)
obj-m        +=  words_rt.o 
words_rt-objs := words.o ../../models/common/nm.o ../../models/int/int.o wordsBuffer.o


.PHONY: all clean modules_install

ifeq ($(KERNEL_24),0)
ifeq ($(KERNELRELEASE),)
all: 	# 2.6 !!!!!!!!!!!!!!!!!
	$(MAKE)  $(EXTRA_CGLAGS) -C  $(KDIR) SUBDIRS=$(PWD)  modules
	cp words_rt.ko ../../modules/neurons/words_rt.ko
clean modules_install:
	$(MAKE) -C $(KDIR) M=$(PWD) SUBDIRS=$(PWD) $@
endif # ifeq ($(KERNELRELEASE),)




else  #################### ifeq ($(KERNEL_24),0)

ifneq ($(KERNELRELEASE),)

export-objs := <export-object-list>

include $(KDIR)/Rules.make
<module-name>.o: $(<module-name>-objs)
	$(Q)$(LD) $(LD_RFLAG) -r -o $@ $(<module-name>-objs)
else  # ifneq ($(KERNELRELEASE),)
all:
	
#	$(MAKE) -C $(KDIR) SUBDIRS=$(PWD) modules
#	$(MAKE)   $(EXTRA_CGLAGS) -C  $(KDIR) SUBDIRS=$(PWD)  modules
#	ld  -r  -m elf_i386 hh.ko  ../../common/nm.ko -o hh.ko
#	mv hh.ko ../../../modules/neurons
clean:
	rm -f *.o
	rm -f *.ko
	rm -f *.cmd
	rm -f *.o.flags
	rm -f *.mod.c
	rm -f *.*~
	rm -f ../modules/hh.ko

endif # ifneq ($(KERNELRELEASE),)


endif #################### ifeq ($(KERNEL_24),0)

