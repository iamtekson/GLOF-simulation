
MODULE_TOPDIR = ../..

SUBDIRS1 = \
	R_code \
	r.avaflow \
	r.avaflow.main \
	r.avaflow.mult

SUBDIRS2 = \
	r.avaflow \
	r.avaflow.main \
	r.avaflow.mult

SUBDIRS = $(SUBDIRS1)

include $(MODULE_TOPDIR)/include/Make/Dir.make

default: parsubdirs

SUBDIRS = $(SUBDIRS2)

install: installsubdirs
