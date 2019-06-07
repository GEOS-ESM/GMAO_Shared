# include the top level include file
include ../Make.inc
#
DIRS=diag psp aztec blocksolve
#
.DEFAULT:
	@$(ECHO) 'Unknown target $@'
	@$(ECHO) 'The only defined target is: clean'
#
clean:
	for d in $(DIRS); do ($(CD) $$d; $(MAKE) clean); done
#
