CC = gcc -Wall -Wextra -pedantic
LM = -lm

deps = GNUmakefile
cprogs = mkspx


subdirs =

clean:
	$(RM) -f *~ *.m.py */*~ */*/*~ */a.out *.tmp $(cprogs) out*.pdb *.pyc
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

# create self-contained python scripts
mpy:
	python pyembed.py -a -o simpdb.m.py simpdb.py

$(cprogs) : % : %.c $(deps)
	$(CC) -o $@ $< $(LM)

.PHONY: clean mpy

