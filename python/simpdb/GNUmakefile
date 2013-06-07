subdirs =

clean:
	$(RM) -f *~ *.m.py */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

# create self-contained python scripts
mpy:
	python pyembed.py -a -o simpdb.m.py simpdb.py

.PHONY: clean mpy

