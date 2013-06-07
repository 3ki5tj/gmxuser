
clean:
	rm -f *.m.py *.pyc *~
	rstrip.py -lRv

subdirs =
#-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done

mpy:
	python pyembed.py \
	  --mod="cspacer,cindent,ccmdx,ccgmx,cc,ccutil,gmxcom" \
	  -o gmxspdr.m.py gmxspdr.py


# backup the directory to usb, delete useless files
usb::
	$(MAKE) -C .. usb

.PHONY: clean usb mpy

