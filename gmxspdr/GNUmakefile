
clean:
	find . -name "*~" | xargs rm -f
	rstrip.py -lRv

subdirs =
#-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done


# backup the directory to usb, delete useless files
usb::
	$(MAKE) -C .. usb


