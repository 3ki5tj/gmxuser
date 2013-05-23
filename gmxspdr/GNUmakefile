
clean:
	find . -name "*~" | xargs rm -f
	rstrip.py -R *.txt *.tex *.htm* *.c *.h *.py *.ma *.gp \
	  *akefile README*

subdirs =
#-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done


# backup the directory to usb, delete useless files
usb::
	$(MAKE) -C .. usb


