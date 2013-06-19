zipfile = gmxcode.zip
usbdir = /media/C3/code

#subdirs = at4 deploy gmx4 gmx45 codenotes buildnotes spider util
subdirs =

clean:
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done
	-rstrip.py -Rv


zip: $(zipfile)

$(zipfile)::
	git archive --format=zip -9 HEAD > $@

usb: $(zipfile)
	mv $< $(usbdir)

.PHONY: clean zip usb

