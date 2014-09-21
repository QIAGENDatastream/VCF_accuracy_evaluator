default:

	$(MAKE) -C vt
	cmake bam-readcount
	$(MAKE) -C bam-readcount deps
	$(MAKE) -C bam-readcount
