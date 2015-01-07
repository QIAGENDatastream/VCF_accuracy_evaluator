default:

	$(MAKE) -C vt
	cmake bam-readcount
	$(MAKE) -C bam-readcount deps
	$(MAKE) -C bam-readcount
	$(MAKE) -C bedtools
	$(MAKE) -C tabix-0.2.6
