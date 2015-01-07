default:

	$(MAKE) -C vt
	cd bam-readcount/build; cmake ../
	$(MAKE) -C bam-readcount/build deps
	$(MAKE) -C bam-readcount/build 
	$(MAKE) -C bedtools
	$(MAKE) -C tabix-0.2.6
