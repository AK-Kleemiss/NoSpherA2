RM := rm

all: sucrose_SF sucrose_IAM fractal disorder_THPP grown_water

sucrose_SF:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx -acc 1\
		&& cp NoSpherA2.log $@.log
	@echo 'Finished running: $@'

sucrose_IAM:
	@echo 'Running test: $@'
	cd sucrose_IAM_SF && ../../NoSpherA2 \
		-cif sucrose.cif \
		-hkl sucrose.hkl \
		-xyz sucrose.xyz -IAM\
		&& cp NoSpherA2.log $@.log
	@echo 'Finished running: $@'

disorder_THPP:
	@echo 'Running test: $@'
	cd disorder && ../../NoSpherA2 \
		-cif thpp.cif \
		-hkl thpp.hkl \
		-mtc olex2/Wfn_job/Part_1/thpp.wfx 0,1\
		olex2/Wfn_job/Part_2/thpp.wfx 0,2 -acc 1\
		&& cp NoSpherA2.log $@.log
	@echo 'Finished running: $@'

grown_water:
	@echo 'Running test: $@'
	cd grown && ../../NoSpherA2 \
		-cif water.cif \
		-hkl water.hkl \
		-wfn water.wfx \
		-acc 1 \
		&& cp NoSpherA2.log $@.log
	@echo 'Finished running: $@'

fractal:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-fractal sucrose_diff.cube \
		&& cp NoSpherA2.log $@.log
	@echo 'Finished running: $@'

clean:
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


