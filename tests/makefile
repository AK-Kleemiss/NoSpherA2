RM := rm

DIFF := diff -q -i

all: sucrose_SF sucrose_IAM sucrose_twin fractal disorder_THPP grown_water properties

sucrose_SF:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 1\
		-no_date\
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

malbac_SF_ECP:
	@echo 'Running test: $@'
	cd ECP_SF && ../../NoSpherA2 \
		-cif malbac.cif \
		-hkl malbac.hkl \
		-wfn malbac_ECP.molden \
		-acc 1\
		-ECP 1\
		-no_date\
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

sucrose_IAM:
	@echo 'Running test: $@'
	cd sucrose_IAM_SF && ../../NoSpherA2 \
		-cif sucrose.cif \
		-hkl sucrose.hkl \
		-xyz sucrose.xyz -IAM\
		-no_date\
		&& mv NoSpherA2.log $@.log\
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
sucrose_twin:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 1 \
		-twin 1 0 0 0 -1 0 0 1 -2\
		-no_date\
		&& mv NoSpherA2.log $@.log
	@echo 'Finished running: $@'

disorder_THPP:
	@echo 'Running test: $@'
	cd disorder && ../../NoSpherA2 \
		-cif thpp.cif \
		-hkl thpp.hkl \
		-mtc olex2/Wfn_job/Part_1/thpp.wfx 0,1\
		olex2/Wfn_job/Part_2/thpp.wfx 0,2\
		-acc 1\
		-no_date\
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

grown_water:
	@echo 'Running test: $@'
	cd grown && ../../NoSpherA2 \
		-cif water.cif \
		-hkl water.hkl \
		-wfn water.wfx \
		-acc 1 \
		-no_date\
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

fractal:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-fractal sucrose_diff.cube \
		-no_date\
		&& mv sucrose_diff.cube_fractal_plot $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
properties:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && ../../NoSpherA2 \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-cif sucrose.cif -lap -eli\
		-resolution 0.5 -test \
		-no_date\
		&& mv NoSpherA2_cube.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

clean:
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


