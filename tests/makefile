RM := rm
DIFF := diff -q -i -b

# Define the base path
BASE_PATH = ../../NoSpherA2
ALL_TARGETS = sucrose_SF sucrose_IAM sucrose_twin fractal disorder_THPP grown_water properties rubredoxin_cmtc wfn_reading malbac_SF_ECP sucrose_ptb Hybrid_mode openBLAS SALTED fourier_transform

# Convert the path to an OS-specific format if necessary
# For Unix-like systems, this will be a no-op
# For Windows, you might need to convert '/' to '\'
ifeq ($(OS),Windows_NT)
  BASE_PATH := $(shell echo $(BASE_PATH) | sed 's/\//\\/g')
	BASE_PATH +=.exe
endif

all: $(ALL_TARGETS)

sucrose_SF:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

openBLAS:
	@echo 'Running test: $@'
	cd OpenBLAS && $(BASE_PATH) \
		-blastest \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

malbac_SF_ECP:
	@echo 'Running test: $@'
	cd ECP_SF && $(BASE_PATH) \
		-cif malbac.cif \
		-hkl malbac.hkl \
		-wfn malbac.gbw \
		-acc 0 \
		-ECP 1 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

SALTED:
	@echo 'Running test: $@'
	cd SALTED && $(BASE_PATH) \
		-SALTED Model \
		-cif test_cysteine.cif \
		-wfn test_cysteine.xyz \
		-dmin 0.73 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

sucrose_IAM:
	@echo 'Running test: $@'
	cd sucrose_IAM_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl sucrose.hkl \
		-xyz sucrose.xyz -IAM\
		-no-date\
		&& mv NoSpherA2.log $@.log\
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
sucrose_twin:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-twin 1 0 0 0 -1 0 0 1 -2 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

sucrose_ptb:
	@echo 'Running test: $@'
	cd sucrose_IAM_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-dmin 0.8 \
		-wfn wfn.xtb \
		-mult 0 \
		-charge 0 \
		-acc 0 \
		-ECP 3 \
		-no-date \
		&& mv NoSpherA2.log $@.log\
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

disorder_THPP:
	@echo 'Running test: $@'
	cd disorder && $(BASE_PATH) \
		-cif thpp.cif \
		-hkl thpp.hkl \
		-mtc olex2/Wfn_job/Part_1/thpp.wfx 0,1 \
		olex2/Wfn_job/Part_2/thpp.wfx 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-acc 0 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

grown_water:
	@echo 'Running test: $@'
	cd grown && $(BASE_PATH) \
		-cif water.cif \
		-hkl water.hkl \
		-wfn water.wfx \
		-acc 0 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

fourier_transform_full:
	@echo 'Running test: $@'
	cd SALTED && $(BASE_PATH) -test_analytical full -no-date && mv NoSpherA2.log $@.log
	@echo 'Finished running: $@'

fourier_transform:
	@echo 'Running test: $@'
	cd SALTED && $(BASE_PATH) -test_analytical -no-date && mv NoSpherA2.log $@.log
	@echo 'Finished running: $@'

fractal:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && $(BASE_PATH) \
		-fractal sucrose_diff.cube \
		-no_date \
		&& mv sucrose_diff.cube_fractal_plot $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
properties:
	@echo 'Running test: $@'
	cd sucrose_fchk_SF && $(BASE_PATH) \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-cif sucrose.cif -lap -eli \
		-resolution 0.5 -test \
		-no-date \
		&& mv NoSpherA2_cube.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
rubredoxin_cmtc:
	@echo 'Running test: $@'
	cd rubredoxin_cmtc && $(BASE_PATH) \
		-hkl 1yk4_h.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

rubredoxin_cmtc_full:
	@echo 'Running test: $@'
	cd rubredoxin_cmtc && $(BASE_PATH) \
		-hkl 1yk4_h_full.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-cpus 24
	@echo 'Finished running: $@'

Hybrid_mode:
	@echo 'Running test: $@'
	cd Hybrid && $(BASE_PATH) \
		-cif ZP2.cif \
		-dmin 0.9 \
		-acc 0 \
		-mtc ZP2_part1.gbw 0,1 ZP2_part2.gbw 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'
	
wfn_reading:
	@echo 'Running test: $@'
	cd wfn_reading && $(BASE_PATH) \
		-hkl test.hkl -acc 0 \
		-cif test.cif -wfn test.wfn \
		-no-date \
		&& mv NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	@echo 'Finished running: $@'

fchk_conversion:
	@echo 'Running test: $@'
	cd NiP3_fchk && $(BASE_PATH) \
                -b def2-TZVP -d ./ \
		-wfn in.ffn -no-date \
		&& ${DIFF} in.fchk good.fchk
	@echo 'Finished running: $@'

.PHONY: all clean dependents
.SECONDARY:


