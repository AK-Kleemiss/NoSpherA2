RM := rm

ALL_TARGETS = sucrose_SF sucrose_IAM sucrose_twin fractal disorder_THPP grown_water properties rubredoxin_cmtc wfn_reading malbac_SF_ECP sucrose_ptb Hybrid_mode openBLAS ri_fit

# Convert the path to an OS-specific format if necessary
# For Unix-like systems, this will be a no-op
ifeq ($(OS),Windows_NT)
  BASE_PATH := ..\..\NoSpherA2.exe
  MV := move
  DIFF := fc /c
else
  BASE_PATH := ../../NoSpherA2
  MV := mv
  DIFF := diff -q -i -b
endif

all: $(ALL_TARGETS) fourier_transform

all_tests: $(ALL_TARGETS) reading_SALTED SALTED fourier_transform_full rubredoxin_cmtc_full fchk_conversion

sucrose_SF:
	@echo 'Running test: $@'
	@cd sucrose_fchk_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

openBLAS:
	@echo 'Running test: $@'
	@cd OpenBLAS && $(BASE_PATH) \
		-blastest \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

malbac_SF_ECP:
	@echo 'Running test: $@'
	@cd ECP_SF && $(BASE_PATH) \
		-cif malbac.cif \
		-hkl malbac.hkl \
		-wfn malbac.gbw \
		-acc 0 \
		-ECP 1 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

SALTED:
	@echo 'Running test: $@'
	@cd SALTED && $(BASE_PATH) \
		-SALTED Model \
		-cif test_cysteine.cif \
		-wfn test_cysteine.xyz \
		-dmin 0.73 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

reading_SALTED:
	@echo 'Running test: $@'
	@cd SALTED && $(BASE_PATH) \
		-test_reading_SALTED_binary \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

sucrose_IAM:
	@echo 'Running test: $@'
	@cd sucrose_IAM_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl sucrose.hkl \
		-xyz sucrose.xyz -IAM\
		-no-date\
		&& ${MV} NoSpherA2.log $@.log\
		&& ${DIFF} $@.log $@.good
	
sucrose_twin:
	@echo 'Running test: $@'
	@cd sucrose_fchk_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-twin 1 0 0 0 -1 0 0 1 -2 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

sucrose_ptb:
	@echo 'Running test: $@'
	@cd sucrose_IAM_SF && $(BASE_PATH) \
		-cif sucrose.cif \
		-dmin 0.8 \
		-wfn wfn.xtb \
		-mult 0 \
		-charge 0 \
		-acc 0 \
		-ECP 3 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log\
		&& ${DIFF} $@.log $@.good

disorder_THPP:
	@echo 'Running test: $@'
	@cd disorder && $(BASE_PATH) \
		-cif thpp.cif \
		-hkl thpp.hkl \
		-mtc olex2/Wfn_job/Part_1/thpp.wfx 0,1 \
		olex2/Wfn_job/Part_2/thpp.wfx 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-mtc_ECP 0 0 \
		-acc 0 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

grown_water:
	@echo 'Running test: $@'
	@cd grown && $(BASE_PATH) \
		-cif water.cif \
		-hkl water.hkl \
		-wfn water.wfx \
		-acc 0 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

fourier_transform_full:
	@echo 'Running test: $@'
	@cd SALTED && $(BASE_PATH) -test_analytical full -no-date && ${MV} NoSpherA2.log $@.log

fourier_transform:
	@echo 'Running test: $@'
	@cd SALTED && $(BASE_PATH) -test_analytical -no-date && ${MV} NoSpherA2.log $@.log

fractal:
	@echo 'Running test: $@'
	@cd sucrose_fchk_SF && $(BASE_PATH) \
		-fractal sucrose_diff.cube \
		-no_date \
		&& ${MV} sucrose_diff.cube_fractal_plot $@.log \
		&& ${DIFF} $@.log $@.good
	
properties:
	@echo 'Running test: $@'
	@cd sucrose_fchk_SF && $(BASE_PATH) \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-cif sucrose.cif -lap -eli \
		-resolution 0.5 -test \
		-no-date \
		&& ${MV} NoSpherA2_cube.log $@.log \
		&& ${DIFF} $@.log $@.good
	
rubredoxin_cmtc:
	@echo 'Running test: $@'
	@cd rubredoxin_cmtc && $(BASE_PATH) \
		-hkl 1yk4_h.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-mtc_ECP 0 0 0 0 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

rubredoxin_cmtc_full:
	@echo 'Running test: $@'
	@cd rubredoxin_cmtc && $(BASE_PATH) \
		-hkl 1yk4_h_full.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-cpus 24


ri_fit:
	@echo 'Running test: $@'
	@cd epoxide_gbw && $(BASE_PATH) \
	    -wfn epoxide.gbw \
	    -cif epoxide.cif \
	    -dmin 0.4 \
	    -ri_fit \
	    -no-date \
	    && ${MV} NoSpherA2.log $@.log \
	    && ${DIFF} $@.log $@.good

Hybrid_mode:
	@echo 'Running test: $@'
	@cd Hybrid && $(BASE_PATH) \
		-cif ZP2.cif \
		-dmin 0.9 \
		-acc 0 \
		-mtc ZP2_part1.gbw 0,1 ZP2_part2.gbw 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-mtc_ECP 0 0 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good
	
wfn_reading:
	@echo 'Running test: $@'
	@cd wfn_reading && $(BASE_PATH) \
		-hkl test.hkl -acc 0 \
		-cif test.cif -wfn test.wfn \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

fchk_conversion:
	@echo 'Running test: $@'
	@cd NiP3_fchk && $(BASE_PATH) \
                -b def2-TZVP -d ./ \
		-wfn in.ffn -no-date \
		&& ${DIFF} in.fchk good.fchk

.PHONY: all clean dependents
.SECONDARY:


