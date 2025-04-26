RM := rm

ALL_TARGETS = sucrose_SF sucrose_IAM sucrose_twin fractal disorder_THPP grown_water properties rubredoxin_cmtc wfn_reading malbac_SF_ECP sucrose_ptb Hybrid_mode openBLAS ri_fit fourier_transform reading_SALTED SALTED alanine_occ

# Convert the path to an OS-specific format if necessary
# For Unix-like systems, this will be a no-op
ifeq ($(OS),Windows_NT)
  BASE_PATH := ..\..\NoSpherA2.exe
  MV := move
  DIFF := fc /c
  RM := del /q /f

  define RUN_TEST
    @echo 'Running test: $1'
    @echo @echo off > run_test.bat
    @echo cd $2 >> run_test.bat
    @echo "$(BASE_PATH)" $3 >> run_test.bat
    @call run_test.bat
    @del run_test.bat
    @if exist $2\NoSpherA2.log (cd $2 && $(MV) NoSpherA2.log $1.log && $(DIFF) $1.log $1.good)
  endef
else
  RM := rm
  BASE_PATH := ../../NoSpherA2
  MV := mv
  DIFF := diff -q -i -b
  define TIME_START
    start_time=$$(date +%s.%N)
  endef
  
  define TIME_END
    end_time=$$(date +%s.%N); \
    elapsed=$$(echo "$$end_time - $$start_time" | bc); \
    printf "Time elapsed: %.3f seconds\n" $$elapsed
  endef
  define RUN_TEST
    @echo 'Running test: $1'
    @eval $(TIME_START); cd $2 && $(BASE_PATH) $3 && \
    if [ -f NoSpherA2.log ]; then $(MV) NoSpherA2.log $1.log; fi && \
    $(DIFF) $1.log $1.good; eval $(TIME_END)
  endef
endif


all: $(ALL_TARGETS)

all_tests: $(ALL_TARGETS) fourier_transform_full rubredoxin_cmtc_full fchk_conversion

alanine_occ:
	@echo 'Running test: $@'
	@cd alanine_occ && $(BASE_PATH) \
		-cif alanine.cif \
		-dmin 0.5 \
		-wfn alanine.owf.fchk \
		-acc 1 \
		-no-date \
		&& ${MV} NoSpherA2.log $@.log \
		&& ${DIFF} $@.log $@.good

sucrose_SF:
	$(call RUN_TEST,$@,sucrose_fchk_SF, \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-no-date)

openBLAS:
	$(call RUN_TEST,$@,OpenBLAS, \
		-blastest \
		-no-date)

malbac_SF_ECP:
	$(call RUN_TEST,$@,ECP_SF, \
		-cif malbac.cif \
		-hkl malbac.hkl \
		-wfn malbac.gbw \
		-acc 0 \
		-ECP 1 \
		-no-date)

SALTED:
	$(call RUN_TEST,$@,SALTED, \
		-SALTED Model \
		-cif test_cysteine.cif \
		-wfn test_cysteine.xyz \
		-dmin 0.73 \
		-no-date)

reading_SALTED:
	$(call RUN_TEST,$@,SALTED, \
		-test_reading_SALTED_binary)

sucrose_IAM:
	$(call RUN_TEST,$@,sucrose_IAM_SF, \
		-cif sucrose.cif \
		-hkl sucrose.hkl \
		-xyz sucrose.xyz -IAM\
		-no-date)
	
sucrose_twin:
	$(call RUN_TEST,$@,sucrose_fchk_SF, \
		-cif sucrose.cif \
		-hkl olex2/Wfn_job/sucrose.hkl \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-acc 0 \
		-twin 1 0 0 0 -1 0 0 1 -2 \
		-no-date)

sucrose_ptb:
	$(call RUN_TEST,$@,sucrose_IAM_SF, \
		-cif sucrose.cif \
		-dmin 0.8 \
		-wfn wfn.xtb \
		-mult 0 \
		-charge 0 \
		-acc 0 \
		-ECP 3 \
		-no-date)

disorder_THPP:
	$(call RUN_TEST,$@,disorder, \
		-cif thpp.cif \
		-hkl thpp.hkl \
		-mtc olex2/Wfn_job/Part_1/thpp.wfx 0,1 \
		olex2/Wfn_job/Part_2/thpp.wfx 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-mtc_ECP 0 0 \
		-acc 0 \
		-no-date)

grown_water:
	$(call RUN_TEST,$@,grown, \
		-cif water.cif \
		-hkl water.hkl \
		-wfn water.wfx \
		-acc 0 \
		-no-date)

fourier_transform_full:
	$(call RUN_TEST,$@,SALTED, \
		-test_analytical full -no-date)

fourier_transform:
	$(call RUN_TEST,$@,SALTED, \
		-test_analytical -no-date)

fractal:
	$(call RUN_TEST,$@,sucrose_fchk_SF, \
		-fractal sucrose_diff.cube \
		-no_date \
		&& ${RM} NoSpherA2.log && ${MV} sucrose_diff.cube_fractal_plot $@.log)
	
properties:
	$(call RUN_TEST,$@,sucrose_fchk_SF, \
		-wfn olex2/Wfn_job/sucrose.wfx \
		-cif sucrose.cif -lap -eli \
		-resolution 0.5 -test \
		-no-date \
		&& ${MV} NoSpherA2_cube.log $@.log)
	
rubredoxin_cmtc:
	$(call RUN_TEST,$@,rubredoxin_cmtc, \
		-hkl 1yk4_h.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-mtc_ECP 0 0 0 0 \
		-no-date)

rubredoxin_cmtc_full:
	$(call RUN_TEST,$@,rubredoxin_cmtc, \
		-hkl 1yk4_h_full.hkl -acc 0 \
		-cmtc residues/1.gbw residues/1.cif 0 \
		residues/2.gbw residues/2.cif 0,1 \
		residues/3.gbw residues/3.cif 0,2 \
		residues/4.gbw residues/4.cif 0 \
		-mtc_mult 1 1 1 1 \
		-mtc_charge 0 0 0 0 \
		-mtc_ECP 0 0 0 0 \
		-cpus 24 \
		-no-date)

ri_fit:
	$(call RUN_TEST,$@,epoxide_gbw, \
	    -wfn epoxide.gbw \
	    -cif epoxide.cif \
	    -dmin 0.4 \
	    -ri_fit \
	    -no-date)

Hybrid_mode:
	$(call RUN_TEST,$@,Hybrid, \
		-cif ZP2.cif \
		-dmin 0.9 \
		-acc 0 \
		-mtc ZP2_part1.gbw 0,1 ZP2_part2.gbw 0,2 \
		-mtc_mult 1 1 \
		-mtc_charge 0 0 \
		-mtc_ECP 0 0 \
		-no-date)
	
wfn_reading:
	$(call RUN_TEST,$@,wfn_reading, \
		-hkl test.hkl -acc 0 \
		-cif test.cif -wfn test.wfn \
		-no-date)

fchk_conversion:
	$(call RUN_TEST,$@,NiP3_fchk, \
                -b def2-TZVP -d ./ \
		-wfn in.ffn -no-date)

.PHONY: all clean dependents SALTED reading_SALTED
.SECONDARY:
