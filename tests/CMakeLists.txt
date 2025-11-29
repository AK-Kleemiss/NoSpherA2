include(CTest)
find_package(Python COMPONENTS Interpreter REQUIRED)

FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runTest.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_nos_test(alanine_occ "-cif alanine.cif -dmin 0.5 -wfn alanine.owf.fchk -acc 1 -no-date")
#add_nos_test(alanine_integrated_occ "-occ alanine.toml -cif alanine.cif -dmin 0.5 -acc 1 -no-date")
add_nos_test(disorder_THPP "-cif thpp.cif -hkl thpp.hkl -mtc olex2/Wfn_job/Part_1/thpp.wfx 0.1 olex2/Wfn_job/Part_2/thpp.wfx 0.2 -mtc_mult 1 1 -mtc_ECP 0 0 -mtc_charge 0 0 -acc 0 -no-date"
    COMPARE_GOOD "disorder_THPP.good" DIRECTORY "disorder")
add_nos_test(fourier_transform "-test_analytical -no-date" DIRECTORY "SALTED")
add_nos_test(fractal "-fractal sucrose_diff.cube -no_date" COMPARE_GOOD "fractal.good" COMPARE_NAME "sucrose_diff.cube_fractal_plot" DIRECTORY "sucrose_fchk_SF")
add_nos_test(grown_water "-cif water.cif -hkl water.hkl -wfn water.wfx -acc 0 -no-date" DIRECTORY "grown")
add_nos_test(Hybrid_mode "-cif ZP2.cif -dmin 0.9 -acc 0 -mtc ZP2_part1.gbw 0.1 ZP2_part2.gbw 0.2 -mtc_mult 1 1 -mtc_charge 0 0 -mtc_ECP 0 0 -no-date" DIRECTORY "Hybrid")
add_nos_test(malbac_SF_ECP "-cif malbac.cif -hkl malbac.hkl -wfn malbac.gbw -acc 0 -ECP 1 -no-date" DIRECTORY "ECP_SF")
add_nos_test(openBLAS "-blastest -no-date" DIRECTORY "OpenBLAS")
# properties is failing, but it seems like the output changed. It doesn't look too bad
add_nos_test(properties "-wfn olex2/Wfn_job/sucrose.wfx -cif sucrose.cif -lap -eli -resolution 0.5 -no-date" COMPARE_NAME "NoSpherA2.log" DIRECTORY "sucrose_fchk_SF")
add_nos_test(reading_SALTED "-test_reading_SALTED_binary" DIRECTORY "SALTED")
add_nos_test(ri_fit "-wfn epoxide.gbw -cif epoxide.cif -dmin 0.4 -ri_fit combo_basis_fit -no-date" DIRECTORY "epoxide_gbw" COMPARE_GOOD "ri_fit.good")
add_nos_test(rubredoxin_cmtc "-hkl 1yk4_h.hkl -acc 0 -cmtc residues/1.gbw residues/1.cif 0 residues/2.gbw residues/2.cif 0.1 residues/3.gbw residues/3.cif 0.2 residues/4.gbw residues/4.cif 0 -mtc_mult 1 1 1 1 -mtc_charge 0 0 0 0 -mtc_ECP 0 0 0 0 -no-date" DIRECTORY "rubredoxin_cmtc")
add_nos_test(SALTED "-SALTED Model -cif test_cysteine.cif -wfn test_cysteine.xyz -dmin 0.73 -no-date" DIRECTORY "SALTED")
add_nos_test(sucrose_IAM "-cif sucrose.cif -hkl sucrose.hkl -xyz sucrose.xyz -IAM -no-date" DIRECTORY "sucrose_IAM_SF")
add_nos_test(sucrose_ptb "-cif sucrose.cif -dmin 0.8 -wfn wfn.xtb -mult 0 -charge 0 -acc 0 -ECP 3 -no-date" DIRECTORY "sucrose_IAM_SF")
add_nos_test(sucrose_SF "-cif sucrose.cif -hkl olex2/Wfn_job/sucrose.hkl -wfn olex2/Wfn_job/sucrose.wfx -acc 0 -no-date" DIRECTORY "sucrose_fchk_SF")
add_nos_test(sucrose_twin "-cif sucrose.cif -hkl olex2/Wfn_job/sucrose.hkl -wfn olex2/Wfn_job/sucrose.wfx -acc 0 -twin 1 0 0 0 -1 0 0 1 -2 -no-date" DIRECTORY "sucrose_fchk_SF")
add_nos_test(wfn_reading "-hkl test.hkl -acc 0 -cif test.cif -wfn test.wfn -no-date" DIRECTORY "wfn_reading")

if (RUN_FULL_TEST)
    add_nos_test(fchk_conversion "-b dev2-TZVP -d ./ -wfn in.ffn -no-date" COMPARE_NAME "good.fchk" DIRECTORY "NiP3_fchk") # This is disabled on the actions, it seems to be missing a config file
    add_nos_test(fourier_transform_full "-test_analytical full -no-date" DIRECTORY "SALTED")
    # This doesn't seem to have a .good file
#    add_nos_test(rubredoxin_cmtc_full "-hkl 1yk4_h_full.hkl -acc 0 -cmtc residues/1.gbw residues/1.cif 0 residues/2.gbw residues/2.cif 0.1 residues/3.gbw residues/3.cif 0.2 residues/4.gbw residues/4.cif 0 -mtc_mult 1 1 1 1 -mtc_charge 0 0 0 0 -mtc_ECP 0 0 0 0 -cpus 24 -no-date" DIRECTORY "rubredoxin_cmtc")
endif ()
