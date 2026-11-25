#include "pch.h"
#include "convenience.h"
#include "help_text.h"


static const char T_getting[] = "GETTING STARTED";
static const char T_core[] = "CORE INPUT AND SCATTERING-FACTOR CALCULATION";
static const char T_partitioning[] = "PARTITIONING, BASIS, AND PERFORMANCE";
static const char T_gpu[] = "GPU ACCELERATION";
static const char T_tsctscb[] = "TSC/TSCB TABLE UTILITIES";
static const char T_multifragment[] = "MULTI-FRAGMENT WORKFLOWS";
static const char T_density[] = "DENSITY, PROPERTIES, AND ANALYSIS";
static const char T_bonding[] = "BONDING ANALYSIS";
static const char T_conversion[] = "CONVERSION, ML, AND SPECIALISED TOOLS";
static const char T_diagnostics[] = "DIAGNOSTICS AND INTERNAL TOOLS";
static const char T_experimental[] = "EXPERIMENTAL AND DEVELOPER COMMANDS";
static const char T_help[] = "HELP";
static const char T_examples[] = "EXAMPLES";

std::vector<help_item>& help_items()
{
	static std::vector<help_item> items = {

	{"", "", "", "", ""},
	{"", "", "", "", "============================================================================\n" " NoSpherA2 -- non-spherical scattering factors and density analysis\n" "============================================================================\n" "Syntax: NoSpherA2 [options]\n" "Values in [brackets] are optional.  Repeatable options may be supplied more\n" "than once.  Input files and options may be given in any order unless noted.\n" "Every option accepts either separator: -rgbi_groups and -rgbi-groups are the\n" "same flag, and the spelling below is the canonical one."},
	{"", "", "", "", ""},
	{T_getting, "", "", "", "GETTING STARTED\n" "  HAR/TSC:  -cif model.cif -hkl data.hkl -wfn wavefunction.wfx -acc 2\n" "  IAM TSC:  -cif model.cif -hkl data.hkl -xyz model.xyz -IAM -acc 2"},
	{T_getting, "", "", "", ""},
	{T_core, "", "", "", "CORE INPUT AND SCATTERING-FACTOR CALCULATION"},
	{T_core, "", "-cif", "<file.cif>", "Structure and atom labels."},
	{T_core, "", "-wfn", "<file>", "Wavefunction input (.wfn, .wfx, .ffn,\n" ".molden, .gbw, .xtb, or fchk/fch);\n" "ORCA natural orbitals .mdci.nat (CCSD)\n" "and .mp2nat (relaxed MP2), or an\n" "orca_2json .json with an AutoCI density\n" "(CCSD(T): %autoci density unrelaxed).\n" "Correlated input is natural orbitals with\n" "fractional occupations: these ORCA files,\n" "a CASSCF .gbw (ActOrbs NatOrbs; full CI\n" "with job.densities beside it) or a molden\n" "with Occup= (PySCF). A Gaussian fchk is\n" "read as its SCF. Engines, keywords, what\n" "is exact or approximate: wiki page\n" "Correlated-Wavefunctions."},
	{T_core, "", "-cube_density", "<file.cube>", "Experimental/computed electron-density\n" "cube used for partitioned SF generation\n" "(requires -cif and -hkl/-dmin)."},
	{T_core, "", "-occ", "<file.toml>", "Run an OCC wavefunction calculation.\n" "method = \"mp2\": relaxed MP2 density;\n" "\"ccsd\" / \"ccsd(t)\": unrelaxed CCSD /\n" "CCSD(T) density; \"mp3\", \"mp4(sdq)\",\n" "\"mp4\": unrelaxed MPn as ORCA's AutoCI\n" "(closed shell). [scf]\n" "keys: ccsd_backend = exact | df (or a\n" "df-/ri- method prefix, ri-basis = aux\n" "basis), ccsd_frozen_core (-1 core, 0\n" "all electrons), ccsd_tol (1e-5, ORCA's\n" "STol), ccsd_max_memory_gb, ccsd_local.\n" "Its natural orbitals feed the scattering\n" "factors only, not -fba/-nbo; wiki\n" "Correlated-Wavefunctions."},
	{T_core, "", "-xyz", "<file.xyz>", "Atomic positions for IAM or SALTED."},
	{T_core, "", "-hkl", "<file.hkl>", "Reflection list."},
	{T_core, "", "-dmin", "<angstrom>", "Generate the resolution sphere to this\n" "d-spacing instead of reading -hkl. With\n" "-hkl_min_max only the symmetry images of\n" "that box are kept from the sphere (the set\n" "a refinement of the box needs); with -ED\n" "the box is ignored."},
	{T_core, "", "-hkl_min_max", "hmin hmax kmin kmax lmin lmax", "\n" "Measured index box; alone it generates the\n" "box and its symmetry images instead of -hkl."},
	{T_core, "", "-hkl_scale", "<n>", "The -hkl file lists h k l in units of 1/n,\n" "whitespace separated, and is taken as it is\n" "(no -dmin, box or twin images): factors at\n" "S = (h,k,l)/n, e.g. the satellites h + m q\n" "of a modulated structure. The table keeps\n" "the scaled integers."},
	{T_core, "", "-IAM", "", "Use Thakkar independent-atom factors."},
	{T_core, "", "-tsc_block", "<n>", "Reflections per block when writing the\n" "tsc [1000]. The table is streamed to disk\n" "a block at a time instead of being held\n" "whole, which is what keeps a protein\n" "inside a small machine. 0 holds it all.\n" "Derived from -mem when that is given\n" "and this is not."},
	{T_core, "", "-acc", "<0..4>", "Numerical grid accuracy [2]; 4 is the\n" "practical maximum."},
	{T_core, "", "-group", "<n ...>", "CIF disorder groups for the asymmetric\n" "unit (prefix a number with + to invert)."},
	{T_core, "", "-charge", "<n> -mult <n>", "Override charge/multiplicity read from\n" "the wavefunction."},
	{T_core, "", "-method", "<RKS|RHF>", "Select DFT or Hartree-Fock treatment."},
	{T_core, "", "-ECP", "[1|2|3]", "Apply ECP correction: 1=def2, 2=xTB,\n" "3=pTB.  Alias: -ecp."},
	{T_core, "", "-ED", "", "Electron-diffraction scattering factors."},
	{T_core, "", "-twin", "<9 values>", "Add a 3x3 twin law; may be repeated."},
	{T_core, "", "-old_tsc", "", "Request the legacy TSC format."},
	{T_core, "", "-tshb", "", "With -ri_fit or -SALTED, write the fitted\n" "atoms as experimental.tshb instead of the\n" "table, one file for all -mtc/-cmtc parts\n" "and the Thakkar fill. The table is kept\n" "for -IAM, grid partitions, -cif, XCW and\n" "fills on tabulated ion data."},
	{T_core, "", "", "", ""},
	{T_partitioning, "", "", "", "PARTITIONING, BASIS, AND PERFORMANCE"},
	{T_partitioning, "", "-b", "<basis> -d <directory>", "Basis name/file and Tonto-style basis\n" "directory (also needed for -fchk)."},
	{T_partitioning, "", "-Becke", "| -TFVC | -mbis | -embis", "Select partitioning (default Hirshfeld)."},
	{T_partitioning, "", "-partition", "<wfn>", "Partition the wavefunction with the scheme\n" "selected above and quit: prints the atomic\n" "charges, dipoles, quadrupoles and molecular\n" "dipole, and writes <stem>.xyzq (label, x y z\n" "in Angstrom, charge in e, dipole in e bohr,\n" "traceless quadrupole xx yy zz xy xz yz in\n" "e bohr^2). -acc sets the grid."},
	{T_partitioning, "", "-ri_fit", "[basis ...]", "RI partitioning; omit a basis or use\n" "auto_aux to generate one automatically, or\n" "auto_aux <elements> for those elements only\n" "(the rest from the basis sets that follow).\n" "With -rho, -esp, -lap or -eli the cubes\n" "come from the fitted density (ELI as\n" "the orbital-free PC07 estimate). Given\n" "before -eli_analysis or -qtaim_eli the\n" "QTAIM basins follow the fitted density."},
	{T_partitioning, "", "-multipole_moments", "<scheme> [N]", "Restrain RI atomic moments through order N\n" "(default 0). Hirshfeld, TFVC, MBIS and EMBIS\n" "support 0-8; Nuclear, Mulliken and Sanderson\n" "use atom-centred charge-only order 0.\n" "Charge and dipole rows use atom-centred\n" "coefficients; higher rows use grid weights.\n" "Does not fix the total electrons exactly.\n" "Implies -ri_fit."},
	{T_partitioning, "", "-multipole_strength", "<x>", "Weight of the restraint rows against the\n" "density-fit metric, default 1."},
	{T_partitioning, "", "-cpus", "<n>", "Maximum worker threads [all available]."},
	{T_partitioning, "", "-mem", "<MB>", "Memory budget for everything sliceable\n" "[unset]. When given, the tsc block size\n" "and XCW's stored integrals are chosen\n" "to fit it, and anything that fits whole\n" "is held whole, that being the fastest\n" "arrangement. An explicit -tsc_block\n" "overrides it."},
	{T_partitioning, "", "-pbc", "<n>", "Periodic-boundary setting."},
	{T_partitioning, "", "", "", ""},
	{T_gpu, "", "", "", "GPU ACCELERATION\n" "  A CUDA or ROCm device is found and used on its own; nothing needs to be\n" "  installed alongside this binary and no BLAS library is required. Every path\n" "  falls back to the CPU when no usable device is present, so these flags are\n" "  safe to leave on a machine that has none. If a GPU is present but this build\n" "  carries no code for it, that is reported rather than silently ignored."},
	{T_gpu, "", "", "", ""},
	{T_gpu, "", "", "", "  On by default:\n" "  (scattering factors)               The Fourier transform runs on the device\n" "                                    whenever one is present and the problem\n" "                                    fits. Roughly 5-30x the CPU loop."},
	{T_gpu, "", "-gpu_itensor", "/ -no_gpu_itensor", "XCW's grid block products and its\n" "stored-integral Fock build. Single\n" "precision by default, which moves the total\n" "energy in the tenth significant figure."},
	{T_gpu, "", "-cpu_itensor_fp32", "/ -no_cpu_itensor_fp32", "\n" "XCW's grid block products on the CPU in\n" "double; single precision, as the device\n" "path runs, is the default."},
	{T_gpu, "", "-gpu_salted", "/ -no_gpu_salted", "The SALTED descriptor combination\n" "(equicomb). It falls back to the CPU when\n" "no usable device is present."},
	{T_gpu, "", "-gpu_grid", "/ -no_gpu_grid", "Atomic integration grid weights (Becke\n" "and TFVC). It falls back to the CPU when\n" "the grid is incompatible or does not fit\n" "on the device."},
	{T_gpu, "", "-gpu_density", "/ -no_gpu_density", "The fitted density on the Gordon-Kim\n" "repulsion grid of -interaction_energy and\n" "the spherical-atom grids (Hirshfeld\n" "surfaces, -hdef, -hirsh)."},
	{T_gpu, "", "", "", ""},
	{T_gpu, "", "", "", "  Off unless asked:"},
	{T_gpu, "", "-gpu_blas", "", "Offer large dense matrix products to the\n" "device. Small ones stay on the CPU, where\n" "they beat what the transfers would allow."},
	{T_gpu, "", "-gpu_cublas", "/ -no_gpu_cublas", "Keep the XCW grid GEMMs off cuBLAS. By\n" "default cuBLAS is used when the machine has\n" "it, being the fastest of the three, and the\n" "built-in kernels otherwise - CUTLASS in\n" "single precision, our own in double.\n" "Nothing is shipped or linked either way:\n" "cuBLAS is opened by name and ignored when\n" "absent. They differ in the last digits, so\n" "a run says which one it used; pin this when\n" "comparing logs across machines."},
	{T_gpu, "", "", "", ""},
	{T_gpu, "", "", "", "  Precision:"},
	{T_gpu, "", "-gpu_fp64", "", "Keep the device in double throughout the\n" "sincos of the transform. Costs almost\n" "nothing on a card with\n" "real double-precision units and a great deal\n" "on one without, which is why it is asked for\n" "rather than detected."},
	{T_gpu, "", "-gpu_fp32", "", "Force single precision even on a card whose\n" "fp64 is fast enough to keep. Mainly for the\n" "tests, which otherwise cannot pin which\n" "kernel a run exercises.\n" "Left alone, the transform picks from the\n" "card's reported fp32:fp64 ratio: double where\n" "that is 2, single where it is 32 or 64. The\n" "single-precision path keeps the phase in\n" "double and agrees with the CPU to ~2e-8."},
	{T_gpu, "", "", "", ""},
	{T_gpu, "", "", "", "  Turning it off and measuring it:"},
	{T_gpu, "", "-no_gpu", "", "Keep everything on the CPU, including the\n" "paths that are on by default. This is what\n" "to use for a CPU baseline."},
	{T_gpu, "", "-gflops", "", "Report achieved GFLOP/s per stage, CPU and\n" "GPU, at the end of the run, with the slowest\n" "single call per stage. Use it to re-derive\n" "the offload thresholds on this machine rather\n" "than inheriting another one's."},
	{T_gpu, "", "", "", ""},
	{T_tsctscb, "", "", "", "TSC/TSCB TABLE UTILITIES"},
	{T_tsctscb, "", "-tscb", "<table.tsc|table.tscb>", "Convert between text .tsc and binary\n" ".tscb, preserving the basename."},
	{T_tsctscb, "", "-tsc_labels", "[<table> <model.cif> [output.tsc]]", "\n" "With table+cif: resolve SCATTERER_IDS\n" "against CIF atoms and write labels .tsc\n" "(default: <table>.labels.tsc). Without\n" "operands: force label-based scatterers in\n" "binary .tscb output (default behavior)."},
	{T_tsctscb, "", "-merge", "<table ...>", "Merge .tsc/.tscb tables with checks, or\n" ".tshb files into combined.tshb (first\n" "file wins for an atom in several);\n" "-cif <cif ...> after the .tshb files\n" "names the run's atoms, and only then can\n" "the result be marked complete."},
	{T_tsctscb, "", "-merge_nocheck", "<table ...>", "Merge tables only when you know their\n" "HKL indices are identical."},
	{T_tsctscb, "", "", "", ""},
	{T_multifragment, "", "", "", "MULTI-FRAGMENT WORKFLOWS"},
	{T_multifragment, "", "-mtc", "<wfn groups ...>", "Multi-TSC calculation: alternating\n" "wavefunction and comma-separated groups.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_multifragment, "", "-cmtc", "<wfn cif groups ...>", "CIF-based multi-TSC: alternating\n" "wavefunction, CIF, and groups.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_multifragment, "", "-mtc_mult", "<n ...>", "Per-input multiplicities for -mtc/-cmtc."},
	{T_multifragment, "", "-mtc_charge", "<n ...>", "Per-input charges (negative: n<value>)."},
	{T_multifragment, "", "-mtc_ECP", "<mode ...>", "Per-input ECP modes."},
	{T_multifragment, "", "-Cation", "<label ...> -Anion <label ...>", "\n" "Explicit ionic fragment labels."},
	{T_multifragment, "", "", "", ""},
	{T_density, "", "", "", "DENSITY, PROPERTIES, AND ANALYSIS\n" "  -rho_cube, the cubes of -rho -elf -eli -lap -rdg -esp -eli_spin -s_rho -def\n" "  -HDEF -fukui, and -fukui_analysis, -eqc, -rgbi, -npa, -npa_summary, -ibo\n" "  take the wavefunction either way: <flag> <wfn>, or -wfn <wfn> <flag>,\n" "  -wfn before or after the flag. A file after a cube flag, -npa or -ibo\n" "  that is not the one -wfn names is refused."},
	{T_density, "", "-rho_cube", "<wfn>", "Write an electron-density cube.  Uses\n" "-radius/-resolution when they precede it\n" "on the command line, and reports the\n" "integrated electron count.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_density, "", "-elf", "-eli -lap -rdg -esp", "Request ELF, ELI, Laplacian, RDG, or\n" "electrostatic-potential properties."},
	{T_density, "", "-eli_spin", "", "Spin-resolved ELI-D cubes: <wfn>_eli_aa,\n" "_eli_bb and _eli_triplet (aa only for a\n" "closed shell), from one pass."},
	{T_density, "", "-esp_isosurface", "[rho]", "Electron-density isosurface [0.002 au]\n" "coloured by the ESP, written as\n" "<wfn>_rho_esp.obj (red negative, blue\n" "positive). Radius is raised to 2.5 A.\n" "-rho, -esp, -lap, -eli and\n" "-esp_isosurface also work from -xyz\n" "<file> with -SALTED <model>: density,\n" "analytic potential, Laplacian and the\n" "PC07 orbital-free ELI come from the ML\n" "prediction, no wavefunction needed."},
	{T_density, "", "-s_rho", "", "Spin density cube."},
	{T_density, "", "-rho_at_points", "<points>", "Density at the points of the file (a\n" "count, then x y z in bohr per line) to\n" "<points>.rho. Needs -wfn before it."},
	{T_density, "", "-cubeb", "", "Write cubes in the binary .cubeb format."},
	{T_density, "", "-cube_convert", "<in> <out>", "Text .cube <-> binary .cubeb; the format\n" "follows the extension of <out>."},
	{T_density, "", "-def", "-HDEF", "Request deformation density or HDEF."},
	{T_density, "", "-MO", "<number|all>", "Generate a molecular-orbital property."},
	{T_density, "fukui", "-fukui", "", "Fukui functions f+/f-/f0 and the dual\n" "descriptor, in the frontier-orbital\n" "approximation. Writes _fukui_plus,\n" "_fukui_minus, _fukui_zero and\n" "_dual_descriptor cubes, plus condensed\n" "per-atom values for all five partitions."},
	{T_density, "fukui", "-fukui_analysis", "<wfn>", "Reactivity analysis of one wavefunction:\n" "frontier orbitals, HOMO-LUMO gap and the\n" "condensed Fukui functions under Hirshfeld,\n" "Becke, TFVC, MBIS and EMBIS. No cubes, so\n" "no grid, radius or CIF is needed.\n" "Built for SCF orbitals and energies; for\n" "natural orbitals see wiki\n" "Correlated-Wavefunctions."},
	{T_density, "eqc", "-eqc", "<gbw> -eqc_frag <atoms> <charge> <mult> <atoms> <charge> <mult> ...", "\n" "Rahm-Hoffmann EQC energy decomposition of\n" "the bonds between the fragments: Delta E\n" "into Delta(nX) and Delta(Vnn-Eee), Q and\n" "the covalency index. occ re-converges the\n" "parent from the gbw, then each fragment\n" "from the parent's density. Atoms are\n" "0-based, e.g. 0,2-4."},
	{T_density, "eqc", "-eqc_method", "<hf|functional>", "occ method, default hf (match the gbw)."},
	{T_density, "eqc", "-eqc_basis", "<name>", "Basis by name instead of from the gbw;\n" "needed when the gbw uses ECPs."},
	{T_density, "eqc", "-eqc_cold", "", "Also converge from occ's own guess."},
	{T_density, "eqc", "-eqc_wfn", "<parent> <frag> <frag> ...", "\n" "EQC from given wavefunctions, no SCF:\n" "wfx or fchk; wfn, gbw or molden with the\n" "ORCA .out of the same name beside them.\n" "Built for SCF orbitals and energies; for\n" "natural orbitals see wiki\n" "Correlated-Wavefunctions."},
	{T_density, "eqc", "", "", "    Writes <stem>.eqc_log. Example, ethane into two CH3 radicals:\n" "      -eqc ethane.gbw -eqc_frag 0,2-4 0 2 1,5-7 0 2\n" "    Theory, output and caveats: docs/wiki/EQC-Energy-Decomposition.md"},
	{T_density, "", "-radius", "<angstrom> -resolution <angstrom>", "\n" "Grid settings for property calculations."},
	{T_density, "", "-hirsh", "<atom-index>", "Hirshfeld analysis for one atom."},
	{T_density, "", "-hirshfeld_surface", "<wfn1> <wfn2>", "Hirshfeld-surface analysis: d_i and d_e\n" "coloured surfaces, shape index and\n" "curvedness in the .dat, fingerprint, and the\n" "ESP-coloured surface when wfn1 has MOs\n" "or, for an xyz, with -SALTED <model>.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_density, "rgbi", "-rgbi", "<wfn>", "Roby-Gould bond-index analysis. Without\n" "<wfn> it reads -wfn. With\n" "-do_XCW it runs on each refined\n" "wavefunction and writes\n" "NA2_<lambda>_RGBI.txt; the .wfn and\n" ".fchk written there cannot carry it."},
	{T_density, "npa", "-npa", "", "Natural atomic orbitals and natural\n" "population analysis (NAO/NPA). Needs a\n" "contracted spherical basis, so a .gbw,\n" ".molden or .fchk, not a .wfn/.wfx."},
	{T_density, "npa", "-npa_summary", "", "-npa without the per-orbital NAO table."},
	{T_density, "ibo", "-ibo", "", "IAO partial charges and intrinsic bond\n" "orbitals (Knizia 2013, MINAO reference),\n" "classified as core, lone pair, bond-like.\n" "Closed shell, .gbw or spherical .molden."},
	{T_density, "ibo", "-ibo_cube", "<list>", "-ibo, then cubes <wfn>_IBO_<n>.cube of the\n" "IBOs in <list>: numbers as -ibo prints\n" "them, all, or C2:C3 for every bond-like\n" "IBO between two atoms, comma-separated\n" "(e.g. 6,C2:C3). Grid as for -MO."},
	{T_density, "ibo", "-nbo_cube", "<list>", "In-house NBO of -wfn, then cubes\n" "<wfn>_NBO_<n>.cube of the NBOs in <list>:\n" "numbers as the NBO table in\n" "NoSpherA2_cube.log prints them, lewis, or\n" "labels with ':' for the blank (LP:Pd1,\n" "BD*:Pd1-P6, LV:Zn2), comma-separated.\n" "Needs .gbw/.molden with all MOs. Open\n" "shell: the table numbers alpha then beta\n" "in one sequence, so a number picks its\n" "own spin and a label takes both; the\n" "cube's comment names the spin."},
	{T_density, "ibo", "-nlmo_cube", "<list>", "The same for the NLMOs of Lewis NBOs:\n" "<wfn>_NLMO_<n>.cube; the log gives each\n" "one's parent-NBO % and atom weights."},
	{T_density, "ibo", "-nbo_lewis", "<file.json>", "NBOs of this Lewis structure for\n" "-nbo_cube/-nlmo_cube (as -nrt_parent; a\n" ".native.nbo.json gives its NRT leader)."},
	{T_density, "rgbi", "-rgbi_no_sym", "", "RGBI without atomic O_h symmetrization. A\n" "free-atom ANO reference stays spherical."},
	{T_density, "rgbi", "-rgbi_basis", "<nao|ano>", "RGBI basis: occupied NAO or ANO [ano]."},
	{T_density, "rgbi", "-rgbi_groups", "<range ...>", "RGBI groups, e.g. 0-5,7; repeat option\n" "for multiple group sets."},
	{T_density, "nci", "-promol_nci", "<a.xyz> <b.xyz> [c.xyz ...] [rcut1 rcut2 rho_max rdg_max colour_max]", "\n" "Promolecular NCI/RDG outputs. Defaults:\n" "rcut1=0.95 and rcut2=0.75."},
	{T_density, "nci", "-promol_nci_single_thread", "", "Disable NCI parallel processing."},
	{T_density, "nci", "", "", ""},
	{T_bonding, "", "", "", "BONDING ANALYSIS\n" "  -topology, -laplacian_bonds, -laplacian_bond_paths, -fba, -eli_analysis,\n" "  -qtaim, -eli_family, -convert_to_47, -nbo_external and -nbo take the\n" "  wavefunction either way: <flag> <wfn> [values], or -wfn <wfn> <flag>\n" "  [values] with -wfn before or after the flag (-wfn w.gbw -eli_analysis\n" "  0.1 3). -qtaim_eli and -nbo_parse take their files only after the flag."},
	{T_bonding, "topology", "-topology", "<wfn>", "All critical points of rho (nuclear,\n" "bond, ring, cage) by Newton-Raphson on\n" "the analytic Hessian, with the\n" "Poincare-Hopf completeness check and any\n" "non-nuclear attractors.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_bonding, "topology", "-laplacian_bonds", "<wfn>", "Plot density Laplacian along bonds.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_bonding, "topology", "-laplacian_bond_paths", "<wfn>", "rho and its Laplacian along the QTAIM\n" "bond paths, traced uphill from each bond\n" "critical point to both nuclei; prints\n" "bond path (BPL) and geometric (GPL)\n" "lengths."},
	{T_bonding, "fba", "-fba", "<wfn>", "Full bonding analysis in one run: RGBI,\n" "the Laplacian along every bond\n" "(-laplacian_bonds), the QTAIM and ELI-D\n" "basins of -eli_analysis, and beside them\n" "on a thread of its own native NBO with\n" "NPA, E2 and NRT (JSON as -nbo\n" "writes it, log block printed last), and\n" "the -eli_contributions table.\n" "Natural orbitals: Mueller DI, 1-RDM ELI-D,\n" "E2(D); wiki Correlated-Wavefunctions."},
	{T_bonding, "qtaim", "-eli_analysis", "<wfn> [resolution radius]", "\n" "QTAIM and ELI-D basins: critical points\n" "as -topology finds them, ELI-D maxima by\n" "gradient ascent, each basin's electrons\n" "on the atomic quadrature grids, the\n" "boundary followed along the analytic\n" "field to rho = 0.0001 e/bohr^3; no cube.\n" "-acc 4 before it tightens the grid sum\n" "from 0.005 to 0.002 e. With -ECP\n" "the core an ECP removed is filled from\n" "Thakkar densities for the QTAIM basins.\n" "Resolution and radius (Angstrom, default\n" "0.05 and 3) set the cube that\n" "-basin_cube and a fitted density still\n" "use, 0.05 A or finer for a real run;\n" "-ri_fit <aux basis> before it (or\n" "-SALTED <model> with an xyz as <wfn>)\n" "takes rho and its gradient from the\n" "fitted density, one pass over the\n" "auxiliary functions instead of the\n" "orbital sum, so the QTAIM basins are\n" "those of the fitted density. ELI-D\n" "stays on the orbitals (a density-only\n" "ELI-D is flat over bonds and lone\n" "pairs) and is skipped without them.\n" "The atomic overlap matrices fall out of\n" "the same point loop, so the QTAIM set\n" "also reports localisation and\n" "delocalization indices; they are left\n" "out, with a note saying so, when the\n" "matrices would need over 512 MB.\n" "Natural orbitals: Mueller DI, 1-RDM ELI-D;\n" "wiki Correlated-Wavefunctions."},
	{T_bonding, "qtaim", "-qtaim", "<wfn> [resolution radius]", "-eli_analysis without the ELI-D basins:\n" "critical points, QTAIM basins and the\n" "delocalization indices; <stem>.qtaim_log."},
	{T_bonding, "qtaim", "-basin_grid", "<n>", "Pull the basin quadrature into the core:\n" "tightest exponent sharpened n^2-fold, the\n" "radial step divided by n, the angular\n" "order up n-1 steps. For heavy atoms."},
	{T_bonding, "qtaim", "-basin_cube", "", "Find the QTAIM and ELI-D basins on the\n" "cube, as before. The default takes the\n" "density's attractors from the analytic\n" "critical-point search, so the grid can no\n" "longer invent one, and walks the analytic\n" "field from every quadrature point for its\n" "basin; use this for the older behaviour."},
	{T_bonding, "qtaim", "-basin_analytic", "", "Accepted and ignored: the analytic ELI-D\n" "boundaries are the default at every\n" "electron count now."},
	{T_bonding, "qtaim", "-basin_step", "<f>", "Scale the trajectory step of the streaming\n" "basin quadrature by f. The step is a fraction\n" "of the cube's voxel; f says how much of it the\n" "basins need. Same basins, f times fewer steps."},
	{T_bonding, "qtaim", "-basin_timing", "", "Print the wall clock of every stage of the\n" "basin analysis: the cube, the critical\n" "points, the beta spheres, the quadrature."},
	{T_bonding, "qtaim", "-basin_gpu", "/ -no_basin_gpu", "Always / never evaluate the streaming QTAIM\n" "and ELI-D trajectories' field and the ELI-D\n" "maxima climbs on the GPU. Default: the GPU\n" "when one is found, a host evaluation costs\n" "at least 15 us (-tune NOS_BASIN_GPU_MIN_US)\n" "and the device outruns all host threads;\n" "the log says which ran and why."},
	{T_bonding, "qtaim", "-adaptive_step", "", "Let a straight run of field earn a longer\n" "ascent step than the validated one: sucrose\n" "363 -> 296 s, UH6 26 -> 13 s, at up to 8e-4\n" "electrons per basin. Off by default."},
	{T_bonding, "qtaim", "-no_adaptive_step", "", "Every ascent step at the validated step\n" "length. The default, and the reference the\n" "equivalence check compares against."},
	{T_bonding, "qtaim", "-beta_spheres", "", "The default; undoes -no_beta_spheres."},
	{T_bonding, "qtaim", "-no_beta_spheres", "", "Climb every quadrature point the whole way\n" "instead of stopping it at the radius around\n" "an attractor that no ascent path can leave.\n" "Same basins, many times the time; it is\n" "there to check that they are the same."},
	{T_bonding, "eli", "-eli_only", "", "With -eli_analysis: skip the density critical\n" "points and QTAIM basins, ELI-D only. With\n" "-eli_basin_cube: only that cube, no basin integration."},
	{T_bonding, "eli", "-eli_contributions", "", "With -eli_analysis: per ELI-D basin, where its electrons\n" "come from: QTAIM atoms (Raub-Jansen), Pipek-Mezey IBOs,\n" "Boys IBOs and NLMOs. Closed shell, .gbw/.molden.\n" "Always on in -fba."},
	{T_bonding, "eli", "-eli_basin_cube", "<A> <B>", "With -eli_analysis: a cube of ELI-D inside the\n" "A-B bond basins (element or label, e.g. Pt Zn),\n" "zero outside, at the -eli_analysis resolution."},
	{T_bonding, "eli", "-eli_mask_cube", "<region.cube>", "With -eli_analysis: only ELI-D on the grid of that\n" "cube where it is nonzero, zero elsewhere (no basins)."},
	{T_bonding, "eli", "-no_spin_eli", "", "Skip the spin-resolved ELI-D basins\n" "(alpha-alpha, beta-beta, triplet) that\n" "-fba and -eli_analysis add for an\n" "unrestricted, spin-polarised wavefunction."},
	{T_bonding, "eli", "-eli_family", "<wfn> [points_file]", "\n" "The rest of Kohout's ELI family as point\n" "values: ELI-D for alpha-alpha, beta-beta and\n" "triplet-coupled pairs and ELI-q, in DGrid 5.2\n" "conventions. Without a points file only the\n" "status of each member is printed; with one\n" "('x y z' per line, bohr) every field at every\n" "point. ELIA needs a correlated 2-matrix and is\n" "not available - the report says why. With\n" "-wfn the points file must not carry a\n" "wavefunction extension (.xyz, .gbw, ...).\n" "Natural orbitals: ELI-D/ELI-q from the\n" "1-RDM; wiki Correlated-Wavefunctions."},
	{T_bonding, "qtaim", "-qtaim_eli", "<rho.cube> <eli.cube> <atoms> [background]", ""},
	{T_bonding, "qtaim", "-qtaim_eli", "<wfn> <atoms> [resolution radius background]", "\n" "Keep ELI only in QTAIM basins of 0-based,\n" "comma-separated atom indices. Runs at\n" "once, so -ri_fit or -SALTED must come\n" "before it (see -eli_analysis).\n" "Natural orbitals: ELI-D/ELI-q from the\n" "1-RDM; wiki Correlated-Wavefunctions."},
	{T_bonding, "nbo", "-convert_to_47", "<wfn>", "Write an NBO File47 (.47). Add\n" "-nbo_keywords \"NRT NRTE2=5\" to fill the\n" "$NBO keylist.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_bonding, "nbo", "-nbo_external", "<wfn>", "Write the .47, run the external NBO on it\n" "and parse the output into <stem>.nbo.json:\n" "NPA, NAO, NBOs with hybridisations, the E2\n" "donor-acceptor table and NRT weights and\n" "bond orders. Options: -nbo_keywords,\n" "-nbo_exe <gennbo>, -nbo_json <file>,\n" "-nbo_dir <dir>. Without -nbo_exe it calls\n" "~/nbo7/gennbo through WSL on Windows and\n" "gennbo from PATH elsewhere; -nbo_exe\n" "wsl:<path> names a binary inside WSL.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_bonding, "nbo", "-nbo_parse", "<out.nbo>", "Parse an NBO output that already exists into\n" "<stem>.nbo.json. Options: -nbo_json <file>."},
	{T_bonding, "nbo", "-nbo", "<wfn>", "The same analysis without the external\n" "program: NAOs and NPA, the NBO search with\n" "hybridisations, E2, and with -nrt natural\n" "resonance theory. Writes\n" "<stem>.native.nbo.json in the format\n" "-nbo_external writes, so\n" "tests/nbo_reference/compare_nbo.py\n" "measures the two against each other.\n" "Options: -nbo_json <file>, -nbo_e2min <kcal>\n" "(default 0.5, what enters the E2 table),\n" "-nbo_threads <n> (0: the OpenMP default),\n" "-nbo_47 <file> to read an existing FILE47\n" "instead of building one from the\n" "wavefunction, -nbo_keep47 to keep the one it\n" "builds, -nbo_ref <scf> for the SCF whose Fock\n" "prices natural-orbital input (E2(D); default\n" "<base>.gbw beside it; wiki\n" "Correlated-Wavefunctions)."},
	{T_bonding, "nrt", "-nrt", "", "Natural resonance theory after the search\n" "(with -nbo): the candidate resonance\n" "structures, the convex quadratic program\n" "that weights them, and the bond orders,\n" "valencies and topology matrices those\n" "weights imply, covalent and ionic apart."},
	{T_bonding, "nrt", "-nrt_e2", "<kcal>", "Donor-acceptor gate on what enters the\n" "resonance search, NBO's NRTE2; default 2.0,\n" "which is between ethane's 2.80 kcal\n" "hyperconjugation and water's strongest\n" "1.01 kcal, so both get the structure count\n" "NBO gives them."},
	{T_bonding, "nrt", "-nrt_arrows", "<n>", "Depth of the arrow-driven candidate\n" "generation, default 2."},
	{T_bonding, "nrt", "-nrt_exhaustive", "", "Enumerate every feasible topology instead\n" "of following arrows. For small systems, and\n" "for checking the arrow search missed none."},
	{T_bonding, "nrt", "-nrt_bond_scale", "<f>", "Longest bond a candidate may place,\n" "f * (r_cov,A + r_cov,B); default 1.75,\n" "looser than the search's own 1.3 because a\n" "resonance structure may bond atoms the\n" "parent does not - ozone's ring structure at\n" "23.6 % spans 2.24 A against 1.32 A."},
	{T_bonding, "nrt", "-nrt_max", "<n>", "Candidate limit, default 4000 capped by a size guard on big molecules."},
	{T_bonding, "nrt", "-nrt_target", "<D(w)>", "Keep searching until the fit error D(w)\n" "reaches this value: arrows are walked again\n" "from every structure the fit weights, then\n" "refitted. Stops early when no new structure\n" "takes weight or the -nrt_max pool is full;\n" "0 = until then. Implies -nrt."},
	{T_bonding, "nrt", "-nrt_rotate", "", "A structure that is a symmetry image of one\n" "whose orbitals converged starts from them,\n" "rotated, instead of the greedy fill: SF6\n" "needs half the sweeps; bond orders move by\n" "up to 2e-4. Default off."},
	{T_bonding, "nrt", "-nrt_screen", "<r|off>", "With -nrt_target or -nrt_ref, a walk of 200+ new\n" "structures builds in full only those whose\n" "reduced cost, from approximate orbitals, is\n" "below r, 16 at a time with their symmetry\n" "images, refitting in between. Default off:\n" "it can miss structures (SF6 at 0.2: D(w)\n" "0.23 % higher); about 2.5 times faster."},
	{T_bonding, "nrt", "-nrt_mirror", "", "The default since 2026-10-14; accepted for\n" "older command lines."},
	{T_bonding, "nrt", "-nrt_no_mirror", "", "Treat a single move that adds no formal-\n" "charge separation (the allyl cation's bond\n" "shift) as half an arrow, complete symmetry\n" "images against every walked state and walk\n" "symmetry references only through new ones:\n" "the behaviour before 2026-10-14, which can\n" "give mirror bonds unequal orders."},
	{T_bonding, "nrt", "-nrt_stall", "<n>", "A structure's orbital sweeps stop once\n" "Anderson mixing has not lowered the largest\n" "vector change by 0.1 % for n sweeps, instead\n" "of running to the cap of 200. Default 0 =\n" "off."},
	{T_bonding, "nrt", "-nrt_atoms", "<a,b,...>", "Subspace NRT: only these 1-based atoms may\n" "change, the rest keep the parent structure's\n" "bonding. NBO's 'NRT <atoms>'."},
	{T_bonding, "nrt", "-nrt_no_symmetry", "", "Keep candidates related by an automorphism\n" "of the molecule instead of collapsing them\n" "onto one weight."},
	{T_bonding, "nrt", "-nrt_no_components", "", "Weight the molecule as one problem instead\n" "of splitting it into connected components of\n" "the delocalisation graph."},
	{T_bonding, "nrt", "-nrt_no_ion", "", "Forbid candidates that move a bond onto a\n" "lone pair (NBO's NRTION off)."},
	{T_bonding, "nrt", "-nrt_parent", "<file.json>", "Start NRT from this Lewis structure instead\n" "of the search's own (NBO's $CHOOSE; implies\n" "-nrt): the NBOs and occupancies are those of\n" "the given structure, the resonance search\n" "runs on top. Closed shell only. The file:\n" "{\"format\": \"NoSpherA2 Lewis structure\",\n" " \"version\": 1, \"atoms\": 3,\n" " \"lone_pairs\": {\"O1\": 2},\n" " \"bonds\": [[\"O1\", \"H2\", 1], [\"O1\", 3, 1]]}\n" "Atoms are element + 1-based index or a bare\n" "index; valence lone pairs only, the cores\n" "are added; bond orders 1-3. A\n" ".native.nbo.json is read through its\n" "\"parent_lewis\". The log compares its rho(NL)\n" "with the search's own best structure.\n" "\"multi_centre\": [[\"B1\", \"H3\", \"B2\"]] adds\n" "one 3- or 4-centre pair per unit."},
	{T_bonding, "nbo", "-nbo_tm_pairs", "<n>", "Lewis pairs a d-block metal may hold (lone\n" "pairs + bonds). 9 (default) = the 18-electron\n" "space, so M-L back-donation shows; 0 = its\n" "valence NAO count, 6, NBO 7's footing. Never lifted\n" "silently (the log says when it must be)."},
	{T_bonding, "nbo", "-nbo_centres", "<k>", "Also place 3..k-centre bonds (k = 3 or 4;\n" "NBO's 3CBOND) over chains of bondable atoms,\n" "when one- and two-centre blocks fail. Kept if\n" "each centre holds 10 % of the orbital and\n" "rho(NL) drops 0.05 e per unit against the\n" "two-centre search (logged); NRT then moves\n" "pairs between bonds and the units."},
	{T_bonding, "nrt", "-nrt_e2_avg", "<percent>", "E(2) per NRT structure whose weight is within\n" "<percent> of the leading one (degenerate\n" "set), each placed under $CHOOSE, averaged\n" "with the NRT weights renormalised over the\n" "set (implies -nrt; closed shell only).\n" "Channels sum orbitals of one type on the same\n" "centres; a channel a structure lacks counts 0."},
	{T_bonding, "nrt", "-nrt_e2_avg_min", "<percent>", "As -nrt_e2_avg, but over every NRT structure\n" "of at least this weight: reaches structures\n" "that change a metal core when the near-\n" "degenerate ones differ only far away. Check\n" "rho(NL) per structure: a high one gives E(2)\n" "that means little."},
	{T_bonding, "nrt", "-nrt_tie", "<fraction>", "The fit fixes the residual, not the weights\n" "(parent + double change = two singles).\n" "Among weights within this fraction of the\n" "best D(w)^2, take the most parent weight,\n" "then single changes, then doubles, ...\n" "Default 1e-6, 0 = off."},
	{T_bonding, "nrt", "", "", ""},
	{T_bonding, "nrt", "-nrt_ref", "<fraction>", "Multi-reference NRT: every structure of at\n" "least this fraction of the leading weight\n" "walks its own arrows into the pool and the\n" "fit repeats (at most 8 cycles); -nrt_tie\n" "ranks by moves to the nearest reference.\n" "The structure nearest the resonance hybrid\n" "re-anchors the walk each cycle; the leading\n" "structure is reported as the parent.\n" "Default 0.5, 0 = off."},
	{T_bonding, "nrt", "", "", ""},
	{T_conversion, "", "", "", "CONVERSION, ML, AND SPECIALISED TOOLS"},
	{T_conversion, "", "-gbw2wfn", "-wfn <file.gbw>", "Convert GBW input to .wfn.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_conversion, "", "-fchk", "<output.fchk>", "Write FCHK output (requires -b and -d)."},
	{T_conversion, "", "-SALTED", "<model> [<model> ...]", "Predict density with a SALTED model.\n" "A model is a directory or a .salted file.\n" "With several, each element is predicted by\n" "the first model trained on it and the\n" "per-atom blocks are stitched together.\n" "With -xyz instead of -wfn the prediction\n" "is the structure's density for -rho,\n" "-esp, -lap, -eli, -esp_isosurface and the Hirshfeld\n" "surface ESP. "},
	{T_conversion, "", "-salted_fold", "<in.salted> <out.salted>", "\n" "Store the model with its weights and\n" "features folded together: smaller and\n" "faster to predict, same prediction to\n" "rounding. Older builds refuse the file."},
	{T_conversion, "", "-salted_ecp", "<in.salted> <out.salted> <mode>", "\n" "Mark a model trained on ECP densities\n" "(mode as -ECP: 1 def2, 2 xTB, 3 pTB);\n" "its atoms get the core back at prediction."},
	{T_conversion, "", "-SALTED_distance", "<model> <xyz> [eigcut]", "\n" "How far each atom lies outside the model's\n" "training span: lambda = 0 kernel novelty\n" "(0 inside, 1 outside), largest kernel and\n" "cosine to a reference environment; eigcut\n" "default 1e-10. For integer zeta also the\n" "share of each neighbour pair channel (and\n" "element) in nov0, ranked by how far it\n" "exceeds the channel's share of the atom's\n" "descriptor, the neighbours the references\n" "have and the atom lacks, the channels the\n" "model's references lack among its kept\n" "features. Also the atoms' channels the\n" "model keeps no features for. Writes\n" "<xyz stem>_salted_distance.csv and .json in\n" "the working directory; -no_gpu_salted, if\n" "wanted, must come before it."},
	{T_conversion, "", "-SALTED_COEFS", "<model-dir>", "Write SALTED_COEFS.npy (requires -wfn)."},
	{T_conversion, "", "-RI_CUBE", "<coefficients.npy>", "Write an RI density cube; use -wfn and\n" "-ri_fit first."},
	{T_conversion, "", "-write_ri_coefs", "", "Write RI_COEFS.npy; use -wfn, -ri_fit\n" "and -multipole_moments first."},
	{T_conversion, "", "-interaction_energy", "<A> <B>", "Electrostatic interaction energy of two\n" "fitted densities: two wavefunctions with\n" "-ri_fit, or two .xyz with -SALTED. With\n" "<A> <A.npy> <B> <B.npy> the coefficient\n" "files are read instead; use -ri_fit.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_conversion, "", "-interaction_energies", "<job>", "Interaction energies of every molecule\n" "pair in contact in a crystal; the job file\n" "lists cif, cutoff, output and one\n" "molecule <xyz> [<coef.npy>] per line.\n" "Use -ri_fit or -SALTED first."},
	{T_conversion, "", "-repulsion_overlap", "<K>", "Exchange-repulsion of -interaction_energy\n" "as K * Int rhoA rhoB, K in Eh bohr^3/e^2,\n" "instead of the Gordon-Kim default."},
	{T_conversion, "", "-repulsion_exchange", "<dirac|pbe|b88|r2scan> Exchange functional of the", "\n" "Gordon-Kim repulsion, default dirac;\n" "r2scan is the deorbitalised r2SCAN-L."},
	{T_conversion, "", "-combine_mos", "<wfn1> <wfn2>", "Combine molecular orbitals.\n" "Built for SCF orbitals and energies; for\n" "natural orbitals see wiki\n" "Correlated-Wavefunctions."},
	{T_conversion, "", "-cmos1", "<MO ...> -cmos2 <MO ...>", "MO selections for -combine_mos."},
	{T_conversion, "", "-QCT", "| -qct", "Enter the legacy QCT workflow."},
	{T_conversion, "", "", "", ""},
	{T_diagnostics, "", "", "", "DIAGNOSTICS AND INTERNAL TOOLS"},
	{T_diagnostics, "", "-v", "| -v2 | -debug", "Verbose diagnostic output."},
	{T_diagnostics, "", "-out", "<file>", "Write the log to <file> instead of\n" "NoSpherA2.log."},
	{T_diagnostics, "", "-rkpts", "/ -skpts", "Read / save the k-point set."},
	{T_diagnostics, "", "-profiling", "[tests-root]", "Run the internal profiling suite\n" "[./tests]. Alias: -profile."},
	{T_diagnostics, "", "-tune", "NAME[=VALUE]", "Set a developer knob, repeatable; a bare NAME\n" "is 1. Diagnostics and cut-offs such as\n" "NOSPHERA2_GPU_BACKEND=hip or NOS_RGBI_DEBUG,\n" "formerly environment variables."},
	{T_diagnostics, "", "-tune", "list", "Print every knob with its default and stop."},
	{T_diagnostics, "", "-no_date", "", "Suppress date information and the GPU notes, so\n" "output does not depend on the machine it ran on."},
	{T_diagnostics, "", "-no_date_but_gpu", "", "As -no_date, but keeps the GPU notes. For the\n" "tests whose reference has to show that the device\n" "did the work, a silent fallback being otherwise\n" "indistinguishable from the CPU result."},
	{T_diagnostics, "", "-draw_orbits", "l,m[,resolution,radius]", "\n" "Draw a spherical-harmonic orbital."},
	{T_diagnostics, "", "-ewal_sum", "<cube> [kmax] [accuracy] Ewald sum of a cube.", ""},
	{T_diagnostics, "", "-atom_dens", "<wfn> [alpha-MOs beta-MOs]", ""},
	{T_diagnostics, "", "-atom_dens_diff", "<gbw1> <gbw2>", "Difference density from two GBW files."},
	{T_diagnostics, "", "-spherical_aver_fukui", "<wfn1> <wfn2>", ""},
	{T_diagnostics, "", "-spherical_aver_hirsh", "<wfn>", ""},
	{T_diagnostics, "", "-dipole_moments", "-polarizabilities <7 wfn files>", "\n" "The five above, correlated input: see -wfn\n" "and wiki Correlated-Wavefunctions."},
	{T_diagnostics, "", "-geometry_aid_cutoff", "<r>", "SOAP cutoff for the geometry-aid descriptor:\n" "  3.5 (default) matches the c_only models,\n" "  3.0 matches the dirty models."},
	{T_diagnostics, "", "-geometry_aid_metals", "", "Add Zn as the metal stand-in species (row-2\n" "  specialist models, 49,686 long)."},
	{T_diagnostics, "", "-geometry_aid_center_weight", "<w>", "SOAP centre-atom weight, 0 masks the centre\n" "  (label-aware models, default 1)."},
	{T_diagnostics, "", "-calc_featomic_descriptor", "", "Write descriptor.npy (requires -wfn)."},
	{T_diagnostics, "", "-calc_featomic_descriptors", "<list>", "Same, for many structures in one run.\n" "  <list> holds one structure path per\n" "  line; each descriptor is written as\n" "  <path>.npy. The radial-basis splines\n" "  cost a fixed 0.7 s and are built once,\n" "  so a batch is far cheaper than N runs."},
	{T_diagnostics, "", "-classify_atoms", "<model.bin> [out]", "Element probabilities per atom from the\n" "  geometry-aid model (requires -wfn).\n" "  Writes probabilities.npy, one row per\n" "  atom, instead of the 42,042-wide\n" "  descriptor. Make the .bin with\n" "  make_geometry_aid_bin.py."},
	{T_diagnostics, "", "-classify_atoms_list", "<list> <model.bin>", "\n" "  Same, for many structures in one run."},
	{T_diagnostics, "", "-SALTED_Training", "", "Create SALTED training data."},
	{T_diagnostics, "", "-all_charges", "", "Print all calculated atomic charges."},
	{T_diagnostics, "", "-atom_sfac", "<wfn1> <wfn2>", "Compare atom-centred scattering factors.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_diagnostics, "", "-density_difference", "<wfn>", "Use a second wavefunction for a\n" "density-difference calculation.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_diagnostics, "", "-e_field", "<value>", "Apply an external electric field."},
	{T_diagnostics, "", "-fractal", "<name>", "Enable the legacy fractal workflow."},
	{T_diagnostics, "", "-get_g", "", "Enable reciprocal-space g calculation."},
	{T_diagnostics, "", "-refine", "[accuracy]", "Set refinement integral accuracy [0.1]."},
	{T_diagnostics, "", "-rgbi_EVs", "", "Include RGBI eigenvectors."},
	{T_diagnostics, "", "-rgbi_theta", "", "Include per-bond Roby-Gould theta-subspace populations and indices."},
	{T_diagnostics, "", "-rgbi_legacy_cutoff", "", "RGBI atomic subspaces from occupation-number thresholds (pre-fix, jumps with geometry)."},
	{T_diagnostics, "", "-sfac_diffuse", "x y z cif wfn dmin", "Calculate diffuse scattering factors.\n" "Correlated (natural-orbital) input: see\n" "-wfn and wiki Correlated-Wavefunctions."},
	{T_diagnostics, "", "", "", ""},
	{T_experimental, "", "", "", "EXPERIMENTAL AND DEVELOPER COMMANDS"},
	{T_experimental, "", "-coef", "<file>", "Use externally supplied SALTED\n" "coefficients."},
	{T_experimental, "", "-convert_XCW", "<stdout> <lambda-step> Convert Tonto XCW lambda-step output.", ""},
	{T_experimental, "", "-do_XCW", "", "XCW mode; anomalous dispersion is read from the CIF."},
	{T_experimental, "", "-mask", "<file.fab>", "Solvent-mask contribution for an XCW\n" "settings file with the mask keyword."},
	{T_experimental, "", "-print_F", "", "Write the F_calc table of every converged\n" "lambda step, NA2_<lambda>_Fcalc.txt."},
	{T_experimental, "", "-xcw_extrapolate", "/ -no_xcw_extrapolate", "\n" "Seed each lambda step from the density\n" "extrapolated through the two previous steps\n" "(the default) / from the last one alone."},
	{T_experimental, "", "-xcw_incremental", "/ -no_xcw_incremental", "\n" "Two-electron Fock build from the change of\n" "the density each iteration, from the stored\n" "integrals when they fit in memory; the\n" "direct build needs -xcw_int_precision 1e-12."},
	{T_experimental, "", "-xcw_int_precision", "<p>", "Integral screening threshold of the XCW\n" "Fock build, default 1e-10 (OCC's own 1e-12)."},
	{T_experimental, "", "-xcw_strong_cutoff", "<z>", "z-score above which a reflection counts as\n" "strong in the XCW diagnostics, default 3;\n" "a settings-file keyword overrides it."},
	{T_experimental, "", "-xcw_weight", "[robust] <scheme> [values] [c]", "\n" "Crystallographic weight of the XCW fit:\n" "sigma (1/sigma^2, default), unit, stl a\n" "((sin th/lambda)^a), shelx a b (WGHT),\n" "shelx_opt a b (WGHT refit every step),\n" "chebyshev A0 .. (1/Sum A_r T_r(t), t from\n" "Fc), chebyshev_opt n (A_r refit every\n" "step); robust adds Tukey's biweight, c\n" "default 6. The settings keyword `weight`."},
	{T_experimental, "", "-xcw_weight_rescale", "0|1", "unit, stl and chebyshev weights times the k\n" "giving S = 1 at the scan's start, held for\n" "the scan (default 1, CRYSTALS REWEIGHT).\n" "The settings keyword `weight_rescale off`."},
	{T_experimental, "", "-xcw_partition", "grid|ri", "How XCW shares the density among the atoms:\n" "grid (default) by the -hirsh, -becke,\n" "-tfvc, -mbis or -embis weights, ri by an RI fit in the\n" "-ri_fit basis (AutoAux without one), its\n" "moments restrained to that grid partition\n" "by -multipole_moments."},
	{T_experimental, "", "-XCW_settings", "<file>", "Keywords for -do_XCW. Besides the\n" "refinement settings:\n" "  save_eri <path> write the stored\n" "                  two-electron integrals\n" "  read_eri <path> reuse them: they depend\n" "                  on structure and basis\n" "                  only; another basis or\n" "                  geometry is refused and\n" "                  they are built\n" "  basis_overrides <file>\n" "                  per-element all-electron\n" "                  BSE JSON orbital bases;\n" "                  other elements keep basis_set\n" "  df_basis <name> density-fit the Fock\n" "                  build with that\n" "                  auxiliary basis from\n" "                  the built-in library,\n" "                  e.g. def2-universal-jkfit\n" "  i_double        the grid's AO values\n" "                  and block products\n" "                  in double precision\n" "  extinction <model> [iso|aniso] [fixed] [start]\n" "                  refine shelx, bc_gaussian\n" "                  or bc_lorentzian extinction;\n" "                  aniso starts from one x or\n" "                  X11 X22 X33 X12 X13 X23\n" "Prefer a GENEROUS budget: holding 12%\n" "of the tensor cost 8% more time than\n" "holding all of it, while shrinking to\n" "0.4% saved a further 7% of memory and\n" "cost another 9% of time."},
	{T_experimental, "", "-partitioning_test", "-NNLS_TEST", "Run internal partitioning/NNLS tests."},
	{T_experimental, "", "-test_RI", "-RI_WFN_DIFF", "RI fitting diagnostics (after -wfn and\n" "-ri_fit)."},
	{T_experimental, "", "-wfn_cif", "", "Write a CIF from the wavefunction."},
	{T_experimental, "", "-lukas_test", "", "SALTED density diagnostic."},
	{T_experimental, "", "-calc_dens_1D", "[atom1 atom2 points padding]", "\n" "One-dimensional density between atoms."},
	{T_experimental, "", "-spherical_harmonic", "-spherical_atoms", "\n" "Internal spherical-atom diagnostics."},
	{T_experimental, "", "-test", "", "Enable test mode."},
	{T_experimental, "", "", "", ""},
	{T_help, "", "", "", "HELP"},
	{T_help, "", "-h", "", "The whole text (also --h, -help, --help)."},
	{T_help, "", "-h", "topics", "The topics and their sub-topics."},
	{T_help, "", "-h", "<word>", "A sub-topic (-h nbo), else the topics whose\n" "title holds <word> (-h gpu)."},
	{T_help, "", "-h", "<-flag>", "The entries of one option, e.g. -h -eqc."},
	{T_help, "", "", "", ""},
	{T_examples, "", "", "", "EXAMPLES\n" "  Normal HAR\n" "    NoSpherA2 -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n" "  Independent-atom (Thakkar) TSC\n" "    NoSpherA2 -cif A.cif -hkl A.hkl -xyz A.xyz -IAM -acc 1 -cpus 7\n" "  Disorder / multi-TSC\n" "    NoSpherA2 -cif A.cif -hkl A.hkl -acc 1 -cpus 7 \\\n" "      -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3 \\\n" "      -mtc_charge 0 0 0 -mtc_mult 1 1 1 -mtc_ECP 0 0 0\n" "  Fragment HAR\n" "    NoSpherA2 -cif A.cif -hkl A.hkl -acc 1 -cpus 7 \\\n" "      -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 \\\n" "      3_2.wfn 3_2.cif 0,2\n" "  Merge tables\n" "    NoSpherA2 -merge A.tsc B.tsc C.tsc\n" "    NoSpherA2 -merge_nocheck A.tsc B.tsc C.tsc  # only identical HKLs\n" "  Convert, label, and GBW conversion\n" "    NoSpherA2 -tscb A.tscb\n" "    NoSpherA2 -tsc_labels A.tscb A.cif A.labels.tsc\n" "    NoSpherA2 -gbw2wfn -wfn A.gbw\n" "  Twin law\n" "    NoSpherA2 -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7 \\\n" "      -twin -1 0 0 0 -1 0 0 0 -1"}
	};
	return items;
}

//an option puts its first text line at column 36 (two spaces after a longer flag column), the rest under it
static std::string render(const help_item& it)
{
	if (it.flag.empty())
		return it.text + "\n";
	const std::string left = "  " + it.flag + (it.args.empty() ? "" : " " + it.args);
	std::string out = left, line;
	std::istringstream in(it.text);
	for (int n = 0; std::getline(in, line); n++)
		if (n > 0)
			out += "\n" + std::string(36, ' ') + line;
		else if (!line.empty())
			out += (left.size() <= 34 ? std::string(36 - left.size(), ' ') : "  ") + line;
	return out + "\n";
}

const std::string& help_text()
{
	static const std::string text = [] {
		const std::vector<help_item>& items = help_items();
		std::string s;
		for (int i = 0; i < items.size(); i++)
			s += render(items[i]);
		return s;
	}();
	return text;
}

//lower case, and a '-' inside a word as '_', so -rgbi-groups finds -rgbi_groups
static std::string norm(std::string s)
{
	for (int i = 0; i < s.size(); i++)
		s[i] = s[i] == '-' && i > 0 && isalnum(static_cast<unsigned char>(s[i - 1])) ? '_' : tolower(static_cast<unsigned char>(s[i]));
	return s;
}

//f as a whole word of s
static bool mentions(const std::string& s, const std::string& f)
{
	auto word = [&](size_t q) { return q < s.size() && (isalnum(static_cast<unsigned char>(s[q])) || s[q] == '_'); };
	for (size_t q = s.find(f); q != std::string::npos; q = s.find(f, q + 1))
		if (!(q > 0 && word(q - 1)) && !word(q + f.size()))
			return true;
	return false;
}

static bool is_header(const help_item& it)
{
	return it.flag.empty() && !it.topic.empty() && it.text.compare(0, it.topic.size(), it.topic) == 0;
}

static std::string help_topics()
{
	const std::vector<help_item>& items = help_items();
	std::string out;
	for (int i = 0; i < items.size(); i++)
	{
		if (!is_header(items[i]))
			continue;
		std::string subs;
		for (int j = i + 1; j < items.size() && items[j].topic == items[i].topic; j++)
			if (!items[j].sub.empty() && !mentions(subs, items[j].sub))
				subs += " " + items[j].sub;
		out += "  " + items[i].topic + (subs.empty() ? "" : ":" + subs) + "\n";
	}
	return out;
}

//-h <word>: the first that matches of a sub-topic, the topics whose title holds the word, the option
//itself (or the line listing it), the options whose text names it; -h alone gives everything
std::string help_lookup(const std::string& want)
{
	const std::vector<help_item>& items = help_items();
	const std::string w = norm(want), f = w[0] == '-' ? w : "-" + w;
	if (w.empty())
		return help_text();
	if (w == "topics")
		return "\n" + help_topics();
	for (int t = 0; t < 4; t++)
	{
		std::string out, topic;
		for (int i = 0; i < items.size(); i++)
		{
			const help_item& it = items[i];
			const bool keep = t == 0 ? w[0] != '-' && it.sub == w
				: t == 1 ? w[0] != '-' && norm(it.topic).find(w) != std::string::npos
				: t == 2 ? !it.flag.empty() && (norm(it.flag) == f || mentions(norm(it.args), f))
				: !it.flag.empty() && mentions(norm(it.text), f);
			if (!keep)
				continue;
			if (it.topic != topic)
				out += "\n" + (is_header(it) ? "" : it.topic + "\n"), topic = it.topic;
			out += render(it);
		}
		if (!out.empty())
			return out;
	}
	return "No help section matches \"" + want + "\"; the topics are:\n" + help_topics();
}
