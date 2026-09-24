"""Copy the benchmark .nbo.json files into tests/nbo_reference/ and write index.json.

    python collect_reference.py <bench_dir> <dest_dir>

Everything in index.json is read back from what the programs actually produced: the level
of theory, charge, multiplicity and convergence come from each ORCA output, the thresholds
from what NBO echoed in its own output (already in the .nbo.json). Nothing is copied from
an input template, so a changed default shows up as drift instead of being absorbed.
"""
import datetime
import json
import os
import re
import shutil
import sys

# molecule -> motif tags; charge, multiplicity and level of theory come from the ORCA output
MOTIFS = {
    "ethane":       ["nonpolar single bond", "C-H"],
    "ethene":       ["double bond", "pi"],
    "acetylene":    ["triple bond", "pi"],
    "n2":           ["triple bond", "lone pair", "homonuclear"],
    "hcn":          ["triple bond", "lone pair", "polar"],
    "water":        ["polar bond", "lone pair rich"],
    "ammonia":      ["polar bond", "lone pair"],
    "formaldehyde": ["polar double bond", "lone pair rich", "pi"],
    "lif":          ["ionic", "very polar"],
    "benzene":      ["aromatic", "resonance", "Kekule NRT"],
    "pyridine":     ["aromatic", "heteroaromatic", "lone pair", "resonance"],
    "ozone":        ["resonance", "delocalised", "hypervalent-looking"],
    "so2":          ["resonance", "hypervalency", "polar"],
    "nitromethane": ["resonance", "delocalised", "polar"],
    "formate":      ["anion", "resonance", "symmetric delocalisation"],
    "sf6":          ["hypervalency", "octahedral"],
    "pf5":          ["hypervalency", "trigonal bipyramidal"],
    "ch3":          ["open shell", "doublet radical"],
    "o2":           ["open shell", "triplet", "pi*"],
    "no":           ["open shell", "doublet", "polar multiple bond"],
    "ticl4":        ["transition metal", "d0", "polar M-Cl"],
    "ni_co_4":      ["transition metal", "d10", "backbonding", "CO ligand"],
}

# The two molecules whose archive deliberately differs from the uniform one. Both are NBO-side
# limits, not wavefunction problems, and a reproduction attempt needs to know about them.
NOTES = {
    "sf6": "NRTSYM=off is required. With NRT's symmetry detection on, NBO aborts the whole "
           "analysis after 0.08 s with 'SYMOPS: generated 48 symmetry operator(s) for Th but "
           "expected 24' and prints no NAO/NBO/E2/NRT at all. With it off the run completes in "
           "21 s CPU. This is the only molecule whose $NBO keylist differs from the uniform one.",
    "so2": "The NRT search cannot be started from NBO's own Lewis structure. NBO localises SO2 "
           "as S-O single plus S=O triple (the $CHOOSE it writes back is 'LONE 1 1 2 3 3 1 / "
           "BOND S 1 2 T 1 3'), and NRT rejects both symmetry-related parents with 'Unable to "
           "find reasonable resonance structure', ending in 'NRTDRV: NRT has no resonance "
           "structures'. Plain NRT fails the same way; NRTDTL only makes the reason visible, and "
           "NRTION, NRTSYM=off and NRTCYC=5 do not help. The archive therefore carries a "
           "hand-written $NRTSTR with the four 9-pair structures (O=S=O, the two O(-)-S(+)=O "
           "forms and the doubly ionic single/single form), appended after the data blocks. With "
           "it the search converges normally: 4 of 10 structures retained, D(w)=0.06896429, "
           "weights 27.58/27.58/22.42/22.42 and symmetric S-O bond orders of 1.5000. NBO checks "
           "the pair count of every supplied structure ('Structure 4 of the $NRTSTR keylist has "
           "1 too many electron pair(s)'), so a wrong $NRTSTR fails loudly rather than quietly.",
    "ni_co_4": "NRTE2=5 is required. At the default delocalisation threshold the search was "
               "still enumerating candidates after 90 minutes on one core; with NRTE2=5 it "
               "converges in 693 s wall / 153 s CPU to 16 of 136 candidates retained, "
               "D(w)=0.06910797. The threshold is not a convergence knob - NRTE2=10 and "
               "NRTE2=20 both truncate to 4 candidates and a worse D(w)=0.07174346 - so this "
               "molecule's NRT numbers may only be compared against a run at NRTE2=5. "
               "NPA, NAO, NBO and E2 are unaffected: they are computed before NRT.",
    "ticl4": "NRTE2=20 is required; at the default threshold the search had 1696 distinct "
             "candidates behind it and was still in its Gram-matrix phase after 90 minutes on "
             "one core. NRTE2=10 was run afterwards as a check and reproduces this reference "
             "exactly where it matters: same 16 of 45 retained, same D(0)=0.08511395 and "
             "D(w)=0.06741150, and all 45 bond orders identical to the printed 4 decimals - but "
             "individual resonance weights differ by up to 7.1 percentage points, and still by "
             "3.1 when matched by rank. That is the sharpest statement in the set of which NRT "
             "quantities survive a changed search: D(w), the bond orders and the valencies do, "
             "the weight of any one structure does not, because a different threshold changes "
             "which structures are in the pool that the weights are distributed over. NRTE2=50 "
             "is a different answer altogether (15 of 62, D(w)=0.08344777).",
}

# key -> (regex, converter). All anchored on what ORCA prints, not on what was passed in.
ORCA_FIELDS = {
    "version":      (r"Program Version (\S+)", str),
    "keyword_line": (r"^\|\s*1>\s*(!.*?)\s*$", str),
    "basis":        (r"Your calculation utilizes the basis:\s*(\S+)", str),
    "exchange":     (r"Exchange Functional\s+Exchange\s+\.+\s*(\S+)", str),
    "correlation":  (r"Correlation Functional\s+Correlation\s+\.+\s*(\S+)", str),
    "hf_type":      (r"Hartree-Fock type\s+HFTyp\s+\.+\s*(\S+)", str),
    "charge":       (r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", int),
    "multiplicity": (r"Multiplicity\s+Mult\s+\.+\s*(\d+)", int),
    "basis_functions": (r"Basis Dimension\s+Dim\s+\.+\s*(\d+)", int),
}


def read_orca(path):
    """First match per field, plus the last final energy and whether the opt converged."""
    out = {}
    pats = {k: re.compile(p, re.M) for k, (p, _) in ORCA_FIELDS.items()}
    energy, converged = None, False
    with open(path, errors="replace") as f:
        for line in f:
            for k, pat in pats.items():
                if k not in out:
                    m = pat.search(line)
                    if m:
                        out[k] = ORCA_FIELDS[k][1](m.group(1))
            if line.startswith("FINAL SINGLE POINT ENERGY"):
                energy = float(line.split()[-1])
            elif "THE OPTIMIZATION HAS CONVERGED" in line:
                converged = True
    out["final_energy_hartree"] = energy
    out["optimisation_converged"] = converged
    return out


def nrtstr_supplied(bench, name):
    """True if the archive carries a hand-written $NRTSTR keylist (so2 does, nothing else)."""
    path = os.path.join(bench, name, name + ".47")
    if not os.path.exists(path):
        return False
    with open(path, errors="replace") as f:
        return "$NRTSTR" in f.read()


def collect_spread(dest):
    """Copy the start-dependence summaries next to the dataset, one file per molecule."""
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(dest, "start_dependence")
    found = {}
    for group, sub in (("keywords", "spread"), ("atom_order", "spread_reorder")):
        root = os.path.join(here, sub)
        if not os.path.isdir(root):
            continue
        for mol in sorted(os.listdir(root)):
            src = os.path.join(root, mol, "spread.json")
            if not os.path.exists(src):
                continue
            os.makedirs(out, exist_ok=True)
            shutil.copyfile(src, os.path.join(out, "%s_%s.json" % (mol, group)))
            #spread_starts.json is the same comparison with NRTCYC=1 left out. That lever caps
            #the number of search cycles, so it truncates the search rather than starting it
            #somewhere else, and including it inflates every range.
            starts = os.path.join(root, mol, "spread_starts.json")
            if os.path.exists(starts):
                shutil.copyfile(starts, os.path.join(out, "%s_%s_starts.json" % (mol, group)))
            found.setdefault(group, []).append(mol)
    return found


def main(argv):
    bench, dest = argv[1], argv[2]
    os.makedirs(dest, exist_ok=True)
    entries, levels = [], set()
    pending = []
    for name in sorted(MOTIFS):
        src = os.path.join(bench, name, name + ".nbo.json")
        out = os.path.join(bench, name, name + ".out")
        if not os.path.exists(src):
            print("MISSING", name)
            pending.append(name)
            continue
        with open(src) as f:
            d = json.load(f)
        #A run still in progress has an NRT section that parses but is not an answer: candidate
        #structures printed, no converged weights or bond orders. Half an NRT section in a
        #reference is worse than none, so such a molecule is listed as pending instead.
        if d["nrt"].get("present") and not d["nrt"].get("bond_orders"):
            print("INCOMPLETE NRT, not collected:", name)
            pending.append(name)
            continue
        orca = read_orca(out) if os.path.exists(out) else {}
        if not orca.get("optimisation_converged", True):
            print("WARNING: optimisation did not converge:", name)
        levels.add("%s / %s" % (orca.get("keyword_line", "?"), orca.get("basis", "?")))
        shutil.copyfile(src, os.path.join(dest, name + ".nbo.json"))
        entries.append({
            "name": name,
            "file": name + ".nbo.json",
            "motifs": MOTIFS[name],
            "open_shell": d["open_shell"],
            "orca": orca,
            "atoms": len(d["npa"]),
            "naos": len(d["nao"]),
            "nbos": len(d["nbos"]),
            "e2_entries": len(d["e2"]),
            "nrt_structures_used": d["nrt"].get("structures_used", 0),
            #The problem size the NRT search actually solved: what it examined, what it kept,
            #how many QP steps that took. This is the cost a screening scheme has to beat.
            "nrt_structures_found": d["nrt"].get("structures_found", 0),
            "nrt_candidates_printed": len(d["nrt"].get("candidates", [])),
            "nrt_weights_listed": len(d["nrt"].get("weights", [])),
            "nrt_zero_weight_structures": sum(1 for w in d["nrt"].get("weights", [])
                                              if w.get("weight_percent", 0) == 0.0),
            "nrt_search_cycles": len(d["nrt"].get("cycles", [])),
            "nrt_qp_iterations": len(d["nrt"].get("qp_iterations", [])),
            "nrt_d_w": d["nrt"].get("d_w", 0.0),
            "nrt_d_0": d["nrt"].get("d_0", 0.0),
            "nrt_symmetry": d["nrt"].get("symmetry", ""),
            "nrt_bond_orders": len(d["nrt"].get("bond_orders", [])),
            "nbo_keywords": d["keywords"],
            #What NBO echoed back as recognised, which is the only keyword evidence a json
            #re-derived with -nbo_parse carries (nbo_keywords is empty for those).
            "nbo_keywords_reported": d.get("keywords_reported", ""),
            "note": NOTES.get(name, ""),
            "nrt_nrtstr_supplied": nrtstr_supplied(bench, name),
            "nbo_version": d["nbo_version"],
            "timings_seconds": d["timings"],
            "thresholds_reported_by_nbo": d["thresholds"],
        })
    index = {
        "description": "NBO 7 reference results for the in-house NBO implementation in "
                       "NoSpherA2. One <name>.nbo.json per molecule, written by "
                       "`NoSpherA2 -nbo <name>.gbw -nbo_keywords \"NRT E2PERT\"`; the "
                       "structures behind the keys are in Src/core/nbo_run.h. Compare a "
                       "candidate against these with compare_nbo.py.",
        "generated": datetime.date.today().isoformat(),
        "provenance": "Every field here is read back from the outputs the programs wrote: "
                      "the level of theory, charge, multiplicity, basis size and final "
                      "energy from each ORCA output, the thresholds from what NBO echoed. "
                      "If a future NBO or ORCA build changes a default, regenerating shows "
                      "the drift instead of absorbing it.",
        "levels_of_theory": sorted(levels),
        "keyword_notes": {
            "NRT": "natural resonance theory: resonance weights and bond orders",
            "E2PERT": "second-order perturbative donor-acceptor (E2) table",
            "NRTLST=0.1": "print the $NRTSTR keylist down to 0.1% weight, so the reference "
                          "carries the bond topology of the small structures too",
            "NRTDTL": "per-candidate TOPO matrices and rhoNL, the QP iteration path, the "
                      "5-decimal weight vector including the zero tail, the ARROWS lines and "
                      "the symmetry-equivalence map. NRTDTL=EXTRA produced a byte-identical "
                      "output on this build; NRTDTL=EXCESS only adds NAO-basis density "
                      "matrices and nothing about the search.",
            "NRTE2": "not passed for any main-group molecule, so their NRT delocalisation-list "
                     "threshold is whatever the NBO default is; the value NBO reported for each "
                     "run is in thresholds_reported_by_nbo.nrt_deloc_kcal. It is passed for the "
                     "transition-metal cases, where the default makes the search intractable - "
                     "see molecules[].note. The threshold is not a convergence parameter: on "
                     "Ni(CO)4, NRTE2=5 retains 16 of 136 candidates at D(w)=0.06910797 while "
                     "NRTE2=10 and NRTE2=20 both collapse to the same 4 candidates at "
                     "D(w)=0.07174346. Two runs at different NRTE2 are not comparable entry by "
                     "entry, and TiCl4 says which entries: 10 against 20 gives identical D(w) "
                     "and identical bond orders with individual weights up to 7.1 percentage "
                     "points apart. Compare D(w), bond orders and valencies across thresholds; "
                     "do not compare the weight of one structure.",
            "NRTSUB": "does not exist in NBO 7.0.9. A subspace analysis is asked for with the "
                      "bracket form 'NRT <atom list>' or by supplying $NRTSTR.",
            "NRTSYM=off": "only sf6, and not a choice: NRT's symmetry detection aborts the "
                          "whole analysis on it. See molecules[].note.",
            "$NRTSTR": "only so2, and not a choice either: NRT cannot start from the Lewis "
                       "structure NBO itself localises. See molecules[].note. A $NRTSTR or "
                       "$CHOOSE keylist must be appended AFTER the data blocks of the .47; "
                       "placed directly behind the $NBO line NBO reports \"$CHOOSE error: "
                       "'END' is not an acceptable orbital type\".",
        },
        "execution": {
            "machine": "FLOWOFFICE, 48 cores, held to about 30% by standing instruction",
            "nbo_threads": 1,
            "nbo_threaded": False,
            "threading_evidence": "the gennbo/nbo7.i8 binaries are statically linked ELF with "
                                  "no GOMP_ strings and no pthread_create symbol, and top "
                                  "shows a single process at ~70-100% of one core for the "
                                  "whole NRT run. NRT timings here are single-thread times.",
            "contention": "ORCA ran at 14 processes for the geometry step of the same "
                          "molecule set, and three sibling NoSpherA2 builds shared the "
                          "machine, so wall times carry noise that CPU seconds do not.",
        },
        "timing_note": "three different numbers, all recorded: timings.nbo_seconds is wall "
                       "time measured around the gennbo process by the wrapper (includes "
                       "process start and WSL crossing), nbo_reported_cpu_seconds and "
                       "nbo_reported_wall_seconds are NBO's own closing line, and the NRT "
                       "section's own Timing(sec) block splits search/Gram/minimise. Use the "
                       "CPU seconds for a speedup claim.",
        "start_dependence": {
            "what": "NRT was re-run from several starting points per molecule and the range "
                    "and standard deviation of every weight, bond order and D(w) recorded. "
                    "Files are start_dependence/<molecule>_<group>.json for all levers and "
                    "start_dependence/<molecule>_<group>_starts.json for the starts alone "
                    "(NRTCYC=1 removed, because it truncates the search instead of moving "
                    "its start). Molecules: benzene, formate, hcn, formaldehyde, n2 and the "
                    "open-shell nh3li.",
            "levers_keywords": ["NRTNBI=off", "NRTSYM=off", "NRTION", "NRTCYC=1",
                                "$CHOOSE forcing resonance structure 2", "$CHOOSE structure 3"],
            "levers_atom_order": ["reversed atom order", "shuffled atom order (fixed seed), "
                                  "single point on the optimised geometry"],
            "result": "Every lever that does not change the candidate-generation rules "
                      "reproduces D(w) to all 8 printed decimals, the same retained/candidate "
                      "counts, zero range in every bond order and zero range in the weights "
                      "compared by rank. Only the labelling moves: a different start can make "
                      "an equivalent structure the reference, which renames every structure's "
                      "Added(Removed) description and shows up as an apparent 22.9 (formate) "
                      "to 37.8 (benzene) percentage-point range if structures are matched by "
                      "description. Measured on six molecules: D(w) range is exactly 0 for "
                      "benzene, formate, hcn, formaldehyde and n2, and every bond order range "
                      "is 0.0000 on all five. NRTCYC=1 truncates the search itself and is not "
                      "a start; with it included the ranges become non-zero (formaldehyde "
                      "17.5 points on the leading weight, hcn 0.77), which is why the _starts "
                      "files leave it out. "
                      "NRTION is the one lever that changes the answer, and only on the "
                      "open-shell case: nh3li 5/11 vs 4/7 candidates, D(w) 0.03564276 vs "
                      "0.03605685, rank-1 weight moved 20.05 points, one bond order 0.2005.",
            "atom_order_result": "Fully insensitive. After mapping the atom numbering back "
                                 "through the permutation, every bond order is identical to "
                                 "the printed 4 decimals and every weight by rank is "
                                 "identical; D(w) moves by 1.2e-7 (formate) and 7e-8 "
                                 "(benzene), which is the geometry print precision of the "
                                 "reoptimised input, not the search.",
            "not_run": "no atom-reorder run for the open-shell nh3li case: its wavefunction "
                       "was set up before this study and its level of theory is not recorded "
                       "in an input file, so a re-run could not be made comparable.",
        },
        "thresholds_note": "thresholds_reported_by_nbo holds the four thresholds NBO "
                           "printed for that run (E2 printing, E2 intermolecular, NRT "
                           "parent structure, NRT delocalisation list). A result produced "
                           "with different thresholds is not comparable entry by entry: "
                           "below a printing threshold a missing entry and a small one "
                           "look the same.",
        "molecules": entries,
        "pending": {
            "molecules": pending,
            "why": "NRT on a transition-metal complex is where the combinatorics bite. With the "
                   "uniform keylist TiCl4 had enumerated 1696 distinct candidate structures and "
                   "then sat in its Gram-matrix phase, and Ni(CO)4 was still generating, both "
                   "past 90 minutes on one core. A second pass at NRTE2=5 (the delocalisation "
                   "threshold that decides which E2 interactions enter the search; the feasibility "
                   "measurement put 1 -> 5 at 15x cheaper for 0.0007 in D(w)) is running. Their "
                   "NPA, NAO, NBO and E2 data exist and are correct; only the NRT section is "
                   "unfinished, and a half-converged NRT section is worse in a reference than an "
                   "absent one. They are the two molecules that should have gone to the AKL "
                   "cluster.",
        } if pending else {},
    }
    index["start_dependence"]["files"] = collect_spread(dest)
    index["validation"] = {
        "open_shell_archive": "ch3_orca_reference.47 is the FILE.47 ORCA 6.1.1 wrote itself "
                              "for an unrestricted r2SCAN0/def2-TZVP CH3 doublet (! r2scan0 "
                              "def2-tzvp d4 NBO). Feeding the same .gbw to NoSpherA2 -nbo and "
                              "running gennbo on both archives gives byte-identical NPA "
                              "summaries in all three blocks (alpha, beta, total), identical "
                              "occupied NBO occupancies and hybridisations and an identical "
                              "E2 table. compare_47.py on the two archives: $OVERLAP agrees "
                              "to 6.6e-16, $DENSITY (both spin blocks, 2450 values) to "
                              "2.1e-10, $LCAOMO to 1.2e-9, $FOCK to 1.5e-5 - the Fock matrix "
                              "is rebuilt from the orbital energies and coefficients rather "
                              "than read, so it carries the SCF residual. The only visible "
                              "consequence is the energies and compositions of the "
                              "zero-occupancy Rydberg orbitals, which no analysed quantity "
                              "depends on.",
        "scripts": ["compare_nbo.py (parsed results)", "compare_47.py (archives)"],
    }
    with open(os.path.join(dest, "index.json"), "w") as f:
        json.dump(index, f, indent=1)
        f.write("\n")
    print("wrote", len(entries), "molecules to", dest)
    print("levels:", sorted(levels))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
