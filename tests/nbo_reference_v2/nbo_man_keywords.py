"""What can NBO 7 be ASKED to write?  Read off the shipped manual, because the binary has no HELP.

`$NBO HELP $END` is not a keyword in NBO 7.0.9: the binary answers "Unrecognized $NBO keyword:
HELP" and then runs the default job, so a run that looks successful has told you nothing.  The
authority is therefore nbo7_man.pdf, which ships beside the binary at
/work/software/bin/NBO/man/nbo7_man.pdf.  That is documentation - nothing here reads, disassembles
or decompiles an executable, and no output of this script is redistributed.

    py -3.12 nbo_man_keywords.py matrix <nbo7_man.pdf>     # the MATKEY table and its default LFNs
    py -3.12 nbo_man_keywords.py cascade <nbo7_man.pdf>    # does ANY keyword name a stage INSIDE
                                                           # the AO -> NAO step?
    py -3.12 nbo_man_keywords.py pages <nbo7_man.pdf> 47 51

The question it serves: the NAO gap sits between gennbo's unit 32 (AO -> pre-NAO) and unit 33
(AO -> NAO), and the two codes are known to agree at 32 and to disagree at 33.  Anything gennbo
could be asked to write in between would settle where.
"""
import re
import sys

from pypdf import PdfReader

# a keyword of the form NAME = [W]nn, which is how the manual prints a default LFN
MAT = re.compile(r"\b([A-Z][A-Z0-9]{2,9})\s*=?\s*W?\s*(\d\d)\b")
# the vocabulary of the published cascade's inner stages.  If the manual exposed any of them as a
# keyword, it would appear beside one of these words.
CASCADE = ["schmidt", "occupancy-weighted", "occupancy weighted", "owso", "symmetry-average",
           "m-average", "rydberg", "natural minimal basis", "pre-orthogonal", "print level",
           "details", "direct access"]


def pages(path):
    r = PdfReader(path)
    out = []
    for i, p in enumerate(r.pages, 1):
        try:
            out.append((i, p.extract_text() or ""))
        except Exception as exc:               # reported, never silently dropped
            out.append((i, "(EXTRACT FAILED: %s)" % exc))
    return out


def cmd_matrix(path):
    kw = {}
    n_lines = 0
    for pno, t in pages(path):
        for line in t.splitlines():
            s = " ".join(line.split())
            if not s:
                continue
            n_lines += 1
            for m in MAT.finditer(s):
                kw.setdefault(m.group(1), set()).add((int(m.group(2)), pno))
    print("scanned %d non-empty lines" % n_lines)
    print("NAME = unit, as the manual prints it:")
    for k in sorted(kw):
        print("  %-12s %s" % (k, " ".join("%d@p%d" % u for u in sorted(kw[k]))))


def cmd_cascade(path):
    ps = pages(path)
    for term in CASCADE:
        pat = re.compile(re.escape(term), re.I)
        hits = [(pno, " ".join(l.split())) for pno, t in ps for l in t.splitlines()
                if pat.search(l)]
        print("\n### %-22s %d hits" % (term, len(hits)))
        for pno, line in hits[:10]:
            print("  p%-4d %s" % (pno, line[:170]))
        if len(hits) > 10:
            print("  ... %d more" % (len(hits) - 10))


def cmd_pages(path, a, b):
    for pno, t in pages(path):
        if int(a) <= pno <= int(b):
            print("\n========== p%d ==========\n%s" % (pno, t))


def main(argv):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if len(argv) < 3:
        print(__doc__)
        return 2
    {"matrix": cmd_matrix, "cascade": cmd_cascade, "pages": cmd_pages}[argv[1]](*argv[2:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
