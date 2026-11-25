"""Byte-level check that two NoSpherA2 builds produce the same results.

run:      run every tests/tests.toml case with one executable and record a sha256 of every file each case
          writes or changes, plus its masked console output and exit code.
compare:  compare two such manifests; any difference is reported and the exit code is 1.

Cases that share a directory run one after another in toml order inside one copy of that directory (later
cases may read what earlier ones wrote, as under ctest); directories run in parallel. Inputs are the files git
tracks under tests/, copied fresh for every run, so leftovers of earlier runs cannot leak in.
Only clock times, dates and durations are masked, and only in log-type text (*.log, *_log, console);
every other file is compared byte for byte.
"""
import argparse, concurrent.futures as cf, hashlib, json, os, re, shutil, subprocess, sys, time
from pathlib import Path

LOG_SUFFIXES = ('.log', '_log')
MASKS = [
	(re.compile(rb'\d{1,2}:\d{2}:\d{2,3}(\.\d+)?'), b'<clock>'),
	(re.compile(rb'\b(Mon|Tue|Wed|Thu|Fri|Sat|Sun)\b[^\n]*\d{4}'), b'<date>'),
	(re.compile(rb'\d{4}-\d{2}-\d{2}'), b'<date>'),
	(re.compile(rb'((seconds|per_structure)\s*=\s*)[-+0-9.eE]+'), rb'\1<t>'),
	(re.compile(rb'\d+(\.\d+)?(e[-+]?\d+)?\s*(ms|s|sec|secs|seconds|min|mins|minutes|h|hours|us)\b(?!=)'), b'<dur>'),
	# derived from timings: host/device speed ratio, share of the run time next to a clock time
	(re.compile(rb'(costs <dur> a point[^\n]*?the device )[0-9.]+'), rb'\1<r>'),
	(re.compile(rb'(<clock>\s+\()[0-9.]+%%\)'), rb'\1<p>)'),
]

# occ's SCF iteration table: iter, energy, dE, commutator, seconds - only the seconds column is masked
OCC_ITER = re.compile(rb'(?m)^(\s*\d+\s+-?\d+\.\d{6,}\s+-?[\d.]+e[-+]\d+\s+-?[\d.]+e[-+]\d+\s+)[\d.]+e[-+]\d+(\s*)$')

# timing fields of the NBO json report ("nao_seconds": 0.12) - only the value is masked
JSON_T = re.compile(rb'("[A-Za-z0-9_]*seconds"\s*:\s*)[-+0-9.eE]+')
JSON_TO = re.compile(rb'("seconds"\s*:\s*\{)[^{}]*\}')

def mask(data):
	data = OCC_ITER.sub(rb'\1<t>\2', data)
	for rx, rep in MASKS:
		data = rx.sub(rep, data)
	return data

def sha(data):
	return hashlib.sha256(data).hexdigest()

def unpath(data, places):
	# outputs that name their own run directory (vmd scripts, logs) must not differ by where the run happened
	for path, tag in places:
		for v in {str(path), str(path).replace('\\', '/')}:
			data = data.replace(v.encode(), tag)
	return data

def file_digest(p, places):
	data = unpath(p.read_bytes(), places)
	if p.name.endswith('.json'):
		data = JSON_TO.sub(rb'\1<t>}', JSON_T.sub(rb'\1<t>', data))
	return sha(mask(data)) if p.name.endswith(LOG_SUFFIXES) else sha(data)

EXTRA = Path(__file__).resolve().parent / 'result_equivalence_extra_cases.toml'

def load_cases(repo, want_full, only, extra=None):
	extra = Path(extra) if extra else EXTRA
	import tomllib
	d = tomllib.loads((repo / 'tests' / 'tests.toml').read_text(encoding='utf-8'))
	defaults = d.pop('defaults', {})
	items = list(d.items())
	if extra.exists():
		x = tomllib.loads(extra.read_text(encoding='utf-8'))
		x.pop('defaults', None)
		items += [('x_' + k, v) for k, v in x.items()]
	cases = []
	for name, sec in items:
		if only and name not in only:
			continue
		if sec.get('full') and not want_full:
			continue
		args = []
		for k, v in list(defaults.items()) + list(sec.get('args', {}).items()):
			vals = v if isinstance(v, list) else [v]
			args.append('-' + k)
			for x in vals:
				if x is True or x in ('', 'true', 'True'):
					continue
				if x is False or x in ('false', 'False'):
					continue
				args.append(repr(x) if isinstance(x, float) else str(x))
		cases.append({'name': name, 'dir': sec['directory'], 'args': args})
	return cases

PLAN = None  # {'cases': [...], 'files': {dir: [paths]}} written by 'plan' for hosts without git / Python 3.11

def tracked(repo, sub):
	if PLAN is not None:
		return PLAN['files'][sub]
	out = subprocess.run(['git', '-C', str(repo), 'ls-files', '-z', 'tests/' + sub], capture_output=True, check=True).stdout
	return [x for x in out.decode().split('\0') if x]

def snapshot(root):
	return {str(p.relative_to(root)).replace('\\', '/'): (st.st_size, st.st_mtime_ns)
		for p in root.rglob('*') if p.is_file() for st in [p.stat()]}

def run_group(repo, exe, out, dirname, group, cpus, timeout, occ_threads):
	root = out / 'work' / dirname.replace('/', '_')
	if root.exists():
		shutil.rmtree(root)
	needed = {dirname}
	for c in group:
		for a in c['args']:
			m = re.match(r'\.\./([^/\\]+)', a)
			if m:
				needed.add(m.group(1))
	for sub in sorted(needed):
		for f in tracked(repo, sub):
			dst = root / f
			dst.parent.mkdir(parents=True, exist_ok=True)
			shutil.copy2(repo / f, dst)
	if occ_threads:
		# occ takes its thread count from the input file, not from -cpus
		for t in (root / 'tests').rglob('*.toml'):
			txt = t.read_text(encoding='utf-8')
			new = re.sub(r'(?m)^threads\s*=\s*\d+', f'threads = {occ_threads}', txt)
			if new != txt:
				t.write_text(new, encoding='utf-8')
	places = [(root, b'<root>'), (Path(exe).parent, b'<exe>')]
	env = dict(os.environ, OCC_DATA_PATH=str(repo / 'occ' / 'share'))
	res = {}
	cwd = root / 'tests' / dirname
	for c in group:
		before = snapshot(root)
		argv = [str(exe)] + (['-cpus', str(cpus)] if cpus else []) + c['args']
		t0 = time.time()
		try:
			p = subprocess.run(argv, cwd=cwd, env=env, capture_output=True, timeout=timeout, stdin=subprocess.DEVNULL)
			rc, console = p.returncode, p.stdout + b'\n--stderr--\n' + p.stderr
		except subprocess.TimeoutExpired:
			rc, console = 'timeout', b''
		after = snapshot(root)
		files = {rel: file_digest(root / rel, places) for rel, st in after.items() if before.get(rel) != st}
		files['<console>'] = sha(mask(unpath(console, places)))
		(out / 'console').mkdir(exist_ok=True)
		(out / 'console' / (c['name'] + '.txt')).write_bytes(console)
		res[c['name']] = {'rc': rc, 'files': files, 'secs': round(time.time() - t0, 1), 'argv': argv[1:]}
	return res

def cmd_run(a):
	repo, exe, out = Path(a.repo).resolve(), Path(a.exe).resolve(), Path(a.out).resolve()
	out.mkdir(parents=True, exist_ok=True)
	global PLAN
	if a.plan:
		PLAN = json.loads(Path(a.plan).read_text())
		only = set(a.cases.split(',')) if a.cases else None
		cases = [c for c in PLAN['cases'] if not only or c['name'] in only]
	else:
		cases = load_cases(repo, a.full, set(a.cases.split(',')) if a.cases else None, a.extra_file)
	groups = {}
	for c in cases:
		groups.setdefault(c['dir'], []).append(c)
	manifest = {}
	order = sorted(groups, key=lambda d: -len(groups[d]))
	with cf.ThreadPoolExecutor(a.jobs) as ex:
		futs = {ex.submit(run_group, repo, exe, out, d, groups[d], a.cpus, a.timeout, a.occ_threads): d for d in order}
		for f in cf.as_completed(futs):
			r = f.result()
			manifest.update(r)
			for n, v in r.items():
				print(f'{n}: rc={v["rc"]} files={len(v["files"])} {v["secs"]}s', flush=True)
	(out / 'manifest.json').write_text(json.dumps({'exe': str(exe), 'cpus': a.cpus, 'cases': manifest}, indent=1, sort_keys=True))
	if not a.keep:
		shutil.rmtree(out / 'work', ignore_errors=True)

def cmd_plan(a):
	repo = Path(a.repo).resolve()
	cases = load_cases(repo, a.full, None, a.extra_file)
	dirs = {c['dir'] for c in cases}
	for c in cases:
		for x in c['args']:
			m = re.match(r'\.\./([^/\\]+)', x)
			if m:
				dirs.add(m.group(1))
	Path(a.out).write_text(json.dumps({'cases': cases, 'files': {d: tracked(repo, d) for d in sorted(dirs)}}, indent=1))
	print(len(cases), 'cases,', len(dirs), 'directories')

def diff(A, B):
	bad = []
	for n in sorted(set(A) | set(B)):
		if n not in A or n not in B:
			bad.append((n, 'case missing on one side', ''))
			continue
		a, b = A[n], B[n]
		if a['rc'] != b['rc']:
			bad.append((n, 'exit code', f'{a["rc"]} vs {b["rc"]}'))
		for f in sorted(set(a['files']) | set(b['files'])):
			if a['files'].get(f) != b['files'].get(f):
				bad.append((n, f, 'missing' if f not in a['files'] or f not in b['files'] else 'bytes differ'))
	return bad

def load(p):
	return json.loads(Path(p).read_text())['cases']

def cmd_gate(a):
	"""Every file must match master at -cpus 8; a file master itself does not reproduce at -cpus 8 (two master
	runs differ) must match master at -cpus 1 instead, where master must reproduce itself."""
	b = Path(a.baseline)
	mt1, mt2, st1, st2 = (load(b / x / 'manifest.json') for x in ('mt1', 'mt2', 'st1', 'st2'))
	cmt, cst = load(Path(a.candidate) / 'mt' / 'manifest.json'), load(Path(a.candidate) / 'st' / 'manifest.json')
	noisy_mt = {(n, f) for n, f, _ in diff(mt1, mt2)}
	noisy_st = {(n, f) for n, f, _ in diff(st1, st2)}
	fails, covered_by_st, unverifiable = [], [], []
	for n, f, why in diff(mt1, cmt):
		if (n, f) in noisy_mt:
			covered_by_st.append((n, f))
		else:
			fails.append(('cpus8', n, f, why))
	for n, f, why in diff(st1, cst):
		if (n, f) in noisy_st:
			unverifiable.append((n, f))
		else:
			fails.append(('cpus1', n, f, why))
	for x in fails:
		print('FAIL', *x)
	for n, f in sorted(noisy_st):
		print('UNVERIFIABLE (master differs from itself at -cpus 1)', n, f)
	print(f'{len(mt1)} cases; {len(noisy_mt)} files master does not reproduce at -cpus 8 (checked at -cpus 1 instead); '
		f'{len(noisy_st)} not reproducible even at -cpus 1; {len(fails)} FAILURES')
	sys.exit(1 if fails else 0)

def cmd_compare(a):
	A = json.loads(Path(a.reference).read_text())['cases']
	B = json.loads(Path(a.candidate).read_text())['cases']
	bad = diff(A, B)
	for n, f, why in bad:
		print(f'DIFF {n}: {f}: {why}')
	print(f'{len(A)} reference cases, {len(B)} candidate cases, {len(bad)} differences')
	sys.exit(1 if bad else 0)

if __name__ == '__main__':
	ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	sp = ap.add_subparsers(dest='cmd', required=True)
	r = sp.add_parser('run')
	r.add_argument('--exe', required=True)
	r.add_argument('--repo', required=True, help='checkout whose tests/ and occ/share are the inputs')
	r.add_argument('--out', required=True)
	r.add_argument('--cpus', type=int, default=0, help='prepend -cpus N (0: leave the default)')
	r.add_argument('--jobs', type=int, default=4)
	r.add_argument('--full', action='store_true', help='include the full = true cases')
	r.add_argument('--cases', default='', help='comma-separated subset')
	r.add_argument('--timeout', type=int, default=7200)
	r.add_argument('--occ-threads', type=int, default=0, help='rewrite threads = N in the copied occ input files')
	r.add_argument('--extra-file', default='', help='extra cases toml instead of result_equivalence_extra_cases.toml beside this script')
	r.add_argument('--plan', default='', help='case plan from the plan command (no git, no tomllib needed)')
	r.add_argument('--keep', action='store_true', help='keep the work directories')
	r.set_defaults(f=cmd_run)
	pl = sp.add_parser('plan', help='resolve the cases and tracked input files into JSON for a run without git/tomllib')
	pl.add_argument('--repo', required=True)
	pl.add_argument('--out', required=True)
	pl.add_argument('--full', action='store_true')
	pl.add_argument('--extra-file', default='')
	pl.set_defaults(f=cmd_plan)
	c = sp.add_parser('compare')
	c.add_argument('reference')
	c.add_argument('candidate')
	c.set_defaults(f=cmd_compare)
	g = sp.add_parser('gate', help='compare a candidate (mt/, st/) against a baseline (mt1/ mt2/ st1/ st2/)')
	g.add_argument('baseline')
	g.add_argument('candidate')
	g.set_defaults(f=cmd_gate)
	a = ap.parse_args()
	a.f(a)
