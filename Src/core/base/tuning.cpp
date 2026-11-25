#include "pch.h"
#include "tuning.h"

static std::map<std::string, std::string> &tuning_knobs()
{
	static std::map<std::string, std::string> knobs;
	return knobs;
}
const char *tuning(const char *name)
{
	const auto &knobs = tuning_knobs();
	const auto it = knobs.find(name);
	return it == knobs.end() ? nullptr : it->second.c_str();
}
void set_tuning(const std::string &name, const char *value)
{
	if (value)
		tuning_knobs()[name] = value;
	else
		tuning_knobs().erase(name);
}
//The knobs and defaults of docs/wiki/Tuning-Knobs.md, which scripts/check_tuning_knobs.py keeps equal to this table
static const std::vector<std::pair<std::string, std::string>> documented_knobs = {
	{"NOSPHERA2_XCW_RCHECK=<fchk>", "unset"},
	{"NOSPHERA2_XCW_ERI_NO_BUDGET", "off"},
	{"NOSPHERA2_ERI_DYNAMIC", "off"},
	{"NOSPHERA2_XCW_ADJOINT_CHECK", "off"},
	{"NOSPHERA2_XCW_TSHB_FIT", "off"},
	{"NOSPHERA2_XCW_TSHB_TABLE", "off"},
	{"NOSPHERA2_XCW_TSCB_BY_F", "off"},
	{"NOSPHERA2_DEBUG_ASYMFACT", "off"},
	{"NOSPHERA2_NO_INVERSION_SCREENING", "off"},
	{"NOSPHERA2_XCW_GRID_DIRECT", "off"},
	{"NOSPHERA2_XCW_GRID_FB_H", "0.25"},
	{"NOSPHERA2_XCW_GRID_DTOL", "0 (off)"},
	{"NOSPHERA2_XCW_GRID_PTOL", "0 (off)"},
	{"NOSPHERA2_XCW_GRID_FULL", "8"},
	{"NOSPHERA2_XCW_GRID_TIME", "off"},
	{"NOSPHERA2_XCW_GRID_AO_TOL", "-acc cutoff"},
	{"NOSPHERA2_XCW_GRID_NO_MORTON", "off"},
	{"NOSPHERA2_XCW_GRID_CHUNK", "1024"},
	{"NOSPHERA2_XCW_GRID_TSCB_OLD", "off"},
	{"NOSPHERA2_XCW_GRID_FCHECK", "off"},
	{"NOSPHERA2_XCW_WEIGHT_UPDATE=k", "-1"},
	{"NOSPHERA2_XCW_WEIGHT_DTOL", "1e-3"},
	{"NOSPHERA2_XCW_WEIGHT_FTOL", "1e-5"},
	{"NOSPHERA2_XCW_WEIGHT_TOL", "1e-6"},
	{"NOSPHERA2_XCW_WEIGHT_MIX", "6"},
	{"NOSPHERA2_XCW_WEIGHT_CHECK", "off"},
	{"NOSPHERA2_XCW_TFVC_OLD", "off"},
	{"NOSPHERA2_XCW_CHI_FALLBACK", "off"},
	{"NOSPHERA2_XCW_RI_FULL", "off"},
	{"NOSPHERA2_XCW_RI_PACK_CHECK", "off"},
	{"NOS_RI_SCREEN", "3e-7"},
	{"NOS_RI_FAR", "1e-14"},
	{"NOSPHERA2_RI_ROWS_NOCUT", "off"},
	{"NOSPHERA2_RI_PENALTY_LOOP", "off"},
	{"NOS_SF_DIRECT", "off"},
	{"NOS_SF_FB_L=<L>", "D/2"},
	{"NOS_SF_FB_EXACT", "off"},
	{"NOS_SF_FB_FORCE", "off"},
	{"NOS_GRID_OLD_FLOOR", "off"},
	{"NOS_MBIS_BECKE", "off"},
	{"NOS_EXP_CRT", "off"},
	{"NOSPHERA2_CHI_DEBUG", "off"},
	{"NOS_DENSITY_ACCURACY", "5e-5"},
	{"NOSPHERA2_GPU_BACKEND=cuda|hip", "probed"},
	{"NOSPHERA2_GPU_BATCH=<n>", "fits the card"},
	{"NOSPHERA2_GRID_LOCAL=1", "off"},
	{"NOSPHERA2_BLAS_GPU_MIN_FLOP", "4e9 x fp64 ratio / 64"},
	{"NOS_BASIN_GPU_MIN_US", "15"},
	{"NOS_BASIN_GPU_FIBERS", "256"},
	{"NOS_BASIN_GPU_PROFILE", "off"},
	{"NOS_BASIN_SKIP", "5e-5"},
	{"NOS_BASIN_STEP_CAP", "2000"},
	{"NOS_BASIN_STEP_RELAX", "1.5"},
	{"NOS_BASIN_SHRINK", "on"},
	{"NOS_BASIN_SHRINK_FLOOR", "0.0625"},
	{"NOS_ADP_CAP", "8"},
	{"NOS_ADP_GROW", "0.99999"},
	{"NOS_ADP_KEEP", "0.999"},
	{"NOS_ADP_REACH", "0.25"},
	{"NOS_BETA_MARGIN", "0.9"},
	{"NOS_BASIN_MEMO", "on"},
	{"NOS_BASIN_MEMO_MB", "256"},
	{"NOS_BASIN_ROTATE", "on"},
	{"NOS_BASIN_SCAN_REVERSE", "off"},
	{"NOS_ELI_SHELL_DIST", "1.2"},
	{"NOS_ELI_SHELL_TOL", "0.05"},
	{"NOS_ELI_SHELL_SADDLE", "0.99"},
	{"NOS_ELI_SHELL_PERSIST", "5e-4"},
	{"NOS_ELI_SHELL_SADDLE_LOG", "off"},
	{"NOS_RGBI_PINV_CUTOFF", "1e-5"},
	{"NOS_RGBI_PAIR_TOL", "1e-4"},
	{"NOS_RGBI_GROUP_LP", "0.99"},
	{"NOS_RGBI_PAIR_LP", "0.999"},
	{"NOS_RGBI_NONREL_FREEATOM", "off"},
	{"NOS_RGBI_FREE_ATOM_MAXITER=<n>", "occ's default"},
	{"NOS_RGBI_NO_FREEATOM_CACHE", "off"},
	{"NOS_FREEATOM_CACHE_DIR=<dir>", "~/.cache/NoSpherA2/free_atom_cache"},
	{"NOS_RGBI_PARALLEL_FREEATOM", "off"},
	{"NOS_RGBI_NO_PIN", "off"},
	{"NOS_RGBI_DEBUG", "off"},
	{"NBO_BOND_FLOOR", "0.15"},
	{"NBO_JSON_DIGITS", "10"},
	{"NAO_PRENAO_NET", "off"},
	{"NAO_LEGACY_CASCADE", "off"},
	{"NAO_OWSO_OFF", "off"},
	{"NAO_CLASS_SPLIT", "off"},
	{"NAO_CORE_POOLED", "off"},
	{"NAO_DUMP_CPRE", "off"},
	{"NAO_DUMP_STEP3", "off"},
	{"NAO_DUMP_C", "off"},
	{"NOSPHERA2_SPLINE_ACCURACY", "1e-6"},
	{"NOSPHERA2_MAX_RADIAL", "6"},
	{"NOSPHERA2_MAX_ANGULAR", "12"},
	{"NOSPHERA2_TIME_SOAP", "off"},
};
bool tuning_known(const std::string &name)
{
	for (int i = 0; i < documented_knobs.size(); i++)
		if (documented_knobs[i].first.substr(0, documented_knobs[i].first.find('=')) == name)
			return true;
	return false;
}
void tuning_list(std::ostream &out)
{
	out << "-tune NAME[=VALUE]: the knobs and their defaults; docs/wiki/Tuning-Knobs.md says what each does\n";
	for (int i = 0; i < documented_knobs.size(); i++)
		out << "  " << std::left << std::setw(34) << documented_knobs[i].first << " " << documented_knobs[i].second << "\n";
	out << std::right;
}
