#include "pch.h"
#include "convenience.h"

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/sysctl.h>
#endif

namespace {

	bool is_valid_occ_data_path(const std::filesystem::path &base_path)
	{
		return std::filesystem::is_directory(base_path) &&
			std::filesystem::is_directory(base_path / "basis") &&
			std::filesystem::is_directory(base_path / "methods");
	}

	std::filesystem::path resolve_executable_directory(const char *argv0)
	{
		// Prefer argv0 — it comes directly from the program's own argument vector
		// and avoids reading /proc/self (which Flawfinder flags as a potential
		// information-hiding / unterminated-result risk on Linux).
		if (argv0 != nullptr)
		{
			try
			{
				auto p = std::filesystem::absolute(std::filesystem::path(argv0)).parent_path();
				if (!p.empty() && std::filesystem::exists(p))
					return p;
			}
			catch (...)
			{
				// Filesystem operations may throw (e.g. on very unusual path strings);
				// fall through to the platform-specific fallback below.
			}
		}

		// Fall back to platform-specific APIs when argv0 is unavailable or unusable.
#ifdef _WIN32
		char module_path[MAX_PATH];
		DWORD len = GetModuleFileNameA(nullptr, module_path, MAX_PATH);
		if (len > 0 && len < MAX_PATH)
		{
			return std::filesystem::path(module_path).parent_path();
		}
#elif defined(__APPLE__)
		uint32_t size = 0;
		_NSGetExecutablePath(nullptr, &size);
		std::string exe_path(size, '\0');
		if (_NSGetExecutablePath(exe_path.data(), &size) == 0)
		{
			return std::filesystem::path(exe_path.c_str()).parent_path();
		}
#else
		char exe_path[4096];
		ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1); /* Flawfinder: ignore - fixed path, bounded, terminated below */
		if (len > 0)
		{
			exe_path[len] = '\0';
			return std::filesystem::path(exe_path).parent_path();
		}
#endif

		return {};
	}

	bool set_occ_data_path(const std::filesystem::path &path)
	{
#ifdef _WIN32
		return _putenv_s("OCC_DATA_PATH", path.string().c_str()) == 0;
#else
		return setenv("OCC_DATA_PATH", path.string().c_str(), 1) == 0;
#endif
	}

	std::filesystem::path choose_occ_data_path_from_exe_dir(const std::filesystem::path &exe_dir)
	{
		const auto parent = exe_dir.parent_path();
		const auto grandparent = parent.parent_path();

		const std::vector<std::filesystem::path> candidates = {
			exe_dir,
			exe_dir / "occ" / "share",
			exe_dir / "share",
			parent / "occ" / "share",
			parent / "share",
			grandparent / "occ" / "share"
		};

		for (const auto &candidate : candidates)
		{
			if (!candidate.empty() && is_valid_occ_data_path(candidate))
			{
				return candidate;
			}
		}
		return {};
	}
}

bool ensure_occ_data_path(const char *argv0)
{
#ifdef _WIN32
	char* occ_data_path_env = nullptr;
	size_t len = 0;
	errno_t err = _dupenv_s(&occ_data_path_env, &len, "OCC_DATA_PATH");

	if (err == 0 && occ_data_path_env != nullptr)
	{
		std::filesystem::path path(occ_data_path_env);
		free(occ_data_path_env);

		if (is_valid_occ_data_path(path))
			return true;
		else
			std::cerr << "OCC DATA PATH is invalid!" << std::endl;
	}
	else if (occ_data_path_env != nullptr)
	{
		free(occ_data_path_env);
	}
#else
	const char* tmp_occ_data_path_env = std::getenv("OCC_DATA_PATH"); // Flawfinder: ignore
	if (tmp_occ_data_path_env != nullptr)
	{
		std::string occ_data_path_env(tmp_occ_data_path_env);
		if (is_valid_occ_data_path(std::filesystem::path(occ_data_path_env)))
			return true;
		else
			std::cerr << "OCC DATA PATH is invalid!" << std::endl;
	}
#endif

	std::filesystem::path exe_dir = resolve_executable_directory(argv0);
	if (exe_dir.empty())
	{
		std::cerr << "OCC_DATA_PATH not set or invalid and executable directory could not be resolved." << std::endl;
		return false;
	}

	std::filesystem::path selected_path = choose_occ_data_path_from_exe_dir(exe_dir);
	if (selected_path.empty())
	{
		std::cerr << "OCC_DATA_PATH not set or invalid. No valid OCC data directory found near executable directory: "
			<< exe_dir << std::endl;
		return false;
	}

	if (!set_occ_data_path(selected_path))
	{
		std::cerr << "Failed to set OCC_DATA_PATH to resolved directory: " << selected_path << std::endl;
		return false;
	}

	std::cerr << "OCC_DATA_PATH not set or invalid. Using: " << selected_path << std::endl;
	return true;
}

std::filesystem::path get_home_path(void)
{
#ifdef _WIN32
	char *homeDrive = nullptr;
	size_t len = 0;
	std::string temp1, temp2;

	errno_t err = _dupenv_s(&homeDrive, &len, "HOMEDRIVE");
	if (err == 0 && homeDrive != nullptr) {
		temp1 = homeDrive;
		// Free the allocated memory
		free(homeDrive);
	}
	else {
		std::cerr << "Failed to retrieve the environment variable." << std::endl;
	}
	err = _dupenv_s(&homeDrive, &len, "HOMEPATH");
	if (err == 0 && homeDrive != nullptr) {
		temp2 = homeDrive;
		// Free the allocated memory
		free(homeDrive);
	}
	else {
		std::cerr << "Failed to retrieve the environment variable." << std::endl;
	}
	temp1.append(temp2);
	return temp1;
#else
	const char *home_env = getenv("HOME"); // Flawfinder: ignore
	if (home_env == nullptr) {
		std::cerr << "Warning: HOME environment variable not set." << std::endl;
		return std::filesystem::path("/tmp"); // Fallback to /tmp
	}

	std::string home = home_env;
	// Basic validation: check if it's a valid path and not empty
	if (home.empty() || home.find_first_of('\0') != std::string::npos) {
		std::cerr << "Warning: Invalid HOME environment variable." << std::endl;
		return std::filesystem::path("/tmp"); // Fallback to /tmp
	}

	return std::filesystem::path(home);
#endif
}


size_t available_memory_bytes()
{
	size_t avail = 0;
#ifdef _WIN32
	MEMORYSTATUSEX status{};
	status.dwLength = sizeof(status);
	if (GlobalMemoryStatusEx(&status))
		avail = static_cast<size_t>(status.ullAvailPhys);
	//A scheduler's or container's job object limit, the Windows counterpart of the cgroup ceiling;
	//nullptr queries this process's job, a process outside one reports no limit flags.
	JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
	DWORD returned = 0;
	if (QueryInformationJobObject(nullptr, JobObjectExtendedLimitInformation,
		&limits, sizeof(limits), &returned)) {
		const DWORD flags = limits.BasicLimitInformation.LimitFlags;
		if (flags & JOB_OBJECT_LIMIT_PROCESS_MEMORY) {
			const size_t l = static_cast<size_t>(limits.ProcessMemoryLimit);
			if (l > 0 && (avail == 0 || l < avail)) avail = l;
		}
		if (flags & JOB_OBJECT_LIMIT_JOB_MEMORY) {
			const size_t l = static_cast<size_t>(limits.JobMemoryLimit);
			if (l > 0 && (avail == 0 || l < avail)) avail = l;
		}
	}
#elif defined(__APPLE__)
	//macOS keeps most memory occupied: inactive and purgeable pages are returned on demand and count as available.
	vm_size_t page_size = 0;
	mach_port_t host = mach_host_self();
	if (host_page_size(host, &page_size) == KERN_SUCCESS) {
		vm_statistics64_data_t vm{};
		mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
		if (host_statistics64(host, HOST_VM_INFO64,
			reinterpret_cast<host_info64_t>(&vm), &count) == KERN_SUCCESS) {
			const uint64_t pages = static_cast<uint64_t>(vm.free_count)
				+ vm.inactive_count + vm.purgeable_count;
			avail = static_cast<size_t>(pages * static_cast<uint64_t>(page_size));
		}
	}
	if (avail == 0) {
		//Fallback: physical size, an upper bound the caller's fraction keeps under.
		uint64_t total = 0;
		size_t len = sizeof(total);
		if (sysctlbyname("hw.memsize", &total, &len, nullptr, 0) == 0)
			avail = static_cast<size_t>(total);
	}
#else
	if (std::FILE *f = std::fopen("/proc/meminfo", "r")) {
		char line[256];
		while (std::fgets(line, sizeof(line), f)) {
			unsigned long long kb = 0;
			if (std::sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
				avail = static_cast<size_t>(kb) * 1024ULL;
				break;
			}
		}
		std::fclose(f);
	}
	//Under a batch system or container the cgroup is the ceiling, while meminfo reports the whole node.
	//SLURM sets the limit on the job's and step's cgroup, never on the mount root (AKL:
	///system.slice/slurmstepd.scope/job_N/step_0/user/task_0), so walk from this process's own
	//cgroup up to the root and keep the smallest limit minus what that cgroup already uses.
	//Reclaimable page cache (inactive_file) is not usage: a large file just written would otherwise eat the budget.
	auto read_ull = [](const std::string& p, unsigned long long& v) {
		std::FILE* f = std::fopen(p.c_str(), "r");
		if (!f) return false;
		const bool ok = std::fscanf(f, "%llu", &v) == 1; //v2 "max" fails the scan
		std::fclose(f);
		return ok;
	};
	auto stat_value = [](const std::string& p, const char* key) {
		unsigned long long v = 0;
		if (std::FILE* f = std::fopen(p.c_str(), "r")) {
			char line[256], k[128];
			unsigned long long x = 0;
			while (std::fgets(line, sizeof(line), f))
				if (std::sscanf(line, "%127s %llu", k, &x) == 2 && std::strcmp(k, key) == 0) { v = x; break; }
			std::fclose(f);
		}
		return v;
	};
	std::string v2, v1;
	if (std::FILE* f = std::fopen("/proc/self/cgroup", "r")) {
		char line[1024];
		while (std::fgets(line, sizeof(line), f)) {
			std::string l(line);
			while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
			const size_t a = l.find(':'), b = a == std::string::npos ? a : l.find(':', a + 1);
			if (b == std::string::npos) continue;
			const std::string ctl = "," + l.substr(a + 1, b - a - 1) + ",";
			if (l.compare(0, a, "0") == 0 && ctl == ",,") v2 = l.substr(b + 1);
			else if (ctl.find(",memory,") != std::string::npos) v1 = l.substr(b + 1);
		}
		std::fclose(f);
	}
	//Hybrid hosts list both; the memory controller then sits on v1
	const bool use_v1 = !v1.empty() || !std::filesystem::exists("/sys/fs/cgroup/cgroup.controllers");
	const std::string mount = use_v1 ? "/sys/fs/cgroup/memory" : "/sys/fs/cgroup";
	const char* max_f = use_v1 ? "/memory.limit_in_bytes" : "/memory.max";
	const char* cur_f = use_v1 ? "/memory.usage_in_bytes" : "/memory.current";
	const char* cache_key = use_v1 ? "total_inactive_file" : "inactive_file";
	for (std::string dir = use_v1 ? v1 : v2;;) {
		if (!dir.empty() && dir.back() == '/') dir.pop_back();
		unsigned long long lim = 0, cur = 0;
		//v1 writes a sentinel near SIZE_MAX when there is no limit
		if (read_ull(mount + dir + max_f, lim) && lim > 0 && lim < (1ULL << 62)) {
			read_ull(mount + dir + cur_f, cur);
			const unsigned long long cache = stat_value(mount + dir + "/memory.stat", cache_key);
			cur = cur > cache ? cur - cache : 0;
			//1, not 0: a full cgroup is "nothing left", 0 would read as "unknown, no cap"
			const size_t l = lim > cur ? static_cast<size_t>(lim - cur) : 1;
			if (avail == 0 || l < avail) avail = l;
		}
		if (dir.empty()) break;
		dir.erase(dir.find_last_of('/'));
	}
#endif
	return avail;
}

size_t resident_memory_bytes()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS pmc{};
	return K32GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)) ? static_cast<size_t>(pmc.WorkingSetSize) : 0;
#elif defined(__APPLE__)
	mach_task_basic_info_data_t info{};
	mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
	return task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS
		? static_cast<size_t>(info.resident_size) : 0;
#else
	unsigned long long pages = 0, resident = 0;
	if (std::FILE* f = std::fopen("/proc/self/statm", "r")) {
		if (std::fscanf(f, "%llu %llu", &pages, &resident) != 2) resident = 0;
		std::fclose(f);
	}
	return static_cast<size_t>(resident) * static_cast<size_t>(sysconf(_SC_PAGESIZE));
#endif
}
