#include "pch.h"
#include "throughput.h"

#include <algorithm>
#include <iomanip>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

namespace throughput {
namespace {

struct Row {
    std::string stage;
    bool on_device = false;
    double flops = 0.0;
    double ms = 0.0;
    double slowest = 0.0;
    long long calls = 0;
};

//A mutex rather than atomics: two doubles and a counter have to move together or a row can
//report work without the time that produced it. Contended once per kernel launch, which is
//nothing next to the launch itself.
std::mutex g_mutex;
std::vector<Row> g_rows;
bool g_enabled = false;

std::string rate(double flops, double ms)
{
    if (ms <= 0.0)
        return "-";
    std::ostringstream s;
    s << std::fixed << std::setprecision(1) << (flops / (ms * 1.0e6));
    return s.str();
}

} //namespace

void set_enabled(bool on)
{
    std::lock_guard<std::mutex> lock(g_mutex);
    g_enabled = on;
}

bool enabled()
{
    std::lock_guard<std::mutex> lock(g_mutex);
    return g_enabled;
}

void record(const char* stage, const bool on_device, const double flops, const double ms)
{
    if (!stage || flops < 0.0)
        return;
    std::lock_guard<std::mutex> lock(g_mutex);
    for (int i = 0; i < g_rows.size(); i++) {
        Row& r = g_rows[i];
        if (r.on_device == on_device && r.stage == stage) {
            r.flops += flops;
            r.ms += ms;
            if (ms > r.slowest) r.slowest = ms;
            r.calls++;
            return;
        }
    }
    g_rows.push_back(Row{stage, on_device, flops, ms, ms, 1});
}

void record_time(const char* stage, const bool on_device, const double ms)
{
    record(stage, on_device, 0.0, ms);
}

void reset()
{
    std::lock_guard<std::mutex> lock(g_mutex);
    g_rows.clear();
}

void report(std::ostream& out)
{
    std::vector<Row> rows;
    {
        std::lock_guard<std::mutex> lock(g_mutex);
        if (!g_enabled || g_rows.empty())
            return;
        rows = g_rows;
    }

    out << "\n\n------------------------- Arithmetic throughput -------------------------\n";
    //A CPU row that looks like slow hardware is usually a thread count of one, and that is
    //easier to read than to deduce.
#ifdef _OPENMP
    int actual = 1;
#pragma omp parallel
    {
#pragma omp single
        actual = omp_get_num_threads();
    }
    out << "host threads: " << actual << " in a parallel region, omp_get_max_threads() = "
        << omp_get_max_threads() << "\n";
#else
    out << "host threads: built without OpenMP\n";
#endif
    //slowest separates a slow stage from one paying a one-off: context creation lands in
    //whichever call touches the GPU first.
    out << std::left << std::setw(34) << "stage" << std::setw(6) << "where"
        << std::right << std::setw(8) << "calls" << std::setw(12) << "time/ms"
        << std::setw(12) << "GFLOP" << std::setw(12) << "GFLOP/s"
        << std::setw(12) << "slowest" << "\n";

    for (int i = 0; i < rows.size(); i++) {
        const Row& r = rows[i];
        //0.000 GFLOP would read as "does no arithmetic" rather than "nobody counted it"
        const bool counted = r.flops > 0.0;
        std::ostringstream gflop;
        if (counted) gflop << std::fixed << std::setprecision(3) << (r.flops / 1.0e9);
        else gflop << "-";
        out << std::left << std::setw(34) << r.stage.substr(0, 33)
            << std::setw(6) << (r.on_device ? "GPU" : "CPU")
            << std::right << std::setw(8) << r.calls
            << std::setw(12) << std::fixed << std::setprecision(1) << r.ms
            << std::setw(12) << gflop.str()
            << std::setw(12) << (counted ? rate(r.flops, r.ms) : std::string("-"))
            << std::setw(12) << std::fixed << std::setprecision(1) << r.slowest << "\n";
    }

    //Where a stage ran in both places the ratio is the whole answer, so do the division here
    //rather than leaving it to be done by eye against two rows that may be far apart.
    bool header = false;
    for (const Row& a : rows) {
        if (a.on_device)
            continue;
        for (const Row& b : rows) {
            if (!b.on_device || b.stage != a.stage || a.ms <= 0.0 || b.ms <= 0.0)
                continue;
            const double cpu_rate = a.flops / a.ms;
            const double gpu_rate = b.flops / b.ms;
            if (cpu_rate <= 0.0)
                continue;
            if (!header) {
                out << "\n";
                header = true;
            }
            out << "  " << a.stage << ": device is "
                << std::fixed << std::setprecision(2) << (gpu_rate / cpu_rate)
                << "x the host rate\n";
        }
    }

    out << "\n  Counting conventions are in throughput.h. The transform's sincos is counted\n"
        << "  as two operations on both sides, so its absolute rate is a convention while\n"
        << "  the CPU/GPU ratio is not. GEMM rows are 2mnk and need no such caveat.\n";
    out << "-------------------------------------------------------------------------\n";
}

} //namespace throughput
