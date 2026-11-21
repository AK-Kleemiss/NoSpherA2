#include "pch.h"
#include "convenience.h"

bool ProgressBar::report_counts = false;

ProgressBar::~ProgressBar()
{
	progress_ = 100.0f;
	write_progress();
	stream_ << std::endl;
	if (report_counts)
	{
		// updates is how often callers reported an item; writes is how often that
		// needed the omp critical. Every item used to take the lock, so writes
		// being far below updates is the whole point of the change.
		const unsigned long long u = update_calls_.load(), w = bar_writes_.load();
		stream_ << "[progress] " << status_text_ << ": " << u << " updates, "
				<< w << " serialised writes";
		if (u > 0)
			stream_ << "  (" << (w * 100.0 / u) << "% of calls took the lock, "
					<< (u > w ? u / (w ? w : 1) : 1) << "x fewer barriers)";
		stream_ << std::endl;
	}
#ifdef _WIN32
	if (taskbarList_)
	{
		taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
		taskbarList_->Release();
		taskbarList_ = nullptr;
	}
#endif
}

void ProgressBar::write_progress()
{
	// No need to write once progress is 100%
	if (progress_ > 100.0f)
		return;

	// Move cursor to the first position on the same line
	// Check if os is a file stream
	if (dynamic_cast<std::filebuf *>(stream_.rdbuf()))
	{
		// Rewinding would overwrite what the loop printed since the last redraw; if the put position
		// moved, keep that text and re-anchor the bar behind it.
		if (stream_.tellp() != barend_)
			linestart = stream_.tellp();
		stream_.seekp(linestart); // Is a file stream
	}
	else
	{
		stream_ << "\r" << std::flush; // Is not a file stream
	}

	// Start bar
	stream_ << "[";

	const auto completed = static_cast<size_t>(progress_ * static_cast<float>(bar_width_) / 100.0);
	for (size_t i = 0; i <= completed; i++)
	{
		stream_ << fill_;
	}

	// End bar
	if (((progress_ < 100.0f) ? progress_ : 100.0f) == 100)
	{
		stream_ << "] 100% " << std::flush;
		barend_ = stream_.tellp();
#ifdef _WIN32
		if (taskbarList_)
		{
			taskbarList_->SetProgressValue(GetConsoleWindow(), 100, 100);
			taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
		}
#endif
		return;
	}

	stream_ << std::flush;
	barend_ = stream_.tellp();

#ifdef _WIN32
	// Update taskbar progress
	if (taskbarList_)
	{
		taskbarList_->SetProgressValue(GetConsoleWindow(), static_cast<ULONGLONG>(progress_), 100);
	}
#endif
}

#ifdef _WIN32
void ProgressBar::initialize_taskbar_progress()
{
	if (SUCCEEDED(CoInitialize(nullptr)))
	{
		if (SUCCEEDED(CoCreateInstance(CLSID_TaskbarList, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&taskbarList_))))
		{
			taskbarList_->HrInit();
			taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NORMAL);
		}
	}
}
#endif
