#pragma once
#include "string_io.h"

//-------------------------Progress_bar--------------------------------------------------
// LMS: My implementation of a progress bar, I would like it to stay within one line that is compatible with parallel loops
class ProgressBar
{
public:
	~ProgressBar();

	ProgressBar(const unsigned long long& worksize, const int& bar_width = 60, const std::string& fill = "#", const std::string& remainder = " ", const std::string& status_text = "", std::ostream& stream_ = std::cout)
		: worksize_(worksize), bar_width_(bar_width), fill_(fill), remainder_(remainder), status_text_(status_text), workdone(0), progress_(0.0f), workpart_(100.0f / worksize), percent_((worksize / 100 > 1) ? worksize / 100 : 1), stream_(stream_)
	{
		int bw = bar_width_ + 2;
		print_centered_text(status_text_, bw, stream_);
		linestart = stream_.tellp();
		barend_ = linestart;
#ifdef _WIN32
			initialize_taskbar_progress();
#endif
	}

	void set_progress()
	{
		progress_ = (float)workdone * workpart_;
	}

	//Called once per reflection; workdone is atomic, so only the call that crosses a reporting boundary takes the lock instead of serialising every worker on the redraw
	void update(const unsigned long long n = 1)
	{
		update_calls_.fetch_add(1, std::memory_order_relaxed);
		const unsigned long long before = workdone.fetch_add(n);
		const unsigned long long after = before + n;
		if (before / percent_ != after / percent_)
		{
#pragma omp critical
			{
				bar_writes_.fetch_add(1, std::memory_order_relaxed);
				set_progress();
				write_progress();
			}
		}
	}

	// How often callers asked, and how often that actually needed the lock.
	unsigned long long update_calls() const { return update_calls_.load(); }
	unsigned long long bar_writes() const { return bar_writes_.load(); }
	static bool report_counts;

	void write_progress();

private:
	std::ostream& stream_;
	const unsigned long long worksize_;
	const float workpart_;
	const unsigned long long percent_;
	int bar_width_;
	std::string fill_;
	std::string remainder_;
	std::string status_text_;
	std::atomic<unsigned long long> workdone;
	std::atomic<unsigned long long> update_calls_{0};
	std::atomic<unsigned long long> bar_writes_{0};
	float progress_;
	std::streampos linestart;
	//end of the bar's last write; a put position elsewhere means the loop printed in between
	std::streampos barend_{};
	bool finished_ = false;
#ifdef _WIN32
	//Assigned only inside initialize_taskbar_progress()'s SUCCEEDED checks; without the initialiser the destructor calls through stack garbage when COM refuses
	ITaskbarList3* taskbarList_ = nullptr;

	void initialize_taskbar_progress();
#endif
};
