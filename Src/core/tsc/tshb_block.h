#pragma once

#include "tsc_block.h"

#include <array>
#include <cstdint>
#include <set>

//One scatterer of a .tshb: its density at full occupancy as primitive shells n e^(-alpha R^2) R^l sum_m c_m Y_lm (real
//Y_lm of spherical_harmonic.h, angstrom, cartesian frame, n folded into c) and spherical records. The label is for the
//reader's messages; the block keys a record by its ScattererLabel, as tsc_block keys a row
struct tshb_record
{
	struct shell
	{
		std::int32_t l = 0;
		double alpha = 0.0;
		vec c; // 2l + 1 coefficients, e/A^3
		bool operator==(const shell&) const = default;
	};
	std::string label;
	std::int32_t part = 0;
	std::int32_t Z = 0;
	std::vector<shell> shells;
	std::vector<std::array<double, 2>> gaussians; // c e^(-a R^2): { c, a }
	std::vector<std::array<double, 3>> slaters;   // c R^n e^(-z R): { c, n, z }
	std::vector<std::array<double, 2>> f_terms;   // a e^(-b s) added to f: { a, b }
	bool operator==(const tshb_record&) const = default;
};

class tshb_block
{
private:
	std::vector<tshb_record> records_; // [scatterer]
	ScattererLabels scatterers_;
	std::string header_;               // lines beyond the ones write() derives from the data
	bool electron_diffraction_ = false;
	bool complete_ = true;             // every scatterer of the run is in the block; a read file says so in its header

	static constexpr std::uint32_t version_ = 4;

	//as tsc_block::make_scatterers: the first of a key stays, a later one is reported and dropped with its record
	template <ScattererValue Label>
	void make_scatterers(std::vector<tshb_record>&& records, const std::vector<Label>& keys)
	{
		if (records.size() != keys.size())
			throw std::runtime_error("TSHB records and scatterers differ in number");
		std::unordered_set<ScattererLabel> seen;
		seen.reserve(keys.size());
		for (std::size_t i = 0; i < keys.size(); ++i)
		{
			const ScattererLabel key(keys[i]);
			if (!seen.insert(key).second)
			{
				std::cout << "Warning: Duplicate scatterer: ";
				std::visit([](const auto& value) { std::cout << value; }, key);
				std::cout << "\n";
				continue;
			}
			scatterers_.push_back(key);
			records_.push_back(std::move(records[i]));
		}
	}

	[[nodiscard]] bool uses_ids() const
	{
		return !scatterers_.empty() && std::holds_alternative<atomID>(scatterers_.front());
	}

	void validate_uniform_scatterer_type() const
	{
		if (scatterers_.empty())
			return;
		const std::size_t expected_index = scatterers_.front().index();
		for (const ScattererLabel& value : scatterers_)
			if (value.index() != expected_index)
				throw std::runtime_error("Scatterer list contains mixed labels and IDs");
	}

	static void check_stream(const std::ios& stream, const char* message)
	{
		if (!stream)
			throw std::runtime_error(message);
	}

	template <typename T>
	static void write_scalar(std::ostream& out, const T& value)
	{
		out.write(reinterpret_cast<const char*>(&value), sizeof(T));
		check_stream(out, "Failed to write binary TSHB value");
	}

	template <typename T>
	static T read_scalar(std::istream& in)
	{
		T value{};
		in.read(reinterpret_cast<char*>(&value), sizeof(T));
		check_stream(in, "Failed to read binary TSHB value");
		return value;
	}

	static void write_string(std::ostream& out, const std::string& s)
	{
		write_scalar(out, static_cast<std::uint32_t>(s.size()));
		out.write(s.data(), static_cast<std::streamsize>(s.size()));
		check_stream(out, "Failed to write binary TSHB string");
	}

	template <std::size_t N>
	static void write_terms(std::ostream& out, const std::vector<std::array<double, N>>& terms)
	{
		write_scalar(out, static_cast<std::uint32_t>(terms.size()));
		for (const auto& t : terms)
			for (const double v : t)
				write_scalar(out, v);
	}

	//the header write() emits: what the data says, then the free lines the block carries
	[[nodiscard]] std::string file_header() const
	{
		std::set<std::int32_t> parts;
		for (const tshb_record& r : records_)
			parts.insert(r.part);
		std::string h = "FRAME: cartesian\nUNITS: angstrom\nSOURCE: NoSpherA2";
		if (electron_diffraction_)
			h += "\nED: true";
		h += std::string("\nCOVERAGE: ") + (complete_ ? "complete" : "partial") + "\nPARTS:";
		for (const std::int32_t p : parts)
			h += " " + std::to_string(p);
		if (!header_.empty())
			h += "\n" + header_;
		return h;
	}

	//keeps what file_header() does not derive, so a read block writes back the same header
	void parse_header(const std::string& h)
	{
		std::istringstream in(h);
		std::string line;
		//a file that does not say it is complete (versions 1-3 never do) is not taken for complete
		complete_ = false;
		while (getline_universal(in, line))
		{
			if (line == "ED: true")
				electron_diffraction_ = true;
			else if (line == "COVERAGE: complete")
				complete_ = true;
			else if (line.rfind("FRAME:", 0) == 0 || line.rfind("UNITS:", 0) == 0 || line.rfind("SOURCE:", 0) == 0 ||
				line.rfind("COVERAGE:", 0) == 0 || line.rfind("PARTS:", 0) == 0)
				continue;
			else
				header_ += (header_.empty() ? "" : "\n") + line;
		}
	}

	template <typename Block>
	void append_impl(Block&& rhs, std::ostream& log)
	{
		if (rhs.scatterers_.empty())
		{
			complete_ = complete_ && rhs.complete_;
			return;
		}
		if (scatterers_.empty())
		{
			const bool was_complete = complete_;
			*this = std::forward<Block>(rhs);
			complete_ = complete_ && was_complete;
			return;
		}
		err_checkf(electron_diffraction_ == rhs.electron_diffraction_, "Appending an ED .tshb block to an X-ray one!", log);
		complete_ = complete_ && rhs.complete_;

		std::unordered_set<ScattererLabel> existing(scatterers_.begin(), scatterers_.end());
		for (std::size_t i = 0; i < rhs.scatterers_.size(); ++i)
		{
			if (!existing.insert(rhs.scatterers_[i]).second)
			{
				std::cout << "Warning: Duplicate scatterer in append: ";
				std::visit([](const auto& value) { std::cout << value; }, rhs.scatterers_[i]);
				std::cout << "\n";
				continue;
			}

			if constexpr (std::is_rvalue_reference_v<Block&&>)
			{
				scatterers_.push_back(std::move(rhs.scatterers_[i]));
				records_.push_back(std::move(rhs.records_[i]));
			}
			else
			{
				scatterers_.push_back(rhs.scatterers_[i]);
				records_.push_back(rhs.records_[i]);
			}
		}
		validate_uniform_scatterer_type();
	}

public:
	tshb_block() = default;
	~tshb_block() = default;
	tshb_block(const tshb_block&) = default;
	tshb_block(tshb_block&&) noexcept = default;
	tshb_block& operator=(const tshb_block&) = default;
	tshb_block& operator=(tshb_block&&) noexcept = default;

	template <ScattererValue Label>
	tshb_block(
		std::vector<tshb_record> given_records,
		const std::vector<Label>& given_scatterers,
		const bool electron_diffraction = false,
		std::string given_header = {})
		: header_(std::move(given_header)), electron_diffraction_(electron_diffraction)
	{
		make_scatterers(std::move(given_records), given_scatterers);
		validate_uniform_scatterer_type();
	}

	//versions 1-3 carry labels alone and are keyed by label: version 1 the shells, 2 the spherical records, 3 Z
	explicit tshb_block(const std::filesystem::path& file_name)
	{
		std::ifstream in(file_name, std::ios::binary | std::ios::ate);
		check_stream(in, "Failed to open TSHB file");
		const std::uint64_t file_size = static_cast<std::uint64_t>(in.tellg());
		in.seekg(0);
		//every count and length is held against the bytes left before anything is allocated: a corrupt or truncated file
		//is an error, never a crash or a gigabyte vector
		const auto left = [&]() { return file_size - static_cast<std::uint64_t>(in.tellg()); };
		const auto bound = [&](const std::uint64_t n, const std::uint64_t bytes_each, const char* what)
		{
			if (n * bytes_each > left())
				throw std::runtime_error(std::string("Corrupt or truncated TSHB file ") + file_name.string() + ": " + what + " runs past the end");
			return n;
		};
		const auto string_ = [&]()
		{
			const std::uint64_t n = bound(read_scalar<std::uint32_t>(in), 1, "a string");
			std::string s(n, '\0');
			if (n) in.read(s.data(), static_cast<std::streamsize>(n));
			check_stream(in, "Failed to read binary TSHB string");
			return s;
		};
		const auto terms = [&]<std::size_t N>(std::vector<std::array<double, N>>& t)
		{
			t.resize(bound(read_scalar<std::uint32_t>(in), 8 * N, "a term list"));
			for (auto& term : t)
				for (double& v : term)
					v = read_scalar<double>(in);
		};
		if (file_size < 8)
			throw std::runtime_error(file_name.string() + " is too short for a TSHB file");
		char magic[4];
		in.read(magic, 4);
		check_stream(in, "Failed to read TSHB magic");
		if (std::string(magic, 4) != "TSHB")
			throw std::runtime_error(file_name.string() + " is not a TSHB file");
		const std::uint32_t version = read_scalar<std::uint32_t>(in);
		if (version < 1 || version > version_)
			throw std::runtime_error("Unsupported TSHB version " + std::to_string(version) + " in " + file_name.string());
		parse_header(string_());
		//the smallest record: label length and shell count, the tag and PART (4), the term counts (2), Z (3)
		const std::uint64_t min_record = 8 + (version >= 4 ? 5 : 0) + (version >= 2 ? 12 : 0) + (version >= 3 ? 4 : 0);
		const std::uint64_t n = bound(read_scalar<std::uint32_t>(in), min_record, "the record count");
		std::vector<tshb_record> records(n);
		ScattererLabels keys;
		keys.reserve(n);
		for (tshb_record& r : records)
		{
			r.label = string_();
			if (version >= 4)
			{
				const std::uint8_t tag = read_scalar<std::uint8_t>(in);
				if (tag == 1)
					keys.emplace_back(read_scalar<atomID>(in));
				else if (tag == 0)
					keys.emplace_back(r.label);
				else
					throw std::runtime_error("Unknown TSHB key tag " + std::to_string(tag) + " in " + file_name.string());
				r.part = read_scalar<std::int32_t>(in);
			}
			else
				keys.emplace_back(r.label);
			r.shells.resize(bound(read_scalar<std::uint32_t>(in), 20, "the shell count"));
			for (tshb_record::shell& s : r.shells)
			{
				s.l = read_scalar<std::int32_t>(in);
				if (s.l < 0 || s.l > constants::ASSOCIATED_LEGENDRE_MAX_L)
					throw std::runtime_error("Corrupt TSHB file " + file_name.string() + ": shell with l = " + std::to_string(s.l));
				s.alpha = read_scalar<double>(in);
				s.c.resize(bound(2 * s.l + 1, 8, "shell coefficients"));
				for (double& c : s.c)
					c = read_scalar<double>(in);
			}
			if (version >= 2)
			{
				terms(r.gaussians);
				terms(r.slaters);
				terms(r.f_terms);
			}
			if (version >= 3)
				r.Z = read_scalar<std::int32_t>(in);
		}
		if (in.peek() != std::char_traits<char>::eof())
			throw std::runtime_error("Unexpected trailing bytes in TSHB file " + file_name.string());
		make_scatterers(std::move(records), keys);
		validate_uniform_scatterer_type();
	}

	[[nodiscard]] static tshb_block read(const std::filesystem::path& file_name)
	{
		return tshb_block(file_name);
	}

	//version 4: every record carries its key (label or the 16 bytes of the atomID SCATTERER_IDS holds) and its PART
	void write(const std::filesystem::path& name = "experimental.tshb") const
	{
		validate_uniform_scatterer_type();
		std::ofstream out(name, std::ios::binary | std::ios::trunc);
		check_stream(out, "Failed to open TSHB file for writing");
		out.write("TSHB", 4);
		write_scalar(out, version_);
		write_string(out, file_header());
		write_scalar(out, static_cast<std::uint32_t>(records_.size()));
		for (std::size_t i = 0; i < records_.size(); ++i)
		{
			const tshb_record& r = records_[i];
			write_string(out, r.label);
			if (std::holds_alternative<atomID>(scatterers_[i]))
			{
				write_scalar(out, std::uint8_t(1));
				write_scalar(out, std::get<atomID>(scatterers_[i]));
			}
			else
				write_scalar(out, std::uint8_t(0));
			write_scalar(out, r.part);
			write_scalar(out, static_cast<std::uint32_t>(r.shells.size()));
			for (const tshb_record::shell& s : r.shells)
			{
				if (s.c.size() != static_cast<std::size_t>(2 * s.l + 1))
					throw std::runtime_error("TSHB shell of " + r.label + " has the wrong number of coefficients");
				write_scalar(out, s.l);
				write_scalar(out, s.alpha);
				for (const double c : s.c)
					write_scalar(out, c);
			}
			write_terms(out, r.gaussians);
			write_terms(out, r.slaters);
			write_terms(out, r.f_terms);
			write_scalar(out, r.Z);
		}
	}

	void append(const tshb_block& rhs, std::ostream& log)
	{
		append_impl(rhs, log);
	}

	void append(tshb_block&& rhs, std::ostream& log)
	{
		append_impl(std::move(rhs), log);
	}

	tshb_block& operator+=(const tshb_block& rhs)
	{
		append_impl(rhs, std::cout);
		return *this;
	}

	tshb_block& operator+=(tshb_block&& rhs)
	{
		append_impl(std::move(rhs), std::cout);
		return *this;
	}

	[[nodiscard]] const tshb_record& get_record(std::size_t nr) const
	{
		return records_.at(nr);
	}

	[[nodiscard]] const ScattererLabel& get_scatterer(std::size_t nr) const
	{
		return scatterers_.at(nr);
	}

	[[nodiscard]] const ScattererLabels& get_scatterers() const noexcept
	{
		return scatterers_;
	}

	[[nodiscard]] std::vector<std::string> get_scatterers_string() const
	{
		std::vector<std::string> result;
		result.reserve(scatterers_.size());
		for (const auto& scatterer : scatterers_)
			if (!uses_ids())
				result.push_back(std::get<std::string>(scatterer));
			else
				result.push_back(std::get<atomID>(scatterer).to_hex_string());
		return result;
	}

	[[nodiscard]] std::size_t scatterer_size() const noexcept
	{
		return records_.size();
	}

	[[nodiscard]] bool is_empty() const noexcept
	{
		return records_.empty();
	}

	void set_ED(bool value) noexcept { electron_diffraction_ = value; }
	[[nodiscard]] bool get_ED() const noexcept { return electron_diffraction_; }
	void set_complete(bool value) noexcept { complete_ = value; }
	[[nodiscard]] bool is_complete() const noexcept { return complete_; }
	[[nodiscard]] const std::string& get_header() const noexcept { return header_; }
};

//-merge for .tshb: first-wins on the key as tsc_block::append, written once to output. COVERAGE stays complete only when
//every input says so and run_atoms (label and atomID of every scatterer of the run, from its CIF) are all among the merged
//keys; without run_atoms nothing tells whether the inputs make up the run, so the file is marked partial
inline bool merge_tshbs(const pathvec& files, const std::filesystem::path& output = "combined.tshb",
	const std::vector<std::pair<std::string, atomID>>* run_atoms = nullptr)
{
	if (files.empty())
		return false;
	try
	{
		tshb_block combined;
		for (std::size_t i = 0; i < files.size(); ++i)
		{
			std::cout << "Reading file number " << i + 1 << ": " << files[i] << '\n';
			const tshb_block part(files[i]);
			const std::size_t before = combined.scatterer_size();
			if (!part.is_complete())
				std::cout << files[i] << " is not marked complete\n";
			combined.append(part, std::cout);
			std::cout << "Read " << part.scatterer_size() << " atoms, " << combined.scatterer_size() - before << " are new.\n";
		}
		bool covered = false;
		if (run_atoms)
		{
			const ScattererLabels& keys = combined.get_scatterers();
			const std::unordered_set<ScattererLabel> have(keys.begin(), keys.end());
			const bool ids = !keys.empty() && std::holds_alternative<atomID>(keys.front());
			std::size_t missing = 0;
			for (const auto& [label, id] : *run_atoms)
				if (!have.count(ids ? ScattererLabel(id) : ScattererLabel(label)) && missing++ < 10)
					std::cout << "In no input: " << label << '\n';
			if (missing)
				std::cout << missing << " scatterer(s) of the CIF are in no input\n";
			covered = missing == 0;
		}
		else
			std::cout << "No CIF to check the inputs against\n";
		combined.set_complete(combined.is_complete() && covered);
		std::cout << "Writing combined file, COVERAGE: " << (combined.is_complete() ? "complete" : "partial") << "...\n";
		combined.write(output);
		std::cout << "Done!\n";
		return true;
	}
	catch (const std::exception& error)
	{
		std::cerr << "Could not merge TSHB files: " << error.what() << '\n';
		return false;
	}
}
