#pragma once

#include "convenience.h"

#include <unordered_set>
using ScattererLabel = std::variant<std::string, std::uint64_t>;
using ScattererLabels = std::vector<ScattererLabel>;

template <typename T>
concept ScattererValue =
std::same_as<std::remove_cv_t<T>, std::string> ||
std::same_as<std::remove_cv_t<T>, std::uint64_t> ||
std::same_as<std::remove_cv_t<T>, ScattererLabel>;

template <typename numtype_index, typename numtype>
class tsc_block
{
private:
    cvec2 sf_;                                         // [scatterer][reflection]
    ScattererLabels scatterers_;
    std::vector<std::vector<numtype_index>> indices_; // [dimension][reflection]
    std::string header_;
    bool anomalous_dispersion_ = false;

    static constexpr const char* nan_message_ = "NaN in SF!";

    template <ScattererValue Label>
    static ScattererLabels make_scatterers(const std::vector<Label>& labels)
    {
        ScattererLabels result;
        result.reserve(labels.size());
        for (const auto& label : labels)
            result.emplace_back(label);
        return result;
    }

    static cvec2 copy_and_validate_sf(
        const std::vector<std::vector<numtype>>& given_sf)
    {
        cvec2 result(given_sf.size());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(given_sf.size()); ++i)
        {
            result[i].resize(given_sf[i].size());
            for (std::size_t j = 0; j < given_sf[i].size(); ++j)
            {
                err_checkf(!is_nan(given_sf[i][j]), nan_message_, std::cout);
                result[i][j] = given_sf[i][j];
            }
        }
        return result;
    }

    static std::vector<std::vector<numtype_index>> copy_indices(
        const std::vector<std::vector<numtype_index>>& given_indices)
    {
        return given_indices;
    }

    template <typename ReflectionRange>
    static std::vector<std::vector<numtype_index>> convert_indices(
        const ReflectionRange& given_indices)
    {
        std::vector<std::vector<numtype_index>> result(3);
        for (auto& dimension : result)
            dimension.reserve(given_indices.size());

        for (const auto& hkl : given_indices)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
            {
                result[dimension].push_back(
                    static_cast<numtype_index>(hkl[dimension]));
            }
        }
        return result;
    }

    template <ScattererValue Label, typename IndexInput>
    void initialize(
        const std::vector<std::vector<numtype>>& given_sf,
        const std::vector<Label>& given_scatterers,
        const IndexInput& given_indices,
        std::string given_header)
    {
        sf_ = copy_and_validate_sf(given_sf);
        scatterers_ = make_scatterers(given_scatterers);

        if constexpr (std::same_as<
            std::remove_cvref_t<IndexInput>,
            std::vector<std::vector<numtype_index>>>)
        {
            indices_ = copy_indices(given_indices);
        }
        else
        {
            indices_ = convert_indices(given_indices);
        }

        header_ = std::move(given_header);
        anomalous_dispersion_ = false;
    }

    [[nodiscard]] bool uses_ids() const
    {
        return !scatterers_.empty() &&
            std::holds_alternative<std::uint64_t>(scatterers_.front());
    }

    void validate_uniform_scatterer_type() const
    {
        if (scatterers_.empty())
            return;

        const std::size_t expected_index = scatterers_.front().index();
        const bool uniform = std::all_of(
            scatterers_.begin(), scatterers_.end(),
            [expected_index](const ScattererLabel& value)
            {
                return value.index() == expected_index;
            });

        if (!uniform)
            throw std::runtime_error("Scatterer list contains mixed labels and IDs");
    }

    void validate_dimensions() const
    {
        if (indices_.size() != 3)
            throw std::runtime_error("TSC index array must contain exactly three dimensions");

        if (sf_.size() != scatterers_.size())
            throw std::runtime_error("Scatterer and structure-factor counts differ");

        if (indices_[0].size() != indices_[1].size() ||
            indices_[1].size() != indices_[2].size())
        {
            throw std::runtime_error("TSC index dimensions have different lengths");
        }

        const std::size_t reflections = indices_[0].size();
        for (const auto& row : sf_)
        {
            if (row.size() != reflections)
                throw std::runtime_error("Structure-factor rows have inconsistent lengths");
        }
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
        check_stream(out, "Failed to write binary TSC value");
    }

    template <typename T>
    static T read_scalar(std::istream& in)
    {
        T value{};
        in.read(reinterpret_cast<char*>(&value), sizeof(T));
        check_stream(in, "Failed to read binary TSC value");
        return value;
    }

    static void write_bytes(
        std::ostream& out,
        const char* data,
        std::size_t size)
    {
        if (size == 0)
            return;
        out.write(data, static_cast<std::streamsize>(size));
        check_stream(out, "Failed to write binary TSC payload");
    }

    static std::string read_bytes(std::istream& in, std::size_t size)
    {
        std::string result(size, '\0');
        if (size != 0)
        {
            in.read(result.data(), static_cast<std::streamsize>(size));
            check_stream(in, "Failed to read binary TSC payload");
        }
        return result;
    }

    void write_scatterer_text(std::ostream& out) const
    {
        validate_uniform_scatterer_type();

        if (uses_ids())
        {
            out << "\nSCATTERER_IDS:" << std::hex;
            for (const auto& scatterer : scatterers_)
                out << ' ' << std::get<std::uint64_t>(scatterer);
            out << std::dec;
        }
        else
        {
            out << "\nSCATTERERS:";
            for (const auto& scatterer : scatterers_)
                out << ' ' << std::get<std::string>(scatterer);
        }
    }

    void write_binary_scatterers(std::ostream& out) const
    {
        validate_uniform_scatterer_type();

        std::string payload;
        int payload_size;
        if (uses_ids())
        {
            payload.resize(scatterers_.size() * sizeof(std::uint64_t));
            char* destination = payload.data();
            for (const auto& scatterer : scatterers_)
            {
                const std::uint64_t id = std::get<std::uint64_t>(scatterer);
                std::memcpy(destination, &id, sizeof(id));
                destination += sizeof(id);
            }
            payload_size = static_cast<int>(scatterers_.size());
        }
        else
        {
            for (std::size_t i = 0; i < scatterers_.size(); ++i)
            {
                payload += std::get<std::string>(scatterers_[i]);
                if (i + 1 < scatterers_.size()) {
                    payload.push_back(' ');
                }
            }
            payload_size =
                static_cast<int>(payload.size());
        }
        write_scalar(out, payload_size);
        write_bytes(out, payload.data(), payload.size());
    }

    void read_binary_scatterers(std::istream& in)
    {
        const int count = read_scalar<int>(in);
        scatterers_.clear();

        const bool id_payload =
            header_.find("SCATTERER_IDS") != std::string::npos;

        if (id_payload)
        {
            const std::string payload = read_bytes(in, count * sizeof(std::uint64_t));
            scatterers_.reserve(count);
            for (std::size_t i = 0; i < count; ++i)
            {
                std::uint64_t id{};
                std::memcpy(
                    &id,
                    payload.data() + i * sizeof(std::uint64_t),
                    sizeof(id));
                scatterers_.emplace_back(id);
            }
        }
        else
        {
            const std::string payload = read_bytes(in, count * sizeof(char));
            const svec labels = split_string<std::string>(payload, " ");
            scatterers_.insert(scatterers_.end(), labels.begin(), labels.end());
        }
    }

    template <typename IndexWriter>
    void write_text_file(
        const std::filesystem::path& cif,
        const std::filesystem::path& name,
        IndexWriter&& write_index) const
    {
        validate_dimensions();

        std::ofstream out(name);
        check_stream(out, "Failed to open TSC file for writing");

        out << "TITLE: " << cif.stem().string()
            << "\nSYMM: expanded";
        if (anomalous_dispersion_)
            out << "\nAD: TRUE";

        write_scatterer_text(out);
        out << "\nDATA:\n";

        for (std::size_t reflection = 0; reflection < reflection_size(); ++reflection)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
            {
                write_index(out, indices_[dimension][reflection]);
                out << ' ';
            }

            for (const auto& scatterer_sf : sf_)
            {
                out << std::scientific << std::setprecision(8)
                    << std::real(scatterer_sf[reflection]) << ','
                    << std::imag(scatterer_sf[reflection]) << ' ';
            }
            out << '\n';
        }
    }

    void verify_compatible_reflections(
        const tsc_block& rhs,
        std::ostream& log) const
    {
        err_checkf(
            reflection_size() == rhs.reflection_size(),
            "Inconsistent number of reflections!",
            log);

        for (std::size_t reflection = 0;
            reflection < reflection_size();
            ++reflection)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
            {
                err_checkf(
                    indices_[dimension][reflection] ==
                    rhs.indices_[dimension][reflection],
                    "Mismatch in indices in append!",
                    log);
            }
        }
    }

    template <typename Block>
    void append_impl(Block&& rhs, std::ostream& log)
    {
        if (reflection_size() == 0)
        {
            *this = std::forward<Block>(rhs);
            return;
        }

        verify_compatible_reflections(rhs, log);

        std::unordered_set<ScattererLabel> existing(
            scatterers_.begin(), scatterers_.end());

        for (std::size_t i = 0; i < rhs.scatterers_.size(); ++i)
        {
            if (!existing.insert(rhs.scatterers_[i]).second)
                continue;

            if constexpr (std::is_rvalue_reference_v<Block&&>)
            {
                scatterers_.push_back(std::move(rhs.scatterers_[i]));
                sf_.push_back(std::move(rhs.sf_[i]));
            }
            else
            {
                scatterers_.push_back(rhs.scatterers_[i]);
                sf_.push_back(rhs.sf_[i]);
            }
        }

        validate_dimensions();
    }

public:
    tsc_block() = default;
    ~tsc_block() = default;
    tsc_block(const tsc_block&) = default;
    tsc_block(tsc_block&&) noexcept = default;
    tsc_block& operator=(const tsc_block&) = default;
    tsc_block& operator=(tsc_block&&) noexcept = default;

    template <ScattererValue Label>
    tsc_block(
        const std::vector<std::vector<numtype>>& given_sf,
        const std::vector<Label>& given_scatterers,
        const std::vector<std::vector<numtype_index>>& given_indices,
        std::string given_header = {})
    {
        initialize(
            given_sf,
            given_scatterers,
            given_indices,
            std::move(given_header));
    }

    template <ScattererValue Label>
    tsc_block(
        const std::vector<std::vector<numtype>>& given_sf,
        const std::vector<Label>& given_scatterers,
        const hkl_list& given_indices,
        std::string given_header = {})
    {
        initialize(
            given_sf,
            given_scatterers,
            given_indices,
            std::move(given_header));
    }

    template <ScattererValue Label>
    tsc_block(
        const std::vector<std::vector<numtype>>& given_sf,
        const std::vector<Label>& given_scatterers,
        const hkl_list_d& given_indices,
        std::string given_header = {})
    {
        initialize(
            given_sf,
            given_scatterers,
            given_indices,
            std::move(given_header));
    }

    explicit tsc_block(const std::filesystem::path& file_name)
    {
        std::ifstream in(file_name, std::ios::binary);
        check_stream(in, "Failed to open TSCB file");

        const std::uint64_t header_size = read_scalar<std::uint64_t>(in);
        header_ = read_bytes(in, static_cast<std::size_t>(header_size));

        read_binary_scatterers(in);

        const std::uint64_t reflection_count = read_scalar<std::uint64_t>(in);
        indices_.assign(3, {});
        for (auto& dimension : indices_)
            dimension.reserve(static_cast<std::size_t>(reflection_count));

        sf_.assign(scatterers_.size(), {});
        for (auto& row : sf_)
            row.reserve(static_cast<std::size_t>(reflection_count));

        for (std::uint64_t reflection = 0;
            reflection < reflection_count;
            ++reflection)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
                indices_[dimension].push_back(read_scalar<numtype_index>(in));

            for (auto& row : sf_)
                row.push_back(read_scalar<std::complex<double>>(in));
        }

        validate_dimensions();
    }

    [[nodiscard]] const cvec& get_sf_for_scatterer(std::size_t nr) const
    {
        return sf_.at(nr);
    }

    [[nodiscard]] const cvec& get_sf_for_scatterer(
        std::size_t nr,
        std::ostream& log) const
    {
        err_checkf(nr < sf_.size(), "Wrong number in get SF", log);
        return sf_.at(nr);
    }

    [[nodiscard]] const ScattererLabel& get_scatterer(std::size_t nr) const
    {
        return scatterers_.at(nr);
    }

    [[nodiscard]] const ScattererLabel& get_scatterer(
        std::size_t nr,
        std::ostream& log) const
    {
        err_checkf(nr < scatterers_.size(), "Invalid nr of scatterer", log);
        return scatterers_.at(nr);
    }

    [[nodiscard]] const ScattererLabels& get_scatterers() const noexcept
    {
        return scatterers_;
    }

    template <typename ReturnType>
        requires (
    std::same_as<ReturnType, std::string> ||
        std::same_as<ReturnType, std::uint64_t>)
        [[nodiscard]] std::vector<ReturnType> get_scatterers_as() const
    {
        std::vector<ReturnType> result;
        result.reserve(scatterers_.size());
        for (const auto& scatterer : scatterers_)
            result.push_back(std::get<ReturnType>(scatterer));
        return result;
    }


    [[nodiscard]] std::vector<std::string> get_scatterers_string() const
    {
        std::vector<std::string> result;
        result.reserve(scatterers_.size());
        for (const auto& scatterer : scatterers_)
            if (!uses_ids())
                result.push_back(std::get<std::string>(scatterer));
            else {
                result.push_back(std::to_string(std::get<uint64_t>(scatterer)));
            }
        return result;
    }

    void set_AD(bool value) noexcept { anomalous_dispersion_ = value; }
    [[nodiscard]] bool get_AD() const noexcept { return anomalous_dispersion_; }

    [[nodiscard]] std::array<numtype_index, 3> get_indices(
        std::size_t reflection) const
    {
        return {
            indices_.at(0).at(reflection),
            indices_.at(1).at(reflection),
            indices_.at(2).at(reflection) };
    }

    [[nodiscard]] numtype_index get_index(
        std::size_t dimension,
        std::size_t reflection) const
    {
        return indices_.at(dimension).at(reflection);
    }

    [[nodiscard]] bool is_empty() const noexcept
    {
        return sf_.empty() && scatterers_.empty() && indices_.empty();
    }

    [[nodiscard]] std::size_t scatterer_size() const noexcept
    {
        return sf_.size() == scatterers_.size() ? scatterers_.size() : 0;
    }

    [[nodiscard]] std::size_t reflection_size() const noexcept
    {
        if (sf_.empty() || indices_.size() != 3)
            return 0;

        const std::size_t size = indices_[0].size();
        if (indices_[1].size() != size || indices_[2].size() != size)
            return 0;

        return std::all_of(
            sf_.begin(), sf_.end(),
            [size](const auto& row) { return row.size() == size; })
            ? size
            : 0;
    }

    [[nodiscard]] const std::vector<std::vector<numtype_index>>&
        get_index_vector() const noexcept
    {
        return indices_;
    }

    void append(const tsc_block& rhs, std::ostream& log)
    {
        append_impl(rhs, log);
    }

    void append(tsc_block&& rhs, std::ostream& log)
    {
        append_impl(std::move(rhs), log);
    }

    void write_tsc_file(
        const std::filesystem::path& cif,
        const std::filesystem::path& name = "experimental.tsc") const
    {
        write_text_file(
            cif,
            name,
            [](std::ostream& out, numtype_index value)
            {
                out << value;
            });
    }

    void write_tsc_file_non_integer(
        const std::filesystem::path& cif,
        const std::filesystem::path& name = "experimental.tsc") const
    {
        write_text_file(
            cif,
            name,
            [](std::ostream& out, numtype_index value)
            {
                out << std::fixed << std::setprecision(3) << value;
            });
    }

    [[nodiscard]] std::string get_tsc_cif_block(const options& opt) const
    {
        std::string partition = "Unknown";
        switch (opt.partition_type)
        {
        case PartitionType::Becke:      partition = "Becke"; break;
        case PartitionType::TFVC:       partition = "TFVC"; break;
        case PartitionType::Hirshfeld:  partition = "Hirshfeld"; break;
        case PartitionType::RI:         partition = "RI-Fit"; break;
        case PartitionType::MBIS:       partition = "MBIS"; break;
        case PartitionType::EMBIS:      partition = "EMBIS"; break;
        default: break;
        }

        std::ostringstream result;
        result << "_aspheric_ffs_partitioning.name     '" << partition << "'\n"
            << "_aspheric_ffs_partitioning.software 'NoSpherA2'\n"
            << "_aspheric_ffs_partitioning.source   "
            "'partitioned molecular wavefunction calculation'\n\n"
            << "loop_\n"
            << "_aspheric_ff.index_h\n"
            << "_aspheric_ff.index_k\n"
            << "_aspheric_ff.index_l\n"
            << "_aspheric_ff.form_factor_real\n"
            << "_aspheric_ff.form_factor_imag\n ";

        for (std::size_t reflection = 0;
            reflection < reflection_size();
            ++reflection)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
                result << indices_[dimension][reflection] << ' ';

            result << "'[";
            for (const auto& row : sf_)
                result << std::real(row[reflection]) << ' ';

            result << "]' '[";
            for (const auto& row : sf_)
                result << std::imag(row[reflection]) << ' ';

            result << "]'\n ";
        }
        return result.str();
    }

    void write_tscb_file(
        const std::filesystem::path & /*cif_name*/ = "test.cif",
        const std::filesystem::path& name = "experimental.tscb") const
    {
        validate_dimensions();
        validate_uniform_scatterer_type();

        std::ofstream out(name, std::ios::binary | std::ios::trunc);
        check_stream(out, "Failed to open TSCB file for writing");

        const int header_size =
            static_cast<int>(header_.size());
        write_scalar(out, header_size);
        write_bytes(out, header_.data(), header_.size());

        write_binary_scatterers(out);

        const int reflections =
            static_cast<int>(reflection_size());
        write_scalar(out, reflections);

        for (std::size_t reflection = 0;
            reflection < reflection_size();
            ++reflection)
        {
            for (std::size_t dimension = 0; dimension < 3; ++dimension)
                write_scalar(out, indices_[dimension][reflection]);

            for (const auto& row : sf_)
                write_scalar(out, row[reflection]);
        }
    }
};

// Merge implementation shared by merge_tscs() and merge_tscs_without_checks().
namespace tsc_merge_detail
{
    enum class ReflectionPolicy
    {
        MatchIncludingFriedelMates,
        AssumeIdenticalOrder
    };

    struct ParsedHeader
    {
        std::vector<ScattererLabel> file_labels;
        std::vector<std::size_t> global_label_indices;
        std::size_t new_label_count = 0;
        std::string header;
    };

    struct ReflectionKey
    {
        int h{};
        int k{};
        int l{};
    };

    inline bool is_zero(const ReflectionKey& key) noexcept
    {
        return key.h == 0 && key.k == 0 && key.l == 0;
    }

    inline bool equals(const ReflectionKey& lhs, const ReflectionKey& rhs) noexcept
    {
        return lhs.h == rhs.h && lhs.k == rhs.k && lhs.l == rhs.l;
    }

    inline bool equals_friedel_mate(
        const ReflectionKey& lhs,
        const ReflectionKey& rhs) noexcept
    {
        return lhs.h == -rhs.h && lhs.k == -rhs.k && lhs.l == -rhs.l;
    }

    inline std::filesystem::path ensure_text_tsc(
        const std::filesystem::path& input)
    {
        if (input.extension() != ".tscb")
            return input;

        std::filesystem::path output = input;
        output.replace_extension(".tsc");

        std::cout << "Converting to: " << output << '\n';
        tsc_block<int, cdouble> block(input);
        block.write_tsc_file(output, output);
        return output;
    }

    inline ScattererLabel parse_label(
        const std::string& token,
        bool ids)
    {
        if (ids)
            return static_cast<std::uint64_t>(std::stoull(token, nullptr, 16));
        return token;
    }

    inline ParsedHeader read_header(
        std::istream& input,
        ScattererLabels& global_labels,
        bool keep_header)
    {
        ParsedHeader result;
        std::unordered_map<ScattererLabel, std::size_t> known_labels;
        known_labels.reserve(global_labels.size());

        for (std::size_t i = 0; i < global_labels.size(); ++i)
            known_labels.emplace(global_labels[i], i);

        std::string line;
        bool found_data = false;

        while (std::getline(input, line))
        {
            const bool is_ids = line.find("SCATTERER_IDS:") != std::string::npos;
            const bool is_names = line.find("SCATTERERS:") != std::string::npos;

            if (is_ids || is_names)
            {
                std::istringstream labels_input(line.substr(line.find(':') + 1));
                std::string token;

                while (labels_input >> token)
                {
                    ScattererLabel label = parse_label(token, is_ids);
                    result.file_labels.push_back(label);

                    auto [position, inserted] = known_labels.emplace(
                        label,
                        global_labels.size());

                    if (inserted)
                    {
                        global_labels.push_back(label);
                        ++result.new_label_count;
                    }

                    result.global_label_indices.push_back(position->second);
                }
                continue;
            }

            if (line.find("DATA:") != std::string::npos)
            {
                found_data = true;
                break;
            }

            if (keep_header && line.find("AD: ") == std::string::npos)
            {
                result.header += line;
                result.header.push_back('\n');
            }
        }

        if (!found_data)
            throw std::runtime_error("TSC file contains no DATA section");

        if (result.file_labels.empty())
            throw std::runtime_error("TSC file contains no scatterers");

        return result;
    }

    inline bool parse_data_row(
        const std::string& line,
        ReflectionKey& key,
        std::vector<std::complex<double>>& values)
    {
        std::istringstream input(line);
        if (!(input >> key.h >> key.k >> key.l))
            return false;

        values.clear();
        std::string token;
        while (input >> token)
        {
            const auto comma = token.find(',');
            if (comma == std::string::npos)
                throw std::runtime_error("Invalid TSC complex value: " + token);

            values.emplace_back(
                std::stod(token.substr(0, comma)),
                std::stod(token.substr(comma + 1)));
        }

        return true;
    }

    inline std::optional<std::pair<std::size_t, bool>> find_reflection(
        const ivec2& indices,
        const ReflectionKey& key)
    {
        for (std::size_t i = 0; i < indices[0].size(); ++i)
        {
            const ReflectionKey existing{
                indices[0][i], indices[1][i], indices[2][i] };

            if (equals(key, existing))
                return std::pair{ i, false };

            if (equals_friedel_mate(key, existing))
                return std::pair{ i, true };
        }

        return std::nullopt;
    }

    inline void append_reflection(
        ivec2& indices,
        const ReflectionKey& key)
    {
        indices[0].push_back(key.h);
        indices[1].push_back(key.k);
        indices[2].push_back(key.l);
    }

    inline void assign_values(
        cvec2& form_factors,
        const ParsedHeader& parsed,
        const std::vector<std::complex<double>>& values,
        std::size_t reflection,
        bool conjugate)
    {
        if (values.size() != parsed.file_labels.size())
            throw std::runtime_error("Scatterer and form-factor counts differ in TSC row");

        for (std::size_t local = 0; local < values.size(); ++local)
        {
            const std::size_t global = parsed.global_label_indices[local];
            if (form_factors[global].size() <= reflection)
                form_factors[global].resize(reflection + 1);

            form_factors[global][reflection] =
                conjugate ? std::conj(values[local]) : values[local];
        }
    }

    inline void read_checked_data(
        std::istream& input,
        std::size_t file_index,
        const ParsedHeader& parsed,
        ivec2& indices,
        cvec2& form_factors)
    {
        std::string line;
        ReflectionKey key;
        std::vector<std::complex<double>> values;

        while (std::getline(input, line))
        {
            if (!parse_data_row(line, key, values) || is_zero(key))
                continue;

            if (file_index == 0)
            {
                if (find_reflection(indices, key))
                    continue;

                const std::size_t reflection = indices[0].size();
                append_reflection(indices, key);
                for (auto& row : form_factors)
                    row.resize(reflection + 1);
                assign_values(form_factors, parsed, values, reflection, false);
                continue;
            }

            const auto match = find_reflection(indices, key);
            if (!match)
                continue;

            assign_values(
                form_factors,
                parsed,
                values,
                match->first,
                match->second);
        }
    }

    inline void read_unchecked_data(
        std::istream& input,
        std::size_t file_index,
        const ParsedHeader& parsed,
        ivec2& indices,
        cvec2& form_factors)
    {
        std::string line;
        ReflectionKey key;
        std::vector<std::complex<double>> values;
        std::size_t reflection = 0;

        while (std::getline(input, line))
        {
            if (!parse_data_row(line, key, values))
                continue;

            if (file_index == 0)
            {
                append_reflection(indices, key);
                for (auto& row : form_factors)
                    row.resize(reflection + 1);
            }
            else if (reflection >= indices[0].size())
            {
                throw std::runtime_error(
                    "A merged TSC file contains more reflections than the first file");
            }

            assign_values(form_factors, parsed, values, reflection, false);
            ++reflection;
        }

        if (file_index != 0 && reflection != indices[0].size())
        {
            throw std::runtime_error(
                "A merged TSC file contains a different number of reflections");
        }
    }

    inline bool merge_tscs_impl(
        const std::string& mode,
        const pathvec& files,
        bool old_tsc,
        ReflectionPolicy policy)
    {
        if (mode.empty() || files.empty())
            return false;

        ScattererLabels labels;
        ivec2 indices(3);
        cvec2 form_factors;
        std::string header;

        try
        {
            for (std::size_t file_index = 0; file_index < files.size(); ++file_index)
            {
                std::cout << "Reading file number " << file_index + 1
                    << ": " << files[file_index] << '\n';

                const std::filesystem::path input_path =
                    ensure_text_tsc(files[file_index]);

                std::ifstream input(input_path);
                if (!input)
                    throw std::runtime_error("Failed to open TSC file: " + input_path.string());

                const std::size_t old_label_count = labels.size();
                ParsedHeader parsed = read_header(
                    input,
                    labels,
                    file_index == 0);

                if (file_index == 0)
                    header = parsed.header;

                form_factors.resize(labels.size());
                for (std::size_t i = old_label_count; i < form_factors.size(); ++i)
                    form_factors[i].resize(indices[0].size());

                std::cout << "Read " << parsed.file_labels.size()
                    << " atoms, " << parsed.new_label_count
                    << " are new.\nReading data block...\n";

                if (policy == ReflectionPolicy::MatchIncludingFriedelMates)
                {
                    read_checked_data(
                        input, file_index, parsed, indices, form_factors);
                }
                else
                {
                    read_unchecked_data(
                        input, file_index, parsed, indices, form_factors);
                }

                std::cout << "Data for " << indices[0].size()
                    << " indices read.\n";
            }

            std::cout << "Writing combined file...\n";
            tsc_block<int, cdouble> combined(
                form_factors,
                labels,
                indices,
                header);

            if (old_tsc)
                combined.write_tsc_file("combined", "combined.tsc");
            else
                combined.write_tscb_file("combined.cif", "combined.tscb");

            std::cout << "Done!\n";
            return true;
        }
        catch (const std::exception& error)
        {
            std::cerr << "Could not merge TSC files: " << error.what() << '\n';
            return false;
        }
    }
} // namespace tsc_merge_detail

inline bool merge_tscs(
    const std::string& mode,
    const pathvec& files,
    bool old_tsc)
{
    return tsc_merge_detail::merge_tscs_impl(
        mode,
        files,
        old_tsc,
        tsc_merge_detail::ReflectionPolicy::MatchIncludingFriedelMates);
}

inline bool merge_tscs_without_checks(
    const std::string& mode,
    const pathvec& files,
    bool old_tsc)
{
    return tsc_merge_detail::merge_tscs_impl(
        mode,
        files,
        old_tsc,
        tsc_merge_detail::ReflectionPolicy::AssumeIdenticalOrder);
}
