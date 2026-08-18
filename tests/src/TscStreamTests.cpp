// The streamed writer must produce exactly the file the one-shot writer does.
// Everything downstream reads these bytes, so "close enough" is not a category
// that exists here: the test compares the two files byte for byte.
#include "pch.h"
#include <gtest/gtest.h>

#include "tsc_block.h"
#include "tsc_stream.h"

#include <filesystem>
#include <fstream>
#include <random>

namespace
{
    using block = tsc_block<int, cdouble>;

    block make_block(const int n_scatterers, const int n_reflections)
    {
        std::mt19937 rng(12345);          // fixed: the comparison must be repeatable
        std::uniform_real_distribution<double> uniform(-5.0, 5.0);

        std::vector<std::vector<cdouble>> sf(n_scatterers,
            std::vector<cdouble>(n_reflections));
        for (auto& row : sf)
            for (auto& value : row)
                value = cdouble(uniform(rng), uniform(rng));

        std::vector<std::string> labels;
        for (int i = 0; i < n_scatterers; ++i)
            labels.push_back("C" + std::to_string(i + 1));

        std::vector<std::vector<int>> idx(3, std::vector<int>(n_reflections));
        for (int r = 0; r < n_reflections; ++r)
        {
            idx[0][r] = r % 17 - 8;
            idx[1][r] = r % 11;
            idx[2][r] = r / 3;
        }
        return block(sf, labels, idx, "TITLE test\nSYMM x,y,z");
    }

    std::string slurp(const std::filesystem::path& p)
    {
        std::ifstream in(p, std::ios::binary);
        return std::string((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
    }
}

TEST(TscStream, StreamedFileMatchesOneShotFile)
{
    const int n_scatterers = 7;
    const int n_reflections = 1000;
    const block source = make_block(n_scatterers, n_reflections);

    const auto reference = std::filesystem::temp_directory_path() / "tsc_stream_ref.tscb";
    const auto streamed = std::filesystem::temp_directory_path() / "tsc_stream_out.tscb";
    source.write_tscb_file("ignored.cif", reference);

    // a block size that does not divide the reflection count, so the last
    // block is short - the case an even split would never exercise
    const std::size_t block_size = 128;
    {
        tsc_stream_writer<int, cdouble> writer(streamed, source, n_reflections, 2);
        std::size_t id = 0;
        for (std::size_t lo = 0; lo < static_cast<std::size_t>(n_reflections); lo += block_size)
        {
            const std::size_t hi = std::min<std::size_t>(lo + block_size, n_reflections);
            std::vector<std::vector<int>> idx(3, std::vector<int>(hi - lo));
            for (std::size_t d = 0; d < 3; ++d)
                for (std::size_t r = lo; r < hi; ++r)
                    idx[d][r - lo] = source.get_indices(r)[d];

            cvec2 chunk(n_scatterers, cvec(hi - lo));
            for (int s = 0; s < n_scatterers; ++s)
                for (std::size_t r = lo; r < hi; ++r)
                    chunk[s][r - lo] = source.get_sf_for_scatterer(s)[r];

            writer.submit(id++, std::move(idx), std::move(chunk));
        }
        writer.finish();
    }

    EXPECT_EQ(slurp(reference), slurp(streamed));
    std::filesystem::remove(reference);
    std::filesystem::remove(streamed);
}

TEST(TscStream, BlocksSubmittedOutOfOrderAreWrittenInOrder)
{
    const int n_scatterers = 3;
    const int n_reflections = 300;
    const block source = make_block(n_scatterers, n_reflections);

    const auto reference = std::filesystem::temp_directory_path() / "tsc_stream_ref2.tscb";
    const auto streamed = std::filesystem::temp_directory_path() / "tsc_stream_out2.tscb";
    source.write_tscb_file("ignored.cif", reference);

    const std::size_t block_size = 100;
    struct pending { std::size_t id; std::vector<std::vector<int>> idx; cvec2 sf; };
    std::vector<pending> blocks;
    std::size_t id = 0;
    for (std::size_t lo = 0; lo < static_cast<std::size_t>(n_reflections); lo += block_size)
    {
        const std::size_t hi = std::min<std::size_t>(lo + block_size, n_reflections);
        std::vector<std::vector<int>> idx(3, std::vector<int>(hi - lo));
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t r = lo; r < hi; ++r)
                idx[d][r - lo] = source.get_indices(r)[d];
        cvec2 chunk(n_scatterers, cvec(hi - lo));
        for (int s = 0; s < n_scatterers; ++s)
            for (std::size_t r = lo; r < hi; ++r)
                chunk[s][r - lo] = source.get_sf_for_scatterer(s)[r];
        blocks.push_back({id++, std::move(idx), std::move(chunk)});
    }
    std::reverse(blocks.begin(), blocks.end());   // hand them over backwards

    {
        tsc_stream_writer<int, cdouble> writer(streamed, source, n_reflections, blocks.size());
        for (auto& b : blocks)
            writer.submit(b.id, std::move(b.idx), std::move(b.sf));
        writer.finish();
    }

    EXPECT_EQ(slurp(reference), slurp(streamed));
    std::filesystem::remove(reference);
    std::filesystem::remove(streamed);
}
