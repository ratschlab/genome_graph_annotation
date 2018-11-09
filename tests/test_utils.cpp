#include "gtest/gtest.h"

// Disable death tests by default
#ifndef _DEATH_TEST
#ifdef EXPECT_DEATH
#undef EXPECT_DEATH
#define EXPECT_DEATH(a, b) (void)0
#endif
#endif


#include "annotate_column_compressed.hpp"
#include "utils.hpp"
#include "bit_vector.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";

const std::vector<sdsl::bit_vector> vectors {
    { 0, 1, 1, 0, 0, 1, 0 },
    { 0, 1, 1, 1, 1, 1, 0 }
};

const std::vector<sdsl::bit_vector> matrix {
    { 0, 0 },
    { 1, 1 },
    { 1, 1 },
    { 0, 1 },
    { 0, 1 },
    { 1, 1 },
    { 0, 0 }
};

const std::vector<std::vector<uint64_t>> indices {
    { 1, 0 }, { 1, 1 },
    { 2, 0 }, { 2, 1 },
    { 3, 1 },
    { 4, 1 },
    { 5, 0 }, { 5, 1 }
};

utils::RowsFromColumnsTransformer generate_rct_file() {
    annotate::ColumnCompressed<> annotation(6);

    annotation.set_labels(1, { "Label0", "Label1" });
    annotation.set_labels(2, { "Label0", "Label1" });
    annotation.set_labels(3, { "Label1" });
    annotation.set_labels(4, { "Label1" });
    annotation.set_labels(5, { "Label0", "Label1" });

    annotation.dump_columns(test_dump_basename);

    utils::RowsFromColumnsTransformer rct(7, {
        test_dump_basename + ".0.raw.column.annodbg",
        test_dump_basename + ".1.raw.column.annodbg"
    });

    return rct;
}

std::vector<bit_vector_small> generate_rows() {
    std::vector<bit_vector_small> rows;
    std::transform(vectors.begin(), vectors.end(), std::back_inserter(rows),
        [](const auto &a) {
            return bit_vector_small(a);
        });

    return rows;
}

std::vector<bit_vector_small const*>
generate_ptrs(const std::vector<bit_vector_small> &rows) {
    std::vector<bit_vector_small const*> rows_ptr;
    std::transform(rows.begin(), rows.end(), std::back_inserter(rows_ptr),
        [](auto &a) {
            return &a;
        });

    return rows_ptr;
}

void check_rows(utils::RowsFromColumnsTransformer&& rct) {
    ASSERT_EQ(2u, rct.columns());
    ASSERT_EQ(7u, rct.rows());
    ASSERT_EQ(8u, rct.values_left());
    // ASSERT_EQ(std::vector<uint64_t>({ 3, 5 }), rct.num_set_bits());

    uint64_t i = 0;
    utils::call_rows([&](auto&& row_indices) {
        sdsl::bit_vector bv(rct.columns());
        for (auto j : row_indices) {
            bv[j] = 1;
        }
        EXPECT_EQ(matrix[i++], bv) << i;
    }, std::move(rct));

    ASSERT_EQ(7u, i);
    EXPECT_EQ(0u, rct.values_left());
}

void check_indices(utils::RowsFromColumnsTransformer&& rct) {
    ASSERT_EQ(2u, rct.columns());
    ASSERT_EQ(7u, rct.rows());
    ASSERT_EQ(8u, rct.values_left());
    // ASSERT_EQ(std::vector<uint64_t>({ 3, 5 }), rct.num_set_bits());

    uint64_t counter = 0;
    for (uint64_t i = 0; i < 8; ++i) {
        rct.call_next([&](auto row, auto column) {
            counter++;
            EXPECT_EQ(std::vector<uint64_t>({ row, column }),
                      indices[i]) << i;
        });
    }

    EXPECT_EQ(8u, counter);
    EXPECT_EQ(0u, rct.values_left());
}

TEST(Utils, RowsFromColumnsTransformerCallRowsFile) {
    auto rct = generate_rct_file();
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesFile) {
    auto rct = generate_rct_file();
    check_indices(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallRowsColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(generate_ptrs(rows));
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesColumns) {
    auto rows = generate_rows();
    utils::RowsFromColumnsTransformer rct(generate_ptrs(rows));
    check_indices(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallRowsConcat) {
    bit_vector_small vectors_small{ 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0 };
    auto rct = utils::RowsFromColumnsTransformer(vectors_small, 7);
    check_rows(std::move(rct));
}

TEST(Utils, RowsFromColumnsTransformerCallIndicesConcat) {
    bit_vector_small vectors_small{ 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0 };
    auto rct = utils::RowsFromColumnsTransformer(vectors_small, 7);
    check_indices(std::move(rct));
}

TEST(Utils, TempFileInitialize) {
    utils::TempFile tmp;
}

TEST(Utils, TempFileOpenWrite) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
    }
}

TEST(Utils, TempFileOpenWriteRead) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
    }
}

TEST(Utils, TempFileCheckStateFlow) {
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEATH(tmp.ofstream().good(), "Can't write after reading");
    }
    {
        utils::TempFile tmp;
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ofstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        ASSERT_TRUE(tmp.ifstream().good());
        EXPECT_DEATH(tmp.ofstream().good(), "Can't write after reading");
    }
}

TEST(Utils, TempFileReadWritePairs) {
    {
        utils::TempFile tmp[2];

        ASSERT_TRUE(tmp[0].ofstream().good());
        ASSERT_TRUE(tmp[0].ifstream().good());

        ASSERT_TRUE(tmp[1].ofstream().good());
        ASSERT_TRUE(tmp[1].ifstream().good());

        EXPECT_TRUE(tmp[0].ifstream().peek() == std::ifstream::traits_type::eof());
        EXPECT_TRUE(tmp[1].ifstream().peek() == std::ifstream::traits_type::eof());
    }

    utils::TempFile tmp[2];

    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;
    tmp[0].ofstream() << "test0 string0" << std::endl;
    tmp[1].ofstream() << "test1 string1" << std::endl;

    std::ifstream& in0 = tmp[0].ifstream();
    std::ifstream& in1 = tmp[1].ifstream();

    std::string temp;
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    in0 >> temp; EXPECT_EQ("test0", temp);
    in1 >> temp; EXPECT_EQ("test1", temp);
    in0 >> temp; EXPECT_EQ("string0", temp);
    in1 >> temp; EXPECT_EQ("string1", temp);
    EXPECT_FALSE(in0 >> temp);
    EXPECT_FALSE(in1 >> temp);
    EXPECT_TRUE(in0.peek() == std::ifstream::traits_type::eof());
    EXPECT_TRUE(in1.peek() == std::ifstream::traits_type::eof());
}

template <typename T>
std::vector<T> uniqueize(const std::vector<T> &input) {
    auto vector = input;
    vector.erase(std::unique(vector.begin(), vector.end()), vector.end());
    return vector;
}

using utils::sample_indexes;

TEST(Utils, sample_indexes) {
    std::mt19937 gen;
    gen.seed(14);

    ASSERT_EQ(0u, sample_indexes(10, 0, gen).size());
    ASSERT_EQ(5u, sample_indexes(10, 5, gen).size());
    ASSERT_EQ(6u, sample_indexes(10, 6, gen).size());
    ASSERT_EQ(9u, sample_indexes(10, 9, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 10, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 11, gen).size());
    ASSERT_EQ(10u, sample_indexes(10, 1000, gen).size());

    ASSERT_EQ(0u, sample_indexes(0, 0, gen).size());
    ASSERT_EQ(0u, sample_indexes(0, 1000, gen).size());

    auto generated = sample_indexes(10, 0, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 5, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 6, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 9, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 10, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 11, gen);
    ASSERT_EQ(generated, uniqueize(generated));
    generated = sample_indexes(10, 1000, gen);
    ASSERT_EQ(generated, uniqueize(generated));

    generated = sample_indexes(100'000, 1'000, gen);
    EXPECT_TRUE(*std::min_element(generated.begin(), generated.end()) < 100'000u / 5);
    EXPECT_TRUE(*std::max_element(generated.begin(), generated.end()) > 4 * 100'000u / 5)
        << *std::max_element(generated.begin(), generated.end());

    ASSERT_EQ(2'000'000u, sample_indexes(100'000'000, 2'000'000, gen).size());
    ASSERT_EQ(12'000'000u, sample_indexes(100'000'000, 12'000'000, gen).size());
}

TEST(KmerExtractor, encode_decode) {
    utils::KmerExtractor encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    EXPECT_EQ('N', encoder.decode(encoder.encode('N')));
    EXPECT_EQ('N', encoder.decode(encoder.encode('X')));
}

TEST(KmerExtractor, encode_decode_kmer) {
    utils::KmerExtractor encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
}

TEST(KmerExtractor, encode_decode_string) {
    utils::KmerExtractor encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto kmers = encoder.sequence_to_kmers(sequence, k);
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0]);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i])[k - 1]);
        }

        EXPECT_EQ(sequence, reconstructed);
    }
}



TEST(KmerExtractor2Bit, encode_decode) {
    utils::KmerExtractor2Bit encoder;
    EXPECT_EQ('A', encoder.decode(encoder.encode('A')));
    EXPECT_EQ('C', encoder.decode(encoder.encode('C')));
    EXPECT_EQ('G', encoder.decode(encoder.encode('G')));
    EXPECT_EQ('T', encoder.decode(encoder.encode('T')));
    EXPECT_EQ('A', encoder.decode(encoder.encode('N')));
    EXPECT_EQ('A', encoder.decode(encoder.encode('X')));
}

TEST(KmerExtractor2Bit, encode_decode_kmer) {
    utils::KmerExtractor2Bit encoder;
    std::string kmer;

    kmer = "ACGT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "AAAAAAAAA";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "TTTTTTTTT";
    EXPECT_EQ(kmer, encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANANA";
    EXPECT_EQ(std::string("AAAAAAAAAAA"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANATANANA";
    EXPECT_EQ(std::string("AAAAATAAAAA"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
    kmer = "ANANANANGNT";
    EXPECT_EQ(std::string("AAAAAAAAGAT"), encoder.kmer_to_sequence(
                        encoder.sequence_to_kmers(kmer, kmer.size())[0])) << kmer;
}

TEST(KmerExtractor2Bit, encode_decode_string) {
    utils::KmerExtractor2Bit encoder;
    std::string sequence = "AAGGCAGCCTACCCCTCTGTAN";
    for (uint64_t k = 2; k <= sequence.length(); ++k) {
        auto kmers = encoder.sequence_to_kmers(sequence, k);
        ASSERT_LT(0u, kmers.size());

        std::string reconstructed = encoder.kmer_to_sequence(kmers[0]);
        for (uint64_t i = 1; i < kmers.size(); ++i) {
            reconstructed.push_back(encoder.kmer_to_sequence(kmers[i])[k - 1]);
        }

        EXPECT_EQ(std::string("AAGGCAGCCTACCCCTCTGTAA"), reconstructed);
    }
}
