#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <string>
#include <vector>


class Config {
  public:
    Config(int argc, const char *argv[]);

    bool verbose = false;
    bool reverse = false;
    bool filename_anno = false;
    bool fasta_anno = false;
    bool dump_raw_anno = false;
    bool sparse = false;
    bool count_labels = false;
    bool suppress_unlabeled = false;
    bool use_kmc = false;
    bool canonical_mode = false;
    bool complete = false;
    bool greedy_brwt = false;

    unsigned int k = 3;
    unsigned int filter_k = 3;
    unsigned int parallel = 1;
    unsigned int max_unreliable_abundance = 0;
    unsigned int unreliable_kmers_threshold = 0;
    unsigned int num_top_labels = -1;
    unsigned int arity_brwt = 2;
    unsigned int relax_arity_brwt = 10;

    double discovery_fraction = 1.0;

    std::vector<std::string> fname;
    std::vector<std::string> anno_labels;
    std::vector<std::string> infbase_annotators;
    std::string outfbase;
    std::string infbase;
    std::string rename_instructions_file;
    std::string refpath;
    std::string fasta_header_delimiter;
    std::string anno_labels_delimiter = ":";
    std::string annotation_label = "";

    enum IdentityType {
        NO_IDENTITY = -1,
        BUILD = 1,
        STATS,
        ANNOTATE,
        MERGE_ANNOTATIONS,
        TRANSFORM_ANNOTATION,
        RELAX_BRWT,
        CLASSIFY
    };
    IdentityType identity = NO_IDENTITY;

    enum AnnotationType {
        ColumnCompressed = 1,
        RowCompressed,
        BRWT,
        BinRelWT_sdsl,
        BinRelWT,
        RowFlat,
        RBFish,
    };

    AnnotationType anno_type = ColumnCompressed;

    static std::string annotype_to_string(AnnotationType state);
    static AnnotationType string_to_annotype(const std::string &string);

    void print_usage(const std::string &prog_name,
                     IdentityType identity = NO_IDENTITY);
};

#endif // __CONFIG_HPP__
