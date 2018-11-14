#include "config.hpp"

#include <cstring>
#include <iostream>

#include "utils.hpp"


Config::Config(int argc, const char *argv[]) {
    // provide help overview if no identity was given
    if (argc == 1) {
        print_usage(argv[0]);
        exit(-1);
    }

    // parse identity from first command line argument
    if (!strcmp(argv[1], "build")) {
        identity = BUILD;
    } else if (!strcmp(argv[1], "stats")) {
        identity = STATS;
    } else if (!strcmp(argv[1], "annotate")) {
        identity = ANNOTATE;
    } else if (!strcmp(argv[1], "merge_anno")) {
        identity = MERGE_ANNOTATIONS;
    } else if (!strcmp(argv[1], "classify")) {
        identity = CLASSIFY;
    } else if (!strcmp(argv[1], "transform_anno")) {
        identity = TRANSFORM_ANNOTATION;
    } else if (!strcmp(argv[1], "relax_brwt")) {
        identity = RELAX_BRWT;
    } else {
        print_usage(argv[0]);
        exit(-1);
    }

    // provide help screen for chosen identity
    if (argc == 2) {
        print_usage(argv[0], identity);
        exit(-1);
    }

    // parse remaining command line items
    for (int i = 2; i < argc; ++i) {
        if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "--verbose")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-r") || !strcmp(argv[i], "--reverse")) {
            reverse = true;
        } else if (!strcmp(argv[i], "--complete")) {
            complete = true;
        } else if (!strcmp(argv[i], "--anno-filename")) {
            filename_anno = true;
        } else if (!strcmp(argv[i], "--anno-header")) {
            fasta_anno = true;
        } else if (!strcmp(argv[i], "--anno-label")) {
            anno_labels.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "--suppress-unlabeled")) {
            suppress_unlabeled = true;
        } else if (!strcmp(argv[i], "--sparse")) {
            sparse = true;
        } else if (!strcmp(argv[i], "-p") || !strcmp(argv[i], "--parallel")) {
            parallel = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "--kmer-length")) {
            k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-abund")) {
            max_unreliable_abundance = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-thres")) {
            unreliable_kmers_threshold = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--filter-k")) {
            filter_k = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--dump-raw-anno")) {
            dump_raw_anno = true;
        } else if (!strcmp(argv[i], "--kmc")) {
            //TODO: add into some USAGE description
            use_kmc = true;
        } else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--canonical")) {
            canonical_mode = true;
        } else if (!strcmp(argv[i], "--discovery-fraction")) {
            discovery_fraction = std::stof(argv[++i]);
        } else if (!strcmp(argv[i], "--count-labels")) {
            count_labels = true;
        } else if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--outfile-base")) {
            outfbase = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--reference")) {
            refpath = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--header-delimiter")) {
            fasta_header_delimiter = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--labels-delimiter")) {
            anno_labels_delimiter = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "--num-top-labels")) {
            num_top_labels = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--anno-type")) {
            anno_type = string_to_annotype(argv[++i]);
        } else if (!strcmp(argv[i], "--rename-cols")) {
            rename_instructions_file = std::string(argv[++i]);
        } else if (!strcmp(argv[i], "-a") || !strcmp(argv[i], "--annotator")) {
            infbase_annotators.emplace_back(argv[++i]);
        } else if (!strcmp(argv[i], "-i") || !strcmp(argv[i], "--infile-base")) {
            infbase = std::string(argv[++i]);
            infbase = utils::remove_suffix(infbase, ".dbg");
        } else if (!strcmp(argv[i], "--greedy")) {
            greedy_brwt = true;
        } else if (!strcmp(argv[i], "--arity")) {
            arity_brwt = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "--relax-arity")) {
            relax_arity_brwt = atoi(argv[++i]);
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            print_usage(argv[0], identity);
            exit(0);
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "\nERROR: Unknown option %s\n\n", argv[i]);
            print_usage(argv[0], identity);
            exit(-1);
        } else {
            fname.push_back(argv[i]);
        }
    }

    if (!fname.size() && identity != STATS) {
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.size())
                fname.push_back(line);
        }
    }

    bool print_usage_and_exit = false;

    if (identity == CLASSIFY && infbase.empty())
        print_usage_and_exit = true;

    if (identity == ANNOTATE && infbase.empty())
        print_usage_and_exit = true;

    if ((identity == ANNOTATE)
            && !filename_anno && !fasta_anno && !anno_labels.size()) {
        std::cerr << "Error: No annotation to add" << std::endl;
        print_usage_and_exit = true;
    }

    if (identity == ANNOTATE && outfbase.empty())
        outfbase = infbase;

    if (identity == MERGE_ANNOTATIONS && outfbase.empty())
        print_usage_and_exit = true;

    if (identity == CLASSIFY && infbase_annotators.empty())
        infbase_annotators.push_back(infbase);

    if (discovery_fraction < 0 || discovery_fraction > 1)
        print_usage_and_exit = true;

    // if misused, provide help screen for chosen identity and exit
    if (print_usage_and_exit) {
        print_usage(argv[0], identity);
        exit(-1);
    }
}

std::string Config::annotype_to_string(AnnotationType state) {
    switch (state) {
        case ColumnCompressed:
            return "column";
        case RowCompressed:
            return "row";
        case BRWT:
            return "brwt";
        case BinRelWT_sdsl:
            return "bin_rel_wt_sdsl";
        case BinRelWT:
            return "bin_rel_wt";
        case RowFlat:
            return "flat";
        case RBFish:
            return "rbfish";
        default:
            assert(false);
            return "Never happens";
    }
}

Config::AnnotationType Config::string_to_annotype(const std::string &string) {
    if (string == "column") {
        return AnnotationType::ColumnCompressed;
    } else if (string == "row") {
        return AnnotationType::RowCompressed;
    } else if (string == "brwt") {
        return AnnotationType::BRWT;
    } else if (string == "bin_rel_wt_sdsl") {
        return AnnotationType::BinRelWT_sdsl;
    } else if (string == "bin_rel_wt") {
        return AnnotationType::BinRelWT;
    } else if (string == "flat") {
        return AnnotationType::RowFlat;
    } else if (string == "rbfish") {
        return AnnotationType::RBFish;
    } else {
        std::cerr << "Error: unknown annotation representation" << std::endl;
        exit(1);
    }
}

void Config::print_usage(const std::string &prog_name, IdentityType identity) {
    fprintf(stderr, "Genome Graph Annotation Schemes -- Version 0.1\n\n");

    const char annotation_list[] = "('column', 'row', 'bin_rel_wt_sdsl', 'bin_rel_wt', 'flat', 'rbfish', 'brwt')";

    switch (identity) {
        case NO_IDENTITY: {
            fprintf(stderr, "Usage: %s <command> [command specific options]\n\n", prog_name.c_str());

            fprintf(stderr, "Available commands:\n");

            fprintf(stderr, "\tbuild\t\tconstruct a graph object from input sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats or integrate sequence\n");
            fprintf(stderr, "\t\t\tfiles in fast[a|q] formats into a given graph\n\n");

            fprintf(stderr, "\tstats\t\tprint graph statistics for given graph(s)\n\n");

            fprintf(stderr, "\tannotate\tgiven a graph and a fast[a|q] file, annotate\n");
            fprintf(stderr, "\t\t\tthe respective kmers\n\n");

            fprintf(stderr, "\tmerge_anno\tmerge annotation columns\n\n");

            fprintf(stderr, "\ttransform_anno\tchange representation of the graph annotation\n\n");

            fprintf(stderr, "\trelax_brwt\toptimize the tree structure in brwt annotator\n\n");

            fprintf(stderr, "\tclassify\tannotate sequences from fast[a|q] files\n\n");

            return;
        }
        case BUILD: {
            fprintf(stderr, "Usage: %s build [options] FASTQ1 [[FASTQ2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for build:\n");
            fprintf(stderr, "\t   --kmc \t\tparse k-mers from precomputed KMC database\n");
            fprintf(stderr, "\t   --reference [STR] \tbasename of reference sequence []\n");
            fprintf(stderr, "\t-c --canonical \t\tindex only canonical k-mers (e.g. for read sets) [off]\n");
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR]\tbasename of output file []\n");
            fprintf(stderr, "\t-k --kmer-length [INT] \tlength of the k-mer to use [3]\n");
            fprintf(stderr, "\t   --filter-abund [INT] threshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] max allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \tlength of k-mers used for counting and filtering [3]\n");
        } break;
        case STATS: {
            fprintf(stderr, "Usage: %s stats [options] GRAPH1 [[GRAPH2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for stats:\n");
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph [off]\n");
            fprintf(stderr, "\t-a --annotator [STR] \tbasename of annotator to update []\n");
            fprintf(stderr, "\t   --anno-type [STR] \tinternal annotation representation [column]\n");
            fprintf(stderr, "\t                     \t  "); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
        } break;
        case ANNOTATE: {
            fprintf(stderr, "Usage: %s annotate -i <graph_basename> [options] <PATH1> [[PATH2] ...]\n"
                            "\tEach path is given as file in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --kmc \t\t\tparse k-mers from precomputed KMC database\n");
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph [off]\n");
            fprintf(stderr, "\t   --reference [STR] \t\tbasename of reference sequence []\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tbasename of annotator to update []\n");
            fprintf(stderr, "\t   --anno-type [STR] \t\tinternal annotation representation [column]\n");
            fprintf(stderr, "\t                     \t\t  "); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            fprintf(stderr, "\t   --sparse \t\t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file [<graph_basename>]\n");
            fprintf(stderr, "\t-r --reverse \t\t\talso annotate reverse complement reads [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t   --anno-filename \t\tinclude filenames as annotation labels [off]\n");
            fprintf(stderr, "\t   --anno-header \t\textract annotation labels from headers of sequences in files [off]\n");
            fprintf(stderr, "\t   --header-delimiter [STR]\tdelimiter for splitting annotation header into multiple labels [off]\n");
            fprintf(stderr, "\t   --anno-label [STR]\t\tadd label to annotation for all sequences from the files passed []\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case MERGE_ANNOTATIONS: {
            fprintf(stderr, "Usage: %s merge_anno [options] -o <annotator_basename> <ANNOT1> [[ANNOT2] ...]\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for annotate:\n");
            fprintf(stderr, "\t   --anno-type [STR] \tinternal annotation representation [column]\n");
            fprintf(stderr, "\t                     \t  "); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            fprintf(stderr, "\t   --sparse \t\tuse the row-major sparse matrix to annotate graph [off]\n");
            // fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case CLASSIFY: {
            fprintf(stderr, "Usage: %s classify -i <graph_basename> [options] <FILE1> [[FILE2] ...]\n"
                            "\tEach file is given in fasta or fastq format.\n\n", prog_name.c_str());

            fprintf(stderr, "Available options for classify:\n");
            fprintf(stderr, "\t   --complete \t\tconstruct a complete graph [off]\n");
            fprintf(stderr, "\t-r --reverse \t\t\tclassify reverse complement sequences [off]\n");
            fprintf(stderr, "\t   --filter-abund [INT] \tthreshold for the abundance of reliable k-mers [0]\n");
            fprintf(stderr, "\t   --filter-thres [INT] \tmax allowed number of unreliable kmers in reliable reads [0]\n");
            fprintf(stderr, "\t   --filter-k [INT] \t\tlength of k-mers used for counting and filtering [3]\n");
            fprintf(stderr, "\t-o --outfile-base [STR] \tbasename of output file []\n");
            fprintf(stderr, "\t-a --annotator [STR] \t\tbasename of annotator [<graph_basename>]\n");
            fprintf(stderr, "\t   --anno-type [STR] \t\tinternal annotation representation [column]\n");
            fprintf(stderr, "\t                     \t\t  "); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            fprintf(stderr, "\t   --sparse \t\t\tuse the row-major sparse matrix to annotate graph [off]\n");
            fprintf(stderr, "\t   --suppress-unlabeled \tdo not show results for sequences missing in graph [off]\n");
            fprintf(stderr, "\t   --count-labels \t\tcount labels for k-mers from querying sequences [off]\n");
            fprintf(stderr, "\t   --num-top-labels \t\tmaximum number of frequent labels to print [off]\n");
            fprintf(stderr, "\t   --discovery-fraction \tfraction of labeled k-mers required for annotation [1.0]\n");
            fprintf(stderr, "\t   --labels-delimiter [STR]\tdelimiter for annotation labels [\":\"]\n");
            fprintf(stderr, "\t-p --parallel [INT] \t\tuse multiple threads for computation [1]\n");
        } break;
        case TRANSFORM_ANNOTATION: {
            fprintf(stderr, "Usage: %s transform_anno [options] -o <annotator_basename> ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --rename-cols [STR]\tfile with rules for renaming annotation labels []\n");
            fprintf(stderr, "\t                      \texample: 'L_1 L_1_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          L_2 L_2_renamed\n");
            fprintf(stderr, "\t                      \t          ... ...........'\n");
            fprintf(stderr, "\t   --anno-type [STR] \ttransform annotations to specified format [column]\n");
            fprintf(stderr, "\t                     \t  "); fprintf(stderr, annotation_list); fprintf(stderr, "\n");
            fprintf(stderr, "\t   --arity  \t\tarity in the brwt tree [2]\n");
            fprintf(stderr, "\t   --greedy  \t\tuse greedy column partitioning in brwt construction [off]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
        case RELAX_BRWT: {
            fprintf(stderr, "Usage: %s relax_brwt [options] -o <annotator_basename> ANNOTATOR\n\n", prog_name.c_str());

            fprintf(stderr, "\t-o --outfile-base [STR] basename of output file []\n");
            fprintf(stderr, "\t   --relax-arity [INT] \trelax brwt tree to optimize arity limited to this number [10]\n");
            fprintf(stderr, "\t-p --parallel [INT] \tuse multiple threads for computation [1]\n");
        } break;
    }

    fprintf(stderr, "\n\tGeneral options:\n");
    fprintf(stderr, "\t-v --verbose \t\tswitch on verbose output [off]\n");
    fprintf(stderr, "\t-q --quiet \t\tproduce as little log output as posible [off]\n");
    fprintf(stderr, "\t-h --help \t\tprint usage info\n");
    fprintf(stderr, "\n");
}
