#include <iostream>

#include "unix_tools.hpp"
#include "config.hpp"
#include "sequence_io.hpp"
#include "annotated_dbg.hpp"
#include "dbg_hash.hpp"
#include "dbg_hash_ordered.hpp"
#include "dbg_sd.hpp"
#include "annotate_row_compressed.hpp"
#include "annotate_column_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "kmc_parser.hpp"

typedef annotate::MultiLabelAnnotation<uint64_t, std::string> Annotator;

const size_t kNumCachedColumns = 10;


template <class Graph = DBGHashOrdered>
Graph* load_critical_graph_from_file(const std::string &filename) {
    auto *graph = new Graph(0);
    if (!graph->load(filename)) {
        std::cerr << "ERROR: can't load graph from file " << filename << std::endl;
        delete graph;
        exit(1);
    }
    return graph;
}

std::string get_filter_filename(std::string filename,
                                size_t k,
                                size_t max_unreliable_abundance,
                                size_t unreliable_kmers_threshold,
                                bool critical = true) {
    filename = filename + ".filter_k" + std::to_string(k)
                        + "_s" + std::to_string(max_unreliable_abundance)
                        + (max_unreliable_abundance
                             ? std::string("_")
                                + std::to_string(unreliable_kmers_threshold)
                             : "");
    if (!critical || std::ifstream(filename).good()) {
        return filename;
    } else if (max_unreliable_abundance == 0) {
        return "";
    }

    std::cerr << "ERROR: read filter "
              << filename << " does not exist."
              << " Filter reads first." << std::endl;
    exit(1);
}


void annotate_data(const std::vector<std::string> &files,
                   AnnotatedDBG *anno_graph,
                   bool reverse,
                   bool use_kmc,
                   size_t filter_k,
                   size_t max_unreliable_abundance,
                   size_t unreliable_kmers_threshold,
                   bool filename_anno,
                   bool fasta_anno,
                   const std::string &fasta_header_delimiter,
                   const std::vector<std::string> &anno_labels,
                   bool verbose) {
    size_t total_seqs = 0;

    std::unique_ptr<Timer> timer;
    if (verbose) {
        timer.reset(new Timer());
    }

    // iterate over input files
    for (const auto &file : files) {
        if (verbose) {
            std::cout << std::endl << "Parsing " << file << std::endl;
        }
        // read files
        if (use_kmc) {
            std::vector<std::string> labels;

            if (filename_anno) {
                labels.push_back(file);
            }

            for (const auto &label : anno_labels) {
                labels.push_back(label);
            }

            kmc::read_kmers(file, [&](std::string&& sequence) {
                anno_graph->annotate_sequence(std::move(sequence), labels);

                total_seqs += 1;
                if (verbose && total_seqs % 10000 == 0) {
                    std::cout << "processed " << total_seqs << " sequences"
                              << ", trying to annotate as ";
                    for (const auto &label : labels) {
                        std::cout << "<" << label << ">";
                    }
                    std::cout << ", " << timer->elapsed() << "sec" << std::endl;
                }
            }, max_unreliable_abundance + 1);

            if (verbose) {
                std::cout << "Finished extracting sequences from file " << file
                          << " in " << timer->elapsed() << "sec" << std::endl;
            }

        } else  if (utils::get_filetype(file) == "FASTA"
                        || utils::get_filetype(file) == "FASTQ") {
            read_fasta_file_critical(file,
                [&](kseq_t *read_stream) {
                    std::vector<std::string> labels;

                    if (fasta_anno) {
                        labels = utils::split_string(read_stream->name.s,
                                                     fasta_header_delimiter);
                    }
                    if (filename_anno) {
                        labels.push_back(file);
                    }

                    for (const auto &label : anno_labels) {
                        labels.push_back(label);
                    }

                    anno_graph->annotate_sequence(read_stream->seq.s, labels);

                    total_seqs += 1;
                    if (verbose && total_seqs % 10000 == 0) {
                        std::cout << "processed " << total_seqs << " sequences"
                                  << ", last was " << read_stream->name.s
                                  << ", trying to annotate as ";
                        for (const auto &label : labels) {
                            std::cout << "<" << label << ">";
                        }
                        std::cout << ", " << timer->elapsed() << "sec" << std::endl;
                    }
                },
                reverse,
                get_filter_filename(
                    file, filter_k,
                    max_unreliable_abundance,
                    unreliable_kmers_threshold
                )
            );
        } else {
            std::cerr << "ERROR: Filetype unknown for file "
                      << file << std::endl;
            exit(1);
        }

        if (verbose) {
            std::cout << "Finished extracting sequences from file " << file
                      << " in " << timer->elapsed() << "sec" << std::endl;
        }
    }

    // join threads if any were initialized
    anno_graph->join();
}


void execute_query(std::string seq_name,
                   std::string sequence,
                   bool count_labels,
                   bool suppress_unlabeled,
                   size_t num_top_labels,
                   double discovery_fraction,
                   std::string anno_labels_delimiter,
                   const AnnotatedDBG &anno_graph) {
    std::ostringstream oss;

    if (count_labels) {
        auto top_labels = anno_graph.get_top_labels(sequence, num_top_labels);

        if (!top_labels.size() && suppress_unlabeled)
            return;

        oss << seq_name << "\t";

        if (top_labels.size()) {
            oss << "<" << top_labels[0].first << ">:" << top_labels[0].second;
        }
        for (size_t i = 1; i < top_labels.size(); ++i) {
            oss << "\t<" << top_labels[i].first << ">:" << top_labels[i].second;
        }
        oss << "\n";
    } else {
        auto labels_discovered
                = anno_graph.get_labels(sequence, discovery_fraction);

        if (!labels_discovered.size() && suppress_unlabeled)
            return;

        oss << seq_name << "\t"
            << utils::join_strings(labels_discovered,
                                   anno_labels_delimiter) << "\n";
    }

    std::cout << oss.str();
}


std::unique_ptr<Annotator> initialize_annotation(const Config &config) {
    std::unique_ptr<Annotator> annotation;
    switch (config.anno_type) {
        case Config::ColumnCompressed: {
            annotation.reset(
                new annotate::ColumnCompressed<>(
                    0, kNumCachedColumns, config.verbose
                )
            );
            break;
        }
        case Config::RowCompressed: {
            annotation.reset(new annotate::RowCompressed<>(0, config.sparse));
            break;
        }
        case Config::BRWT: {
            annotation.reset(new annotate::BRWTCompressed<>());
            break;
        }
        case Config::BinRelWT_sdsl: {
            annotation.reset(new annotate::BinRelWT_sdslAnnotator());
            break;
        }
        case Config::BinRelWT: {
            annotation.reset(new annotate::BinRelWTAnnotator());
            break;
        }
        case Config::RowFlat: {
            annotation.reset(new annotate::RowFlatAnnotator());
            break;
        }
        case Config::RBFish: {
            annotation.reset(new annotate::RainbowfishAnnotator());
            break;
        }
    }
    return annotation;
}


int main(int argc, const char *argv[]) {
    std::unique_ptr<Config> config { new Config(argc, argv) };

    if (config->verbose) {
        std::cout << "#############################\n"
                  << "### Welcome to AnnoGraph! ###\n"
                  << "#############################\n" << std::endl;
    }

    const auto &files = config->fname;

    switch (config->identity) {
        case Config::BUILD: {
            if (config->complete) {
                if (config->verbose)
                    std::cout << "Build complete De Bruijn graph with k-mer size k="
                              << config->k << std::endl;

                Timer timer;
                std::unique_ptr<DeBruijnGraph> graph {
                    new DBGSD(config->k, config->canonical_mode)
                };

                std::cout << "Graph with " << graph->num_nodes()
                          << " k-mers was built in "
                          << timer.elapsed() << "sec" << std::endl;

                if (config->outfbase.size()) {
                    timer.reset();
                    graph->serialize(config->outfbase);
                    std::cout << timer.elapsed() << "sec" << std::endl;
                }

                return 0;
            }

            std::unique_ptr<DeBruijnGraph> graph {
                new DBGHashOrdered(config->k, config->canonical_mode)
            };

            if (config->verbose)
                std::cout << "Build De Bruijn Graph with k-mer size k="
                          << graph->get_k() << std::endl;

            Timer timer;

            if (config->verbose) {
                std::cout << "Start reading data and extracting k-mers" << std::endl;
            }

            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }

                Timer data_reading_timer;

                if (config->use_kmc) {
                    bool warning_different_k = false;
                    kmc::read_kmers(file, [&](std::string&& sequence) {
                        if (!warning_different_k && sequence.size() != graph->get_k()) {
                            std::cerr << "Warning: k-mers parsed from KMC database "
                                      << file << " have length " << sequence.size()
                                      << " but graph is constructed for k=" << graph->get_k()
                                      << std::endl;
                            warning_different_k = true;
                        }
                        graph->add_sequence(sequence);
                    }, config->max_unreliable_abundance + 1);
                } else if (utils::get_filetype(file) == "FASTA"
                            || utils::get_filetype(file) == "FASTQ") {
                    read_fasta_file_critical(file,
                        [&](kseq_t *read_stream) {
                            graph->add_sequence(read_stream->seq.s);
                        },
                        false,
                        get_filter_filename(file, config->filter_k,
                                            config->max_unreliable_abundance,
                                            config->unreliable_kmers_threshold)
                    );
                } else {
                    std::cerr << "ERROR: Filetype unknown for file "
                              << file << std::endl;
                    exit(1);
                }
                if (config->verbose) {
                    std::cout << "Finished extracting sequences from file " << file
                              << " in " << timer.elapsed() << "sec"
                              << ", number of k-mers in graph: " << graph->num_nodes() << std::endl;
                }
                if (config->verbose) {
                    std::cout << "File processed in "
                              << data_reading_timer.elapsed()
                              << "sec, current mem usage: "
                              << get_curr_mem2() / (1<<20) << " MB"
                              << ", total time: " << timer.elapsed()
                              << "sec" << std::endl;
                }
            }

            std::cout << "Graph with " << graph->num_nodes()
                      << " k-mers was built in "
                      << timer.elapsed() << "sec" << std::endl;

            if (config->outfbase.size()) {
                timer.reset();
                graph->serialize(config->outfbase);
                std::cout << timer.elapsed() << "sec" << std::endl;
            }

            return 0;
        }
        case Config::ANNOTATE: {
            // load graph
            DeBruijnGraph *graph = nullptr;

            if (config->complete) {
                graph = load_critical_graph_from_file<DBGSD>(config->infbase);
            } else {
                graph = load_critical_graph_from_file(config->infbase);
            }

            AnnotatedDBG anno_dbg(
                graph,
                config->parallel
            );

            // initialize empty annotation
            if (config->anno_type == Config::RowCompressed) {
                anno_dbg.set_annotation(
                    new annotate::RowCompressed<>(anno_dbg.num_anno_rows(),
                                                  config->sparse)
                );
            } else if (config->anno_type == Config::ColumnCompressed) {
                anno_dbg.set_annotation(
                    new annotate::ColumnCompressed<>(anno_dbg.num_anno_rows(),
                                                     kNumCachedColumns,
                                                     config->verbose)
                );
            }

            if (config->infbase_annotators.size()
                    && !anno_dbg.get_annotation().merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            if (!anno_dbg.check_compatibility(true)) {
                std::cerr << "ERROR: graph cannot be annotated" << std::endl;
                exit(1);
            }

            annotate_data(files,
                          &anno_dbg,
                          config->reverse,
                          config->use_kmc,
                          config->filter_k,
                          config->max_unreliable_abundance,
                          config->unreliable_kmers_threshold,
                          config->filename_anno,
                          config->fasta_anno,
                          config->fasta_header_delimiter,
                          config->anno_labels,
                          config->verbose);

            anno_dbg.get_annotation().serialize(config->outfbase);

            return 0;
        }
        case Config::MERGE_ANNOTATIONS: {
            std::unique_ptr<Annotator> annotation;
            if (config->anno_type == Config::RowCompressed) {
                throw std::runtime_error("To be implemented");
                annotation.reset(new annotate::RowCompressed<>(0, config->sparse));
            } else {
                annotation.reset(
                    new annotate::ColumnCompressed<>(
                        0, kNumCachedColumns, config->verbose
                    )
                );
            }

            if (!annotation->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }

            annotation->serialize(config->outfbase);

            return 0;
        }
        case Config::CLASSIFY: {
            // load graph
            DeBruijnGraph *graph = nullptr;

            if (config->complete) {
                graph = load_critical_graph_from_file<DBGSD>(config->infbase);
            } else {
                graph = load_critical_graph_from_file(config->infbase);
            }
            
            AnnotatedDBG anno_dbg(
                graph,
                config->parallel
            );

            anno_dbg.set_annotation(initialize_annotation(*config).release());

            if (!anno_dbg.get_annotation().merge_load(config->infbase_annotators)) {
                std::cerr << "ERROR: can't load annotations for graph "
                          << config->infbase + ".dbg"
                          << ", file corrupted" << std::endl;
                exit(1);
            }

            utils::ThreadPool thread_pool(std::max(1u, config->parallel) - 1);

            std::unique_ptr<Timer> timer { config->verbose ? new Timer() : NULL };

            if (!anno_dbg.check_compatibility(true)) {
                std::cerr << "ERROR: graph cannot be annotated" << std::endl;
                exit(1);
            }

            // iterate over input files
            for (const auto &file : files) {
                if (config->verbose) {
                    std::cout << std::endl << "Parsing " << file << std::endl;
                }

                size_t seq_count = 0;

                read_fasta_file_critical(file,
                    [&](kseq_t *read_stream) {
                        thread_pool.enqueue(execute_query,
                            std::to_string(seq_count++) + "\t"
                                + std::string(read_stream->name.s),
                            std::string(read_stream->seq.s),
                            config->count_labels,
                            config->suppress_unlabeled,
                            config->num_top_labels,
                            config->discovery_fraction,
                            config->anno_labels_delimiter,
                            std::ref(anno_dbg)
                        );
                    },
                    config->reverse,
                    get_filter_filename(file, config->filter_k,
                                        config->max_unreliable_abundance,
                                        config->unreliable_kmers_threshold)
                );
                if (config->verbose) {
                    std::cout << "Finished extracting sequences from file " << file
                              << " in " << timer->elapsed() << "sec" << std::endl;
                }

                // wait while all threads finish processing the current file
                thread_pool.join();
            }

            return 0;
        }
        case Config::STATS: {
            for (const auto &file : files) {
                DeBruijnGraph *graph = nullptr;

                if (config->complete) {
                    graph = load_critical_graph_from_file<DBGSD>(file);
                } else {
                    graph = load_critical_graph_from_file(file);
                }

                AnnotatedDBG anno_dbg(
                    graph,
                    config->parallel
                );

                std::cout << "Statistics for graph " << file << std::endl;
                std::cout << "nodes: " << graph->num_nodes() << std::endl;
                std::cout << "k: " << graph->get_k() << std::endl;
                if (dynamic_cast<DBGHashOrdered*>(graph)) {
                    std::cout << "canonical: "
                              << (dynamic_cast<DBGHashOrdered&>(*graph).is_canonical_mode()
                                    ? "yes"
                                    : "no")
                              << std::endl;
                } else if (dynamic_cast<DBGSD*>(graph)) {
                    std::cout << "canonical: "
                              << (dynamic_cast<DBGSD&>(*graph).is_canonical_mode()
                                    ? "yes"
                                    : "no")
                              << std::endl;
                }

            }

            for (const auto &file : config->infbase_annotators) {
                auto annotation = initialize_annotation(*config);
                if (!annotation->load(file)) {
                    std::cerr << "ERROR: can't load annotation from file "
                              << file << std::endl;
                    exit(1);
                }

                std::cout << "Statistics for annotation " << file << std::endl;
                std::cout << "labels: " << annotation->num_labels() << std::endl;
                std::cout << "density: " << std::scientific
                                          << static_cast<double>(annotation->num_relations())
                                                / annotation->num_objects()
                                                / annotation->num_labels() << std::endl;
                std::cout << "representation: "
                          << config->annotype_to_string(config->anno_type) << std::endl;
            }

            return 0;
        }
        case Config::TRANSFORM_ANNOTATION: {
            Timer timer;

            auto annotator = std::make_unique<annotate::ColumnCompressed<>>(
                0, kNumCachedColumns, config->verbose
            );

            if (config->verbose)
                std::cout << "Loading annotator...\t" << std::flush;

            if (!annotator->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->rename_instructions_file.size()) {
                if (config->verbose)
                    std::cout << "Renaming...\t" << std::flush;

                std::map<std::string, std::string> dict;
                std::ifstream instream(config->rename_instructions_file);
                if (!instream.is_open()) {
                    std::cerr << "ERROR: Can't open file "
                              << config->rename_instructions_file << std::endl;
                    exit(1);
                }
                std::string old_name;
                std::string new_name;
                while (instream.good() && !(instream >> old_name).eof()) {
                    instream >> new_name;
                    if (instream.fail() || instream.eof()) {
                        std::cerr << "ERROR: wrong format of the rules for"
                                  << " renaming annotation columns passed in file "
                                  << config->rename_instructions_file << std::endl;
                        exit(1);
                    }
                    dict[old_name] = new_name;
                }
                annotator->rename_columns(dict);

                annotator->serialize(config->outfbase);
                if (config->verbose)
                    std::cout << timer.elapsed() << "sec" << std::endl;
            }

            switch (config->anno_type) {
                case Config::ColumnCompressed:
                    break;
                case Config::RowCompressed: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    annotate::RowCompressed<> row_annotator(0);
                    annotator->convert_to_row_annotator(&row_annotator,
                                                        config->parallel);
                    annotator.reset();

                    row_annotator.serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BRWT: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto brwt_annotator = config->greedy_brwt
                        ? annotate::convert_to_greedy_BRWT<annotate::BRWTCompressed<>>(std::move(*annotator))
                        : annotate::convert_to_simple_BRWT<annotate::BRWTCompressed<>>(
                            std::move(*annotator),
                            config->arity_brwt
                        );

                    annotator.reset();

                    brwt_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BinRelWT_sdsl: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto binrelwt_sdsl_annotator
                            = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    binrelwt_sdsl_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::BinRelWT: {
                    if (config->verbose)
                        std::cout << "Converting...\t" << std::flush;

                    auto binrelwt_annotator = annotate::convert<annotate::BinRelWTAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                       std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    binrelwt_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::RowFlat: {
                    if (config->verbose)
                        std::cout << "Converting to flat annotator...\t" << std::flush;

                    auto flat_annotator = annotate::convert<annotate::RowFlatAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    flat_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
                case Config::RBFish: {
                    if (config->verbose)
                        std::cout << "Converting to rainbowfish annotator...\t" << std::flush;

                    auto flat_annotator = annotate::convert<annotate::RainbowfishAnnotator>(
                        std::move(*annotator)
                    );
                    annotator.reset();

                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;

                    if (config->verbose)
                        std::cout << "Serializing to " << config->outfbase
                                  << "...\t" << std::flush;
                    flat_annotator->serialize(config->outfbase);
                    if (config->verbose)
                        std::cout << timer.elapsed() << "sec" << std::endl;
                    break;
                }
            }
            return 0;
        }
        case Config::RELAX_BRWT: {
            Timer timer;

            auto annotator = std::make_unique<annotate::BRWTCompressed<>>();

            if (config->verbose)
                std::cout << "Loading annotator...\t" << std::flush;

            if (!annotator->merge_load(files)) {
                std::cerr << "ERROR: can't load annotations" << std::endl;
                exit(1);
            }
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            if (config->verbose)
                std::cout << "Relaxing BRWT tree...\t" << std::flush;

            annotate::relax_BRWT<annotate::BRWTCompressed<>>(annotator.get(),
                                                             config->relax_arity_brwt);

            annotator->serialize(config->outfbase);
            if (config->verbose)
                std::cout << timer.elapsed() << "sec" << std::endl;

            return 0;
        }
        default:
            std::cerr << "This mode is not supported in AnnoGraph" << std::endl;
            exit(1);
    }

    return 0;
}
