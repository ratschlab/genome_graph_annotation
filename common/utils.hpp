#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <queue>
#include <atomic>
#include <bitset>

#include <sdsl/sd_vector.hpp>

#if _USE_FOLLY
#include <folly/FBVector.h>
#include <folly/small_vector.h>
template <typename T>
using Vector = folly::fbvector<T>;
typedef folly::small_vector<uint32_t, 2, uint32_t> SmallVector;
#else
template <typename T>
using Vector = std::vector<T>;
typedef std::vector<uint32_t> SmallVector;
#endif


#include "serialization.hpp"
#include "kmer.hpp"
#include "kmer_packed.hpp"


struct SmallVectorHash {
    std::size_t operator()(const SmallVector &vector) const;
};

namespace utils {

    class KmerExtractor {
        static constexpr size_t kLogSigma = 3;

      public:
        // alphabet for k-mer representation
        typedef uint8_t TAlphabet;
        // kmer type
        typedef KMer<uint64_t, kLogSigma> Kmer;

        KmerExtractor();

        // extract k-mers from sequence
        Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                       size_t k) const;

        std::string kmer_to_sequence(const Kmer &kmer) const;

        // map input character to k-mer character
        TAlphabet encode(char s) const;
        // map k-mer character to input character
        char decode(TAlphabet c) const;

      private:
        std::vector<TAlphabet> encode(const std::string &sequence) const;
        std::string decode(const std::vector<TAlphabet> &sequence) const;

        Vector<Kmer>
        sequence_to_kmers(std::vector<TAlphabet>&& seq, size_t k) const;

        static const std::string alphabet;
        static const TAlphabet kCharToNucleotide[128];
    };

    class KmerExtractor2Bit {
        static constexpr size_t kLogSigma = 2;

      public:
        // alphabet for k-mer representation
        typedef uint8_t TAlphabet;
        // kmer type
        typedef KMerPacked<uint64_t, kLogSigma> Kmer;

        KmerExtractor2Bit();

        // extract k-mers from sequence
        Vector<Kmer> sequence_to_kmers(const std::string &sequence,
                                       size_t k) const;

        std::string kmer_to_sequence(const Kmer &kmer) const;

        // map input character to k-mer character
        TAlphabet encode(char s) const;
        // map k-mer character to input character
        char decode(TAlphabet c) const;

      private:
        std::vector<TAlphabet> encode(const std::string &sequence) const;
        std::string decode(const std::vector<TAlphabet> &sequence) const;

        Vector<Kmer>
        sequence_to_kmers(std::vector<TAlphabet>&& seq, size_t k) const;

        static const std::string alphabet;
        static const TAlphabet kCharToNucleotide[128];
    };

    std::string remove_suffix(const std::string &str, const std::string &suffix);

    template <typename... String>
    std::string remove_suffix(const std::string &str, const std::string &suffix,
                                                      const String&... other_suffices) {
        return remove_suffix(remove_suffix(str, suffix), other_suffices...);
    }

    std::string join_strings(const std::vector<std::string> &strings,
                             const std::string &delimiter);

    std::vector<std::string> split_string(const std::string &string,
                                          const std::string &delimiter);


    uint64_t kFromFile(const std::string &infbase);
    /**
     *  This function checks whether two given strings are identical.
     */
    template <class String>
    bool seq_equal(const String &s1, const String &s2, size_t start = 0) {
        if (s1.size() != s2.size())
            return false;

        for (size_t i = start; i < s1.size(); ++i) {
            if (s1.at(i) != s2.at(i))
                return false;
        }
        return true;
    }

    /**
     *  This function checks whether string s1 is co-lexicographically
     *  greater than s2.
     */
    template <class String>
    bool colexicographically_greater(const String &s1, const String &s2) {
        size_t ss1 = s1.size();
        size_t ss2 = s2.size();
        for (size_t i = 1; i <= std::min(ss1, ss2); ++i) {
            if (s1.at(ss1 - i) != s2.at(ss2 - i))
                return (s1.at(ss1 - i) > s2.at(ss2 - i));
        }
        return ss1 > ss2;
    }

    std::string get_filetype(const std::string &fname);

    std::deque<std::string> generate_strings(const std::string &alphabet,
                                             size_t length);

    uint32_t code_length(uint64_t a);


    /**
     * The code was copied and has been modified from:
     * https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
     */
    class ThreadPool {
      public:
        ThreadPool(size_t num_workers, size_t max_num_tasks = -1);

        template <class F, typename... Args>
        auto enqueue(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
            using return_type = typename std::result_of<F(Args...)>::type;
            auto task = std::make_shared<std::packaged_task<return_type()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

            if (!workers.size()) {
                (*task)();
                return task->get_future();
            } else {
                std::unique_lock<std::mutex> lock(queue_mutex);
                full_condition.wait(lock, [this]() {
                    return this->tasks.size() < this->max_num_tasks_;
                });
                tasks.emplace([task](){ (*task)(); });
            }
            empty_condition.notify_one();

            return task->get_future();
        }

        void join();

        ~ThreadPool();

      private:
        void initialize(size_t num_threads);

        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;
        size_t max_num_tasks_;

        std::mutex queue_mutex;
        std::condition_variable empty_condition;
        std::condition_variable full_condition;

        bool joining_;
        bool stop_;
    };


    class AsyncActivity {
      public:
        template <class F, typename... Args>
        auto run_async(F&& f, Args&&... args) -> decltype(f(args...)) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                parallel_jobs_++;
            }
            auto result = f(args...);
            {
                std::unique_lock<std::mutex> lock(mutex_);
                parallel_jobs_--;
            }
            cond_var_.notify_one();
            return std::move(result);
        }

        template <class F, typename... Args>
        auto run_unique(F&& f, Args&&... args) -> decltype(f(args...)) {
            std::unique_lock<std::mutex> lock(mutex_);
            while (parallel_jobs_ > 0) {
                cond_var_.wait(lock);
            }
            return f(args...);
        }

      private:
        size_t parallel_jobs_ = 0;

        std::mutex mutex_;
        std::condition_variable cond_var_;
    };


    /** A faster alternative to std::allocator<T>
     *
     * The code was copied and has been modified from:
     * https://probablydance.com/2014/11/09/plalloc-a-simple-stateful-allocator-for-node-based-containers/
     */
    template <typename T>
    class plalloc {
      public:
        typedef T value_type;

        plalloc() = default;
        template <typename U>
        plalloc(const plalloc<U>&) {}
        plalloc(const plalloc&) {}
        plalloc& operator=(const plalloc&) { return *this; }
        plalloc(plalloc&&) = default;
        plalloc& operator=(plalloc&&) = default;

        typedef std::true_type propagate_on_container_copy_assignment;
        typedef std::true_type propagate_on_container_move_assignment;
        typedef std::true_type propagate_on_container_swap;

        bool operator==(const plalloc &other) const { return this == &other; }
        bool operator!=(const plalloc &other) const { return !(*this == other); }

        T* allocate(size_t num_to_allocate) {
            if (num_to_allocate != 1)
                return static_cast<T*>(::operator new(sizeof(T) * num_to_allocate));

            if (available.size()) {
                T *result = available.back();
                available.pop_back();
                return result;
            }

            // first allocate 8, then double whenever
            // we run out of memory
            size_t to_allocate = 8 << memory.size();
            available.reserve(to_allocate);
            std::unique_ptr<value_holder[]> allocated(new value_holder[to_allocate]);
            value_holder *first_new = allocated.get();
            memory.emplace_back(std::move(allocated));
            size_t to_return = to_allocate - 1;
            for (size_t i = 0; i < to_return; ++i) {
                available.push_back(std::addressof(first_new[i].value));
            }
            return std::addressof(first_new[to_return].value);
        }
        void deallocate(T *ptr, size_t num_to_free) {
            if (num_to_free == 1) {
                available.push_back(ptr);
            } else {
                ::operator delete(ptr);
            }
        }

        // boilerplate that shouldn't be needed, except
        // libstdc++ doesn't use allocator_traits yet
        template<typename U>
        struct rebind {
            typedef plalloc<U> other;
        };

        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

        template<typename U, typename... Args>
        void construct(U *object, Args&&... args) {
            new (object) U(std::forward<Args>(args)...);
        }
        template<typename U, typename... Args>
        void construct(const U *object, Args &&... args) = delete;

        template<typename U>
        void destroy(U *object) { object->~U(); }

      private:
        union value_holder {
            value_holder() {}
            ~value_holder() {}
            T value;
        };

        std::vector<std::unique_ptr<value_holder[]>> memory;
        std::vector<T*> available;
    };

    template <typename T>
    struct Hash {
        size_t operator()(const T &x) const {
            return hasher(reinterpret_cast<const std::bitset<sizeof(T) * 8>&>(x));
        }

        std::hash<std::bitset<sizeof(T) * 8>> hasher;
    };

    void decompress_sd_vector(const sdsl::sd_vector<> &vector,
                              sdsl::bit_vector *out);

    // indexes - positions of inserted elements in the final vector
    template <typename Index, class Vector>
    void insert_default_values(const std::vector<Index> &indexes, Vector *vector);

    template <typename T>
    void erase(std::vector<T> *vector, const std::vector<bool> &erase_mask);


    // Read indices of set bits from a vector of VectorStreams
    class RowsFromColumnsTransformer {
      public:
        // Files store serialized vectors as plain indexes of the set bits
        RowsFromColumnsTransformer(uint64_t num_rows,
                                   const std::vector<std::string> &files);

        template <typename BitVectorPtr>
        RowsFromColumnsTransformer(const std::vector<BitVectorPtr> &columns);

        RowsFromColumnsTransformer(const bit_vector_small &columns_concatenated,
                                   uint64_t column_size);

        using ValueCallback = std::function<void(uint64_t /*row*/,
                                                 uint64_t /*column*/)>;
        void call_next(ValueCallback callback);

        uint64_t rows() const { return num_rows_; }
        uint64_t columns() const { return streams_.size(); }

        // get the number of set bits in the vectors left
        uint64_t values_left() const { return num_set_bits_left_; }

      private:
        void init_heap();

        std::vector<std::unique_ptr<VectorStream>> streams_;
        uint64_t num_set_bits_left_ = 0;
        uint64_t num_rows_;

        // store pair of kmer index and label index
        struct kmer_label_pair {
            uint64_t row_id;
            uint64_t col_id;

            bool operator<(const kmer_label_pair &other) const {
                return row_id > other.row_id || (row_id == other.row_id
                                                    && col_id > other.col_id);
            }
        };

        std::priority_queue<kmer_label_pair,
                            std::vector<kmer_label_pair>> index_heap_;
    };

    using SetBitPositions = std::vector<uint64_t>;
    void call_rows(const std::function<void(const SetBitPositions &)> &callback,
                   RowsFromColumnsTransformer&& transformer);

    template <typename... Args>
    void call_rows(const std::function<void(const SetBitPositions &)> &callback,
                   Args&&... args) {
        call_rows(callback,
                  RowsFromColumnsTransformer(std::forward<Args>(args)...));
    }

    template <class BitVectorType = bit_vector_stat>
    std::vector<std::unique_ptr<bit_vector>>
    transpose(const std::vector<std::unique_ptr<bit_vector>> &matrix);


    class TempFile {
      public:
        // The state flow:
        //    init -> APPEND -> READ -> deinit
        enum State { APPEND, READ };

        TempFile(const std::string &tmp_dir = "");
        ~TempFile();

        std::ofstream& ofstream();
        std::ifstream& ifstream();

      private:
        std::string tmp_file_name_;
        State state_;
        std::unique_ptr<std::ofstream> tmp_ostream_;
        std::unique_ptr<std::ifstream> tmp_istream_;
    };

    // Partitions a range of numbers in [0,n) into groups.
    // The groups must not be empty.
    class RangePartition {
      public:
        typedef uint32_t T;
        typedef uint32_t G;
        typedef uint32_t R;

        RangePartition() {}
        RangePartition(const std::vector<uint64_t> &arrangement,
                       const std::vector<size_t> &group_sizes);
        explicit RangePartition(std::vector<std::vector<uint64_t>>&& partition);

        explicit RangePartition(const RangePartition &) = default;
        RangePartition& operator=(const RangePartition &) = default;
        RangePartition(RangePartition&&) = default;
        RangePartition& operator=(RangePartition&&) = default;

        // get group that contains value
        G group(T value) const;

        // get index of value in its group
        R rank(T value) const;

        // get value given its group and the rank
        T get(G group, R rank) const;

        uint64_t num_groups() const;
        uint64_t size() const;

        bool load(std::istream &in);
        void serialize(std::ostream &out) const;

      private:
        // Based on |partition_|, initializes groups and ranks.
        // Returns false if partition is invalid.
        bool initialize_groups_and_ranks();

        std::vector<std::vector<T>> partition_;
        std::vector<G> groups_;
        std::vector<R> ranks_;
    };

    template <typename T>
    std::vector<T> arange(T first, size_t size) {
        std::vector<T> result(size);
        std::iota(result.begin(), result.end(), first);
        return result;
    }


    // indexes are distinct and sorted
    sdsl::bit_vector subvector(const bit_vector &col,
                               const std::vector<uint64_t> &indexes);

    std::vector<uint64_t> sample_indexes(uint64_t universe_size,
                                         uint64_t sample_size,
                                         std::mt19937 &gen);

} // namespace utils

template <typename T>
std::set<T> convert_to_set(const std::vector<T> &vector);

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

std::set<std::pair<std::string, size_t>>
to_set(const std::vector<std::pair<std::string, size_t>> &vector);

#endif // __UTILS_HPP__
