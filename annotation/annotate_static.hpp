#ifndef __ANNOTATE_STATIC_HPP__
#define __ANNOTATE_STATIC_HPP__

#include <memory>

#include "annotate.hpp"


namespace annotate {

template <class BinaryMatrixType, typename Label = std::string>
class StaticBinRelAnnotator : public MultiLabelEncoded<uint64_t, Label> {
  public:
    using Index = typename MultiLabelEncoded<uint64_t, Label>::Index;
    using VLabels = typename MultiLabelEncoded<uint64_t, Label>::VLabels;

    StaticBinRelAnnotator() : matrix_(new BinaryMatrixType()) {}
    StaticBinRelAnnotator(std::unique_ptr<BinaryMatrixType>&& matrix,
                          const LabelEncoder<Label> &label_encoder);

    bool has_label(Index i, const Label &label) const ;
    bool has_labels(Index i, const VLabels &labels) const ;

    VLabels get_labels(Index i) const;
    // Get labels that occur at least in |presence_ratio| rows.
    // If |presence_ratio| = 0, return all occurring labels.
    VLabels get_labels(const std::vector<Index> &indices,
                       double presence_ratio) const;

    void serialize(const std::string &filename) const;
    bool merge_load(const std::vector<std::string> &filenames);

    uint64_t num_objects() const;
    size_t num_labels() const;
    uint64_t num_relations() const;

    void set_labels(Index, const VLabels &) { except_dyn(); }
    void add_label(Index, const Label &) { except_dyn(); }
    void add_labels(Index, const VLabels &) { except_dyn(); }
    void add_labels(const std::vector<Index> &, const VLabels &) { except_dyn(); }
    void insert_rows(const std::vector<Index> &) { except_dyn(); }

    const BinaryMatrixType& data() const { return *matrix_; }

  private:
    std::vector<uint64_t> count_labels(const std::vector<Index> &indices) const;

    void except_dyn();

    std::unique_ptr<BinaryMatrixType> matrix_;

    LabelEncoder<Label> &label_encoder_ {
        MultiLabelEncoded<uint64_t, Label>::label_encoder_
    };

    static const std::string kFileExtension;
};

} // namespace annotate

#endif // __ANNOTATE_STATIC_HPP__
