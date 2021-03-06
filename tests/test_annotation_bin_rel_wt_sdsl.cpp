#include <random>

#include "gtest/gtest.h"

#include "annotate_column_compressed.hpp"
#include "static_annotators_def.hpp"
#include "annotation_converters.hpp"
#include "utils.hpp"

const std::string test_data_dir = "../tests/data";
const std::string test_dump_basename = test_data_dir + "/dump_test";
const std::string test_dump_basename_vec_bad = test_dump_basename + "_bad_filename";
const std::string test_dump_basename_vec_good = test_dump_basename + "_bin_rel_wt";


class BinRelWT_sdslTest : public ::testing::Test {
 protected:
    static annotate::BinRelWT_sdslAnnotator *annotation;
    virtual void SetUp() {
        annotate::ColumnCompressed<> column_annotator(5);
        column_annotator.add_labels(0, {"Label0", "Label2", "Label8"});
        column_annotator.add_labels(2, {"Label1", "Label2"});
        column_annotator.add_labels(3, {"Label1", "Label2", "Label8"});
        column_annotator.add_labels(4, {"Label2"});
        annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
            std::move(column_annotator)
        ).release();
    }

    virtual void TearDown() {
         delete annotation;
    }

};

annotate::BinRelWT_sdslAnnotator* BinRelWT_sdslTest::annotation = NULL;

TEST(BinRelWT_sdsl, EmptyConstructor) {
    annotate::ColumnCompressed<> empty_column_annotator(5);
    auto empty_annotation = annotate::convert<annotate::BinRelWT_sdslAnnotator>(
        std::move(empty_column_annotator)
    );
    EXPECT_EQ(0u, empty_annotation->num_labels());
    EXPECT_EQ(0u, empty_annotation->num_objects());
    EXPECT_EQ(0u, empty_annotation->num_relations());
}

//TODO: Add tests for copying from other annotators

std::set<std::string> convert_to_set(const std::vector<std::string> &vector);

TEST_F(BinRelWT_sdslTest, GetLabels) {
    EXPECT_EQ(convert_to_set({"Label0", "Label2", "Label8"}),
              convert_to_set(annotation->get_labels(0)));
    EXPECT_EQ(std::vector<std::string>({}),
              annotation->get_labels(1));
    EXPECT_EQ(convert_to_set({"Label1", "Label2"}),
              convert_to_set(annotation->get_labels(2)));
    EXPECT_EQ(convert_to_set({"Label1", "Label2", "Label8"}),
              convert_to_set(annotation->get_labels(3)));
    EXPECT_EQ(convert_to_set({"Label2"}),
              convert_to_set(annotation->get_labels(4)));
}

TEST_F(BinRelWT_sdslTest, NumObjects) {
    EXPECT_EQ(5u, annotation->num_objects());
}

TEST_F(BinRelWT_sdslTest, NumLables) {
    EXPECT_EQ(4u, annotation->num_labels());
}

TEST_F(BinRelWT_sdslTest, HasLabel) {
    EXPECT_TRUE(annotation->has_label(0, "Label0"));
    EXPECT_TRUE(annotation->has_label(0, "Label2"));
    EXPECT_TRUE(annotation->has_label(0, "Label8"));
    EXPECT_FALSE(annotation->has_label(0, "Label1"));

    EXPECT_FALSE(annotation->has_label(1, "Label0"));

    EXPECT_TRUE(annotation->has_label(2, "Label1"));
    EXPECT_TRUE(annotation->has_label(2, "Label2"));
    EXPECT_FALSE(annotation->has_label(2, "Label0"));

    EXPECT_TRUE(annotation->has_label(3, "Label1"));
    EXPECT_TRUE(annotation->has_label(3, "Label2"));
    EXPECT_TRUE(annotation->has_label(3, "Label8"));
    EXPECT_FALSE(annotation->has_label(3, "Label0"));

    EXPECT_TRUE(annotation->has_label(4, "Label2")) << annotation->get_labels(4)[0];
    EXPECT_FALSE(annotation->has_label(4, "Label0"));
    EXPECT_FALSE(annotation->has_label(4, "Label1"))<< annotation->get_labels(4)[0];
    EXPECT_FALSE(annotation->has_label(4, "Label8"));
}


TEST_F(BinRelWT_sdslTest, SerializationAndLoad) {
    {
       annotation->serialize(test_dump_basename + "_bin_rel_wt");
    }
    {
        ASSERT_FALSE(annotation->merge_load({test_dump_basename_vec_bad}));
        ASSERT_TRUE(annotation->merge_load({test_dump_basename_vec_good}));

        EXPECT_EQ(convert_to_set({ "Label0", "Label2", "Label8" }), convert_to_set(annotation->get_labels(0)));
        EXPECT_EQ(convert_to_set({}), convert_to_set(annotation->get_labels(1)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2" }), convert_to_set(annotation->get_labels(2)));
        EXPECT_EQ(convert_to_set({ "Label1", "Label2", "Label8" }), convert_to_set(annotation->get_labels(3)));
        EXPECT_EQ(convert_to_set({ "Label2" }), convert_to_set(annotation->get_labels(4)));
    }
}

TEST_F(BinRelWT_sdslTest, Sparsity){
    EXPECT_EQ(1 - static_cast<double>(annotation->num_relations())
                                          / annotation->num_objects()
                                          / annotation->num_labels(), 0.55);
}

