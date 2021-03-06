#include "bit_vector.hpp"

#include <cassert>
#include <libmaus2/bitio/putBit.hpp>

#include "serialization.hpp"


uint64_t inner_prod(const sdsl::bit_vector &first,
                    const sdsl::bit_vector &second) {
    assert(first.size() == second.size());

    if (first.empty())
        return 0;

    const uint64_t *first_data = first.data();
    const uint64_t *second_data = second.data();

    uint64_t result = sdsl::bits::cnt((*first_data) & (*second_data));

    for (typename sdsl::bit_vector::size_type i = 1; i < (first.capacity() >> 6); ++i) {
        result += sdsl::bits::cnt(*(++first_data) & *(++second_data));
    }
    if (first.bit_size() & 0x3F) {
        result -= sdsl::bits::cnt((*first_data) & (*second_data)
                                    & (~sdsl::bits::lo_set[first.bit_size() & 0x3F]));
    }
    return result;
}


uint64_t bit_vector::rank0(uint64_t id) const {
    return std::min(id + 1, size()) - rank1(id);
}

std::vector<bool> bit_vector::to_vector() const {
    std::vector<bool> result(size(), false);
    call_ones([&result](auto i) { result[i] = true; });
    return result;
}

void bit_vector::add_to(std::vector<bool> *other) const {
    assert(other);
    assert(other->size() == size());
    call_ones([other](auto i) { (*other)[i] = true; });
}

std::ostream& operator<<(std::ostream &os, const bit_vector &bv) {
    for (bool a : bv.to_vector()) {
        os << a;
    }
    return os;
}

template <class Vector>
Vector bit_vector::convert_to() {
    if (dynamic_cast<Vector*>(this)) {
        return Vector(dynamic_cast<Vector&&>(*this));
    } else if (dynamic_cast<bit_vector_stat*>(this)) {
        return Vector(std::move(dynamic_cast<bit_vector_stat*>(this)->vector_));
    } else {
        sdsl::bit_vector bv(size(), 0);
        call_ones([&bv](auto i) { bv[i] = true; });
        return Vector(std::move(bv));
    }
}
template bit_vector_dyn bit_vector::convert_to<bit_vector_dyn>();
template bit_vector_stat bit_vector::convert_to<bit_vector_stat>();
template bit_vector_sd bit_vector::convert_to<bit_vector_sd>();
template bit_vector_rrr<3> bit_vector::convert_to<bit_vector_rrr<3>>();
template bit_vector_rrr<8> bit_vector::convert_to<bit_vector_rrr<8>>();
template bit_vector_rrr<15> bit_vector::convert_to<bit_vector_rrr<15>>();
template bit_vector_rrr<31> bit_vector::convert_to<bit_vector_rrr<31>>();
template bit_vector_rrr<63> bit_vector::convert_to<bit_vector_rrr<63>>();
template bit_vector_rrr<127> bit_vector::convert_to<bit_vector_rrr<127>>();
template bit_vector_rrr<255> bit_vector::convert_to<bit_vector_rrr<255>>();
template bit_vector_small bit_vector::convert_to<bit_vector_small>();
template sdsl::bit_vector bit_vector::convert_to<sdsl::bit_vector>();

template <class Vector>
Vector bit_vector::copy_to() const {
    if (dynamic_cast<const Vector*>(this)) {
        return Vector(dynamic_cast<const Vector&>(*this));
    } else if (dynamic_cast<const bit_vector_stat*>(this)) {
        auto bv = dynamic_cast<const bit_vector_stat*>(this)->vector_;
        return Vector(std::move(bv));
    } else {
        sdsl::bit_vector bv(size(), 0);
        call_ones([&bv](auto i) { bv[i] = true; });
        return Vector(std::move(bv));
    }
}
template bit_vector_dyn bit_vector::copy_to<bit_vector_dyn>() const;
template bit_vector_stat bit_vector::copy_to<bit_vector_stat>() const;
template bit_vector_sd bit_vector::copy_to<bit_vector_sd>() const;
template bit_vector_rrr<3> bit_vector::copy_to<bit_vector_rrr<3>>() const;
template bit_vector_rrr<8> bit_vector::copy_to<bit_vector_rrr<8>>() const;
template bit_vector_rrr<15> bit_vector::copy_to<bit_vector_rrr<15>>() const;
template bit_vector_rrr<31> bit_vector::copy_to<bit_vector_rrr<31>>() const;
template bit_vector_rrr<63> bit_vector::copy_to<bit_vector_rrr<63>>() const;
template bit_vector_rrr<127> bit_vector::copy_to<bit_vector_rrr<127>>() const;
template bit_vector_rrr<255> bit_vector::copy_to<bit_vector_rrr<255>>() const;
template bit_vector_small bit_vector::copy_to<bit_vector_small>() const;
template sdsl::bit_vector bit_vector::copy_to<sdsl::bit_vector>() const;

/////////////////////////////
// bit_vector_dyn, libmaus //
/////////////////////////////

template <class BitVector>
std::vector<uint64_t> pack_bits(const BitVector &v) {
    std::vector<uint64_t> bits((v.size() + 63) / 64);
    for (size_t i = 0; i < v.size(); ++i) {
        libmaus2::bitio::putBit(bits.data(), i, v[i]);
    }
    return bits;
}

bit_vector_dyn::bit_vector_dyn(const std::vector<uint64_t> &bits_packed,
                               size_t num_bits)
      : vector_(num_bits, bits_packed.data()) {}

bit_vector_dyn::bit_vector_dyn(uint64_t size, bool value)
      : vector_(size, value) {}

template <class BitVector>
bit_vector_dyn::bit_vector_dyn(const BitVector &v)
      : bit_vector_dyn(pack_bits(v), v.size()) {}

template bit_vector_dyn::bit_vector_dyn(const std::vector<bool> &);
template bit_vector_dyn::bit_vector_dyn(const sdsl::bit_vector &);

bit_vector_dyn::bit_vector_dyn(std::initializer_list<bool> init)
      : bit_vector_dyn(pack_bits(std::vector<bool>(init)), init.size()) {}

uint64_t bit_vector_dyn::rank1(uint64_t id) const {
    return vector_.rank1(id);
}

uint64_t bit_vector_dyn::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= rank1(size() - 1));
    return vector_.select1(id - 1);
}

void bit_vector_dyn::set(uint64_t id, bool val) {
    vector_.set(id, val);
}

void bit_vector_dyn::setBitQuick(uint64_t id, bool val) {
    vector_.setBitQuick(id, val);
}

bool bit_vector_dyn::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

uint64_t bit_vector_dyn::get_int(uint64_t id, uint32_t width) const {
    assert(id + width <= size());
    uint64_t pos = id + width - 1;
    uint64_t word = 0;
    while (pos >= id) {
        word = (word << 1) + vector_[pos--];
    }
    return word;
}

void bit_vector_dyn::insertBit(uint64_t id, bool val) {
    vector_.insertBit(id, val);
}

void bit_vector_dyn::deleteBit(uint64_t id) {
    assert(size() > 0);
    vector_.deleteBit(id);
}

bool bit_vector_dyn::load(std::istream &in) {
    if (!in.good())
        return false;

    //TODO: catch reading errors
    vector_.deserialise(in);
    return true;
}

uint64_t bit_vector_dyn::serialize(std::ostream &out) const {
    vector_.serialise(out);
    // TODO
    return 0;
}

void bit_vector_dyn::call_ones(const std::function<void(uint64_t)> &callback) const {
    for (uint64_t i = 0; i < vector_.size(); ++i) {
        if (vector_[i])
            callback(i);
    }
}

///////////////////////////////////////////////
// bit_vector_stat, sdsl rank-select support //
///////////////////////////////////////////////

bit_vector_stat::bit_vector_stat(uint64_t size, bool value)
      : vector_(size, value) {
    if (value)
        num_set_bits_ = size;
}

bit_vector_stat::bit_vector_stat(const std::vector<bool> &v)
      : vector_(v.size(), 0) {
    for (uint64_t i = 0; i < v.size(); ++i) {
        if (v[i]) {
            vector_[i] = 1;
            num_set_bits_++;
        }
    }
}

bit_vector_stat::bit_vector_stat(const bit_vector_stat &other) {
    *this = other;
}

uint64_t count_num_set_bits(const sdsl::bit_vector &vector) {
    uint64_t count = 0;
    uint64_t i = 0;
    for (; i + 64 <= vector.size(); i += 64) {
        count += sdsl::bits::cnt(vector.get_int(i));
    }
    for (; i < vector.size(); ++i) {
        if (vector[i])
            count++;
    }
    return count;
}

bit_vector_stat::
bit_vector_stat(const std::function<void(const std::function<void(uint64_t)>&)> &call_ones,
                uint64_t size,
                uint64_t)
      : vector_(size, false) {
    call_ones([&](uint64_t pos) { vector_[pos] = true; });
}

bit_vector_stat::bit_vector_stat(sdsl::bit_vector&& vector) noexcept
      : vector_(std::move(vector)),
        num_set_bits_(count_num_set_bits(vector_)) {}

bit_vector_stat::bit_vector_stat(bit_vector_stat&& other) noexcept {
    *this = std::move(other);
}

bit_vector_stat::bit_vector_stat(std::initializer_list<bool> init)
      : bit_vector_stat(sdsl::bit_vector(init)) {}

bit_vector_stat& bit_vector_stat::operator=(const bit_vector_stat &other) {
    vector_ = other.vector_;
    num_set_bits_ = other.num_set_bits_;
    if (!other.requires_update_) {
        rk_ = other.rk_;
        rk_.set_vector(&vector_);
        slct_ = other.slct_;
        slct_.set_vector(&vector_);
        requires_update_ = false;
    }
    return *this;
}

bit_vector_stat& bit_vector_stat::operator=(bit_vector_stat&& other) noexcept {
    vector_ = std::move(other.vector_);
    num_set_bits_ = other.num_set_bits_;
    if (!other.requires_update_) {
        rk_ = std::move(other.rk_);
        rk_.set_vector(&vector_);
        slct_ = std::move(other.slct_);
        slct_.set_vector(&vector_);
        requires_update_ = false;
    }
    return *this;
}

uint64_t bit_vector_stat::rank1(uint64_t id) const {
    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();
    //the rank method in SDSL does not include id in the count
    return rk_(id >= this->size() ? this->size() : id + 1);
}

uint64_t bit_vector_stat::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits_);
    assert(num_set_bits_ == rank1(size() - 1));

    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();
    return slct_(id);
}

bool bit_vector_stat::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

uint64_t bit_vector_stat::get_int(uint64_t id, uint32_t width) const {
    return vector_.get_int(id, width);
}

void bit_vector_stat::set(uint64_t id, bool val) {
    if (vector_[id] == val)
        return;

    if (val) {
        num_set_bits_++;
    } else {
        num_set_bits_--;
    }

    vector_[id] = val;
    requires_update_ = true;
}

void bit_vector_stat::insertBit(uint64_t id, bool val) {
    assert(id <= size());

    if (val)
        num_set_bits_++;

    vector_.resize(size() + 1);
    std::copy_backward(vector_.begin() + id, vector_.end() - 1, vector_.end());

    vector_[id] = val;
    requires_update_ = true;
}

void bit_vector_stat::deleteBit(uint64_t id) {
    assert(size() > 0);
    assert(id < size());

    if (vector_[id])
        num_set_bits_--;

    std::copy(vector_.begin() + id + 1, vector_.end(), vector_.begin() + id);

    vector_.resize(vector_.size() - 1);
    requires_update_ = true;
}

uint64_t bit_vector_stat::serialize(std::ostream &out) const {
    uint64_t num_written_bytes = 0;
    num_written_bytes += vector_.serialize(out);

    if (requires_update_)
        const_cast<bit_vector_stat*>(this)->init_rs();

    num_written_bytes += serialize_number(out, num_set_bits_);
    num_written_bytes += rk_.serialize(out);
    num_written_bytes += slct_.serialize(out);

    return num_written_bytes;
}

bool bit_vector_stat::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);

        try {
            num_set_bits_ = load_number(in);
            rk_.load(in, &vector_);
            slct_.load(in, &vector_);
            requires_update_ = false;
        } catch (...) {
            std::cerr << "Warning: Loading from file without bit_vector rank"
                      << " and select support dumped. Reserialize to"
                      << " make the loading faster." << std::endl;

            num_set_bits_ = count_num_set_bits(vector_);
            requires_update_ = true;
            init_rs();
        }
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_stat." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

void call_ones(const sdsl::bit_vector &vector,
               const std::function<void(uint64_t)> &callback) {
    uint64_t j = 64;
    uint64_t i = 0;
    uint64_t word;
    for (; j <= vector.size(); j += 64) {
        word = vector.get_int(i);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (vector[i])
                callback(i);
        }
    }
    for (; i < vector.size(); ++i) {
        if (vector[i])
            callback(i);
    }
}

void call_zeros(const sdsl::bit_vector &vector,
                const std::function<void(uint64_t)> &callback) {
    uint64_t j = 64;
    uint64_t i = 0;
    uint64_t word;
    for (; j <= vector.size(); j += 64) {
        word = ~vector.get_int(i);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for (; i < j; ++i) {
            if (!vector[i])
                callback(i);
        }
    }
    for (; i < vector.size(); ++i) {
        if (!vector[i])
            callback(i);
    }
}

void bit_vector_stat::call_ones(const std::function<void(uint64_t)> &callback) const {
    if (num_set_bits_ * 2 < size()) {
        ::call_ones(vector_, callback);
    } else {
        for (uint64_t i = 0; i < vector_.size(); ++i) {
            if (vector_[i])
                callback(i);
        }
    }
}

void bit_vector_stat::init_rs() {
    std::unique_lock<std::mutex> lock(mu_);

    if (!requires_update_)
        return;

    rk_ = sdsl::rank_support_v5<>(&vector_);
    slct_ = sdsl::select_support_mcl<>(&vector_);

    requires_update_ = false;
}

////////////////////////////////////////////////////////////////
// bit_vector_small, sdsl hybrid method                       //
////////////////////////////////////////////////////////////////

// TODO: proper checking for sizes
using default_bit_vector = bit_vector_rrr<>;

bit_vector_small::bit_vector_small(uint64_t size, bool value)
      : vector_(new default_bit_vector(size, value)) {}

template <class BitVector>
bit_vector_small::bit_vector_small(const BitVector &vector)
      : vector_(new default_bit_vector(
            bit_vector_stat(vector).convert_to<sdsl::bit_vector>())) {}

template bit_vector_small::bit_vector_small(const std::vector<bool> &);

bit_vector_small::bit_vector_small(const sdsl::bit_vector &vector)
      : vector_(new default_bit_vector(vector)) {}

bit_vector_small::bit_vector_small(const bit_vector_small &other) {
    *this = other;
}

bit_vector_small
::bit_vector_small(const std::function<void(const std::function<void(uint64_t)>&)> &call_ones,
                 uint64_t size,
                 uint64_t num_set_bits)
      : vector_(new default_bit_vector(
            bit_vector_stat(call_ones, size, num_set_bits).convert_to<sdsl::bit_vector>())) {}

bit_vector_small& bit_vector_small::operator=(const bit_vector_small &other) {
    vector_ = std::make_unique<default_bit_vector>(*dynamic_cast<default_bit_vector*>(other.vector_.get()));
    return *this;
}

bit_vector_small::bit_vector_small(std::initializer_list<bool> init)
      : vector_(new default_bit_vector(init)) {}

uint64_t bit_vector_small::rank1(uint64_t id) const {
    return vector_->rank1(id);
}

uint64_t bit_vector_small::select1(uint64_t id) const {
    return vector_->select1(id);
}

void bit_vector_small::set(uint64_t id, bool val) {
    vector_->set(id, val);
}

bool bit_vector_small::operator[](uint64_t id) const {
    return vector_->operator[](id);
}

uint64_t bit_vector_small::get_int(uint64_t id, uint32_t width) const {
    return vector_->get_int(id, width);
}

void bit_vector_small::insertBit(uint64_t id, bool val) {
    vector_->insertBit(id, val);
}

void bit_vector_small::deleteBit(uint64_t id) {
    vector_->deleteBit(id);
}

bool bit_vector_small::load(std::istream &in) {
    return vector_->load(in);
}

uint64_t bit_vector_small::serialize(std::ostream &out) const {
    return vector_->serialize(out);
}

uint64_t bit_vector_small::size() const {
    return vector_->size();
}

std::vector<bool> bit_vector_small::to_vector() const {
    return vector_->to_vector();
}

void bit_vector_small::call_ones(const std::function<void(uint64_t)> &callback) const {
    vector_->call_ones(callback);
}

bool bit_vector_small::is_inverted() const {
    if (dynamic_cast<bit_vector_sd*>(vector_.get()))
        return dynamic_cast<bit_vector_sd*>(vector_.get())->is_inverted();

    return false;
}

////////////////////////////////////////////////////////////////
// bit_vector_sd, sdsl compressed with rank-select support    //
////////////////////////////////////////////////////////////////

bit_vector_sd::bit_vector_sd(uint64_t size, bool value)
      : inverted_(value && size) {
    sdsl::sd_vector_builder builder(size, 0);
    vector_ = decltype(vector_)(builder);
    rk1_ = decltype(rk1_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    slct0_ = decltype(slct0_)(&vector_);
}

bit_vector_sd::bit_vector_sd(const sdsl::bit_vector &vector) {
    // check if it needs to be inverted
    uint64_t num_set_bits = count_num_set_bits(vector);

    if (num_set_bits <= vector.size() / 2) {
        // vector is sparse, no need to invert
        vector_ = sdsl::sd_vector<>(vector);
        inverted_ = false;
    } else {
        // invert
        sdsl::sd_vector_builder builder(vector.size(), vector.size() - num_set_bits);
        ::call_zeros(vector, [&builder](uint64_t i) { builder.set(i); });
        vector_ = sdsl::sd_vector<>(builder);
        inverted_ = true;
    }
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

template <class BitVector>
bit_vector_sd::bit_vector_sd(const BitVector &vector) {
    // check if it needs to be inverted
    uint64_t num_set_bits = std::count(vector.begin(), vector.end(), true);

    if (num_set_bits <= vector.size() / 2) {
        // vector is sparse, no need to invert
        sdsl::sd_vector_builder builder(vector.size(), num_set_bits);
        for (uint64_t i = 0; i < vector.size(); ++i) {
            if (vector[i])
                builder.set(i);
        }
        vector_ = sdsl::sd_vector<>(builder);
        inverted_ = false;
    } else {
        // invert
        sdsl::sd_vector_builder builder(vector.size(), vector.size() - num_set_bits);
        for (uint64_t i = 0; i < vector.size(); ++i) {
            if (!vector[i])
                builder.set(i);
        }
        vector_ = sdsl::sd_vector<>(builder);
        inverted_ = true;
    }
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

template bit_vector_sd::bit_vector_sd(const std::vector<bool> &);

bit_vector_sd::bit_vector_sd(const bit_vector_sd &other) {
    *this = other;
}

bit_vector_sd::bit_vector_sd(bit_vector_sd&& other) noexcept {
    *this = std::move(other);
}

bit_vector_sd
::bit_vector_sd(const std::function<void(const std::function<void(uint64_t)>&)> &call_ones,
                 uint64_t size,
                 uint64_t num_set_bits)
      : inverted_(num_set_bits > size / 2) {
    sdsl::sd_vector_builder builder(size,
                                    !inverted_
                                    ? num_set_bits
                                    : size - num_set_bits);
    if (inverted_) {
        uint64_t last_pos = 0;
        call_ones([&](uint64_t pos) {
            while (last_pos < pos) {
                builder.set(last_pos++);
            }
            ++last_pos;
        });
        while (last_pos < size) {
            builder.set(last_pos++);
        }
    } else {
        call_ones([&](uint64_t pos) { builder.set(pos); });
    }

    assert(builder.items() == builder.capacity());

    vector_ = sdsl::sd_vector<>(builder);
    slct0_ = decltype(slct0_)(&vector_);
    slct1_ = decltype(slct1_)(&vector_);
    rk1_ = decltype(rk1_)(&vector_);
}

bit_vector_sd::bit_vector_sd(std::initializer_list<bool> init)
      : bit_vector_sd(sdsl::bit_vector(init)) {}

bit_vector_sd& bit_vector_sd::operator=(const bit_vector_sd &other) {
    inverted_ = other.inverted_;
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

bit_vector_sd& bit_vector_sd::operator=(bit_vector_sd&& other) noexcept {
    inverted_ = other.inverted_;
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

uint64_t bit_vector_sd::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    size_t idx = id >= this->size() ? this->size() : id + 1;
    return !inverted_ ? rk1_(idx)
                      : idx - rk1_(idx);
}

uint64_t bit_vector_sd::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return !inverted_ ? slct1_(id) : slct0_(id);
}

void bit_vector_sd::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_sd::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id] != inverted_;
}

uint64_t bit_vector_sd::get_int(uint64_t id, uint32_t width) const {
    if (inverted_)
        return ~vector_.get_int(id, width) & sdsl::bits::lo_set[width];

    return vector_.get_int(id, width);
}

void bit_vector_sd::insertBit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

void bit_vector_sd::deleteBit(uint64_t) {
    throw std::runtime_error("Not supported");
}

bool bit_vector_sd::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        inverted_ = in.get();
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_sd." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

uint64_t bit_vector_sd::serialize(std::ostream &out) const {
    uint64_t num_written_bytes = 0;
    num_written_bytes += vector_.serialize(out);
    if (!out.put(inverted_).good())
        throw std::ofstream::failure("Error when dumping bit_vector_sd");

    return num_written_bytes + 1;
}

std::vector<bool> bit_vector_sd::to_vector() const {
    std::vector<bool> vector(size(), inverted_);
    uint64_t max_rank = rk1_(size());
    for (uint64_t i = 1; i <= max_rank; ++i) {
        vector[slct1_(i)] = !inverted_;
    }
    return vector;
}

void bit_vector_sd::call_ones(const std::function<void(uint64_t)> &callback) const {
    if (!inverted_) {
        // sparse
        uint64_t num_ones = num_set_bits();
        for (uint64_t i = 1; i <= num_ones; ++i) {
            callback(select1(i));
        }
    } else {
        // dense, vector_ keeps positions of zeros
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = rk1_(size());
        for (uint64_t r = 1; r <= num_zeros; ++r) {
            zero_pos = slct1_(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < size()) {
            callback(one_pos++);
        }
    }
}



////////////////////////////////////////////////////////////////
//  bit_vector_rrr, sdsl compressed with rank-select support  //
////////////////////////////////////////////////////////////////

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(uint64_t size, bool value)
      : vector_(sdsl::bit_vector(size, value)),
        rk1_(&vector_),
        slct1_(&vector_),
        slct0_(&vector_) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(const sdsl::bit_vector &vector)
      : vector_(vector),
        rk1_(&vector_),
        slct1_(&vector_),
        slct0_(&vector_) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(const bit_vector_rrr<log_block_size> &other) {
    *this = other;
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(bit_vector_rrr<log_block_size>&& other) noexcept {
    *this = std::move(other);
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>
::bit_vector_rrr(std::initializer_list<bool> init)
      : bit_vector_rrr(sdsl::bit_vector(init)) {}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>&
bit_vector_rrr<log_block_size>::operator=(const bit_vector_rrr<log_block_size> &other) {
    vector_ = other.vector_;
    rk1_ = other.rk1_;
    rk1_.set_vector(&vector_);
    slct1_ = other.slct1_;
    slct1_.set_vector(&vector_);
    slct0_ = other.slct0_;
    slct0_.set_vector(&vector_);
    return *this;
}

template <size_t log_block_size>
bit_vector_rrr<log_block_size>&
bit_vector_rrr<log_block_size>::operator=(bit_vector_rrr<log_block_size>&& other) noexcept {
    vector_ = std::move(other.vector_);
    rk1_ = std::move(other.rk1_);
    rk1_.set_vector(&vector_);
    slct1_ = std::move(other.slct1_);
    slct1_.set_vector(&vector_);
    slct0_ = std::move(other.slct0_);
    slct0_.set_vector(&vector_);
    return *this;
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::rank1(uint64_t id) const {
    //the rank method in SDSL does not include id in the count
    return rk1_(id >= this->size() ? this->size() : id + 1);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::select1(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return slct1_(id);
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::select0(uint64_t id) const {
    assert(id > 0 && size() > 0 && id <= size() - num_set_bits());
    assert(num_set_bits() == rank1(size() - 1));

    return slct0_(id);
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::set(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
bool bit_vector_rrr<log_block_size>::operator[](uint64_t id) const {
    assert(id < size());
    return vector_[id];
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>
::get_int(uint64_t id, uint32_t width) const {
    return vector_.get_int(id, width);
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::insertBit(uint64_t, bool) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>::deleteBit(uint64_t) {
    throw std::runtime_error("Not supported");
}

template <size_t log_block_size>
bool bit_vector_rrr<log_block_size>::load(std::istream &in) {
    if (!in.good())
        return false;

    try {
        vector_.load(in);
        if (!in.good())
            return false;
        rk1_ = decltype(rk1_)(&vector_);
        slct1_ = decltype(slct1_)(&vector_);
        slct0_ = decltype(slct0_)(&vector_);
        return true;
    } catch (const std::bad_alloc &exception) {
        std::cerr << "ERROR: Not enough memory to load bit_vector_rrr." << std::endl;
        return false;
    } catch (...) {
        return false;
    }
}

template <size_t log_block_size>
uint64_t bit_vector_rrr<log_block_size>::serialize(std::ostream &out) const {
    uint64_t num_written_bytes = vector_.serialize(out);

    if (!out.good())
        throw std::ofstream::failure("Error when dumping bit_vector_rrr");

    return num_written_bytes;
}

template <size_t log_block_size>
std::vector<bool> bit_vector_rrr<log_block_size>::to_vector() const {
    if (2 * num_set_bits() < size()) {
        // sparse
        std::vector<bool> vector(size(), 0);
        uint64_t max_rank = rank1(size() - 1);
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct1_(i)] = 1;
        }
        return vector;
    } else {
        // dense
        std::vector<bool> vector(size(), 1);
        uint64_t max_rank = rank0(size() - 1);
        for (uint64_t i = 1; i <= max_rank; ++i) {
            vector[slct0_(i)] = 0;
        }
        return vector;
    }
}

template <size_t log_block_size>
void bit_vector_rrr<log_block_size>
::call_ones(const std::function<void(uint64_t)> &callback) const {
    if (2 * num_set_bits() < size()) {
        // sparse
        uint64_t num_ones = rank1(size() - 1);
        for (uint64_t i = 1; i <= num_ones; ++i) {
            callback(select1(i));
        }
    } else {
        // dense
        uint64_t one_pos = 0;
        uint64_t zero_pos = 0;
        uint64_t num_zeros = rank0(size() - 1);
        for (uint64_t r = 1; r <= num_zeros; ++r) {
            zero_pos = select0(r);
            while (one_pos < zero_pos) {
                callback(one_pos++);
            }
            one_pos++;
        }
        while (one_pos < size()) {
            callback(one_pos++);
        }
    }
}

template class bit_vector_rrr<3>;
template class bit_vector_rrr<8>;
template class bit_vector_rrr<15>;
template class bit_vector_rrr<31>;
template class bit_vector_rrr<63>;
template class bit_vector_rrr<127>;
template class bit_vector_rrr<255>;
