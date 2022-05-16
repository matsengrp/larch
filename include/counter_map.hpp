#pragma once

#include <map>

template <typename T>
class CounterMap {
public:

    void Add(const T& item, size_t count = 1);
    size_t GetCount(const T& item) const;

private:
    std::map<T, size_t> counts_;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
void CounterMap<T>::Add(const T& item, size_t count) {
    counts_[item] += count;
}

template <typename T>
size_t CounterMap<T>::GetCount(const T& item) const {
    auto i = counts_.find(item);
    if (i == counts_.end()) return 0;
    return *i;
}
