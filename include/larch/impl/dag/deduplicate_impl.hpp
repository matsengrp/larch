#include <memory>
#ifndef DAG_DEFINITIONS
#error "Don't include this header"
#endif

template <typename Feature>
template <typename Id>
Deduplicate<Feature>::GlobalData<Id>::GlobalData() {
  //   deduplicated_storage_.insert(Deduplicate<Feature>::Empty);
}

template <typename Feature>
template <typename Id>
const Feature& Deduplicate<Feature>::GlobalData<Id>::Get(
    const Deduplicate<Feature>& self, Id) const {
  return *self.feature_;
}

template <typename Feature>
template <typename Id>
void Deduplicate<Feature>::GlobalData<Id>::Set(Deduplicate<Feature>& self, Id id,
                                               Feature&& feature) {
  auto iter = deduplicated_storage_.insert(std::forward<Feature>(feature));
  self.feature_ = std::addressof(*iter.first);
}
