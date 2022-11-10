#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
class DefaultContainerBase {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<std::vector<Features>...>;

  template <typename Feature>
  static inline constexpr bool contains_feature =
      tuple_contians_v<FeaturesType, std::vector<Feature>>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<DAG::is_mutable,
                                  FeatureWriter<Features, ViewType<DAG>>,
                                  FeatureReader<Features, ViewType<DAG>>>... {};

  template <typename Feature>
  Feature& GetFeatureAt(Id id);

  template <typename Feature>
  const Feature& GetFeatureAt(Id id) const;

  template <typename Feature>
  typename Feature::GlobalData& GetGlobalData();

  template <typename Feature>
  const typename Feature::GlobalData& GetGlobalData() const;

 private:
  DAG_VIEW_FRIENDS;
  [[no_unique_address]] std::tuple<std::vector<Features>...> features_;
  [[no_unique_address]] std::tuple<
      std::conditional_t<has_global_data<Features>, typename Features::GlobalData,
                         NoGlobalData<Features>>...>
      global_data_;
};