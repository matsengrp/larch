#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Id, typename DV>
struct ElementView;

template <typename Id, typename DV>
struct element_view_base {
  using type = std::conditional_t<
      DV::is_mutable,
      typename DV::template MutableElementViewBase<Id, ElementView<Id, DV>>,
      typename DV::template ConstElementViewBase<Id, ElementView<Id, DV>>>;
};
template <typename Id, typename DV>
using element_view_base_t = typename element_view_base<Id, DV>::type;

/**
 * A view into a single node or edge of a DAG.
 */
template <typename Id, typename DV>
struct ElementView : element_view_base_t<Id, std::decay_t<DV>> {
 public:
  /**
   * This operator= is used for setting specailly handles per-element
   * features, like Deduplicate.
   */
  using element_view_base_t<Id, DV>::operator=;

  ElementView(DV dag_view, Id id);

  operator Id() const;

  DV GetDAG() const;
  Id GetId() const;

  template <typename F>
  auto& GetFeatureStorage() const;

  template <typename F>
  auto& GetFeatureExtraStorage() const;

 private:
  std::decay_t<DV> dag_view_;
  Id id_;
};