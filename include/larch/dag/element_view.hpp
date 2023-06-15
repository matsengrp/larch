#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * A view into a single node or edge of a DAG.
 */
template <Component C, typename DAGViewType>
struct ElementView : DAGViewType::BaseType::template ElementViewBase<C> {
 public:
  /**
   * This operator= is used for setting specailly handled per-element
   * features, like Deduplicate.
   */
  using DAGViewType::BaseType::template ElementViewBase<C>::operator=;

  template <typename Feature>
  static const bool contains_feature;

  ElementView(DAGViewType dag_view, Id<C> id);

  auto Const() const;

  operator Id<C>() const;

  DAGViewType GetDAG() const;
  Id<C> GetId() const;

  template <typename Feature>
  auto& GetFeatureStorage() const;

  template <typename Feature>
  auto& GetFeatureExtraStorage() const;

 private:
  std::decay_t<DAGViewType> dag_view_;
  Id<C> id_;
};

template <typename Id, typename DAGViewType>
ElementView(DAGViewType dag_view, Id id) -> ElementView<ComponentOf<Id>, DAGViewType>;
