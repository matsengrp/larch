#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * A view into a single node or edge of a DAG.
 */
template <typename Id, typename DAGViewType>
struct ElementView : DAGViewType::BaseType::template ElementViewBase<Id> {
 public:
  /**
   * This operator= is used for setting specailly handles per-element
   * features, like Deduplicate.
   */
  using DAGViewType::BaseType::template ElementViewBase<Id>::operator=;

  ElementView(DAGViewType dag_view, Id id);

  operator Id() const;

  DAGViewType GetDAG() const;
  Id GetId() const;

  template <typename Feature>
  auto& GetFeatureStorage() const;

  template <typename Feature>
  auto& GetFeatureExtraStorage() const;

 private:
  std::decay_t<DAGViewType> dag_view_;
  Id id_;
};