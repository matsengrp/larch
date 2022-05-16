#include "history_dag_edge_storage.hpp"
#include <cassert>

template <typename T>
EdgeView<T>::EdgeView(T dag, EdgeId id) : dag_{dag}, id_{id} {
    static_assert(std::is_same_v<T, HistoryDAG&> or
        std::is_same_v<T, const HistoryDAG&>);
	assert(id.value != NoId);
    assert(id.value < dag_.edges_.size());
}

template <typename T>
EdgeView<T>::operator Edge() const {
	return {dag_, id_};
}

template <typename T>
T EdgeView<T>::GetDAG() const { return dag_; }

template <typename T>
EdgeId EdgeView<T>::GetId() const { return id_; }

template <typename T>
typename EdgeView<T>::NodeType EdgeView<T>::GetParent() const {
    return {dag_, GetStorage().parent_};
}

template <typename T>
typename EdgeView<T>::NodeType EdgeView<T>::GetChild() const {
    return {dag_, GetStorage().child_};
}

template <typename T>
CladeIdx EdgeView<T>::GetClade() const { return GetStorage().clade_; }

template <typename T>
bool EdgeView<T>::IsRoot() const {
    return GetParent().IsRoot();
}

template <typename T>
bool EdgeView<T>::IsLeaf() const {
    return GetChild().IsLeaf();
}

template <typename T>
double EdgeView<T>::GetProbability() const {
    return GetStorage().probability_;
}

template <typename T>
const auto& EdgeView<T>::GetWeight() const {
    GetStorage().weight_;
}

template <typename T>
std::optional<EdgeView<T>> EdgeView<T>::FindNextSibling() const {
	auto parent = GetParent();
    auto child = GetChild();
    auto children = parent.GetChildren();
    for (auto i = children.begin(); i != children.end(); ++i) {
        if ((*i).GetChild().GetId() == child.GetId()) {
            if (++i != children.end()) {
                return *i;
            }
            break;
        }
    }
    return std::nullopt;
}

template <typename T>
const auto& EdgeView<T>::GetStorage() const {
	return dag_.edges_.at(id_.value);
}

template <typename T>
auto& EdgeView<T>::GetStorage() {
	return dag_.edges_.at(id_.value);
}

template <typename T>
inline bool operator==(EdgeView<T> lhs, EdgeView<T> rhs) {
	return std::addressof(lhs.dag_) == std::addressof(rhs.dag_) &&
		lhs.id_ == rhs.id_;
}
