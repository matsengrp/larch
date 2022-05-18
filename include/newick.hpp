#pragma once

#include <stack>
#include <vector>
#include <string>
#include <optional>

template <typename T, typename N, typename E>
void ParseNewick(const T& source, N&& on_node, E&& on_edge) {
  struct Node {
    size_t id;
    std::string label;
  };
  std::stack<std::vector<Node>> nodes;
  std::string label;
  size_t node_id = 0;

  auto begin_node = [&]() { nodes.push({}); };

  auto end_label = [&]() {
    if (!nodes.empty()) {
      for (auto& i : nodes.top()) {
        on_edge(node_id, i.id);
      }
      nodes.pop();
      if (!nodes.empty()) {
        nodes.top().push_back({node_id, label});
      }
    }
    std::string branch_length;
    bool have_branch_length = false;
    for (char i : label | ranges::views::reverse) {
      if (i == ':') {
        have_branch_length = true;
        break;
      }
      if ((i < '0' || i > '9') && i != '.') break;
      branch_length = i + branch_length;
    }
    if (have_branch_length) {
      label.erase(label.size() - branch_length.size() - 1);
      on_node(node_id++, label, std::stod(branch_length));
    } else {
      on_node(node_id++, label, std::nullopt);
    }
    label = "";
  };

  begin_node();
  for (auto i = source.begin(); i != source.end(); ++i) {
    switch (*i) {
      case '(':
        begin_node();
        break;
      case ',':
        end_label();
        begin_node();
        break;
      case ')':
      case ';':
        end_label();
        break;
      default:
        label += *i;
        break;
    }
  }
}
