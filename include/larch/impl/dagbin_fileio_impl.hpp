
bool DagbinFileIO::IsFileDagbinFormat(std::string_view path) {
  std::ifstream infile;
  infile.open(std::string{path}, std::ios::binary);
  auto is_dagbin_file = CheckMagicNumber(infile, false);
  infile.close();
  return is_dagbin_file;
}

MADAGStorage<> DagbinFileIO::ReadDAG(std::string_view path) {
  std::ifstream infile;
  std::optional<Header> header = std::nullopt;
  MADAGStorage<> dag_storage = MADAGStorage<>::EmptyDefault();
  auto dag = dag_storage.View();

  infile.open(std::string{path}, std::ios::binary);
  auto labeled_linked_list = ReadLabeledLinkedList(infile);
  infile.close();

  auto read_section = [&infile, &dag, &header](SectionId section_id) {
    switch (section_id) {
      case SectionId::Header:
        header = ReadHeader(infile, dag);
        break;
      case SectionId::RefSeq:
        ReadReferenceSequence(infile, dag);
        break;
      case SectionId::Nodes:
        ReadNodes(infile, dag);
        break;
      case SectionId::Edges:
        ReadEdges(infile, dag);
        break;
      default:
        assert(false && "ERROR: Invalid section_id");
    }
  };

  infile.open(std::string{path}, std::ios::binary);
  CheckMagicNumber(infile);
  for (auto [section_begpos, section_id] : labeled_linked_list) {
    [[maybe_unused]] auto true_section_id = ReadData<SectionId>(infile);
    Assert(section_id == true_section_id);
    auto section_endpos = ReadData<std::streampos>(infile);

    read_section(section_id);
    [[maybe_unused]] auto true_section_endpos = infile.tellg();
    Assert(section_endpos == true_section_endpos);
  }
  infile.close();

  dag.BuildConnections();
  dag.AssertUA();
  // dag.RecomputeCompactGenomes();

  for (auto node : dag.GetNodes()) {
    if (node.IsLeaf() and not node.HaveSampleId()) {
      node = SampleId{node.GetCompactGenome().ToString()};
    }
  }
  Assert(header.has_value());
  Assert(dag.GetNodesCount() == header.value().node_count);
  Assert(dag.GetEdgesCount() == header.value().edge_count);
  Assert(dag.GetLeafsCount() == header.value().leaf_count);

  dag.GetRoot().Validate(true, true);
  return dag_storage;
}

template <typename DAG>
void DagbinFileIO::WriteDAG(DAG dag, std::string_view path) {
  std::ofstream outfile(std::string{path}, std::ios::binary);
  std::vector<std::streampos> offsets;

  // write magic number
  WriteMagicNumber(outfile);

  // write header
  offsets.push_back(outfile.tellp());
  WriteData(outfile, SectionId::Header);
  WriteData(outfile, std::streampos{-1});
  WriteHeader(outfile, dag);

  // write reference sequence
  offsets.push_back(outfile.tellp());
  WriteData(outfile, SectionId::RefSeq);
  WriteData(outfile, std::streampos{-1});
  WriteReferenceSequence(outfile, dag);

  // write nodes
  for (size_t i = 0; i < dag.GetNodesCount(); i += batch_size) {
    offsets.push_back(outfile.tellp());
    WriteData(outfile, SectionId::Nodes);
    WriteData(outfile, std::streampos{-1});
    WriteNodes(outfile, dag, i, {i + batch_size});
  }

  // write edges
  for (size_t i = 0; i < dag.GetEdgesCount(); i += batch_size) {
    offsets.push_back(outfile.tellp());
    WriteData(outfile, SectionId::Edges);
    WriteData(outfile, std::streampos{-1});
    WriteEdges(outfile, dag, i, {i + batch_size});
  }

  // build linked list of offsets
  offsets.push_back(outfile.tellp());
  WriteLinkedList(outfile, offsets);

  outfile.close();
}

template <typename DAG>
void DagbinFileIO::AppendDAG(DAG dag, std::string_view path) {
  std::fstream file(std::string{path}, std::ios::in | std::ios::out | std::ios::binary);
  std::vector<std::streampos> offsets;

  // Read old header
  file.seekg(0, std::ios::beg);
  CheckMagicNumber(file);
  auto section_id = ReadData<SectionId>(file);
  auto section_offset = ReadData<std::streampos>(file);
  Assert(section_id == SectionId::Header);
  Header old_header = ReadHeader(file, dag);

  // Overwrite header
  file.seekp(0, std::ios::beg);
  WriteMagicNumber(file);
  WriteData(file, section_id);
  WriteData(file, section_offset);
  WriteHeader(file, dag);

  // Prepare to append
  file.clear();
  file.seekp(0, std::ios::end);

  // Write new nodes
  for (size_t i = old_header.node_count; i < dag.GetNodesCount(); i += batch_size) {
    offsets.push_back(file.tellp());
    WriteData(file, SectionId::Nodes);
    WriteData(file, std::streampos{42});
    WriteNodes(file, dag, i, {i + batch_size});
  }

  // Write new edges
  for (size_t i = old_header.edge_count; i < dag.GetEdgesCount(); i += batch_size) {
    offsets.push_back(file.tellp());
    WriteData(file, SectionId::Edges);
    WriteData(file, std::streampos{42});
    WriteEdges(file, dag, i, {i + batch_size});
  }

  // Update linked list of offsets
  offsets.push_back(file.tellp());
  WriteLinkedList(file, offsets);

  file.close();
}

template <typename T, typename iostream>
T DagbinFileIO::ReadData(iostream &infile) {
  T data;
  infile.read(reinterpret_cast<char *>(&data), sizeof(data));
  return data;
}

template <typename T, typename iostream>
void DagbinFileIO::WriteData(iostream &outfile, const T data) {
  outfile.write(reinterpret_cast<const char *>(&data), sizeof(data));
}

template <typename iostream>
std::string DagbinFileIO::ReadString(iostream &infile) {
  auto str_len = ReadData<size_t>(infile);
  std::string str(str_len, '\0');
  infile.read(&str[0], str_len);
  return str;
}

template <typename iostream>
void DagbinFileIO::WriteString(iostream &outfile, const std::string &str) {
  WriteData(outfile, str.size());
  outfile.write(str.c_str(), str.size());
}

template <typename iostream>
std::vector<std::pair<std::streampos, DagbinFileIO::SectionId>>
DagbinFileIO::ReadLabeledLinkedList(iostream &infile) {
  std::vector<std::pair<std::streampos, SectionId>> offsets;
  auto start_pos = infile.tellg();

  SectionId section_id;
  std::streampos section_offset;
  infile.seekg(0, std::ios::beg);
  CheckMagicNumber(infile);
  section_offset = infile.tellg();
  while (!infile.eof()) {
    infile.seekg(section_offset);
    section_id = ReadData<SectionId>(infile);
    if (infile.eof()) {
      break;
    }
    offsets.push_back({section_offset, section_id});
    section_offset = ReadData<std::streampos>(infile);
  }

  infile.seekg(start_pos);
  return offsets;
}

template <typename iostream>
std::vector<std::streampos> DagbinFileIO::ReadLinkedList(iostream &infile) {
  std::vector<std::streampos> offsets;
  auto start_pos = infile.tellg();

  std::streampos section_offset;
  infile.seekg(0, std::ios::beg);
  CheckMagicNumber(infile);
  section_offset = infile.tellg();
  while (!infile.eof()) {
    offsets.push_back(section_offset);
    infile.seekg(section_offset + static_cast<std::streamoff>(1));
    section_offset = ReadData<std::streampos>(infile);
  }

  infile.seekg(start_pos);
  return offsets;
}

template <typename iostream>
void DagbinFileIO::WriteLinkedList(iostream &outfile,
                                   const std::vector<std::streampos> &offsets) {
  auto start_pos = outfile.tellp();

  for (size_t i = 1; i < offsets.size(); i++) {
    auto cur_offset = offsets[i - 1] + static_cast<std::streamoff>(1);
    auto nxt_offset = offsets[i];
    outfile.seekp(cur_offset);
    WriteData(outfile, nxt_offset);
  }

  outfile.seekp(start_pos);
}

template <typename iostream>
bool DagbinFileIO::CheckMagicNumber(iostream &infile, bool do_assert) {
  std::vector<unsigned char> magic_number(MAGIC_NUMBER.size());
  infile.read(reinterpret_cast<char *>(magic_number.data()), MAGIC_NUMBER.size());
  if (do_assert) Assert(magic_number == MAGIC_NUMBER)
  return (magic_number == MAGIC_NUMBER);
}

template <typename iostream>
void DagbinFileIO::WriteMagicNumber(iostream &outfile) {
  outfile.write(reinterpret_cast<const char *>(MAGIC_NUMBER.data()),
                MAGIC_NUMBER.size());
}

template <typename iostream, typename DAG>
DagbinFileIO::Header DagbinFileIO::ReadHeader(iostream &infile, DAG dag) {
  std::ignore = dag;
  return ReadData<Header>(infile);
}

template <typename iostream, typename DAG>
void DagbinFileIO::WriteHeader(iostream &outfile, const DAG dag) {
  Header header = {dag.GetNodesCount(), dag.GetEdgesCount(), dag.GetLeafsCount()};
  WriteData<Header>(outfile, header);
}

template <typename iostream, typename DAG>
void DagbinFileIO::ReadReferenceSequence(iostream &infile, DAG dag) {
  dag.SetReferenceSequence(ReadString(infile));
}

template <typename iostream, typename DAG>
void DagbinFileIO::WriteReferenceSequence(iostream &outfile, const DAG dag) {
  WriteString(outfile, dag.GetReferenceSequence());
}

template <typename iostream, typename DAG>
void DagbinFileIO::ReadNodes(iostream &infile, DAG dag) {
  // write node count
  auto node_count = ReadData<size_t>(infile);

  // read nodes
  for (size_t i = 0; i < node_count; i++) {
    auto node_id = ReadData<size_t>(infile);
    auto new_node = dag.AddNode({node_id});
    auto node_is_leaf = ReadData<bool>(infile);
    if (node_is_leaf) {
      auto sample_id = ReadString(infile);
      new_node = SampleId{sample_id};
    }
  }
}

template <typename iostream, typename DAG>
void DagbinFileIO::WriteNodes(iostream &outfile, const DAG dag,
                              std::optional<size_t> min_id_opt,
                              std::optional<size_t> max_id_opt) {
  size_t min_id = min_id_opt.value_or(0);
  min_id = std::max(size_t{0}, min_id);
  size_t max_id = max_id_opt.value_or(dag.GetNodesCount());
  max_id = std::min(dag.GetNodesCount(), max_id);
  std::streampos start_pos, end_pos;
  start_pos = outfile.tellp();
  // write node count
  WriteData(outfile, max_id - min_id);

  // write nodes
  size_t node_count = 0;
  for (size_t id = min_id; id < max_id; id++) {
    auto node = dag.Get(NodeId{id});

    WriteData(outfile, node.GetId());
    WriteData(outfile, node.IsLeaf());
    if (node.IsLeaf()) {
      WriteString(outfile, node.GetSampleId().value());
    }
    node_count++;
  }

  // update node count
  end_pos = outfile.tellp();
  outfile.seekp(start_pos);
  WriteData(outfile, node_count);
  outfile.seekp(end_pos);
}

template <typename iostream, typename DAG>
void DagbinFileIO::ReadEdges(iostream &infile, DAG dag) {
  // read edge_count
  auto edge_count = ReadData<size_t>(infile);

  // read edges
  for (size_t i = 0; i < edge_count; i++) {
    auto edge_id = ReadData<size_t>(infile);
    auto parent_id = ReadData<size_t>(infile);
    auto child_id = ReadData<size_t>(infile);
    auto parent_clade = ReadData<size_t>(infile);
    auto edge = dag.AddEdge({edge_id}, {parent_id}, {child_id}, {parent_clade});

    // read mutations
    EdgeMutations mutations;
    auto mutation_count = ReadData<size_t>(infile);
    for (size_t j = 0; j < mutation_count; j++) {
      auto pos = ReadData<size_t>(infile);
      auto nuc_first = ReadData<char>(infile);
      auto nuc_second = ReadData<char>(infile);
      mutations[{pos}] = {{nuc_first}, {nuc_second}};
    }
    edge.SetEdgeMutations(std::move(mutations));
  }
}

template <typename iostream, typename DAG>
void DagbinFileIO::WriteEdges(iostream &outfile, const DAG dag,
                              std::optional<size_t> min_id_opt,
                              std::optional<size_t> max_id_opt) {
  size_t min_id = min_id_opt.value_or(0);
  min_id = std::max(size_t{0}, min_id);
  size_t max_id = max_id_opt.value_or(dag.GetEdgesCount());
  max_id = std::min(dag.GetEdgesCount(), max_id);
  std::streampos start_pos, end_pos;
  start_pos = outfile.tellp();
  // write edge count
  WriteData(outfile, max_id - min_id);

  // write edges
  size_t edge_count = 0;
  for (size_t id = min_id; id < max_id; id++) {
    auto edge = dag.Get(EdgeId{id});
    WriteData(outfile, edge.GetId().value);
    WriteData(outfile, edge.GetParentId().value);
    WriteData(outfile, edge.GetChildId().value);
    WriteData(outfile, edge.GetClade().value);

    // write mutations
    WriteData(outfile, edge.GetEdgeMutations().size());
    for (auto [pos, nucs] : edge.GetEdgeMutations()) {
      WriteData(outfile, pos.value);
      WriteData(outfile, nucs.first.ToChar());
      WriteData(outfile, nucs.second.ToChar());
    }
    edge_count++;
  }

  // update edge count
  end_pos = outfile.tellp();
  outfile.seekp(start_pos);
  WriteData(outfile, edge_count);
  outfile.seekp(end_pos);
}
