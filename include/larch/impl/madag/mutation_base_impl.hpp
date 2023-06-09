
inline const char MutationBase::DNA::ambiguous_char = 'N';
inline const MutationBase MutationBase::DNA::A{{1, 0, 0, 0}};
inline const MutationBase MutationBase::DNA::C{{0, 1, 0, 0}};
inline const MutationBase MutationBase::DNA::G{{0, 0, 1, 0}};
inline const MutationBase MutationBase::DNA::T{{0, 0, 0, 1}};
inline const MutationBase MutationBase::DNA::N{{1, 1, 1, 1}};
inline const std::map<MutationBase, char> MutationBase::DNA::mut_to_char_map = {
    {MutationBase::DNA::A, 'A'},
    {MutationBase::DNA::C, 'C'},
    {MutationBase::DNA::G, 'G'},
    {MutationBase::DNA::T, 'T'},
    {MutationBase::DNA::N, MutationBase::DNA::ambiguous_char}};
inline const std::map<MutationBase, MutationBase> MutationBase::DNA::complement_map = {
    {MutationBase::DNA::A, MutationBase::DNA::T},
    {MutationBase::DNA::C, MutationBase::DNA::G},
    {MutationBase::DNA::G, MutationBase::DNA::C},
    {MutationBase::DNA::T, MutationBase::DNA::A},
    {MutationBase::DNA::N, MutationBase::DNA::N}};

MutationBase::MutationBase(const MutationBase::BitArray m_value_in) {
  value = m_value_in;
};

MutationBase::MutationBase(const char m_char_in) {
  for (const auto &[m_base, m_char] : DNA::mut_to_char_map) {
    if (m_char_in == m_char) {
      value = m_base.value;
      return;
    }
  }
  Fail("Error: Invalid char given for MutationBase constructor.");
}

bool MutationBase::IsAmbiguous() const {
  size_t count = 0;
  for (auto i : value) {
    count += i;
  }
  return count > 1;
}

MutationBase MutationBase::GetComplementaryBase() const {
  MutationBase m_out;
  for (size_t i = 0; i < BitCount; i++) {
    m_out.value[i] = value[BitCount - i - 1];
  }
  return m_out;
}

MutationBase MutationBase::GetFirstBase() const {
  MutationBase m_out;
  for (size_t i = 0; i < BitCount; i++) {
    if (value[i]) {
      m_out.value[i] = true;
      return m_out;
    }
  }
  Fail("Error: Cannot GetFirstBase() of empty MutationBase.");
}

MutationBase MutationBase::GetCommonBases(const MutationBase &rhs) const {
  MutationBase m_out;
  for (size_t i = 0; i < BitCount; i++) {
    if (value[i] & rhs.value[i]) {
      m_out.value[i] = true;
    }
  }
  return m_out;
}

MutationBase MutationBase::GetFirstCommonBase(const MutationBase &rhs) const {
  MutationBase m_out;
  for (size_t i = 0; i < BitCount; i++) {
    if (value[i] & rhs.value[i]) {
      m_out.value[i] = true;
      return m_out;
    }
  }
  Fail("Error: Cannot GetFirstCommonBase() of two disjoint MutationBases.");
}

bool MutationBase::IsCompatible(const MutationBase &rhs) const {
  for (size_t i = 0; i < BitCount; i++) {
    if (value[i] & rhs.value[i]) {
      return true;
    }
  }
  return false;
}

char MutationBase::ToChar() const {
  auto found_char = DNA::mut_to_char_map.find(value);
  if (found_char == DNA::mut_to_char_map.end()) {
    return DNA::ambiguous_char;
  }
  return DNA::mut_to_char_map.find(value)->second;
}

std::string MutationBase::ToString(std::vector<MutationBase> m_in) {
  std::string str_out = "";
  for (size_t i = 0; i < BitCount; i++) {
    str_out[i] += m_in[i].ToChar();
  }
  return str_out;
}

inline std::ostream &operator<<(std::ostream &os, const MutationBase m_in) {
  os << m_in.ToChar();
  return os;
}

inline std::ostream &operator<<(std::ostream &os,
                                const std::vector<MutationBase> m_in) {
  os << MutationBase::ToString(m_in);
  return os;
}

inline std::string &operator+=(std::string &str, const MutationBase m_in) {
  str += m_in.ToChar();
  return str;
}

bool MutationBase::operator==(const MutationBase &rhs) const {
  return value == rhs.value;
}

bool MutationBase::operator!=(const MutationBase &rhs) const {
  return value != rhs.value;
}

bool MutationBase::operator<(const MutationBase &rhs) const {
  return value < rhs.value;
}

bool MutationBase::operator==(const MutationBase::BitArray &rhs) const {
  return value == rhs;
}

bool MutationBase::operator!=(const MutationBase::BitArray &rhs) const {
  return value != rhs;
}

bool MutationBase::operator<(const MutationBase::BitArray &rhs) const {
  return value < rhs;
}

bool MutationBase::operator==(const char &rhs) const { return ToChar() == rhs; }

bool MutationBase::operator!=(const char &rhs) const { return ToChar() != rhs; }

bool MutationBase::operator<(const char &rhs) const { return ToChar() < rhs; }

inline bool operator==(const char &lhs, const MutationBase &rhs) {
  return lhs == rhs.ToChar();
}

inline bool operator!=(const char &lhs, const MutationBase &rhs) {
  return lhs != rhs.ToChar();
}

inline bool operator<(const char &lhs, const MutationBase &rhs) {
  return lhs < rhs.ToChar();
}

inline MutationBase::BitArray operator&(const MutationBase::BitArray &lhs,
                                        const MutationBase::BitArray &rhs) {
  MutationBase::BitArray res;
  size_t i = 0;
  for (auto &res_i : res) {
    res_i = (lhs[i] & rhs[i]);
    i++;
  }
  return res;
}
