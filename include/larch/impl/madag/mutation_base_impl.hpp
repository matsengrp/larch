inline const MutationBase MutationBase::DNA::A{{1, 0, 0, 0}};
inline const MutationBase MutationBase::DNA::C{{0, 1, 0, 0}};
inline const MutationBase MutationBase::DNA::G{{0, 0, 1, 0}};
inline const MutationBase MutationBase::DNA::T{{0, 0, 0, 1}};
inline const char MutationBase::DNA::ambiguous_char = 'N';
inline const std::map<MutationBase, char> MutationBase::DNA::mut_to_char_map = {
    {MutationBase::DNA::A, 'A'},
    {MutationBase::DNA::C, 'C'},
    {MutationBase::DNA::G, 'G'},
    {MutationBase::DNA::T, 'T'}};
inline const std::map<MutationBase, MutationBase>
    MutationBase::DNA::mut_to_complement_map = {
        {MutationBase::DNA::A, MutationBase::DNA::T},
        {MutationBase::DNA::C, MutationBase::DNA::G},
        {MutationBase::DNA::G, MutationBase::DNA::C},
        {MutationBase::DNA::T, MutationBase::DNA::A}};

MutationBase::MutationBase(const MutationBase::BitArray m_value) { value = m_value; };

MutationBase::MutationBase(const char m_char_in) {
  for (const auto &[m_base, m_char] : DNA::mut_to_char_map) {
    if (m_char_in == m_char) {
      value = m_base.value;
      return;
    }
  }
  Fail("ERROR: Invalid char given for MutationBase constructor.");
}

bool MutationBase::IsAmbiguous() const {
  auto count = std::count(value.begin(), value.end(), true);
  return count != 1;
}

bool MutationBase::HasCommonBase(MutationBase other) const {
  BitArray tmp;
  for (size_t i = 0; i < tmp.size(); i++) {
    tmp[i] = value[i] & other.value[i];
  }
  auto count = std::count(tmp.begin(), tmp.end(), true);
  return count > 0;
}

MutationBase MutationBase::GetComplementaryBase() const {
  if (IsAmbiguous()) {
    MutationBase m_out;
    m_out.value[0] = value[3];
    m_out.value[1] = value[2];
    m_out.value[2] = value[1];
    m_out.value[3] = value[0];
    return m_out;
  }
  MutationBase m_out{DNA::mut_to_complement_map.find(*this)->second};
  return m_out;
}

char MutationBase::ToChar() const {
  if (IsAmbiguous()) {
    return DNA::ambiguous_char;
  }
  return DNA::mut_to_char_map.find(value)->second;
}

std::string MutationBase::ToString(std::vector<MutationBase> m_in) {
  std::string str_out = "";
  for (size_t i = 0; i < m_in.size(); i++) {
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
