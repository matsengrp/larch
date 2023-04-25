// ** Static objects

inline const MutationBase MutationBase::DNA::A{{0, 0}};
inline const MutationBase MutationBase::DNA::C{{0, 1}};
inline const MutationBase MutationBase::DNA::G{{1, 0}};
inline const MutationBase MutationBase::DNA::T{{1, 1}};
inline const std::map<MutationBase, char> MutationBase::DNA::mut_to_char_map = {
    {MutationBase::DNA::A, 'A'},
    {MutationBase::DNA::C, 'C'},
    {MutationBase::DNA::G, 'G'},
    {MutationBase::DNA::T, 'T'}};
inline const std::map<MutationBase, MutationBase> MutationBase::DNA::complement_map = {
    {MutationBase::DNA::A, MutationBase::DNA::T},
    {MutationBase::DNA::C, MutationBase::DNA::G},
    {MutationBase::DNA::G, MutationBase::DNA::C},
    {MutationBase::DNA::T, MutationBase::DNA::A}};

// ** Constructors

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

// ** Miscellaneous

MutationBase MutationBase::GetComplementaryBase() const {
  MutationBase m_out{{!value[0], !value[1]}};
  return m_out;
}

// ** I/O

char MutationBase::ToChar() const { return DNA::mut_to_char_map.find(value)->second; }

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

// ** Comparators

bool MutationBase::operator==(const MutationBase rhs) const {
  return value == rhs.value;
}

bool MutationBase::operator<(const MutationBase rhs) const { return value < rhs.value; }

bool MutationBase::operator==(const MutationBase::BitArray rhs) const {
  return value == rhs;
}

bool MutationBase::operator==(const char rhs) const { return ToChar() == rhs; }

bool MutationBase::operator!=(const char rhs) const { return ToChar() != rhs; }
