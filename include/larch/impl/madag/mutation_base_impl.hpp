
bool MutationBase::IsAmbiguous() const { return Count() != 1; }

MutationBase MutationBase::GetFirstBase() const { return GetFirst(); }

MutationBase MutationBase::GetCommonBases(const MutationBase& rhs) const {
  auto result = Common(rhs);
  if (not result.has_value()) {
    Fail("Disjoint MutationBases");
  }
  return *result;
}

MutationBase MutationBase::GetFirstCommonBase(const MutationBase& rhs) const {
  return GetCommonBases(rhs).GetFirst();
}

bool MutationBase::IsCompatible(const MutationBase& rhs) const {
  return Common(rhs).has_value();
}

std::string MutationBase::ToString(const std::vector<MutationBase>& m_in) {
  std::string str_out;
  for (const auto& i : m_in) {
    str_out += i.ToChar();
  }
  return str_out;
}

// inline std::ostream& operator<<(std::ostream& os, const MutationBase& m_in) {
//   os << m_in.ToChar();
//   return os;
// }

// inline std::ostream& operator<<(std::ostream& os,
//                                 const std::vector<MutationBase>& m_in) {
//   os << MutationBase::ToString(m_in);
//   return os;
// }

inline std::string& operator+=(std::string& str, const MutationBase& m_in) {
  str += m_in.ToChar();
  return str;
}

bool MutationBase::operator==(const MutationBase& rhs) const {
  return value == rhs.value;
}

bool MutationBase::operator!=(const MutationBase& rhs) const {
  return value != rhs.value;
}

bool MutationBase::operator<(const MutationBase& rhs) const {
  return value < rhs.value;
}

bool MutationBase::operator==(const char& rhs) const { return ToChar() == rhs; }

bool MutationBase::operator!=(const char& rhs) const { return ToChar() != rhs; }

bool MutationBase::operator<(const char& rhs) const { return ToChar() < rhs; }

inline bool operator==(const char& lhs, const MutationBase& rhs) {
  return lhs == rhs.ToChar();
}

inline bool operator!=(const char& lhs, const MutationBase& rhs) {
  return lhs != rhs.ToChar();
}

inline bool operator<(const char& lhs, const MutationBase& rhs) {
  return lhs < rhs.ToChar();
}
