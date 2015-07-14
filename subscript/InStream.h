#include <iostream>
#include <string>
using std::string;
using std::istream;


class InStream {

private:
  static const char m_delim = '\t';
  static const char m_lineSep = '\n';
  istream& m_stream;
  const string& m_filename;
  int m_num;
  
public:
  InStream(istream& targetStream, const string& filename);
  const string nextValue();
  const string nextValueWithCheck();
  const string nextLine();
  bool eof();
};

