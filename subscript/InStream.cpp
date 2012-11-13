#include "InStream.h"
#include <fstream>
#include <string>
#include <cstdlib>
using std::string;
using std::istream;
using std::cerr;
using std::endl;


// construct & set the pointer of 'ris' to private variable
InStream::InStream(istream& targetStream, const string& filename) :
m_stream(targetStream),
m_filename(filename)
{}

// get separated value
const string InStream::nextValue() {
  string col;
  while ((m_num = m_stream.get()) != m_delim && m_num != m_lineSep && m_num != EOF) {
    col += m_num;
  }
  return col;
}

const string InStream::nextValueWithCheck() {
  string col;
  bool valFlg = false;
  while ((m_num = m_stream.get()) != m_delim && m_num != m_lineSep && m_num != EOF) {
    col += m_num;
    valFlg = true;
  }
  if (!valFlg) {
	  cerr << m_filename << " is broken." << endl;
	  exit(1);
  }
  return col;
}

// get line value
const string InStream::nextLine() {
  string line;
  getline(m_stream, line);
  return line;
}

// is end of file
bool InStream::eof() {
  return m_stream.eof();
}

