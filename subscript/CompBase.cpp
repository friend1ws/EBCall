#include <string>
#include <map>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <string.h>
#include <limits>
#include "InStream.h"
using std::string;
using std::map;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::numeric_limits;

static const char m_delim = '\t';
static const char m_com = ',';

struct LineData {
  LineData():
  chr(""),
  pos(-1),
  ref(""),
  depth(-1),
  aP(-1),
  aM(-1),
  cP(-1),
  cM(-1),
  gP(-1),
  gM(-1),
  tP(-1),
  tM(-1)
  {}
  string chr;
  int pos;
  string ref;
  int depth;
  int aP;
  int aM;
  int cP;
  int cM;
  int gP;
  int gM;
  int tP;
  int tM;
};

struct Base2idx {
  map<string, int> bases;
  Base2idx() {
    bases["A"] = 0;
    bases["C"] = 1;
    bases["G"] = 2;
    bases["T"] = 3;
  }
};

struct Idx2base {
map<int, string> bases;
 Idx2base() {
    bases[0] = "A";
    bases[1] = "C";
    bases[2] = "G";
    bases[3] = "T";
  }
};

inline int getMisIndex(const int iRef, const int iA, const int iC, const int iG, const int iT) {
  int misIdx = 0;
  int iMax = 0;
  if (iMax < iA && iRef != 0) { iMax = iA; misIdx = 0; }
  if (iMax < iC && iRef != 1) { iMax = iC; misIdx = 1; }
  if (iMax < iG && iRef != 2) { iMax = iG; misIdx = 2; }
  if (iMax < iT && iRef != 3) { iMax = iT; misIdx = 3; }
  return misIdx;
}

inline int getIdxValue(const int idx, const int iA, const int iC, const int iG, const int iT) {
  switch (idx) {
	  case 0: return iA;
	  case 1: return iC;
	  case 2: return iG;
	  case 3: return iT;
  }
  return 0;
}

inline double getMisRate(const int misNum, const int depth)
{
  return (((double)misNum) / depth);
}

inline void checkInFileFail(ifstream& fin, const string& fileName)
{
  if(fin.fail()){
    cerr << " cannot open " << fileName << endl;
    exit(1);
  }
}

inline void checkNumString(char *str, int num)
{
  int strSize = strlen(str);
  for(int i = 0; i < strSize; i++){
    if((str[i] < '0') || ('9' < str[i])) {
      cerr << " number format error : " << str << endl;
      cerr << " please, check argv[" << num << "]" << endl;
      exit(1);
    }
  }
}

inline int getDepth(LineData& line)
{
  return (line.aP + line.aM + line.cP + line.cM + line.gP + line.gM + line.tP + line.tM);
}

inline int getMisNum(struct LineData& line, int & misIdx)
{
    int nA = (line.aP + line.aM);
    int nC = (line.cP + line.cM);
    int nG = (line.gP + line.gM);
    int nT = (line.tP + line.tM);
    int nMisNum = getIdxValue(misIdx, nA, nC, nG, nT);
    return nMisNum;
}

inline void outputData(struct LineData& tLine, int misIdx, struct LineData& nLine, Idx2base& m_idx2base)
{
    int tP = (tLine.aP + tLine.cP + tLine.gP + tLine.tP);
    int tMisP = getIdxValue(misIdx, tLine.aP, tLine.cP, tLine.gP, tLine.tP);
    int tM = (tLine.aM + tLine.cM + tLine.gM + tLine.tM);
    int tMisM = getIdxValue(misIdx, tLine.aM, tLine.cM, tLine.gM, tLine.tM);
    int nP = (nLine.aP + nLine.cP + nLine.gP + nLine.tP);
    int nMisP = getIdxValue(misIdx, nLine.aP, nLine.cP, nLine.gP, nLine.tP);
    int nM = (nLine.aM + nLine.cM + nLine.gM + nLine.tM);
    int nMisM = getIdxValue(misIdx, nLine.aM, nLine.cM, nLine.gM, nLine.tM);
    cout << tLine.chr << m_delim << tLine.pos << m_delim << tLine.ref << m_delim
    << m_idx2base.bases[misIdx] << m_delim
    << tP << m_com << tMisP << m_com << tM << m_com << tMisM << m_delim
    << nP << m_com << nMisP << m_com << nM << m_com << nMisM << endl;
}


inline struct LineData getLineData(InStream& inStream)
{
	struct LineData line;
    line.chr = inStream.nextValue();
    if (inStream.eof()) return line;
    line.pos = atoi((inStream.nextValueWithCheck()).c_str());
    line.ref = inStream.nextValueWithCheck();
    line.depth = atoi((inStream.nextValueWithCheck()).c_str());
    line.aP = atoi((inStream.nextValueWithCheck()).c_str());
    line.aM = atoi((inStream.nextValueWithCheck()).c_str());
    line.cP = atoi((inStream.nextValueWithCheck()).c_str());
    line.cM = atoi((inStream.nextValueWithCheck()).c_str());
    line.gP = atoi((inStream.nextValueWithCheck()).c_str());
    line.gM = atoi((inStream.nextValueWithCheck()).c_str());
    line.tP = atoi((inStream.nextValueWithCheck()).c_str());
    line.tM = atoi((inStream.nextValueWithCheck()).c_str());
    return line;
}


int main(int argc, char** argv) {

  if (argc <= 7) {
    cerr << "wrong number of arguments. argc : " << argc << endl;
    return 1;
  }
  const string inTumorFile = (argv[1]);
  ifstream inTumorFileStream(inTumorFile.c_str());
  checkInFileFail(inTumorFileStream, inTumorFile);
  InStream isTumor(inTumorFileStream, inTumorFile);

  const string inNormalFile = (argv[2]);
  ifstream inNormalFileStream(inNormalFile.c_str());
  checkInFileFail(inNormalFileStream, inNormalFile);
  InStream isNormal(inNormalFileStream, inNormalFile);

  checkNumString(argv[3],3);
  const int minTumorDepth = atoi(argv[3]);

  checkNumString(argv[4],4);
  const int minNormalDepth = atoi(argv[4]);

  checkNumString(argv[5],5);
  const int minTumorVariantRead = atoi(argv[5]);

  const double minTumorAlleleFreq = atof(argv[6]);
  const double maxNormalAlleleFreq = atof(argv[7]);

  Base2idx m_base2idx;
  Idx2base m_idx2base;

  struct LineData tLine;
  struct LineData nLine;

  // skip header line
  isTumor.nextLine();
  isNormal.nextLine();

  while(!isTumor.eof() || !isNormal.eof()){

    while(!isTumor.eof()){
      tLine = getLineData(isTumor);

      if(tLine.pos > nLine.pos || tLine.chr != nLine.chr || isTumor.eof()) {
        break;
      } else if (tLine.pos == nLine.pos && tLine.chr == nLine.chr){
        int tDepth = getDepth(tLine);
        if (tDepth < minTumorDepth) continue;

        int nDepth = getDepth(nLine);
        if (nDepth < minNormalDepth) continue;

        char iRef = m_base2idx.bases[tLine.ref];
        int tA = (tLine.aP + tLine.aM);
        int tC = (tLine.cP + tLine.cM);
        int tG = (tLine.gP + tLine.gM);
        int tT = (tLine.tP + tLine.tM);
        int misIdx = getMisIndex(iRef, tA, tC, tG, tT);
        int tMisNum = getIdxValue(misIdx, tA, tC, tG, tT);
        double tMisRate = getMisRate(tMisNum, tDepth);
        if (tMisRate < minTumorAlleleFreq || tMisNum < minTumorVariantRead) continue;

        int nMisNum = getMisNum(nLine, misIdx);
        double nMisRate = getMisRate(nMisNum, nDepth);
        if (nMisRate >= maxNormalAlleleFreq) continue;

        outputData(tLine, misIdx, nLine, m_idx2base);
      }
    }

    while (!isNormal.eof()) {
      nLine = getLineData(isNormal);

      if(tLine.pos < nLine.pos || tLine.chr != nLine.chr || isNormal.eof() ) {
        break;
      } else if (tLine.pos == nLine.pos && tLine.chr == nLine.chr){
        int tDepth = getDepth(tLine);
        if (tDepth < minTumorDepth) continue;

        int nDepth = getDepth(nLine);
        if (nDepth < minNormalDepth) continue;

        char iRef = m_base2idx.bases[tLine.ref];
        int tA = (tLine.aP + tLine.aM);
        int tC = (tLine.cP + tLine.cM);
        int tG = (tLine.gP + tLine.gM);
        int tT = (tLine.tP + tLine.tM);
        int misIdx = getMisIndex(iRef, tA, tC, tG, tT);
        int tMisNum = getIdxValue(misIdx, tA, tC, tG, tT);
        double tMisRate = getMisRate(tMisNum, tDepth);
        if (tMisRate < minTumorAlleleFreq || tMisNum < minTumorVariantRead) continue;

        int nMisNum = getMisNum(nLine, misIdx);
        double nMisRate = getMisRate(nMisNum, nDepth);
        if (nMisRate >= maxNormalAlleleFreq) continue;

        outputData(tLine, misIdx, nLine, m_idx2base);
      }
    }
  }
  return 0;
}


