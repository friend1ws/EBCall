#include <string>
#include <map>
#include <fstream>
#include <ostream>
#include <algorithm>
#include <string.h>
#include "InStream.h"
using std::string;
using std::map;
using std::ifstream;
using std::ofstream;
using std::transform;
using std::cerr;
using std::cout;
using std::endl;

static const char m_delim = '\t';

struct LineData {
  LineData():
  chr(""),
  pos(-1),
  ref(""),
  depth(-1),
  bases(""),
  quals("")
  {}
  string chr;
  int pos;
  string ref;
  int depth;
  string bases;
  string quals;
};

struct Str2lowc {
  map<string, char> bases;
  Str2lowc() {
    bases["A"] = 'a';
    bases["C"] = 'c';
    bases["G"] = 'g';
    bases["T"] = 't';
    bases["N"] = 'n';
  }
};

struct Str2upc {
  map<string, char> bases;
  Str2upc() {
    bases["A"] = 'A';
    bases["C"] = 'C';
    bases["G"] = 'G';
    bases["T"] = 'T';
    bases["N"] = 'N';
  }
};

struct Base2idx {
  map<char, int> bases;
  Base2idx() {
    bases['A'] = 0;
    bases['C'] = 1;
    bases['G'] = 2;
    bases['T'] = 3;
    bases['N'] = 4;
    bases['a'] = 5;
    bases['c'] = 6;
    bases['g'] = 7;
    bases['t'] = 8;
    bases['n'] = 9;
  }
};

struct ToUpRef {
  map<string, string> bases;
  ToUpRef() {
    bases["A"] = "A";
    bases["C"] = "C";
    bases["G"] = "G";
    bases["T"] = "T";
    bases["N"] = "N";
    bases["a"] = "A";
    bases["c"] = "C";
    bases["g"] = "G";
    bases["t"] = "T";
    bases["n"] = "N";
  }
};

struct ToUpper{
   char operator()(char c) { return toupper(c); }
};

inline void checkInFileFail(ifstream& fin, const string& fileName)
{
  if(fin.fail()){
    cerr << " cannot open " << fileName << endl;
    exit(1);
  }
}

inline void checkOutFileFail(ofstream& fon, const string& fileName)
{
  if(fon.fail()){
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

int main(int argc, char** argv) {

  if (argc <= 7) {
    cerr << "wrong number of arguments. argc : " << argc << endl;
    return 1;
  }
  checkNumString(argv[1],1);
  int thres = atoi(argv[1]);
  if ( thres < 0 || 40 < thres){
    cerr << "input qualith threshould is not appropriate : " << thres << endl;
    return 1;
  }
  const int minQual = 32 + thres;

  checkNumString(argv[2],2);
  int minDepth = atoi(argv[2]);

  const string inFile = (argv[3]);
  ifstream fin(inFile.c_str());
  checkInFileFail(fin, inFile);
  InStream inStream(fin, inFile);

  const string baseOutFile = (argv[4]);
  ofstream baseofs(baseOutFile.c_str());
  checkOutFileFail(baseofs, baseOutFile);

  const string insOutFile = (argv[5]);
  ofstream insofs(insOutFile.c_str());
  checkOutFileFail(insofs, insOutFile);

  const string delOutFile = (argv[6]);
  ofstream delofs(delOutFile.c_str());
  checkOutFileFail(delofs, delOutFile);

  const string depthOutFile = (argv[7]);
  ofstream depthofs(depthOutFile.c_str());
  checkOutFileFail(depthofs, depthOutFile);

  Base2idx base2idx;
  Str2upc str2upc;
  Str2lowc str2lowc;
  ToUpRef toupref;

  map<string,int> mapInsStrandP;
  map<string,int> mapInsStrandM;
  map<string,int> mapDelStrandP;
  map<string,int> mapDelStrandM;

  baseofs << "chr\tpos\tref\tdepth\tA\ta\tC\tc\tG\tg\tT\tt" << endl;

  while (!inStream.eof()) {
    int iBaseArray[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int iDepthArray[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    struct LineData line;
    line.chr = inStream.nextValue();
    if (inStream.eof()) break;
    line.pos = atoi(inStream.nextValueWithCheck().c_str());
    line.ref = inStream.nextValueWithCheck();
    line.depth = atoi(inStream.nextValueWithCheck().c_str());
    line.bases = inStream.nextValueWithCheck();
    line.quals = inStream.nextValueWithCheck();

    if(line.depth <  minDepth) continue;

    string upRefStr = toupref.bases[line.ref]; 
    char upRef = str2upc.bases[upRefStr];
    char lowRef = str2lowc.bases[upRefStr];
    string readBases = line.bases;
    string baseQuals = line.quals;
    int maxSize = readBases.size();
    int j = 0; // j is  base qualities pointer
    for (int i = 0;i < maxSize;i++) {
      // i is read bases pointer
      // skip the start of a read segment and the following character
      if (readBases[i] == '^') {
        i++;
        continue;

      // skip the end of a read segment
      } else if(readBases[i] == '$'){
          continue;

      // if base is insertion
      } else if(readBases[i] == '+'){
        i++; // skip a symbole '+'
        string numStr;
        for(;'0' <= readBases[i] && readBases[i] <= '9';i++){
          numStr += readBases[i];
        }
        int len = atoi(numStr.c_str());
        string key = readBases.substr(i, len);
        i += (len - 1);
        bool lowerFlg = false;
        bool upperFlg = false;

        if(isupper(key[0]) != 0) {
          upperFlg = true;
        } else if(islower(key[0]) != 0) {
          transform(key.begin(), key.end(), key.begin(), ToUpper());
          lowerFlg = true;
        } else {
          cerr << "something is wrong! insertion." << line.chr << " " << line.pos << endl;
          exit(1);
        }

        if (mapInsStrandP.find(key) == mapInsStrandP.end() || mapInsStrandM.find(key) == mapInsStrandM.end()) {
          if (upperFlg){
            mapInsStrandP[key] = 1; mapInsStrandM[key] = 0;
          } else if(lowerFlg){
            mapInsStrandP[key] = 0; mapInsStrandM[key] = 1;
          }
        } else {
          if (upperFlg){
            mapInsStrandP[key] += 1;
          } else if(lowerFlg){
            mapInsStrandM[key] += 1;
          }
        }
      
      // if base is deletion
      } else if (readBases[i] == '-') {
        i++; // skip a symbole '-'
        string numStr;
        for (;'0' <= readBases[i] && readBases[i] <= '9'; i++) {
          numStr += readBases[i];
        }

        int len = atoi(numStr.c_str());
        string key = readBases.substr(i, len);
        i += (len-1);

        bool lowerFlg = false;
        bool upperFlg = false;
        if (isupper(key[0]) != 0) {
          upperFlg = true;
        } else if (islower(key[0]) != 0) {
          transform(key.begin(), key.end(), key.begin(), ToUpper());
          lowerFlg = true;
        } else {
          cerr << "something is wrong! deletion." << line.chr << " " << line.pos << endl;
        }

        if (mapDelStrandP.find(key) == mapDelStrandP.end() || mapDelStrandM.find(key) == mapDelStrandM.end()) {
          if      (upperFlg) { mapDelStrandP[key] = 1; mapDelStrandM[key] = 0; }
          else if (lowerFlg) { mapDelStrandP[key] = 0; mapDelStrandM[key] = 1; }
        } else {
          if      (upperFlg) { mapDelStrandP[key] += 1; }
          else if (lowerFlg) { mapDelStrandM[key] += 1; }
        }

      // if SNV 
      } else {
        if (readBases[i] == '.') {
          int idx = base2idx.bases[upRef];
          if (baseQuals[j] > minQual) {
            iBaseArray[idx] += 1;
          }
          iDepthArray[idx] += 1;
        } else if (readBases[i] == ',') {
          int idx = base2idx.bases[lowRef];
          if (baseQuals[j] > minQual) {
            iBaseArray[idx] += 1;
          }
          iDepthArray[idx] += 1;
        } else {
          if (base2idx.bases.find(readBases[i]) != base2idx.bases.end()){
            int idx = base2idx.bases[readBases[i]];
            if (baseQuals[j] > minQual) {
              iBaseArray[idx] += 1;
            }
            iDepthArray[idx] += 1;
          }
        }
        j++; // the pointer moves next base qualities
      }
    }

    int baseQualLen = baseQuals.length();
    if (j != baseQualLen) {
      cerr << "something is wrong!!! " << j << " : " << baseQuals.length() << endl;
    }

    if (!mapInsStrandP.empty()) {
      map<string, int>::iterator it = mapInsStrandP.begin();
      for(;it != mapInsStrandP.end(); it++) {
        insofs << line.chr << m_delim << line.pos << m_delim << upRefStr << m_delim << line.depth << m_delim
        << it->first << m_delim << it->second << m_delim << mapInsStrandM[it->first] << endl;
      }
      mapInsStrandP.clear();
      mapInsStrandM.clear();
    }

    if (!mapDelStrandP.empty()) {
      map<string, int>::iterator it = mapDelStrandP.begin();
      for(;it != mapDelStrandP.end(); it++) {
        delofs << line.chr << m_delim << line.pos << m_delim << upRefStr << m_delim << line.depth << m_delim
        << it->first << m_delim << it->second << m_delim << mapDelStrandM[it->first] << endl;
      }
      mapDelStrandP.clear();
      mapDelStrandM.clear();
    }
    
    baseofs << line.chr << m_delim << line.pos << m_delim << upRefStr << m_delim << line.depth << m_delim
    << iBaseArray[0] << m_delim
    << iBaseArray[5] << m_delim
    << iBaseArray[1] << m_delim
    << iBaseArray[6] << m_delim
    << iBaseArray[2] << m_delim
    << iBaseArray[7] << m_delim
    << iBaseArray[3] << m_delim
    << iBaseArray[8] << endl;

    depthofs << line.chr << m_delim << line.pos << m_delim << upRefStr << m_delim << line.depth << m_delim
    << (iDepthArray[0] + iDepthArray[1] + iDepthArray[2] + iDepthArray[3] + iDepthArray[4]) << m_delim
    << (iDepthArray[5] + iDepthArray[6] + iDepthArray[7] + iDepthArray[8] + iDepthArray[9]) << endl;
  }
  return 0;
}

