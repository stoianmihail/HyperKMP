#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>

#define MAX_N 2000000
#define MAX_COUNT 1000

using namespace std;

unsigned n, m, count;
char str[MAX_N + 5], pattern[MAX_N + 5];
unsigned pi[MAX_N + 1];
unsigned match[MAX_COUNT + 1];

void preprocessing() {
  unsigned k = 0;
  for (unsigned q = 1; q < m; q++) {
    while ((k > 0) && (pattern[q] != pattern[k])) {
      k = pi[k - 1];
    }
    if (pattern[q] == pattern[k]) {
      k++;
    }
    pi[q] = k;
  }	
}
	
/** Cauta "p" in "s" si retine potrivirile in "match". **/
void kmp() {
  if (m > n) {
    cerr << "pattern bigger than string" << endl;
    exit(0);
  }
  preprocessing();

#if 1
    cerr << "pi vector " << endl;
    for (unsigned index = 0; index < m; ++index)
      cerr << pi[index] << " ";
    cerr << endl;
#endif
  
    cerr << "lens = pattern : " << m << " str : " << n << endl; 
    for (unsigned index = 0; index < m; index++)
      cerr << index << " -> " << pattern[index] << endl;
    //cerr << "pattern : " << endl << pattern << endl;
    
  unsigned q = 0;
  for (unsigned i = 0; i < n; i++) {
    cerr << "analyse " << i << " with " << str[i] << " and state = " << q << endl;
    unsigned inside = 0;
    char tracker = pattern[q];
    while ((q > 0) && (str[i] != pattern[q])) {
      q = pi[q - 1];
#if 1
      ++inside;
      if (q) {
        cerr << pattern[q] << " ";
        if (pattern[q] != tracker) {
          cerr << "jetzt sind wir bei " << i << " mit tracker.." << tracker << "..\n";
          assert(0);
        }
      }
#endif
    }
    if (inside) {
      cerr << "end " << endl;
    }
    if (str[i] == pattern[q]) {
      q++;
    }
    if (q == m) {
      if ((++count) <= MAX_COUNT) {
        match[count - 1] = i - m + 1;
      }
      // goes a step back, because the condition in while will be already set to false
      q = pi[q - 1];
    }
  }
}
	
int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << "file " << endl;
    return 0;
  }
  ifstream in(argv[1]);
  in >> pattern >> str;
  m = strlen(pattern);
  pattern[m] = '#';
  n = strlen(str);

  //cerr << pattern << endl << str << endl;
  
  kmp();
#if 0
  ofstream out("response.out");
  out << count << "\n";
  count = count > MAX_COUNT ? MAX_COUNT : count;
  for (unsigned index = 0; index < count; ++index)
    out << match[index] << " ";
  out << "\n";
  out.close();
#endif
  return 0;
}
