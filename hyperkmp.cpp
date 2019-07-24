#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <cctype>
#include <vector>
#include <algorithm>

#define MAX_N 2000000
#define MAX_COUNT 1000

using namespace std;

unsigned n, m, matchCount;
char str[MAX_N + 5], pattern[MAX_N + 5];
unsigned pi[MAX_N + 1];
std::vector<unsigned> match;
std::vector<std::pair<unsigned, unsigned>> store;
unsigned **jumps;
uint64_t *configurations;

static void binary(uint64_t x) {
    cerr << "bits : ";
    for (unsigned bit = 0; bit < 64; ++bit)
        if ((x >> bit) & 1ULL)
            cerr << bit << " ";
    cerr << endl;
}


static inline int getBit(uint64_t x, unsigned index) {
  return (x >> index) & 1;
}

static void setBit(uint64_t& x, unsigned index) {
  x |= (1ULL << index);
}

#define SIGMA 26
static int encode(char c)
// keep only lower letters
{
    if (islower(c))
        return c - 'a';
    if (isupper(c))
        return c - 'A' + SIGMA;
    return c - '0' + 2 * SIGMA;
}

void preprocessing() {
  // assert all letters are from the lower english alphabet
  for (unsigned index = 0; index < m; ++index)
    assert(isalnum(pattern[index]));

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
#if 1
    cerr << m << endl;
  unsigned maximum = 0;
  jumps = new unsigned*[m];
  configurations = new uint64_t[m];
  for (unsigned index = 0; index < m; ++index) {
  	unsigned q = index;
    uint64_t configuration = 0;

    // store all jumps, inclusive the first comparison
    // with the initial stateunsigned maximum = 0;
    
    // for strings like 'TTTTTTTTTTTTT', this takes quite long
    // TODO: replace with topological sort on DAG
    store.clear();
    while (q > 0) {
        q = pi[q - 1];
        if (q) {
            if (!getBit(configuration, encode(pattern[q]))) {
                setBit(configuration, encode(pattern[q]));
                store.push_back(make_pair(encode(pattern[q]), q));
            }
        }
    }
    
    
    //cerr << "for index = " << index << " conf = ";
    //binary(configuration);
    
    //cerr << maximum << endl;
    
    if (!store.empty()) {
        maximum = std::max(maximum, static_cast<unsigned>(store.size()));
        //cerr << "stored : " << store.size() << endl;
    // sort in ascending order the pushed chars
      sort(store.begin(), store.end());

#if 0
      if (store.size() == 2) {
      
          cerr << "store " << store.begin()->first << " " << store.begin()->second << endl;
      cerr << (--store.end())->first << " " << (--store.end())->second << " end strre" << endl;
      }
#endif
      // get only the indexes where they are stored
      // with other words, exactly the position which
      // is returned by the while-loop in the initial KMP-//  algoritm
      jumps[index] = new unsigned[store.size()];
      //cerr << "bei index = " << index << endl;
      //cerr << "jumps : ";
      for (unsigned pos = 0; pos < store.size(); ++pos) {
        jumps[index][pos] = store[pos].second;
      }
    }
    configurations[index] = configuration;
    //binary(configurations[index]);
  }
  cerr << "maximum stored " << maximum << endl;
#endif    
}

void obtainNextState(unsigned& state, char currChar)
// get the next state, by analysing the stored jumps
{
  if (!isalnum(currChar)) {
    state = 0;
    return;
}
  // Check if inside
  
  //cerr << "come with " << currChar << " mit " << encode(currChar) << endl;
  //binary(configurations[state]);
  
  if (!getBit(configurations[state], encode(currChar))) {
    state = 0;
    return;
    }
  
  // only one bit
  if ((configurations[state] & (configurations[state] - 1)) == 0) {
    state = jumps[state][0];
    return;
    }
  
  // Get its index;
  unsigned code = encode(currChar);
  // count how many bits are set before currChar
  uint64_t MASK = (1ULL << code) - 1;
  uint64_t remainedBits = configurations[state] & MASK;
  unsigned pos = __builtin_popcountll(remainedBits);

#if 0
  binary(MASK);
  binary(remainedBits);
  cerr << "pos = " << pos << " mit " << currChar;
  binary(configurations[state]);
  
  cerr << "get rest " << jumps[state][pos] << endl;
#endif
  // And get the final result
  state = jumps[state][pos];
}

/** Cauta "p" in "s" si retine potrivirile in "match". **/
void kmp() {
  if (m > n) {
    cerr << "pattern bigger than string" << endl;
    exit(0); 
  }
  preprocessing();
    cerr << "end preprocessing" << endl;
#if 0
    cerr << "pi vector " << endl;
    for (unsigned index = 0; index < m; ++index)
      cerr << pi[index] << " ";
    cerr << endl;

    cerr << "lens = pattern : " << m << " str : " << n << endl;
    for (unsigned index = 0; index < m; index++)
      cerr << index << " -> " << pattern[index] << endl;
    //cerr << "pattern : " << endl << pattern << endl;
#endif
  unsigned q = 0;
  for (unsigned i = 0; i < n; i++) {
#if 0
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
#endif  
    //cerr << "before " << q << endl;
    // it shoud remain as it is now!
#if 0
    if ((q > 0) && (str[i] != pattern[q]))
        obtainNextState(q, str[i]);
#else
    while ((q > 0) && (str[i] != pattern[q]))
        q = pi[q - 1];
#endif

    if (str[i] == pattern[q]) {
      q++;
    }
    if (q == m) {
      if ((++matchCount) <= MAX_COUNT) {
        match.push_back(i - m + 1);
      }
      // goes a step back, because the condition in while will be already set to false
        q = pi[q - 1];
        
    }
  }
}

int main(int argc, char** argv) {
    ifstream in;
#if 1
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << "file " << endl;
    return 0;
  }
  in.open(argv[1]);
#else
  in.open("strmatch.in");
#endif
  in >> pattern >> str;
  m = strlen(pattern);
  n = strlen(str);

  if (m > n)
        goto print;
  
  //cerr << pattern << endl << str << endl;

  kmp();

  //cerr << encode('9') << " " << encode('b') << " " << encode('C');
  
  print : {
      ofstream out;
#if 0
    out.open("strmatch.out");
#else
    out.open("response.out");
#endif
    out << matchCount << "\n";
    matchCount = matchCount > MAX_COUNT ? MAX_COUNT : matchCount;
    for (unsigned index = 0; index < matchCount; ++index)
        out << match[index] << " ";
    out << "\n";
    out.close();
  }
  return 0;
}
