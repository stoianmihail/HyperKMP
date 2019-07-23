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

unsigned n, m, count;
char str[MAX_N + 5], pattern[MAX_N + 5];
unsigned pi[MAX_N + 1];
unsigned match[MAX_COUNT + 1];
std::vector<std::pair<char, unsigned>> store;
unsigned **jumps;
unsigned *configurations;

static inline int getBit(unsigned x, unsigned index) {
  return (x >> index) & 1;
}

static void setBit(unsigned& x, unsigned index) {
  x |= (1 << index);
}

static int encode(char c)
// keep only lower letters
{
  return c - 'a';
}

void preprocessing() {
  // assert all letters are from the lower english alphabet
  for (unsigned index = 0; index < m; ++index)
    assert(islower(pattern[index]));

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

  jumps = new unsigned*[m];
  configurations = new unsigned[m];
  for (unsigned index = 0; index < m; ++index) {
  	unsigned q = index, configuration = 0;

    // store all jumps, inclusive the first comparison
    // with the initial state
    while (q > 0) {
      if (!getBit(configuration, encode(pattern[q]))) {
        setBit(configuration, encode(pattern[q]));
        store.push_back(make_pair(pattern[q], q));
      }
      q = pi[q - 1];

    }
    if (!store.empty()) {
      // sort in ascending order the pushed chars
      sort(store.begin(), store.end());

      // get only the indexes where they are stored
      // with other words, exactly the position which
      // is returned by the while-loop in the initial KMP-//  algoritm
      jumps[index] = new unsigned[store.size()];
      for (unsigned pos = 0; pos < store.size(); ++pos) {
        jumps[index][pos] = store[pos].second;
      }
      configurations[index] = configuration;
    }
  }
}

void obtainNextState(unsigned& state, char currChar)
// get the next state, by analysing the stored jumps
{
  if (!islower(currChar))
    goto resetState;

  // Check if inside
  if (!getBit(configurations[state], encode(currChar)))
    goto resetState;

  // Get its index;
  unsigned code = encode(currChar);
  // count how many bits are set before currChar
  unsigned MASK = (1 << code) - 1;
  unsigned remainedBits = configurations[state] & MASK;
  unsigned pos = __builtin_popcount(remainedBits);

  // And get the final result
  state = jumps[state][pos];
  return;

  resetState : {
    state = 0;
    return;
  }

}

/** Cauta "p" in "s" si retine potrivirile in "match". **/
void kmp() {
  if (m > n) {
    cerr << "pattern bigger than string" << endl;
    exit(0);
  }
  preprocessing();

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
    obtainNextState(q, str[i]);
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
  n = strlen(str);

  //cerr << pattern << endl << str << endl;

  kmp();
#if 1
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
