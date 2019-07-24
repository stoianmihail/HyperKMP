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
unsigned match[MAX_COUNT];
#if 0
std::vector<std::pair<unsigned, unsigned>> store;
unsigned **jumps;
uint64_t *configurations;
unsigned shifted;
#endif

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

static uint64_t log2(uint64_t x)
// log2 when x power of 2
{
    uint64_t lg = 1;
    while ((1ULL << lg) < x)
        ++lg;
    return lg;
}

#if 0
void preprocessing() {
  // assert all letters are from the lower english alphabet

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

    cerr << m << endl;
  unsigned maximum = 0;
  jumps = new unsigned*[m];
  configurations = new uint64_t[m];
  for (unsigned index = 0; index < m; ++index) {
  	unsigned q = index;
    uint64_t configuration;

    // store all jumps, inclusive the first comparison
    // with the initial stateunsigned maximum = 0;

    // for strings like 'TTTTTTTTTTTTT', this takes quite long
    // TODO: replace with topological sort on DAG
#if 0
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
#else
    // Look at pi[q - 1]
    if (q == 0)
        continue;
    q = pi[q - 1];
    if (q == 0)
        continue;

    configuration = configurations[q];

    // Connect to the node q
    if (configuration == 0) {
        //cerr << "maaa kann das sein " << ++shifted << endl;
        jumps[index] = new unsigned[1];
        setBit(configurations[index], encode(pattern[q]));
        jumps[index][0] = q;
        //cerr << "am bagat " << q << " si " << pattern[q] << endl;
        continue;
    }

    if (getBit(configuration, encode(pattern[q]))) {
        // don't change anything in jumps;
        // pointer to it?
        //cerr << "viele" << endl;
        configurations[index] = configuration;
#if 1
        // Create only when configuration has more than one bit
        //if (configuration & (configuration - 1)) {
            // Recover the most fresh position!
            unsigned size = __builtin_popcountll(configuration);
            jumps[index] = new unsigned[size];
            uint64_t save = configuration;
            unsigned ptr = 0;
            unsigned toMerge = encode(pattern[q]);

            while (save) {
                unsigned lg = log2(save & -save);
                if (toMerge != lg)
                    jumps[index][ptr] = jumps[q][ptr];
                else
                    jumps[index][ptr] = q;
                ++ptr;
                save &= save - 1;
            }
        //}
#else
    jumps[index] = jumps[q];
#endif
    } else {
        // add only one
        //cerr << "kann das sein mit " << index << " " << q << ++shifted << endl;

        //cerr << configuration << endl;
        unsigned size = __builtin_popcountll(configuration);
        //cerr << configuration << endl;

        unsigned toMerge = encode(pattern[q]);
        //cerr << "toMerge " << toMerge << endl;

        jumps[index] = new unsigned[size + 1];
        uint64_t save = configuration, ptr = 0;
        bool inserted = false;

        //cerr << "reduce " << save << endl;

        while (save) {
            unsigned lg = log2(save & -save);
            //cerr << "between " << lg << endl;
            if (lg < toMerge) {
                //cerr << "1 " << ptr << endl;
                jumps[index][ptr] = jumps[q][ptr];
                ++ptr;
            } else if (!inserted && lg > toMerge) {
                //cerr << "2 " << ptr << endl;
                jumps[index][ptr++] = q;
                inserted = true;

                // And add the one with which you compared
                jumps[index][ptr] = jumps[q][ptr - 1];
                ptr++;
            } else {
                //cerr << "3 " << ptr << endl;
                jumps[index][ptr] = jumps[q][ptr - 1];
                ++ptr;
            }
            save &= save - 1;
        }

        setBit(configuration, encode(pattern[q]));
        configurations[index] = configuration;

#if 0
        binary(configuration);
        for (ptr = 0; ptr < size + 1; ptr++) {
            cerr << jumps[index][ptr] << " ";
        }
        cerr << "end merge" << endl;
#endif
    }
  }
#endif

#if 0
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

  //cerr << "Realyy" << currChar << " " << encode(currChar) << endl;

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
#if 1
//    unsigned save = q;
    if ((q > 0) && (str[i] != pattern[q]))
        obtainNextState(q, str[i]);
#else
    unsigned take = q;

    q = save;
    while ((q > 0) && (str[i] != pattern[q]))
        q = pi[q - 1];
#endif

#if 0
    if (q != take) {
        //cerr << "Bai ce facem " << i << " " << q << " vs " << take << "cu char " << str[i] << " si pa " << pattern[q] << endl;

    }
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
#endif

unsigned pi[MAX_N + 2];
unsigned psi[MAX_N + 2];
//unsigned lastSeen[MAX_N + 2];

void normalPi() {
  unsigned k = 0;
	// Compute pi and psi
//pi = new unsigned[m];
	//psi = new unsigned[m + 1]; // special for the last char
	for (unsigned q = 1; q < m; ++q) {
	  //cerr << "wir sind jetzt bei " << q << endl;
    while ((k > 0) && (pattern[q] != pattern[k])) {
      k = pi[k - 1];
    }
    if (pattern[q] == pattern[k]) {
      k++;
    }
    pi[q] = k;
	}
}

#if 1
void compressPi() {
	unsigned k = 0;
	// Compute pi and psi
//pi = new unsigned[m];
	//psi = new unsigned[m + 1]; // special for the last char
	for (unsigned q = 1; q < m; ++q) {
	  //cerr << "wir sind jetzt bei " << q << endl;
    while ((k > 0) && (pattern[q] != pattern[k])) {
      k = pi[k - 1];
    }
    if (pattern[q] == pattern[k]) {
      k++;
    }
    pi[q] = k;

		// Compres possbile path
		// If no jump, connect directly to 0
		//cerr << q << " trigger " << (q + 1) << " ";
		if (pi[q] == 0) {
        //cerr << " = 0" << endl;
			psi[q + 1] = 0;
		// Try to hang the edge (q, pi[q]) at the root of pi[q]
		} else if ((psi[pi[q]] != 0) && (pattern[pi[q]] == pattern[psi[pi[q]]])) {
			psi[q + 1] = psi[pi[q]];
			//lastSeen[q + 1] = pi[q];
      //cerr << "second if -> " << endl;
      //cerr << psi[q + 1] << endl;
		} else {
			// Put the edge where it was
			psi[q + 1] = pi[q];
			//lastSeen[q + 1] = pi[q];
			//cerr << "third if -> " << endl;
			//cerr << psi[q + 1] << endl;
		}
  }
  //delete[] pi;
}
#endif

void hyperKmp() {
	if (m > n) {
		cerr << "pattern biger than string!" << endl;
		return;
	}
	compressPi();
  // cerr << "go voer pi" << endl;

	unsigned q = 0;
	unsigned maximum = 0;
	for (unsigned index = 0; index < n; ++index) {
		// Search now on psi
  //unsigned loopsCtr = 0;

#if 1
		//unsigned bef = q;
		int last = -1;
		//cerr << "start with " << q << endl;
		while ((q > 0) && (str[index] != pattern[q])) {
			last = pi[q - 1];
			//cerr << "set last = " << last << endl;
			q = psi[q];
      //loopsCtr++;
		}
		//cerr << "end with " << q << " and " << last << endl;

    if (last != -1)
      q = last;
      //q = (last == -1) ? q : last;
		//unsigned after = q;
#else
  //q = bef;
  unsigned ctr = 0;
  while ((q > 0) && (str[index] != pattern[q])) {
    q = pi[q - 1];
    loopsCtr++;
  }
#if 0
  if (q != after) {
    cerr << "assert : at " << index << " " << q << " vs " << after << endl;
    assert(0);
  }
#endif
#endif
  //if (loopsCtr > maximum)
    //maximum = loopsCtr;

		if (str[index] == pattern[q]) {
			++q;
		}
		if (q == m) {
			if ((++matchCount) <= MAX_COUNT) {
        match[matchCount - 1] = index - m + 1;
      }
      // goes a step back, because the condition in while will be already set to false
			//q = pi[q - 1];
		}
	}
	cerr << maximum << endl;
}

int main(int argc, char** argv) {
    ifstream in;
#if 0
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << "file " << endl;
    return 0;
  }
  in.open(argv[1]);
#else
  in.open(0 ? "strmatch.in" : "test/grader_test30.in");
#endif
  in >> pattern >> str;
  m = strlen(pattern);
  n = strlen(str);

  if (m > n)
        goto print;

  //cerr << pattern << endl << str << endl;

#if 0
  compressPi();

  for (unsigned index = 0; index < m; ++index) {
      cerr << index << " -> " << psi[index] << endl;
  }
  cerr << endl;
#endif
  hyperKmp();

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
