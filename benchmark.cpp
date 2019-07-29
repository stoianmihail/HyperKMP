#include <iostream>
#include <fstream>
#include <cstring>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <chrono>
#include "hyperKMP.hpp"

#include <emmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

#include <boost/iostreams/device/mapped_file.hpp>
//---------------------------------------------------------------------------
#define HYPER_MODE 0
#define NORMAL_MODE 1
#define TEST_MODE 2
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
hyperKMP::hyperKMP() {
}
//---------------------------------------------------------------------------
hyperKMP::hyperKMP(char* pattern, unsigned mode) {
  this->mode = mode;
  this->length = strlen(pattern);
  this->pattern = new char[this->length + 2];
  strcpy(this->pattern, pattern); // the last char should be null
  
  // Transform pattern without underscores
  for (unsigned index = 0; index < this->length; ++index)
      if (this->pattern[index] == '_') 
        this->pattern[index] = ' ';
  cerr << "modified pattern: " << this->pattern << endl;
      
  // For benchmarking
  this->hyperSum = this->hyperMaximum = this->normalSum = this->normalMaximum = 0;
  
  this->valid = new bool[this->length + 1]();
  this->psi = new unsigned[this->length + 1]();
  this->omega = new unsigned[this->length + 1]();
  this->degrees = new unsigned[this->length + 1]();
  this->edges = new unsigned*[this->length + 1];
  compressPi();
  
  // Only for benchmark reasons
  this->pi = new unsigned[this->length + 1]();
  normalPi();
}
//---------------------------------------------------------------------------
void hyperKMP::normalPi() 
// compute the normal pi[] table
{
  unsigned k = 0;
  for (unsigned q = 1; q < this->length; ++q) {
    while ((k > 0) && (pattern[q] != pattern[k]))
      k = pi[k - 1];
    if (pattern[q] == pattern[k])
      k++;
    pi[q] = k;
  }
#if 1
  cerr << "Debug pi[]" << endl;
  for (unsigned index = 0; index < this->length; ++index) {
    cerr << index << " -> " << pi[index] << endl;
  }
  cerr << endl;
#endif
}
//---------------------------------------------------------------------------
void hyperKMP::compressPi()
// compress pi[] and create two new arrays: psi[] and omega[]
// description in README
{  
  // Initialize the first values
  // 0 is the first state. At this step it receives 2 sons
  psi[0] = 0;
  psi[1] = 0;
  degrees[0] = 2;
  // Compute pi and psi
  unsigned k = 0;
  for (unsigned q = 1; q < this->length; ++q) {
    // Compute psi
    while ((k > 0) && (pattern[q] != pattern[k])) {
      k = psi[k];
    }
    if (pattern[q] == pattern[k]) {
      k++;
    }
    // Compress possible path
    // For understandability reasons, keep the formula like this
    // TODO: it could be replaced by: (pattern[q + 1] == pattern[k]) ? psi[k] : k
#if 1
    psi[q + 1] = (!k) ? 0 : ((pattern[q + 1] == pattern[k]) ? psi[k] : k);
#else
    // If no jump, connect directly to 0
    pi[q] = k;
    if (pi[q] == 0) {
      psi[q + 1] = 0;
    } else if (pattern[q + 1] == pattern[pi[q]]) {
      // Try to hang the edge (q, pi[q]) at the root of pi[q]
      psi[q + 1] = psi[pi[q]];
    } else {
      // Put the edge back where it was
      psi[q + 1] = pi[q];
    }
    // build the graph of psi, by counting the degree of each node
#endif
    degrees[psi[q + 1]]++;
  }
  
  // Alloc the graph and reset the degrees to 0
  for (unsigned q = 0; q <= this->length; ++q) {
    if (degrees[q])
      edges[q] = new unsigned[degrees[q]];
    degrees[q] = 0;
  }
  // Compute the graph
  for (unsigned q = 1; q <= this->length; ++q) {
    edges[psi[q]][degrees[psi[q]]++] = q; 
  }

#if 0
  cerr << "Debug" << endl;
  for (unsigned index = 0; index <= this->length; ++index) {
    cerr << "Main node : " << index << " ";
    for (unsigned ptr = 0; ptr < degrees[index]; ++ptr) {
      cerr << edges[index][ptr] << " ";
    }
    cerr << endl;
  }
  cerr << endl;
#endif
  
  // Compute omega[] with dfs
  unordered_map<char, unsigned> indexes;
  vector<unsigned> activeNodes;
  dfs(0, indexes, activeNodes, 0);
  assert(indexes.empty());
  assert(activeNodes.empty());  

#if 1
  cerr << "Debug psi[]" << endl;
  for (unsigned index = 0; index <= this->length; ++index) {
    cerr << index << " with " << pattern[index] << " psi = " << psi[index] << " and omega = " << omega[index] << endl;
  }
  cerr << endl;
#endif

  // delete the graph
  for (unsigned index = 0; index <= this->length; ++index) {
    if (degrees[index])
      delete[] edges[index];
  }
  delete[] degrees;
}
//---------------------------------------------------------------------------
void hyperKMP::dfs(unsigned state, unordered_map<char, unsigned>& indexes, vector<unsigned>& activeNodes, unsigned minimumHalt)
// compute the omega[] in linear time
{
  // Save into the modified recursion stack
  activeNodes.push_back(state);
  valid[state] = true;
    
  // Additional variables for reset
  bool replaced = false;
  unsigned replacer;
  if (state) {
    // Get the last index of the curren state
    auto iter = indexes.find(state);
    
    // Is it the first time we see it?
    if (iter == indexes.end()) {
      // Add its index
      indexes[pattern[state]] = state;
    } else {
      // Replace the last with the new position
      replaced = true;
      replacer = iter->second;
      
      // Unmark the last index and update with the current index
      valid[replacer] = false;
      indexes[pattern[state]] = state;
    }
    
    // Compute the minimum halt
    while (!valid[activeNodes[minimumHalt]])
      ++minimumHalt;
  }
  omega[state] = activeNodes[minimumHalt];
  
  // Explore the neighbours
  for (unsigned index = 0; index < degrees[state]; ++index)
    dfs(edges[state][index], indexes, activeNodes, minimumHalt);
  
  // Reset the stack
  activeNodes.pop_back();
  valid[state] = false;
  if (state) {
    if (!replaced) {
      // Erase the state
      indexes.erase(pattern[state]);
    } else {
      // Mark the last index and update it
      valid[replacer] = true;
      indexes[pattern[state]] = replacer;
    }
  }
}
//---------------------------------------------------------------------------
bool hyperKMP::search(const char* str, unsigned n)
// checks if the pattern can be found in str
{
  if (this->length > n)
    return false;
  
  switch (this->mode) {
    case HYPER_MODE : {
      unsigned q = 0; // current state
      for (unsigned index = 0; index < n; ++index) {
        // Search only on psi now
#if 1
        while ((q > 0) && (str[index] != pattern[q]))
          q = psi[q];
#else
        // lastHalt represents the last state which should be checked (starting from q)
        unsigned lastHalt = omega[q];
        while ((q > 0) && (str[index] != pattern[q]))
          q = psi[q];
#endif
        if (str[index] == pattern[q]) {
          ++q;
        }
        // Found?
        if (q == this->length)
          return true;
      }
      return false;
    }
    case NORMAL_MODE : {
      unsigned q = 0; // current state
      for (unsigned index = 0; index < n; ++index) {
        while ((q > 0) && (str[index] != pattern[q]))
          q = pi[q - 1];
        if (str[index] == pattern[q]) {
          ++q;
        }
        // Found?
        if (q == this->length)
          return true;
      }
      return false;
    }
    case TEST_MODE : {
      unsigned q = 0;
      unsigned hyperLoopsCtr, normalLoopsCtr;
      for (unsigned index = 0; index < n; ++index) {
        hyperLoopsCtr = 0, normalLoopsCtr = 0;
        // lastHalt represents the last state which should be discovered (starting from q)
        unsigned before = q;
        unsigned lastHalt = omega[q];
        while ((q > lastHalt) && (str[index] != pattern[q])) {
          q = psi[q];
          hyperLoopsCtr++;
        }
        unsigned after = q;
    
        q = before;
        //if (q >= this->length / 2 && str[index] != pattern[q] && pi[q - 1] != 0) cerr << '#' << endl;
        while ((q > 0) && (str[index] != pattern[q])) {
          q = pi[q - 1];
          normalLoopsCtr++;
        }
        assert(q == after);
#if 0
        //cerr << q << " vs " << after << endl;
        if (q != after) {
          cerr << "assert : at " << index << " " << q << " vs " << after << endl;
          assert(0);
        }
#endif
        // Compute the average and the maximum loop steps
        hyperSum += hyperLoopsCtr;
        normalSum += normalLoopsCtr;
        hyperMaximum = std::max(hyperMaximum, hyperLoopsCtr);
        normalMaximum = std::max(normalMaximum, normalLoopsCtr);
    
        // Continue after testing
        if (str[index] == pattern[q]) {
          ++q;
        }
        // Found?
        if (q == this->length) {
          return true;
        }
      }
      return false;
    }
  }
  cerr << "Something went wrong!" << endl;
  return false;
}
//---------------------------------------------------------------------------
void hyperKMP::benchmark() {
  cerr << "Show benchmark" << endl;
  cerr << "Hyper: sum = " << hyperSum << ", max = " << hyperMaximum << endl; 
  cerr << "Normal: sum = " << normalSum << ", max = " << normalMaximum << endl;
  cerr << "Win " << (normalSum - hyperSum) << " comparisons saved. Relative saving = " << 100 * ((double)(normalSum - hyperSum) / normalSum) << "%" << endl;
}
//---------------------------------------------------------------------------
static inline const char* findNl(const char* reader, const char* readerLimit)
{
#if 0
   auto bars=_mm_set1_epi8('\n');
   auto limit=readerLimit-16;
   for (;reader<=limit;reader+=16) {
      auto nextChars=_mm_lddqu_si128(reinterpret_cast<const __m128i*>(reader));
      auto cmp=_mm_cmpeq_epi8(nextChars,bars);
      unsigned mask=_mm_movemask_epi8(cmp);
      if (mask)
         return reader=reader+__builtin_ctzl(mask);
   }
   for (;reader<=readerLimit;++reader)
      if ((*reader)=='\n')
         return reader;
   return nullptr;
#else
   return static_cast<const char*>(memchr(reader, '\n', readerLimit-reader));
#endif
}
//---------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 4) {
    cerr << "Usage: " << argv[0] << "\n<pattern> (spaces are also allowed, simply an underscore instead: 'aa bb' -> 'aa_bb')\n<file>\n<mode>(0 -> hyperKMP, 1 -> normalKMP, 2 -> benchmark to see comparisons)" << endl;
    return 0;
  }
  // ifstream in;
  // in.open(argv[2]);
  
  unsigned mode = atoi(argv[3]);
  if (mode > TEST_MODE) {
    cerr << "Mode not available!" << endl;
    exit(0);
  }
  hyperKMP hyper(argv[1], mode);
  
  const char* str;
  unsigned testCases = 0;
  unsigned countMatches = 0;
  
  boost::iostreams::mapped_file_source f(argv[2]);
  auto fileBegin = f.data();
  auto fileEnd = fileBegin + f.size();
  auto reader = fileBegin;
  
  auto start = std::chrono::high_resolution_clock::now();
  while (true) {
    if (!reader) 
      break;
    ++testCases;
    
    str = reader;
    unsigned n;
    {
      auto eol = findNl(reader, fileEnd);
      n = eol ? (eol - reader) : (fileEnd - reader);
      reader = eol? (eol + 1) : nullptr;
    }
    countMatches += hyper.search(str, n);
  }
  cout << "Matches = " << countMatches << " out of " << testCases << endl;
  if (mode == TEST_MODE) {
    hyper.benchmark();
  } else {
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = stop - start;
    cerr << "Mode (0 -> hyperKMP, 1 -> kmp, 2 -> benchmark) " << mode << " took: " << duration.count() * 1e-6 << "ms" << endl;
  }
  return 0;
}
