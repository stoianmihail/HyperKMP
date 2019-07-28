#include <vector>
#include <unordered_map>
//---------------------------------------------------------------------------
class hyperKMP;
//---------------------------------------------------------------------------
class hyperKMP 
{
  private:
  unsigned length, mode;
  char* pattern;
  unsigned* pi;
  unsigned* psi;
  unsigned* omega;
  unsigned* degrees;
  unsigned** edges;
  bool* valid;
  
  void compressPi();
  void dfs(unsigned state, std::unordered_map<char, unsigned>& indexes, std::vector<unsigned>& activeNodes, unsigned minimumHalt);
  
  // For benchmarking
  unsigned testCases, hyperSum, normalSum, hyperMaximum, normalMaximum;
      
  public:
  hyperKMP();
  hyperKMP(char* pattern, unsigned mode);
  // TODO: After benchmark, place the const back!
  bool search(char* str);
  void benchmark();
};
//---------------------------------------------------------------------------
