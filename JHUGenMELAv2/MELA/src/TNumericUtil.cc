#include "TNumericUtil.hh"
#include <numeric>
#include <algorithm>

using namespace std;

void TNumericUtil::PermutationGenerator(int n, int k, std::vector<std::vector<int>>& perms, int valbegin, int increment){
  if (n<=0 || k<=0 || k>n) return;
  vector<int> d(n);
  iota(d.begin(), d.end(), 1);
  do{
    vector<int> aperm;
    for (int i=0; i<k; i++) aperm.push_back(d[i]);
    perms.push_back(aperm);
    reverse(d.begin()+k, d.end());
  } while (next_permutation(d.begin(), d.end()));
  if (valbegin!=1 || increment!=1){ for (auto& p:perms){ for (auto& pp:p) pp = increment*(pp-1)+valbegin; } }
}
void TNumericUtil::CombinationGenerator(int n, int k, std::vector<std::vector<int>>& perms, int valbegin, int increment){
  if (n==k && k>0){
    vector<int> d(n);
    iota(d.begin(), d.end(), 1);
    perms.push_back(d);
  }
  else{
    PermutationGenerator(n, k, perms);
    for (auto& perm:perms) sort(perm.begin(), perm.end());
    for (auto ip=perms.begin(); ip!=perms.end(); ip++){
      for (auto jp=perms.rbegin(); (jp.base()-1)!=ip; jp++){
        if (*(ip)==*(jp)) perms.erase(jp.base()-1);
      }
    }
  }
  if (valbegin!=1 || increment!=1){ for (auto& p:perms){ for (auto& pp:p) pp = increment*(pp-1)+valbegin; } }
}
