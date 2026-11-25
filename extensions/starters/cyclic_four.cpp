// One T(6) orbit on 15 parts of size four, including a short orbit of four points.
#include <algorithm>
#include <array>
#include <bit>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
struct V{int e,s,x;};
int diff28[28][28],sameClass[28],M[28],mode,mult2,mu7;vector<V>tile;long long nodes=0,finals=0,colours=0;
int rot(int e){int t=0;int p[]={0,1,2,4,5,3};for(int i=0;i<6;i++)if(e>>i&1)t+=1<<p[i];return t;}
int bit(int a,int b){return 1<<a|1<<b;}
int edgeid(V a,V b){
 if(a.s>b.s)swap(a,b);
 if(a.s==2)return -1;
 if(a.s==b.s){
   int d=diff28[a.x][b.x];if(d/4==0)return -1;
   return 12*a.s+sameClass[d];
 }
 if(b.s==2)return 52+4*a.s+((b.x-a.x+4)%4);
 return 24+diff28[a.x][b.x];
}
bool put(V a,uint64_t&m){
 for(auto b:tile){if(a.s==b.s&&a.x==b.x)return false;if(!(a.e&b.e))continue;int i=edgeid(a,b);if(i<0||(m>>i&1))return false;m|=1ULL<<i;}
 tile.push_back(a);return true;
}
vector<pair<int,int>> freeOrb;
bool rec(int i,uint64_t m){
 if(i==(int)freeOrb.size()){
  ++finals;if(popcount(m)!=60)return false;
  cout<<mode;for(auto a:tile)cout<<' '<<a.e<<' '<<a.s<<' '<<a.x;cout<<'\n';return false;
 }++nodes;auto [s,e]=freeOrb[i];
 for(int x=0;x<28;x++){
  int period=e==rot(e)?1:3;
  if(period==1&&x!=M[x])continue;
  if(period==3&&x==M[x])continue;
  int cur=x,edge=e;uint64_t z=m;bool ok=true;auto sz=tile.size();
  for(int j=0;j<period;j++){if(!put({edge,s,cur},z)){ok=false;break;}cur=M[cur];edge=rot(edge);}
  if(ok&&rec(i+1,z))return true;tile.resize(sz);
 }return false;
}
int main(){
 for(mode=3;mode<=4;mode++){
  for(int x=0;x<28;x++)for(int y=0;y<28;y++){
   int sign=(mode==4 && (x&1))?-1:1;
   diff28[x][y]=4*((sign*(y/4-x/4)+7)%7)+(y-x+28)%4;
  }
  std::fill(sameClass,sameClass+28,-1);int num=0;
  for(int x=4;x<28;x++)if(sameClass[x]<0){sameClass[x]=sameClass[diff28[x][0]]=num++;}
  if(num!=12)return 1;
  for(int x=0;x<28;x++)M[x]=4*(2*(x/4)%7)+x%4;
  vector<int>E;set<int>seen;for(int a=0;a<6;a++)for(int b=0;b<a;b++){int e=bit(a,b);if(seen.count(e))continue;E.push_back(e);for(int t=0;t<3;t++,e=rot(e))seen.insert(e);}
  int root=bit(0,1);vector<int>oth;for(int e:E)if(e!=root)oth.push_back(e);
  for(unsigned mask=0;mask<1U<<oth.size();mask++){
   vector<V>cols{{root,2,0}};for(unsigned i=0;i<oth.size();i++){int e=oth[i],v=e;do{cols.push_back({v,int(mask>>i&1),0});v=rot(v);}while(v!=e);}
   array<int,3>ct{};int same=0,short0=0;
   for(auto a:cols){ct[a.s]++;for(auto b:cols)if(a.e>b.e&&(a.e&b.e)){same+=a.s==0&&b.s==0;short0+=(a.s==2&&b.s==0)||(a.s==0&&b.s==2);}}
   if(ct[0]!=7||ct[1]!=7||same!=12||short0!=4)continue;
   ++colours;freeOrb.clear();tile.clear();uint64_t m=0;
   put({root,2,0},m);
   // Centralizer translations of each full sheet need not preserve the short-to-sheet labels;
   // they permute the four legitimate difference classes. They normalize its fixed entry to zero.
   bool ok=true;for(unsigned i=0;i<oth.size();i++){int e=oth[i],s=mask>>i&1;
    if(e==rot(e))ok&=put({e,s,0},m);else freeOrb.push_back({s,e});}
   if(ok&&rec(0,m))return 0;
  }
 }cerr<<"colour patterns "<<colours<<" nodes "<<nodes<<" full "<<finals<<"\n";
}
