/* Exact GF(2^d) arithmetic, O(B^2) storage and O(B^2 log B) sort checks.
   Representation: bit i = coefficient of X^i, modulus includes X^d.
   No logarithms or full q-1 cycle are formed. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
static unsigned d; static uint64_t top,lo,mask;
static uint64_t mul(uint64_t a,uint64_t b){uint64_t z=0;while(b){z^=a&-(b&1);b>>=1;uint64_t carry=a>>(d-1);a=((a<<1)&mask)^(lo&-carry);}return z;}
static uint64_t powgf(uint64_t a,uint64_t e){uint64_t z=1;while(e){if(e&1)z=mul(z,a);a=mul(a,a);e>>=1;}return z;}
static int cmp(const void*a,const void*b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return(x>y)-(x<y);}
static uint64_t rng(uint64_t *s){uint64_t z=(*s+=UINT64_C(0x9e3779b97f4a7c15));z=(z^(z>>30))*UINT64_C(0xbf58476d1ce4e5b9);z=(z^(z>>27))*UINT64_C(0x94d049bb133111eb);return z^(z>>31);}
static unsigned rankgf(uint64_t *basis,unsigned b){uint64_t p[64]={0};unsigned r=0;for(unsigned i=0;i<b;i++){uint64_t x=basis[i];for(int j=d-1;j>=0;j--)if((x>>j)&1){if(p[j])x^=p[j];else{p[j]=x;r++;break;}}}return r;}
static uint64_t duplicates(uint64_t*a,size_t n,uint64_t*first){uint64_t c=0;for(size_t i=1;i<n;i++)if(a[i]==a[i-1]){if(!c)*first=a[i];c++;}return c;}
static uint64_t hash(uint64_t*a,size_t n){uint64_t h=UINT64_C(14695981039346656037);for(size_t i=0;i<n;i++)for(int j=0;j<8;j++){h^=(a[i]>>(8*j))&255;h*=UINT64_C(1099511628211);}return h;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec/1e9;}
int main(int argc,char**argv){if(argc<7){fprintf(stderr,"usage: gf_check d modulus b generator seed attempts [basis...]\n");return 2;}d=atoi(argv[1]);uint64_t poly=strtoull(argv[2],0,0);unsigned b=atoi(argv[3]);uint64_t g=strtoull(argv[4],0,0),seed=strtoull(argv[5],0,0);unsigned attempts=atoi(argv[6]);if(d>=63||b>=24||!b||d<b)return 2;top=UINT64_C(1)<<d;mask=top-1;lo=poly^top;if(poly<top||poly>=2*top)return 2;
 size_t B=(size_t)1<<b,N=B-1,T=N*(N+1)/2,R=N*(N-1);uint64_t *basis=calloc(b,8),*W=calloc(B,8),*inv=calloc(B,8),*a=malloc((R>T?R:T)*8),*shift=malloc(R*8);if(!basis||!W||!inv||!a||!shift)return 3;
 for(unsigned attempt=0;attempt<attempts;attempt++){double start=now();if(argc==7+(int)b){for(unsigned i=0;i<b;i++)basis[i]=strtoull(argv[7+i],0,0);}else for(unsigned i=0;i<b;i++)basis[i]=rng(&seed)&mask;
 unsigned rank=rankgf(basis,b);W[0]=0;for(size_t i=1;i<B;i++){unsigned k=__builtin_ctzll(i);W[i]=W[i^((size_t)1<<k)]^basis[k];}uint64_t pc=0,rc=0,ac=0,shortc=0,pfirst=0,rfirst=0,afirst=0,ph=0,rh=0;size_t t=0;
 if(rank==b){for(size_t i=1;i<B;i++)for(size_t j=i;j<B;j++)a[t++]=mul(W[i],W[j]);qsort(a,T,8,cmp);pc=duplicates(a,T,&pfirst);ph=hash(a,T);
 for(size_t i=1;i<B;i++)inv[i]=powgf(W[i],mask-1);t=0;for(size_t i=1;i<B;i++)for(size_t j=1;j<B;j++)if(i!=j)a[t++]=mul(W[i],inv[j]);qsort(a,R,8,cmp);rc=duplicates(a,R,&rfirst);rh=hash(a,R);
 for(size_t i=0;i<R;i++)shift[i]=mul(g,a[i]);qsort(shift,R,8,cmp);size_t i=0,j=0;while(i<R&&j<R){if(a[i]==shift[j]){if(!ac)afirst=a[i];ac++;i++;j++;}else if(a[i]<shift[j])i++;else j++;}
 uint64_t gi=powgf(g,mask-1),shorts[4]={g,gi,mul(g,g),mul(gi,gi)};for(unsigned k=0;k<4;k++)if(bsearch(shorts+k,a,R,8,cmp))shortc++;
 }
 int ok=rank==b&&!pc&&!rc&&!ac&&!shortc;printf("{\"attempt\":%u,\"basis\":[",attempt);for(unsigned i=0;i<b;i++)printf("%s%"PRIu64,i?",":"",basis[i]);printf("],\"rank\":%u,\"product_count\":%zu,\"product_duplicate_adjacent_count\":%"PRIu64",\"product_first_collision\":%"PRIu64",\"product_sorted_fnv1a64\":\"%016"PRIx64"\",\"ratio_count\":%zu,\"ratio_duplicate_adjacent_count\":%"PRIu64",\"ratio_first_collision\":%"PRIu64",\"ratio_sorted_fnv1a64\":\"%016"PRIx64"\",\"S_intersect_gS_count\":%"PRIu64",\"adjacency_first_value\":%"PRIu64",\"short_forbidden_hits\":%"PRIu64",\"passed\":%s,\"seconds\":%.6f}\n",rank,T,pc,pfirst,ph,R,rc,rfirst,rh,ac,afirst,shortc,ok?"true":"false",now()-start);fflush(stdout);if(ok)break;if(argc==7+(int)b)break;
 }free(basis);free(W);free(inv);free(a);free(shift);return 0;}
