/* A small cross-implementation bridge: no extra GF implementation. */
#define main search_main_not_called
#include "gf_check.c"
#undef main
int main(int argc,char**argv){if(argc!=4)return 2;d=atoi(argv[1]);uint64_t poly=strtoull(argv[2],0,0),g=strtoull(argv[3],0,0);top=UINT64_C(1)<<d;mask=top-1;lo=poly^top;uint64_t a,b;while(scanf("%"SCNu64" %"SCNu64,&a,&b)==2){uint64_t inv=powgf(b,mask-1);printf("%"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n",mul(a,b),inv,mul(a,inv),mul(g,a));}return 0;}
