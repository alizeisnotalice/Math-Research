#include <vector>
#include <cmath>
#include <algorithm>
extern "C" void profiles(int ns, int n, int nr, const double* rr,
                        const double* bits, const double* gg,
                        int shift, double central_mass, double shifted_mass,
                        double band_mass, double* out) {
  std::vector<double> a(n+1),b(n+1);
  int mid=n/2;
  for(int x=0;x<ns;++x) {
    for(int ir=0;ir<nr;++ir) {
      double vals[3];
      for(int per=0;per<3;++per) {
        std::fill(a.begin(),a.end(),0.0);
        std::fill(b.begin(),b.end(),0.0); a[0]=1;
        for(int i=0;i<n;++i) {
          int id=(per*ns+x)*n+i;
          double p=(1-rr[ir])*bits[id]+rr[ir]*gg[id];
          b[0]=a[0]*(1-p);
          for(int k=1;k<=i+1;++k) b[k]=a[k]*(1-p)+a[k-1]*p;
          a.swap(b);
        }
        if(per==0) vals[0]=(a[mid-1]+a[mid]+a[mid+1])/band_mass;
        if(per==1) vals[1]=(a[mid-shift]+a[mid+shift])/shifted_mass;
        if(per==2) {
          vals[2]=a[mid]/central_mass;
          out[(0*ns+x)*nr+ir]=vals[2];
          out[(1*ns+x)*nr+ir]=(a[mid-shift]+a[mid+shift])/shifted_mass;
          out[(2*ns+x)*nr+ir]=(a[mid-1]+a[mid]+a[mid+1])/band_mass;
        }
      }
      out[(3*ns+x)*nr+ir]=(vals[0]+vals[1]+vals[2])/3;
    }
  }
}
