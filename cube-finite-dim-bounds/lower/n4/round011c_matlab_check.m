% Fixed registered endpoints only; numerical crosscheck, no search.
n4r11cP0=[1+sqrt(2)/2,(2+2^(3/4)*3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2];
n4r11cs=[1,-1,-1,1];n4r11ccos=[0,0;0,1;1,0;1,1];
n4r11crows=zeros(6,5);n4r11cidx=0;
for n4r11ck=0:2
 for n4r11ct=[1/20,1]
  n4r11cP=n4r11cP0;n4r11ch=[n4r11ct/50,0,0,0];n4r11ceps=-n4r11ct/5;
  if n4r11ck==0;n4r11ceps=n4r11ct/10;n4r11ch=[n4r11ct/50,n4r11ct/50,n4r11ct/100,0];end
  if n4r11ck==2;n4r11cP(1)=n4r11cP0(1)/(1+n4r11ct/25);end
  n4r11clpts=[];n4r11chpts=[];
  for n4r11ci=1:4
   n4r11co=n4r11ccos(n4r11ci,:);n4r11cp=[];
   for n4r11ck1=0:1
    for n4r11ck2=0:2
     for n4r11ck3=0:1
      for n4r11ck4=0:2
       if mod(n4r11ck1,2)==n4r11co(1) && mod(n4r11ck2,2)==n4r11co(2)
        n4r11cp(end+1,:)=( [n4r11ck1,n4r11ck2,n4r11ck3,n4r11ck4]+n4r11cs(n4r11ci)*n4r11ch).*n4r11cP;
       end
      end
     end
    end
   end
   n4r11chpts=[n4r11chpts;n4r11cp];
   if n4r11ci~=4;n4r11clpts=[n4r11clpts;n4r11cp];end
  end
  n4r11crL=(30-6*n4r11ceps)^(1/4)/2;n4r11crH=36^(1/4)/2;
  n4r11cF=2*n4r11cP(1)*n4r11ch(1)+n4r11crL-n4r11crH;
  n4r11clowL=max(n4r11clpts,[],1)-n4r11crL;n4r11clowH=max(n4r11chpts,[],1)-n4r11crH;
  assert(abs(n4r11clowH(1)-n4r11clowL(1)-n4r11cF)<1e-12);
  assert(size(n4r11clpts,1)==30 && size(n4r11chpts,1)==36);
  assert(abs(sum(1+n4r11cs*n4r11ceps)-4)<1e-14);
  assert((n4r11ct==1 && n4r11cF>0)||(n4r11ct<1 && n4r11cF<0));
  n4r11cidx=n4r11cidx+1;n4r11crows(n4r11cidx,:)=[n4r11ck,n4r11ct,30-6*n4r11ceps,n4r11cF,4];
 end
end
format long g;disp(n4r11crows)
