% Four registered shears only; exact-index identities/numerical geometry crosscheck.
n4shP=[1+sqrt(2)/2,(2+2^(3/4)*3^(1/4))/4,(1+2^(1/4))/2,(1+2^(1/4))/2];
n4shnum=[-1,-1,1,1];n4shden=[16,32,32,16];n4shout=zeros(4,6);
for n4shi=1:4
 n4shp=n4shnum(n4shi);n4shq=n4shden(n4shi);n4sht=n4shp/n4shq;n4shS=eye(4);n4shS(1,2)=n4sht*n4shP(1)/n4shP(2);assert(abs(det(n4shS)-1)<1e-15);
 n4shcount=0;
 for n4shk1=-6:6
  for n4shk2=-6:6
   n4shm=floor(n4shk2/n4shq);n4shj=n4shk2-n4shq*n4shm;n4shc=floor(n4shp*n4shj/n4shq);n4shphi=n4shp*n4shj/n4shq-n4shc;n4shn1=n4shk1+n4shp*n4shm+n4shc;
   assert(n4shj>=0 && n4shj<n4shq);assert(abs(n4shn1+n4shphi-(n4shk1+n4sht*n4shk2))<1e-14);assert(n4shk1==n4shn1-n4shp*n4shm-n4shc);n4shcount=n4shcount+1;
  end
 end
 n4shatoms=[];
 for n4shk1=0:1
  for n4shk2=0:2
   for n4shk3=0:1
    for n4shk4=0:1
     n4shatoms(end+1,:)=[(n4shk1+n4sht*n4shk2)*n4shP(1),n4shk2*n4shP(2),n4shk3*n4shP(3),n4shk4*n4shP(4)];
    end
   end
  end
 end
 n4shspan=max(n4shatoms,[],1)-min(n4shatoms,[],1);assert(abs(n4shspan(1)-n4shP(1)*(1+2*abs(n4sht)))<1e-14);assert(size(n4shatoms,1)==24);
 assert(abs(det(diag([n4shP(1),n4shq*n4shP(2),n4shP(3),n4shP(4)]))-n4shq*prod(n4shP))<1e-12);
 n4shout(n4shi,:)=[n4sht,n4shcount,n4shq,24,24^(1/4)/2,n4shspan(1)];
end
format long g;disp(n4shout)
