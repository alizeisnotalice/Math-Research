ClearAll[k,g,c,v];
Print["Wolfram ",$Version];
Print["Assumed R16 model constants I/u^(15/2)=k³, J/u^(15/2)=k^5/64; k>0."];
ii=k^3;jj=k^5/64;
Print["I/J: ",FullSimplify[ii/jj,k>0]];
Print["k=6: I,J,I/J,I-J = ",{ii,jj,ii/jj,ii-jj}/.k->6];
Print["Young completion residual: ",FullSimplify[-g+c v Sqrt[g]-(c^2 v^2/4-(Sqrt[g]-c v/2)^2),g>=0&&v>=0&&c>=0]];
Print["Young bound violation: ",Reduce[g>=0&&c>=0&&v>=0&&-g+c v Sqrt[g]>c^2 v^2/4,{g,c,v},Reals]];
