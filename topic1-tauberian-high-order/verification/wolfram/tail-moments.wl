Print["Version: ",$Version];
ass=Element[J,Integers] && J>=1 && v>0 && theta>0;
s0=Sum[theta v/2^j,{j,1,J},Assumptions->ass];s1=Sum[j theta v/2^j,{j,1,J},Assumptions->ass];
Print["sum0=",InputForm[s0]];Print["sum1=",InputForm[s1]];
Print["checks=",InputForm[{FullSimplify[s0==theta v(1-2^-J),ass],FullSimplify[s1==theta v(2-(J+2)2^-J),ass]}]];
