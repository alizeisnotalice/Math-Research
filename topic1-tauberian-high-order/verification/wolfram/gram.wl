Print["Version: ", $Version];
lam = 1/(1 + l/2); ass = r >= s > 0 && d >= 0 && l >= 0 && r-s <= l d;
checks = {"noncontainment" -> Resolve[ForAll[{r,s,d,l}, Implies[ass && d>0, lam (r-s)/2 < d]], Reals], "overlap_implies_d_lt_s" -> Resolve[ForAll[{r,s,d,l}, Implies[ass && d < lam (r+s)/2, d < s]], Reals], "first_inequality" -> Resolve[ForAll[{r,s,d,l}, Implies[ass, 1+(r-s)/(2 s)-d/(lam s) <= 1-d/s]], Reals], "second_inequality" -> Resolve[ForAll[{r,s,d,l}, Implies[ass, 1-d/s <= 1-d/r]], Reals]};
Print[InputForm[checks]];
Print["counterexample set: ", InputForm[Reduce[ass && (1+(r-s)/(2 s)-d/(lam s) > 1-d/s || 1-d/s > 1-d/r), {r,s,d,l}, Reals]]];
