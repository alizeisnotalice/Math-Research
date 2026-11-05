(* Exact compact direction at fixed h. No numerical optimization claim. *)
Print["Version: ", $Version];
b[x_] := Piecewise[{{(4-6 x^2+3 Abs[x]^3)/6, Abs[x] <= 1}, {(2-Abs[x])^3/6, 1<Abs[x]<2}}, 0];
oneD = Integrate[b[x], {x, -2, 2}];
square = Integrate[b[x]^2, {x, -2, 2}];
h=1/32;
Print["Integral B = ", oneD];
Print["Integral B^2 = ", square];
Print["Origin tensor halfmass = ", h^2 oneD^2/2];
Print["Perturbation halfmass = ", -epsilon h^2 oneD^2/2];
Print["Support = [-1/16,1/16]^2; symmetric nonproportional perturbation."];
If[oneD =!= 1 || h^2 oneD^2/2 =!= 1/2048, Exit[2]];
Exit[0];
