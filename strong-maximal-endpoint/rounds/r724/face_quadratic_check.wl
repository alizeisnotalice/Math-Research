mm = 1-aa z^2/2+uu z^3;
ff = 1-bb z^2/2+vv z^3;
ex = 1-(mm/.z->-z)ff-mm(ff/.z->-z)+mm(mm/.z->-z);
Print[Coefficient[Expand[ex],z,2]];
Print[Coefficient[Expand[ex],z,3]];
