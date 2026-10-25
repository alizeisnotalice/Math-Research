% cb5 exact-binary-rational cross-check; no clear and no path change.
% All arithmetic below uses small integers divided by 8, hence is exact in double.
cb5_axis = -2:2;
cb5_num = zeros(5,5,5);
cb5_orbit_numerators = [22,5,2,1];
for cb5_i=1:5
    for cb5_j=1:5
        for cb5_k=1:5
            cb5_vec = [cb5_axis(cb5_i),cb5_axis(cb5_j),cb5_axis(cb5_k)];
            if max(abs(cb5_vec))<2
                cb5_num(cb5_i,cb5_j,cb5_k)=cb5_orbit_numerators(1+nnz(cb5_vec));
            end
        end
    end
end
cb5_pad=zeros(7,7,7);
cb5_pad(2:6,2:6,2:6)=cb5_num;
cb5_lap_num=zeros(5,5,5);
for cb5_i=2:6
    for cb5_j=2:6
        for cb5_k=2:6
            cb5_lap_num(cb5_i-1,cb5_j-1,cb5_k-1)=...
                cb5_pad(cb5_i-1,cb5_j,cb5_k)+cb5_pad(cb5_i+1,cb5_j,cb5_k)+...
                cb5_pad(cb5_i,cb5_j-1,cb5_k)+cb5_pad(cb5_i,cb5_j+1,cb5_k)+...
                cb5_pad(cb5_i,cb5_j,cb5_k-1)+cb5_pad(cb5_i,cb5_j,cb5_k+1)-...
                6*cb5_pad(cb5_i,cb5_j,cb5_k);
        end
    end
end
cb5_off_mask=true(5,5,5); cb5_off_mask(3,3,3)=false;
assert(all(cb5_lap_num(cb5_off_mask)>=0));
assert(cb5_lap_num(3,3,3)==-102);
assert(sum(cb5_lap_num(:))==0);
assert(sum(cb5_num(:))==84);
cb5_face_num=zeros(2,2);
for cb5_j=1:2
    for cb5_k=1:2
        cb5_face_num(cb5_j,cb5_k)=cb5_num(2,cb5_j+2,cb5_k+2)+...
            cb5_num(4,cb5_j+2,cb5_k+2)-2*cb5_num(3,cb5_j+2,cb5_k+2);
    end
end
assert(isequal(cb5_face_num,[-34,-6;-6,-2]));
cb5_face_midpoint=sum(cb5_face_num(:))/32;
cb5_test_pairing=cb5_face_midpoint/16;
assert(cb5_face_midpoint==-3/2 && cb5_test_pairing==-3/32);
fprintf('cb5 PASS: 125 nodes; origin Delta = %.8g; Q1 mass = %.8g; face midpoint = %.8g; test pairing = %.8g\n',...
    cb5_lap_num(3,3,3)/8,sum(cb5_num(:))/8,cb5_face_midpoint,cb5_test_pairing);
