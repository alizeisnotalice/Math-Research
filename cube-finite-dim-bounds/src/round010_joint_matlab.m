joint10_root=fileparts(fileparts(mfilename('fullpath')));
joint10_in=jsondecode(fileread(fullfile(joint10_root,'results','round010_matlab_input.json')));
joint10_options=optimoptions('linprog','Display','none');
[joint10_x,joint10_f,joint10_flag,joint10_output]=linprog(joint10_in.objective,joint10_in.A,joint10_in.rhs,[],[],[0;0],[],joint10_options);
assert(joint10_flag>0);
assert(max(joint10_in.A*joint10_x-joint10_in.rhs)<1e-8);
assert(abs(-joint10_f-joint10_in.python_gain)<1e-8);
fprintf('JOINT10 weights %.17g %.17g gain %.17g flag %d max_violation %.17g\n',joint10_x(1),joint10_x(2),-joint10_f,joint10_flag,max(joint10_in.A*joint10_x-joint10_in.rhs));
disp(joint10_output.algorithm);
