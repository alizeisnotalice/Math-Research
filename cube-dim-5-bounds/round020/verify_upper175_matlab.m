% New175 candidate scalar cross-check. Private MATLAB batch, no clear/path reset.
math175dir = fileparts(mfilename('fullpath'));
math175out = fullfile(math175dir,'matlab175_result.json');
assert(~isfile(math175out));
math175clock = tic;
math175data = jsondecode(fileread(fullfile(math175dir,'matlab175_inputs.json')));
math175new = str2sym(math175data.candidate_UH);
math175old = str2sym(math175data.old_UH);
math175gain = str2sym(math175data.gain);
math175g = str2sym(math175data.g);
math175J = str2sym(math175data.J);
math175I = math175new + math175J;
math175V = sym(3125)*math175I^4/math175g^5;
math175oldV = str2sym(math175data.old_fourth_power);
math175portable = str2sym(math175data.old_portable);
math175flags = [isAlways(math175new>0),isAlways(math175old>math175new), ...
 isAlways(math175I>0),isAlways(math175V<math175oldV), ...
 isAlways(math175V<sym(32)^4),isAlways(math175V<math175portable^4), ...
 isAlways(math175old-math175new==math175gain),isAlways(math175V>0)];
assert(all(math175flags));
math175receipt = struct('status','PASS_8_new175_exact_scalar_checks', ...
 'checks',math175flags,'seconds',toc(math175clock),'pid',feature('getpid'), ...
 'version',version,'candidate_input_sha256',math175data.candidate_sha256, ...
 'new_real_node_roots_gains_H_log_calls',0, ...
 'qualification','candidate scalar cross-check;full175 audit and main proof required');
math175write(math175out,jsonencode(math175receipt));
disp(math175receipt.status);
function math175write(math175path,math175text)
 math175fid=fopen(math175path,'w','n','UTF-8');assert(math175fid>=0);
 math175closer=onCleanup(@()fclose(math175fid));
 fwrite(math175fid,math175text,'char');fwrite(math175fid,newline,'char');
end
