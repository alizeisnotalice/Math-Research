% New geometry/event interface examples only; no old union-volume rerun.
% Main agent executes serially, all variables prefixed n3r11.
n3r11root=fileparts(mfilename('fullpath'));n3r11seed=jsondecode(fileread(fullfile(n3r11root,'round008_certificate_input.json')));n3r11ex=jsondecode(fileread(fullfile(n3r11root,'round011_examples.json')));n3r11p=zeros(1,3);n3r11ph=zeros(3,3);n3r11w=zeros(1,3);
for n3r11j=1:3,n3r11p(n3r11j)=n3r11parse(n3r11seed.periods{n3r11j});n3r11w(n3r11j)=n3r11parse(n3r11seed.weights{n3r11j});end
for n3r11b=1:3,for n3r11j=1:3,n3r11ph(n3r11b,n3r11j)=n3r11parse(n3r11seed.phases{n3r11b}{n3r11j});end,end
n3r11tail=sum(n3r11w)*prod(1./n3r11p+1/6);assert(n3r11tail<1);fprintf('N3 round011 input-specific r>=3 density tail bound=%.17g <1\n',n3r11tail);
for n3r11i=1:numel(n3r11ex.stable_center_probes)
 n3r11row=n3r11ex.stable_center_probes(n3r11i);n3r11x=zeros(1,3);for n3r11j=1:3,n3r11x(n3r11j)=n3r11parse(n3r11row.point{n3r11j});end;n3r11r=n3r11parse(n3r11row.radius);n3r11lo=ceil((n3r11x-n3r11r)./n3r11p-n3r11ph);n3r11hi=floor((n3r11x+n3r11r)./n3r11p-n3r11ph);n3r11counts=prod(max(0,n3r11hi-n3r11lo+1),2)';assert(isequal(n3r11counts,double(n3r11row.orbit_counts(:)')));n3r11m=sum(n3r11counts.*n3r11w);fprintf('N3 round011 center%d r=%.17g counts=[%g %g %g] mass=%.17g\n',n3r11row.center_id,n3r11r,n3r11counts(1),n3r11counts(2),n3r11counts(3),n3r11m);
end
for n3r11i=1:numel(n3r11ex.strict_empty_examples)
 n3r11row=n3r11ex.strict_empty_examples(n3r11i);n3r11key=n3r11row.key(:)';n3r11lo=reshape(n3r11key(1:9),3,3)';n3r11hi=reshape(n3r11key(10:18),3,3)';n3r11counts=prod(max(0,n3r11hi-n3r11lo+1),2);n3r11active=n3r11counts>0;n3r11m=sum(n3r11counts'.*n3r11w);n3r11LC=max((n3r11hi(n3r11active,:)+n3r11ph(n3r11active,:)).*n3r11p,[],1);n3r11HC=min((n3r11lo(n3r11active,:)+n3r11ph(n3r11active,:)).*n3r11p,[],1);n3r11axis=n3r11row.axis+1;
 if n3r11row.kind==0,n3r11rhs=(n3r11LC(n3r11axis)-n3r11HC(n3r11axis))^3;elseif n3r11row.kind==1,n3r11rhs=8*(n3r11LC(n3r11axis)-n3r11p(n3r11axis))^3;else,n3r11rhs=8*(-n3r11HC(n3r11axis))^3;end
 assert(n3r11rhs>n3r11m);fprintf('N3 round011 rawkey%d kind%d axis%d cubic gap=%.17g >0\n',n3r11row.raw_key_index,n3r11row.kind,n3r11axis,n3r11rhs-n3r11m);
end
function n3r11v=n3r11parse(n3r11s)
 n3r11parts=split(string(n3r11s),'/');n3r11v=str2double(n3r11parts(1));if numel(n3r11parts)==2,n3r11v=n3r11v/str2double(n3r11parts(2));end
end
