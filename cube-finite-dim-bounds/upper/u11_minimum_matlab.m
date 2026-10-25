% 11K independent finite certificate sign/denominator calibration.
% Full exact rational reconstruction is the Python/FLINT checker, not this script.
% No G sampling, I/cost/tau computation, clear, or path changes.
upU11MinFile=fullfile(fileparts(mfilename('fullpath')),'u11_minimum_certificate_results.json');
upU11MinObj=jsondecode(fileread(upU11MinFile));
upU11MinExpected=[1739,1702,1702];upU11MinTotal=0;
assert(pi<22/7);
for upU11MinI=1:3
    upU11MinCase=upU11MinObj.physical_cases(upU11MinI);
    upU11MinCoeff=upU11MinCase.P_coefficients;
    assert(numel(upU11MinCoeff)==upU11MinExpected(upU11MinI));
    for upU11MinJ=1:numel(upU11MinCoeff)
        upU11MinText=char(upU11MinCoeff(upU11MinJ).value);
        % Canonical nonnegative integer/rational: sign tested as text, no overflow.
        assert(~isempty(regexp(upU11MinText,'^[0-9]+(/[0-9]+)?$','once')));
        assert(sum(upU11MinCoeff(upU11MinJ).powers)==20);
    end
    upU11MinR=double(upU11MinCase.cone_rays).';
    upU11MinFactors=upU11MinCase.denominator_factor_sign_proof;
    for upU11MinJ=1:numel(upU11MinFactors)
        upU11MinExpr=char(upU11MinFactors(upU11MinJ).factor);
        assert(~isempty(regexp(upU11MinExpr,'^[a-d0-9+\- *]+$','once')));
        upU11MinFun=str2func(['@(a,b,c,d)',upU11MinExpr]);
        upU11MinValues=zeros(1,4);
        for upU11MinK=1:4
            upU11MinX=upU11MinR(:,upU11MinK);
            upU11MinValues(upU11MinK)=upU11MinFun(upU11MinX(1),upU11MinX(2),upU11MinX(3),upU11MinX(4));
        end
        assert(all(upU11MinValues>=0) && any(upU11MinValues>0));
    end
    upU11MinTotal=upU11MinTotal+numel(upU11MinCoeff);
    fprintf('11K case%d nonnegative_coefficients %d degree20 denominator_linear_factors %d PASS\n',upU11MinI,numel(upU11MinCoeff),numel(upU11MinFactors));
end
fprintf('11K total_coefficients %d target %.17g pi_upper %.17g PASS; no standalone weak-type claim\n',upU11MinTotal,21/25000,22/7);
