#!/usr/bin/env python3
"""Materialize a minimal analytic-aware checker copy using stdlib only.

No imports of flint/checkers and no mathematical execution. The only allowed
AST changes are four log keywords, eighteen pw arguments and one pw parameter.
"""
import ast
from copy import deepcopy
import difflib
from hashlib import sha256
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ORIGINAL_SHA='4bd7e823bdcb336f8e20c188d68d6d8c1e4aa22ec52353885059c24445317db5'
DRIVER_SHA='079906db354febdd57d537f0d747428bc90bd37fa88297449cd1d46c7446310d'
BASELINE_SHA='21d39978590d3659021c093b3100e4e4d5e8dafe4834d46aa55d6c86c6a5e172'


class EraseOnlyForwarding(ast.NodeTransformer):
    def __init__(self):self.counts={'log_keywords':0,'pw_arguments':0,'pw_parameters':0}
    def visit_Call(self,node):
        self.generic_visit(node)
        if isinstance(node.func,ast.Attribute) and node.func.attr=='log' and node.keywords:
            if len(node.keywords)!=1 or node.keywords[0].arg!='analytic' or not isinstance(node.keywords[0].value,ast.Name) or node.keywords[0].value.id not in ('an','analytic'):
                raise ValueError('Unexpected log patch AST.')
            node.keywords=[];self.counts['log_keywords']+=1
        if isinstance(node.func,ast.Name) and node.func.id=='pw' and len(node.args)==3 and isinstance(node.args[-1],ast.Name) and node.args[-1].id=='an':
            node.args.pop();self.counts['pw_arguments']+=1
        return node
    def visit_Lambda(self,node):
        self.generic_visit(node)
        if [arg.arg for arg in node.args.args]==['z','e','an']:
            node.args.args.pop();self.counts['pw_parameters']+=1
        return node


def main():
    target=ROOT/'fixed_checker'
    if (ROOT/'patch_manifest.json').exists() or (ROOT/'hlcert.patch').exists() or (target.exists() and any(target.iterdir())):
        raise SystemExit('Existing patch artifacts; no overwrite or regenerated provenance.')
    source=ROOT/'original_checker'
    if not (source/'hlcert.py').is_file():source=ROOT.parent/'upstream-v4/checker'
    original_bytes=(source/'hlcert.py').read_bytes();driver=(source/'cert_run.py').read_bytes()
    if sha256(original_bytes).hexdigest()!=ORIGINAL_SHA or sha256(driver).hexdigest()!=DRIVER_SHA:
        raise ValueError('Pinned original checker/driver SHA mismatch.')
    original=original_bytes.decode('utf-8')
    h0=original.index('    def _Hint(');h1=original.index('    def _outside_taylor(',h0)
    hint=original[h0:h1]
    old_pw='        pw = lambda z, e: ((z.log()) * e).exp()'
    new_pw='        pw = lambda z, e, an: ((z.log(analytic=an)) * e).exp()'
    if hint.count(old_pw)!=1:raise ValueError('Original _Hint pw definition changed.')
    hint=hint.replace(old_pw,new_pw)
    patterns=[('pw(t + bb, -alc)',3),('pw(t + bb, -alc - 1)',2),('pw(t + bb, -alc - 2)',1),
        ('pw(aa - t, -1 - alc)',3),('pw(aa + t, -1 - alc)',3),
        ('pw(aa - t, -2 - alc)',2),('pw(aa + t, -2 - alc)',2),
        ('pw(aa - t, -3 - alc)',1),('pw(aa + t, -3 - alc)',1)]
    replacements=[]
    for text,count in patterns:
        if hint.count(text)!=count:raise ValueError('Original _Hint callback occurrence count changed: '+text)
        replacement=text[:-1]+', an)';hint=hint.replace(text,replacement)
        replacements.append({'original':text,'fixed':replacement,'occurrences':count})
    patched=original[:h0]+hint+original[h1:]
    p0=patched.index('    def Apsi_outside_point(')
    p1=patched.index('    # ---------------------------------------------------------------- K',p0)
    exterior=patched[p0:p1]
    if exterior.count('.log()')!=3:raise ValueError('Original exterior point logarithms changed.')
    exterior=exterior.replace('.log()','.log(analytic=analytic)')
    patched=patched[:p0]+exterior+patched[p1:]
    old_ast=ast.parse(original);new_ast=ast.parse(patched);eraser=EraseOnlyForwarding()
    restored=eraser.visit(deepcopy(new_ast))
    expected={'log_keywords':4,'pw_arguments':18,'pw_parameters':1}
    if eraser.counts!=expected:raise ValueError('Forwarding-only AST count mismatch.')
    if ast.dump(old_ast,include_attributes=False)!=ast.dump(restored,include_attributes=False):
        raise ValueError('Non-forwarding formula/control AST change detected.')
    ast.parse(driver.decode('utf-8'))
    fixed_bytes=patched.encode('utf-8')
    diff=''.join(difflib.unified_diff(original.splitlines(keepends=True),patched.splitlines(keepends=True),
        fromfile='upstream/checker/hlcert.py',tofile='fixed_checker/hlcert.py'))
    if not diff:raise ValueError('No actual patch produced.')
    target.mkdir(exist_ok=True)
    (target/'hlcert.py').write_bytes(fixed_bytes);(target/'cert_run.py').write_bytes(driver)
    (ROOT/'hlcert.patch').write_text(diff)
    manifest={'status':'static_patch_materialized_no_math_execution','original_hlcert_sha256':ORIGINAL_SHA,
        'fixed_hlcert_sha256':sha256(fixed_bytes).hexdigest(),'original_cert_run_sha256':DRIVER_SHA,
        'fixed_cert_run_sha256':sha256(driver).hexdigest(),'patch_sha256':sha256(diff.encode()).hexdigest(),
        'make_patch_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'baseline_sha256':BASELINE_SHA,'forwarding_only_AST_normalization_passed':True,
        'AST_forwarding_counts':eraser.counts,'hint_pw_replacements':replacements,
        'driver_byte_identical':driver==(target/'cert_run.py').read_bytes(),
        'original_files_modified':False,'flint_imports':0,'mathematical_executions':0,
        'fixed_files':{'hlcert.py':'fixed_checker/hlcert.py','cert_run.py':'fixed_checker/cert_run.py'},
        'planned_unique_math_run':{'kernel':'B1 only','hard_seconds':180,'Q':8,'maxdepth':7,
            'Kdepth':6,'procs':6,'only':'MKA','precision':'Original 96bit unchanged'}}
    (ROOT/'patch_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':manifest['status'],'fixed_hlcert_sha256':manifest['fixed_hlcert_sha256'],
        'driver_sha256':DRIVER_SHA,'AST_counts':eraser.counts,'mathematical_executions':0},indent=2))


if __name__=='__main__':main()
