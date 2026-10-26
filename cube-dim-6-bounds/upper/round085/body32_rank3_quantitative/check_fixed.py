#!/usr/bin/env python3
"""A085: one frozen seven-configuration local body32-U batch, no 6D kernel."""
from fractions import Fraction as F
from itertools import product
from math import factorial, comb, prod, pi
from pathlib import Path
import hashlib
import json
import os
import resource
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
resource.setrlimit(resource.RLIMIT_CPU, (27, 27))
START, CPU_START = time.monotonic(), time.process_time()
OUT = {"schema": "A085-fixed-results-v1", "pid": os.getpid(), "status": "RUNNING",
       "configurations_started": 0, "configurations_completed": 0, "configurations": [],
       "actual_work": {"column_transforms": 0, "scalar_two_cluster_assignments": 0,
                       "three_dimensional_density_paths": 0, "upper_halfplane_residues": 0,
                       "amplitude_point_partials": 0, "amplitude_box_partial_bounds": 0,
                       "remainder_box_partial_bounds": 0, "quotient_matrix_inversions": 0,
                       "dual_sign_norm_evaluations": 0, "outside_label_basis_mappings": 0,
                       "full6D_kernel_calls": 0, "body32_full_kernel_calls": 0,
                       "old_baseline_reruns": 0, "native_calls": 0, "matlab_calls": 0,
                       "recursive_agents": 0, "automatic_retries": 0}}


def save():
    OUT["wall_seconds"] = time.monotonic()-START
    OUT["cpu_seconds"] = time.process_time()-CPU_START
    (ROOT/"results.json").write_text(json.dumps(OUT, indent=2)+"\n")


def begin(identifier):
    OUT["configurations_started"] += 1
    OUT["active_configuration"] = identifier
    save()


def finish(identifier, detail):
    OUT["configurations"].append({"id": identifier, "status": "PASS", **detail})
    OUT["configurations_completed"] += 1
    OUT["active_configuration"] = None
    save()


def indices(order):
    return [(a,b,k-a-b) for k in range(order+1) for a in range(k+1) for b in range(k-a+1)]


def afact(alpha):
    return prod(factorial(a) for a in alpha)


def multiply(p, q, order):
    out = {}
    for a, va in p.items():
        for b, vb in q.items():
            c = tuple(a[i]+b[i] for i in range(3))
            if sum(c) <= order:
                out[c] = out.get(c,F(0))+va*vb
    return {k:v for k,v in out.items() if v}


def inverse_linear(d, gradient, power, order):
    """Exact Taylor coefficients of (d+gradient.h)^(-power)."""
    out = {}
    for alpha in indices(order):
        k = sum(alpha)
        val = F((-1)**k * comb(power+k-1,k)*factorial(k), afact(alpha))
        val *= prod(g**a for g,a in zip(gradient,alpha))*d**(-power-k)
        if val:
            out[alpha] = val
    return out


def add_scaled(out, p, scale):
    for a,v in p.items():
        out[a] = out.get(a,F(0))+scale*v


def derivative_table(poly, order):
    return {a:poly.get(a,F(0))*afact(a) for a in indices(order)}


def tensor_l1_bounds(table, order):
    return [sum(F(factorial(k),afact(a))*abs(table[a]) for a in table if sum(a)==k)
            for k in range(order+1)]


def rows(table):
    return [{"multiindex":list(a),"value":str(v)} for a,v in table.items()]


try:
    reg = json.loads((ROOT/"registration.json").read_text())
    for name, expected in reg["frozen_local_sha256"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected, name
    inp = json.loads((ROOT/"input_manifest.json").read_text())
    assert len(inp["configurations"])==7
    import sympy as s
    OUT["runtime"]={"python":sys.version,"executable":sys.executable,"sympy":s.__version__}
    OUT["registration_sha256"]=hashlib.sha256((ROOT/"registration.json").read_bytes()).hexdigest()
    h,T,rho,R0=map(F,[inp["tangential_box_halfwidth"],inp["time_upper_bound"],inp["local_ball_radius"],inp["scalar_cutoff_radius"]])
    center=list(map(F,inp["tangential_center"]))
    body=[(1,)+v for v in product((1,-1),repeat=5)]
    rot=s.Matrix([[F(1,2)]*4+[0,0],[0,0,0,0,1,0],[0,0,0,0,0,1],
                  [F(1,2),F(1,2),F(-1,2),F(-1,2),0,0],
                  [F(1,2),F(-1,2),F(1,2),F(-1,2),0,0],
                  [F(1,2),F(-1,2),F(-1,2),F(1,2),0,0]])
    begin("Q001-geometry")
    assert rot*rot.T==s.eye(6)
    point=s.Matrix(list(map(F,inp["point"])))
    assert list(rot*point)==list(map(s.Rational,center))+[0,0,0]
    mapped=[]; inside=[]; outside=[]; groups={}
    for index,col in enumerate(body):
        q=tuple(map(F,rot*s.Matrix(col)))
        record={"body_index_zero_based":index,"column":list(col),"tangent":list(map(str,q[:3])),"normal":list(map(str,q[3:]))}
        mapped.append(record)
        if q[3:]==(0,0,0): inside.append(record)
        else:
            outside.append(record)
            groups.setdefault(q[3:],[]).append(record)
        OUT["actual_work"]["column_transforms"]+=1
    assert len(inside)==4 and len(outside)==28 and len(groups)==7
    assert all(len(v)==4 for v in groups.values())
    vin=s.Matrix.hstack(*[s.Matrix(list(map(F,r["tangent"]))) for r in inside])
    assert {tuple(vin[:,i]) for i in range(4)}=={(2,a,b) for a,b in product((1,-1),repeat=2)}
    assert abs(vin[:,:3].det())==8
    cout=s.Matrix.hstack(*[s.Matrix(list(map(F,r["tangent"]))) for r in outside])
    assert cout*cout.T==s.diag(16,28,28)
    values=list(map(F,inp["point"]))
    best=None; minimizers=[]
    for bits in product((0,1),repeat=5):
        labels=(0,)+bits
        loss=F(0)
        for group in (0,1):
            vals=[v for v,label in zip(values,labels) if label==group]
            if vals:
                mean=sum(vals,F(0))/len(vals)
                loss+=sum((v-mean)**2 for v in vals)
        OUT["actual_work"]["scalar_two_cluster_assignments"]+=1
        if best is None or loss<best: best,minimizers=loss,[labels]
        elif loss==best: minimizers.append(labels)
    assert best==F(1,8)
    assert F(1,4*6**3)>F(1,32**2)
    assert 16*rho<=F(1,32) and R0==rho/81 and 4*rho<=h
    assert 32**2*6<=80**2
    affine=[list(map(F,r)) for r in inp["inplane_scalar_centers_affine_coefficients_y0_y1_y2"]]
    centers=[sum(a*y for a,y in zip(row,center)) for row in affine]
    assert centers==[0,F(-1,4),F(1,2),F(1,4)]
    # The eliminated first three variables are s+a, b-s, c-s.
    yy0,yy1,yy2,ww=s.symbols("y0 y1 y2 w",real=True)
    solved=s.Matrix([ww+(yy1+yy2)/2,yy0/4-yy2/2-ww,yy0/4-yy1/2-ww,ww])
    assert s.simplify(vin*solved-s.Matrix([yy0,yy1,yy2]))==s.zeros(3,1)
    finish("Q001-geometry",{
        "rotation_rows":[list(map(str,row)) for row in rot.tolist()],
        "inside_columns":inside,"outside_columns":outside,
        "quotient_groups":[{"normal":list(map(str,k)),"body_indices":[r["body_index_zero_based"] for r in v],"multiplicity":len(v)} for k,v in groups.items()],
        "inplane_basis_absolute_determinant":8,"fiber_jacobian":"1/8",
        "outside_tangent_CCT":[[str(v) for v in row] for row in (cout*cout.T).tolist()],
        "true_S2_distance_squared":"1/8",
        "two_cluster_minimizing_assignments":[list(v) for v in minimizers],
        "nonincident_rank3_distance_squared_lower":"1/864",
        "usable_nonincident_distance_lower":"1/32",
        "rho":str(rho),"R0":str(R0),"unique_incident_U_status":"depends on parent's independent primary_unique_U appendix; not computed by this batch",
        "scalar_centers":list(map(str,centers))
    })

    begin("Q002-inplane-density")
    t=s.symbols("t",real=True)
    residue_terms=[]
    for j in range(4):
        term=s.Integer(1)
        for ell in range(4):
            if ell!=j:
                d=s.Rational(centers[j]-centers[ell])
                term/=d*(d+2*s.I*t)
        residue_terms.append(term/8)
        OUT["actual_work"]["upper_halfplane_residues"]+=1
    density=s.factor(s.cancel(sum(residue_terms)))
    OUT["actual_work"]["three_dimensional_density_paths"]+=1
    assert s.cancel(density-s.conjugate(density))==0
    assert s.cancel(density-density.subs(t,-t))==0
    leading=F(density.subs(t,0))
    correction=s.factor(s.cancel((density-s.Rational(leading))/t**2))
    correction0=F(correction.subs(t,0))
    assert leading==F(2560,9)
    finish("Q002-inplane-density",{"pi3_density_divided_by_t3":str(density),
        "pi3_amplitude_at_center":str(leading),"pi3_t5_correction":str(correction),
        "pi3_t5_correction_at_time_zero":str(correction0),
        "residue_count":4,"normalization":"p_U=t^3/pi^3 times displayed expression; Jacobian1/8 included"})

    begin("Q003-amplitude-C4")
    a_poly={}; upper_poly={}; lower=F(0); differences=[]
    for j in range(4):
        exact={(0,0,0):F(1)}; upper={(0,0,0):F(1)}; lowterm=F(1)
        row=[]
        for ell in range(4):
            if ell==j: continue
            d=centers[j]-centers[ell]
            grad=[affine[j][k]-affine[ell][k] for k in range(3)]
            dmin=abs(d)-h*sum(map(abs,grad))
            dmax=abs(d)+h*sum(map(abs,grad))
            assert dmin>0
            exact=multiply(exact,inverse_linear(d,grad,2,4),4)
            upper=multiply(upper,inverse_linear(dmin,[-abs(g) for g in grad],2,4),4)
            lowterm*=dmax**-2
            row.append({"ell":ell,"center_difference":str(d),"gradient":list(map(str,grad)),"abs_lower":str(dmin),"abs_upper":str(dmax)})
        differences.append(row)
        add_scaled(a_poly,exact,F(1,8)); add_scaled(upper_poly,upper,F(1,8))
        lower+=lowterm/8
    ajets=derivative_table(a_poly,4); abound=derivative_table(upper_poly,4)
    assert len(ajets)==len(abound)==35 and ajets[(0,0,0)]==leading
    assert all(abs(ajets[a])<=abound[a] for a in ajets)
    Aop=tensor_l1_bounds(abound,4)
    h4=[1,0,0,0,0,-126,420,-540,315,-70]
    H=[F(1)]+[sum(F(abs(h4[n])*factorial(n),factorial(n-k)) for n in range(k,10)) for k in range(1,5)]
    H[1]=F(315,128) # exact maximum 630[v(1-v)]^4.
    d1,d2=F(8,7)/rho,F(2,7)/rho**2
    bump=[F(1),H[1]*d1,H[2]*d1**2+H[1]*d2,
          H[3]*d1**3+3*H[2]*d1*d2,
          H[4]*d1**4+6*H[3]*d1**2*d2+3*H[2]*d2**2]
    Aext=[sum(F(comb(k,r))*bump[r]*Aop[k-r] for r in range(k+1)) for k in range(5)]
    OUT["actual_work"]["amplitude_point_partials"]=35
    OUT["actual_work"]["amplitude_box_partial_bounds"]=35
    finish("Q003-amplitude-C4",{"units":"all entries are pi^3 times the indicated amplitude derivative/bound",
        "point_partials":rows(ajets),"uniform_box_partial_bounds":rows(abound),
        "uniform_amplitude_lower":str(lower),"uniform_amplitude_upper":str(Aop[0]),
        "operator_norm_majorants_order0to4":list(map(str,Aop)),
        "separation_data_by_residue":differences,
        "h4_derivative_majorants_order0to4":list(map(str,H)),
        "radial_bump_operator_majorants_order0to4":list(map(str,bump)),
        "extended_amplitude_operator_majorants_order0to4":list(map(str,Aext)),
        "extension_agreement_radius":str(3*rho),"extension_support_radius":str(4*rho)})

    begin("Q004-time-remainder-C2")
    rbound_poly={}
    for j in range(4):
        polys=[]
        for tau in (T,-T,F(0)):
            poly={(0,0,0):F(1)}
            for item in differences[j]:
                dmin=F(item["abs_lower"]); grad=[-abs(F(g)) for g in item["gradient"]]
                assert dmin-2*abs(tau)>0
                poly=multiply(poly,inverse_linear(dmin,grad,1,2),2)
                poly=multiply(poly,inverse_linear(dmin-2*tau,grad,1,2),2)
            polys.append(poly)
        for alpha in indices(2):
            value=(polys[0].get(alpha,0)+polys[1].get(alpha,0)-2*polys[2].get(alpha,0))/(16*T**2)
            assert value>=0
            rbound_poly[alpha]=rbound_poly.get(alpha,F(0))+value
    rbounds=derivative_table(rbound_poly,2); Rop=tensor_l1_bounds(rbounds,2)
    Rext=[sum(F(comb(k,r))*bump[r]*Rop[k-r] for r in range(k+1)) for k in range(3)]
    assert len(rbounds)==10 and abs(correction0)<=rbounds[(0,0,0)]
    OUT["actual_work"]["remainder_box_partial_bounds"]=10
    finish("Q004-time-remainder-C2",{"units":"pi^3 times r_U derivative bound in p_U=t^3 a_U+t^5 r_U",
        "time_interval":["0",str(T)],"uniform_box_partial_bounds":rows(rbounds),
        "operator_norm_majorants_order0to2":list(map(str,Rop)),
        "extended_remainder_operator_majorants_order0to2":list(map(str,Rext)),
        "method":"positive two-variable power-series majorant; even time coefficients and exact +/-T averaging; no time grid"})

    begin("Q005-labelled-quotient-bases")
    fixed_cols=[tuple(v) for v in inp["fixed_quotient_basis_body_columns"]]
    bycol={tuple(r["column"]):r for r in outside}
    fixed_records=[bycol[col] for col in fixed_cols]
    fixed_B=s.Matrix.hstack(*[s.Matrix(list(map(F,r["normal"]))) for r in fixed_records])
    fixed_C=s.Matrix.hstack(*[s.Matrix(list(map(F,r["tangent"]))) for r in fixed_records])
    assert fixed_B==2*s.eye(3)
    assert (fixed_C/2)*(fixed_C/2).T==s.Matrix([[0,0,0],[0,F(3,4),F(-1,4)],[0,F(-1,4),F(3,4)]])
    basis_types={}; label_maps=[]; sum_D=F(0)
    for normal,records in groups.items():
        target=records[0]
        if sum(v!=0 for v in normal)==1:
            axis=next(i for i,v in enumerate(normal) if v)
            extras=[fixed_records[k] for k in range(3) if k!=axis]
        else: extras=[fixed_records[1],fixed_records[2]]
        matrix=s.Matrix.hstack(s.Matrix(normal),*[s.Matrix(list(map(F,r["normal"]))) for r in extras])
        det=abs(F(matrix.det())); inv=matrix.inv()
        OUT["actual_work"]["quotient_matrix_inversions"]+=1
        norms=[]
        for sign in product((1,-1),repeat=3):
            dual=inv.T*s.Matrix(sign)
            norms.append(F((dual.T*dual)[0]))
            OUT["actual_work"]["dual_sign_norm_evaluations"]+=1
        normsq=max(norms)
        expected=F(3,4) if sum(v!=0 for v in normal)==1 else F(9,2)
        assert normsq==expected
        key=",".join(map(str,normal))
        basis_types[key]={"normal":list(map(str,normal)),"absolute_determinant":str(det),
                          "inverse":[list(map(str,row)) for row in inv.tolist()],
                          "operator_2to1_norm_squared":str(normsq),
                          "other_basis_body_indices":[r["body_index_zero_based"] for r in extras]}
        for record in records:
            indices_now=[record["body_index_zero_based"]]+[r["body_index_zero_based"] for r in extras]
            assert len(set(indices_now))==3
            label_maps.append({"weighted_label":record["body_index_zero_based"],"basis_body_indices":indices_now,"basis_type":key,
                               "pi3_W0":str(1/det),"pi3_W1":str(1/(2*det)),"pi3_W2":str(1/det)})
            sum_D+=1/det
            OUT["actual_work"]["outside_label_basis_mappings"]+=1
    assert sum_D==F(11,2) and len(label_maps)==28
    finish("Q005-labelled-quotient-bases",{"basis_types":basis_types,"all_28_label_mappings":label_maps,
        "fixed_differentiation_basis_body_indices":[r["body_index_zero_based"] for r in fixed_records],
        "fixed_quotient_basis":"2I3","fixed_tangent_transport_operator_norm_upper":"1",
        "outside_C_operator_norm_squared":"28",
        "pi3_M0_majorant":"1/8","pi3_M1_sum_majorant":"11/4","pi3_M2_sum_majorant":"11/2",
        "weighted_or_tangent_laws_merged":False})

    begin("Q006-quotient-coercivity")
    expected_groups={(F(2),0,0),(0,F(2),0),(0,0,F(2)),(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)}
    assert set(groups)==expected_groups
    aa,bb,cc,rr=s.symbols("a b c r", nonnegative=True)
    assert s.expand((3*aa+rr)**2-8*(aa**2+rr**2)-(aa-rr)*(aa+7*rr))==0
    assert 8**2*11 < 16**2*3
    assert (16**2*2)==512 and (16**2*3)==768
    Nlower=F(1,786432); Nupper=F(3,1048576)
    assert Nlower==F(3,4*(16**2*3)**2)
    assert Nupper==F(3,4*(16**2*2)**2)
    finish("Q006-quotient-coercivity",{
        "psi_formula":"8*l1+max(16*linfinity,8*l1)",
        "coercivity_exact":"16*sqrt(2)","upper_norm_exact":"16*sqrt(3)",
        "lower_equality_direction":"(1,1,0)/sqrt(2)","upper_equality_direction":"(1,1,1)/sqrt(3)",
        "N_norm_lower_coefficient_pi1":str(Nlower),"N_Lipschitz_upper_coefficient_pi1":str(Nupper),
        "signed_renormalized_layer_Hessian_bound":"3/(2097152*sqrt(2)*pi)",
        "proof_scope":"Only the marginal quotient norm law uses multiplicity4. All signed/tangent fields retain28 labels."})

    begin("Q007-single-cusp-budget")
    cusp_lower=leading*Nlower; cusp_upper=leading*Nupper
    assert cusp_lower==F(5,13824) and cusp_upper==F(5,6144)
    cusp_hessian_pi4=Nupper*(rho*Aop[2]+2*Aop[1])
    cusp_lipschitz_pi4=Nupper*(Aop[0]+rho*Aop[1])
    finish("Q007-single-cusp-budget",{
        "center_amplitude":"2560/(9*pi^3)",
        "center_aN_lower_coefficient_pi4":str(cusp_lower),
        "center_aN_upper_coefficient_pi4":str(cusp_upper),
        "single_negative_cusp_semiconcavity_bound_coefficient_pi4":str(cusp_hessian_pi4),
        "single_negative_cusp_semiconcavity_rational_upper_using_pi_ge_3":str(cusp_hessian_pi4/81),
        "single_negative_cusp_lipschitz_bound_coefficient_pi4":str(cusp_lipschitz_pi4),
        "orientation_only_decimal_semiconcavity":float(cusp_hessian_pi4)/pi**4,
        "orientation_only_decimal_lipschitz":float(cusp_lipschitz_pi4)/pi**4,
        "complete_G_semiconcavity_constant_obtained":False,
        "unique_incident_U_status":"depends on parent's independent primary_unique_U appendix; final source hash/review to be bound in report",
        "unpaid":["parent unique-U appendix final hash and acceptance", "chosen-U cutoff removal and quadratic/signed remainder aggregation", "rank>=4 aggregate contribution owned by086", "large-time full-G tail constant", "matched true scalar vertices"]})
    assert OUT["configurations_completed"]==7
    OUT["status"]="PASS_ALL_SEVEN_FIXED_CONFIGURATIONS"
    OUT["complete_body32_Gamma_obtained"]=False
    OUT["global_g_certificate"]=False
    save()
    print(json.dumps({"status":OUT["status"],"pid":os.getpid(),"wall_seconds":OUT["wall_seconds"],"cpu_seconds":OUT["cpu_seconds"],"actual_work":OUT["actual_work"]}))
except BaseException as exc:
    OUT["status"]="FAILED_FIRST_STOP_NO_AUTOMATIC_RETRY"
    OUT["failure"]={"type":type(exc).__name__,"message":str(exc),"traceback":traceback.format_exc()}
    save()
    print(json.dumps({"status":OUT["status"],"pid":os.getpid(),"failure":str(exc)}))
    raise
