from optvl import OVLSolver
import numpy as np


def finite_diff(ovl, var, step=1e-6):
    forces = ovl.get_total_forces()
    val = ovl.get_variable(var)
    ovl.set_variable(var, val + step)
    ovl.execute_run()
    dforces = ovl.get_total_forces()

    for key in dforces:
        dforces[key] = (dforces[key] - forces[key]) / step

    return dforces


for geom in ["../geom_files/rectangle.avl", "../geom_files/aircraft.avl", "../geom_files/supra.avl"]:
    ovl = OVLSolver(geo_file=geom, debug=False)

    alpha = 10.0
    ovl.set_variable('alpha', alpha)
    ovl.execute_run()

    forces = ovl.get_total_forces()
    R = np.array([[np.sin(np.pi/180*alpha), 0, -1*np.cos(np.pi/180*alpha)]])

    CF = np.array([[forces['CX'], forces['CY'], forces['CZ']]]).T

    # Compute dCL/dalpha using chain rule:
    # CL(alpha) = R(alpha) @ CF(alpha)
    # dCL/dalpha = dR/dalpha @ CF + R @ dCF/dalpha

    dR_dalpha = np.array([[np.cos(np.pi/180*alpha)*np.pi/180, 0, np.sin(np.pi/180*alpha)*np.pi/180]])

    sens = ovl.execute_run_sensitivities(['CL', 'CX', 'CY', 'CZ'])
    dCF_dalpha = np.array([[sens['CX']['alpha'], sens['CY']['alpha'], sens['CZ']['alpha']]]).T

    body_derivs = ovl.get_body_axis_derivs()
    dCF_dV = np.zeros((3, 6))
    for idx_vel, vel_comp in enumerate(["u", "v", "w", "p", "q", "r"]):
        for idx_f, f_comp in enumerate(["X", "Y", "Z"]):
            dCF_dV[idx_f, idx_vel] = body_derivs[f"dC{f_comp}/d{vel_comp}"]

    dV_dalpha = np.array([[-np.sin(np.pi/180*alpha)*np.pi/180, 0, np.cos(np.pi/180*alpha)*np.pi/180, 0, 0, 0]]).T

    stab_derivs = ovl.get_stab_derivs()
    dforces = finite_diff(ovl, 'alpha', step=1e-8)
    dCL_dalpha_avl = stab_derivs['dCL/dalpha'] * np.pi/180  # stab derivs use radians, others use degrees
    dCL_dalpha_p_ad = (dR_dalpha @ CF + R @ dCF_dalpha)[0, 0]
    dCL_dalpha_analytic = (dR_dalpha @ CF + R @ dCF_dV @ dV_dalpha)[0, 0]

    print(f'--- comparison for {geom} ----')
    print('AVL    ', dCL_dalpha_avl)
    print('full AD', sens['CL']['alpha'], (sens['CL']['alpha'] - dCL_dalpha_avl) / dCL_dalpha_avl)
    print('part AD', dCL_dalpha_p_ad, (dCL_dalpha_p_ad - dCL_dalpha_avl) / dCL_dalpha_avl)
    print('hand   ', dCL_dalpha_analytic, (dCL_dalpha_analytic - dCL_dalpha_avl) / dCL_dalpha_avl)
    print('FD     ', dforces['CL'], (dforces['CL'] - dCL_dalpha_avl) / dCL_dalpha_avl)