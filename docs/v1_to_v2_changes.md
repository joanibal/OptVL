# Migration Guide: OptVL v1 to v2

This guide outlines the breaking changes introduced in OptVL v2 and provides instructions for updating your code.
OptVL v2 includes improvements to AVL consistency, API clarity, and additional functionality. 

## Summary of Changes

1. **Renamed force and derivative coefficients** for consistency with AVL and to distinguish between body axis and stability axis
2. **Introduced body axis derivatives** via a new method `get_body_axis_derivs()`
3. **Updated API** for setting variables, control deflections, and constraints
4. **Upgraded to AVL v3.52 source code** with updated physics models
<!-- 5. **Improved testing infrastructure** for generating reference data directly from AVL -->

## Force Coefficient Changes

### Changes to `get_total_forces()`

The moment coefficient naming has been updated to be more consistent with AVL conventions and to clearly distinguish between body axis and stability axis frames.

| Old Key | New Key | Description                         |
| :-----: | :-----: | :---------------------------------: |
| CR BA   | Cl      | roll moment coef. (body frame)      |
| CM      | Cm      | pitch moment coef. (body frame)     |
| CN BA   | Cn      | yaw moment coef. (body frame)       |
| CR SA   | Cl'     | roll moment coef. (stability frame) |
| CN SA   | Cn'     | yaw moment coef. (stability frame)  |


**example:**

old code
```python
forces = ovl.get_total_forces()
roll_moment_body = forces["CR BA"]
pitch_moment = forces["CM"]
yaw_moment_body = forces["CN BA"]
```
new code
```python
# v2 code
forces = ovl.get_total_forces()
roll_moment_body = forces["Cl"]
pitch_moment = forces["Cm"]
yaw_moment_body = forces["Cn"]
```

## Stability Derivative Changes

### Changes to `get_stab_derivs()`

The stability derivatives now use the AVL notation for rotations around the stability axis (p', q', r'). 

**v1 stability derivative keys:**

|     | alpha      | beta      | roll rate      | pitch rate      | yaw rate      |
| --- | ---------- | --------- | -------------- | --------------- | ------------- |
| CL  | dCL/dalpha | dCL/dbeta | dCL/droll rate | dCL/dpitch rate | dCL/dyaw rate |
| CD  | dCD/dalpha | dCD/dbeta | dCD/droll rate | dCD/dpitch rate | dCD/dyaw rate |
| CY  | dCY/dalpha | dCY/dbeta | dCY/droll rate | dCY/dpitch rate | dCY/dyaw rate |
| CR  | dCR/dalpha | dCR/dbeta | dCR/droll rate | dCR/dpitch rate | dCR/dyaw rate |
| CM  | dCM/dalpha | dCM/dbeta | dCM/droll rate | dCM/dpitch rate | dCM/dyaw rate |
| CN  | dCN/dalpha | dCN/dbeta | dCN/droll rate | dCN/dpitch rate | dCN/dyaw rate |

**v2 stability derivative keys:**

|     | alpha       | beta       | p'       | q'       | r'       |
| --- | ----------- | ---------- | -------- | -------- | -------- |
| CL  | dCL/dalpha  | dCL/dbeta  | dCL/dp'  | dCL/dq'  | dCL/dr'  |
| CD  | dCD/dalpha  | dCD/dbeta  | dCD/dp'  | dCD/dq'  | dCD/dr'  |
| CY  | dCY/dalpha  | dCY/dbeta  | dCY/dp'  | dCY/dq'  | dCY/dr'  |
| Cl' | dCl'/dalpha | dCl'/dbeta | dCl'/dp' | dCl'/dq' | dCl'/dr' |
| Cm  | dCm/dalpha  | dCm/dbeta  | dCm/dp'  | dCm/dq'  | dCm/dr'  |
| Cn' | dCn'/dalpha | dCn'/dbeta | dCn'/dp' | dCn'/dq' | dCn'/dr' |

**Migration example:**
```python
# v1 code
derivs = ovl.get_stab_derivs()
roll_damping = derivs["CR/roll rate"]

# v2 code
derivs = ovl.get_stab_derivs()
roll_damping = derivs["Cl'/p'"]
```

## New Body Axis Derivatives

Version 2 introduces a new method `get_body_axis_derivs()` that provides derivatives with respect to body axis velocities and angular rates.

```python
body_derivs = ovl.get_body_axis_derivs()
```

This method returns derivatives in the following format:

|     | p      | q      | r      | u      | v      | w      |
| --- | ------ | ------ | ------ | ------ | ------ | ------ |
| CX  | dCX/dp | dCX/dq | dCX/dr | dCX/du | dCX/dv | dCX/dw |
| CY  | dCY/dp | dCY/dq | dCY/dr | dCY/du | dCY/dv | dCY/dw |
| CZ  | dCZ/dp | dCZ/dq | dCZ/dr | dCZ/du | dCZ/dv | dCZ/dw |
| Cl  | dCl/dp | dCl/dq | dCl/dr | dCl/du | dCl/dv | dCl/dw |
| Cm  | dCm/dp | dCm/dq | dCm/dr | dCm/du | dCm/dv | dCm/dw |
| Cn  | dCn/dp | dCn/dq | dCn/dr | dCn/du | dCn/dv | dCn/dw |

## API Changes for Variables and Constraints

The API for setting variables, control surface deflections, and constraints has been updated for clarity and consistency.

### Setting Flight Condition Variables

Use `set_variable()` to set angle of attack, sideslip angle, and angular rates:

```python
# v2 syntax
ovl.set_variable("alpha", 5.0)
ovl.set_variable("beta", 0.0)
ovl.set_variable("roll rate", 0.0)
ovl.set_variable("pitch rate", 0.0)
ovl.set_variable("yaw rate", 0.0)
```

### Setting Control Surface Deflections

Use the dedicated `set_control_deflection()` method for control surfaces:

```python
# v1 code - control surfaces mixed with flight variables
ovl.set_variable("Elevator", -2.0)

# v2 code - dedicated method for control surfaces
ovl.set_control_deflection("Elevator", -2.0)
ovl.set_control_deflection("Aileron", 5.0)
ovl.set_control_deflection("Rudder", 0.0)
```

### Setting Constraints

The `set_constraint()` method now requires all three parameters explicitly, the variable to change, the constraint variable and constraint value.

```python
# v2 syntax - constraint variable is required
ovl.set_constraint("Elevator", "Cm", 0.0)  # Trim elevator to achieve Cm = 0
ovl.set_constraint("alpha", "CL", 0.5)     # Adjust alpha to achieve CL = 0.5
```

**API changes example:**
```python
# v1 code
ovl.set_constraint("alpha", 5.0)
ovl.set_constraint("Elevator", 0.00)
ovl.set_constraint("alpha", 0.5, con_var="CL")  

# v2 code
ovl.set_variable("alpha", 5.0)
ovl.set_control_deflection("Elevator", -2.0)
ovl.set_constraint("alpha", "CL", 0.5)  # All parameters explicitly required
```

## AVL v3.52 Source Code Updates

OptVL v2 is based on AVL version 3.52, which multiple changes including:

- Updated physics models with different default vortex core radius calculations
- New Core geometry file input to explicitly adjust the vortex cores of segments.

!!! warning "Slight Numerical Differences"
    Due to the updated physics models in AVL v3.52, you may observe small differences in analysis results compared to v1. These differences are generally minor but reflect improvements in the underlying AVL code.

## Quick Update Checklist

When migrating from v1 to v2, update your code as follows:

- Replace `"CR BA"` with `"Cl"` in `get_total_forces()` calls
- Replace `"CM"` with `"Cm"` in `get_total_forces()` calls
- Replace `"CN BA"` with `"Cn"` in `get_total_forces()` calls
- Replace `"CR SA"` with `"Cl'"` in `get_total_forces()` calls
- Replace `"CN SA"` with `"Cn'"` in `get_total_forces()` calls
- Update stability derivative keys: `"roll rate"` → `"p'"`, `"pitch rate"` → `"q'"`, `"yaw rate"` → `"r'"`
- Replace moment coefficient keys in `get_stab_derivs()`: `"CR"` → `"Cl'"`, `"CM"` → `"Cm"`, `"CN"` → `"Cn'"`
- Use `set_control_deflection()` instead of `set_constraint()` for control surfaces
- Use `set_variable()` instead of `set_constraint()` for flight condition variables
- Ensure all three parameters are provided to `set_constraint()`
- Test your code and verify results are as expected (minor numerical differences may occur)
