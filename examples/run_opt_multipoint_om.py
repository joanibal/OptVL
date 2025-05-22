"""
Multi-point optimization using OptVL within OpenMDAO.
Minimizes a weighted average of drag from two flight conditions,
subject to achieving a target CL in each.
"""

import numpy as np
import openmdao.api as om

# Attempt to import OVLGroup, use placeholder if not found
try:
    # Assuming OVLGroup is in optvl.openmdao_comp or similar
    from optvl.openmdao_comp import OVLGroup
    print("Successfully imported OVLGroup from optvl.openmdao_comp")
    # A way to get num_twist_variables, assuming OVLGroup has a class/static method or property
    # This is highly dependent on OVLGroup's API.
    # For now, we'll instantiate a temporary one later if needed.
    # NUM_TWIST_VARIABLES = OVLGroup.get_default_num_twist_vars() # Example API
except ImportError:
    print("Failed to import OVLGroup. Using placeholder OVLGroupPlaceholder.")

    class OVLGroupPlaceholder(om.ExplicitComponent):
        """Placeholder for OptVL's OpenMDAO Group/Component."""
        def initialize(self):
            self.options.declare('geom_file', types=str)
            self.options.declare('flight_condition_name', types=str, default='default')
            # Based on typical OVLGroup usage from docs/optimization_setup_om.md
            self.options.declare('num_twist_variables', types=int, default=5) # Default for placeholder

        def setup(self):
            num_twist = self.options['num_twist_variables']

            # Inputs
            self.add_input('Mach', val=0.1)
            self.add_input('Elevator', val=0.0) # Elevator deflection in degrees
            self.add_input('Wing:aincs', val=np.zeros(num_twist)) # Wing incidence angles (twist)

            # Outputs
            self.add_output('CL', val=0.0)
            self.add_output('CD', val=0.01)
            self.add_output('CM', val=0.0)
            # Add other outputs if OVLGroup provides them, e.g., structural metrics

            # For placeholder, declare partials for all inputs against all outputs as finite difference
            self.declare_partials('*', '*', method='fd')
            print(f"Placeholder OVLGroup '{self.pathname}' (flight condition: {self.options['flight_condition_name']}) setup with {num_twist} twist variables.")
        
        def compute(self, inputs, outputs):
            Mach = inputs['Mach']
            Elevator = inputs['Elevator']
            aincs = inputs['Wing:aincs']
            
            # Dummy physics for placeholder
            outputs['CL'] = 0.1 + Mach * 0.5 + np.sum(aincs) * 0.01 + Elevator * 0.05
            outputs['CD'] = 0.01 + outputs['CL']**2 / (np.pi * 8) + np.sum(np.abs(aincs)) * 0.001 + np.abs(Elevator) * 0.002
            outputs['CM'] = -0.05 + Elevator * -0.02 # Elevator for trim
            
            # print(f"Placeholder OVLGroup '{self.pathname}' compute: Mach={Mach:.2f}, Elev={Elevator:.2f}, CL={outputs['CL']:.3f}, CD={outputs['CD']:.4f}")

    OVLGroup = OVLGroupPlaceholder # Use the placeholder

# --- Constants for the Optimization Problem ---
GEOM_FILE = "aircraft.avl" # Make sure this file exists or is handled by OVLGroup

# Condition 1: Cruise
MACH_CRUISE = 0.8
TARGET_CL_CRUISE = 0.5
WEIGHT_CRUISE = 0.7 # Weight for this condition in the objective

# Condition 2: Loiter
MACH_LOITER = 0.4
TARGET_CL_LOITER = 0.7
WEIGHT_LOITER = 0.3 # Weight for this condition in the objective

# Dynamically determine num_twist_variables by instantiating a temporary OVLGroup
# This assumes OVLGroup either has default num_twist_variables or can deduce from geom_file
try:
    # If the real OVLGroup is used, it might need geom_file to determine this
    temp_ovl_group = OVLGroup(geom_file=GEOM_FILE)
    # Accessing options if it's an ExplicitComponent or Group with options
    # This depends on OVLGroup's specific API
    if hasattr(temp_ovl_group, 'options') and 'num_twist_variables' in temp_ovl_group.options:
         NUM_TWIST_VARIABLES = temp_ovl_group.options['num_twist_variables']
    elif hasattr(OVLGroup, 'DEFAULT_NUM_TWIST_VARS'): # Hypothetical class attribute
         NUM_TWIST_VARIABLES = OVLGroup.DEFAULT_NUM_TWIST_VARS
    else: # Fallback if not easily determined
        print("Warning: Could not dynamically determine NUM_TWIST_VARIABLES. Using default of 5.")
        NUM_TWIST_VARIABLES = 5
except Exception as e:
    print(f"Error instantiating temporary OVLGroup to find NUM_TWIST_VARIABLES: {e}. Using default of 5.")
    NUM_TWIST_VARIABLES = 5

print(f"Using NUM_TWIST_VARIABLES = {NUM_TWIST_VARIABLES}")


class WeightedDragObjective(om.ExplicitComponent):
    """
    Calculates the weighted average of drag coefficients from two flight conditions.
    """
    def setup(self):
        self.add_input('CD_cruise', val=0.0)
        self.add_input('CD_loiter', val=0.0)
        self.add_output('avg_CD', val=0.0)

        self.declare_partials('avg_CD', 'CD_cruise', val=WEIGHT_CRUISE)
        self.declare_partials('avg_CD', 'CD_loiter', val=WEIGHT_LOITER)

    def compute(self, inputs, outputs):
        outputs['avg_CD'] = WEIGHT_CRUISE * inputs['CD_cruise'] + WEIGHT_LOITER * inputs['CD_loiter']

    # compute_partials can be omitted if vals are provided in declare_partials and don't change.

class MultiPointOptimizationGroup(om.Group):
    """
    Main OpenMDAO Group for multi-point optimization.
    """
    def setup(self):
        # Instantiate OVLGroup for each flight condition
        cruise_point = self.add_subsystem('cruise_point',
                                          OVLGroup(geom_file=GEOM_FILE,
                                                   flight_condition_name='cruise',
                                                   num_twist_variables=NUM_TWIST_VARIABLES),
                                          promotes_inputs=['Wing:aincs']) # Promote shared DVs

        loiter_point = self.add_subsystem('loiter_point',
                                          OVLGroup(geom_file=GEOM_FILE,
                                                   flight_condition_name='loiter',
                                                   num_twist_variables=NUM_TWIST_VARIABLES),
                                          promotes_inputs=['Wing:aincs']) # Promote shared DVs

        # Independent variables for point-specific settings (Mach, Elevator)
        # These will be set by prob.set_val before running
        point_ivc = om.IndepVarComp()
        point_ivc.add_output('Mach_cruise', val=MACH_CRUISE)
        point_ivc.add_output('Elevator_cruise', val=0.0) # Initial guess for elevator
        point_ivc.add_output('Mach_loiter', val=MACH_LOITER)
        point_ivc.add_output('Elevator_loiter', val=0.0) # Initial guess for elevator
        self.add_subsystem('point_ivc', point_ivc, promotes_outputs=['*'])

        # Independent variables for shared design variables (wing twist)
        # This component's output `wing_aincs` will be connected to the promoted `Wing:aincs`
        shared_ivc = om.IndepVarComp()
        shared_ivc.add_output('wing_aincs_source', val=np.zeros(NUM_TWIST_VARIABLES))
        self.add_subsystem('shared_ivc', shared_ivc, promotes_outputs=['wing_aincs_source'])


        # Objective component
        objective_comp = self.add_subsystem('objective_comp', WeightedDragObjective())

        # Connections
        # Shared DVs
        self.connect('wing_aincs_source', 'Wing:aincs') # Connect source to the promoted input

        # Point-specific DVs to OVLGroups
        self.connect('Mach_cruise', 'cruise_point.Mach')
        self.connect('Elevator_cruise', 'cruise_point.Elevator')
        self.connect('Mach_loiter', 'loiter_point.Mach')
        self.connect('Elevator_loiter', 'loiter_point.Elevator')

        # OVLGroup outputs to Objective component
        self.connect('cruise_point.CD', 'objective_comp.CD_cruise')
        self.connect('loiter_point.CD', 'objective_comp.CD_loiter')


# --- Main Execution Block ---
if __name__ == "__main__":
    prob = om.Problem(model=MultiPointOptimizationGroup())

    # Setup Driver
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['tol'] = 1e-6
    prob.driver.options['disp'] = True # Display optimizer progress
    # prob.driver.options['debug_print'] = ['desvars', 'nl_cons', 'objs']


    # Add Design Variables
    # Shared design variable: wing twist (aincs)
    # The source of this DV is 'shared_ivc.wing_aincs_source'
    prob.model.add_design_var('wing_aincs_source', lower=-5.0, upper=5.0, scaler=0.1)

    # Point-specific design variables: elevator deflections
    # The source of these DVs are 'point_ivc.Elevator_cruise' and 'point_ivc.Elevator_loiter'
    prob.model.add_design_var('Elevator_cruise', lower=-15.0, upper=15.0, scaler=0.1)
    prob.model.add_design_var('Elevator_loiter', lower=-15.0, upper=15.0, scaler=0.1)

    # Add Objective
    prob.model.add_objective('objective_comp.avg_CD', scaler=10.0) # Scaler can help optimizer

    # Add Constraints
    prob.model.add_constraint('cruise_point.CL', equals=TARGET_CL_CRUISE, scaler=1.0)
    prob.model.add_constraint('loiter_point.CL', equals=TARGET_CL_LOITER, scaler=1.0)
    # Optional: Add trim constraints if Elevator is not a DV or for CM !=0 targets
    # prob.model.add_constraint('cruise_point.CM', equals=0.0)
    # prob.model.add_constraint('loiter_point.CM', equals=0.0)

    # Setup Problem
    # For OptVL, 'rev' mode is often required or beneficial for efficiency
    try:
        prob.setup(check=True, mode='rev')
        print("OpenMDAO problem setup successful.")
    except Exception as e:
        print(f"Error during OpenMDAO problem setup: {e}")
        exit()

    # Set initial values for IndepVarComps that are NOT design variables (Mach numbers)
    # Design variable initial values are handled by their 'val' in IndepVarComp
    # or can be set here too if needed, but driver usually starts from IVC 'val'.
    prob.set_val('Mach_cruise', MACH_CRUISE) # From point_ivc
    prob.set_val('Mach_loiter', MACH_LOITER) # From point_ivc
    
    # Set initial values for design variables (optional if 'val' in IVC is desired start)
    prob.set_val('wing_aincs_source', np.full(NUM_TWIST_VARIABLES, 0.0)) # Start with no twist
    prob.set_val('Elevator_cruise', 0.0) # Start with zero elevator
    prob.set_val('Elevator_loiter', 0.0) # Start with zero elevator


    # Run Optimization
    print("\nStarting OpenMDAO optimization driver...")
    try:
        prob.run_driver()
        print("OpenMDAO optimization driver completed.")
    except Exception as e:
        print(f"Error during OpenMDAO run_driver: {e}")
        # prob.cleanup() # Might be useful if OptVL instances need explicit closing
        exit()

    # Print Results
    print("\n--- Optimization Results ---")
    if prob.driver.result.success:
        print("Optimization Successful!")
    else:
        print("Optimization Failed.")
        print(f"  Message: {prob.driver.result.message}")


    optimal_wing_aincs = prob.get_val('wing_aincs_source')
    optimal_elevator_cruise = prob.get_val('Elevator_cruise')
    optimal_elevator_loiter = prob.get_val('Elevator_loiter')
    final_avg_CD = prob.get_val('objective_comp.avg_CD')

    print(f"\nOptimal Shared Wing Twist (aincs): {optimal_wing_aincs}")
    print(f"Optimal Cruise Elevator: {optimal_elevator_cruise:.4f} deg")
    print(f"Optimal Loiter Elevator: {optimal_elevator_loiter:.4f} deg")
    print(f"\nFinal Weighted Average CD: {final_avg_CD:.6f}")

    print("\n--- Final Values for Cruise Point ---")
    print(f"  Mach: {prob.get_val('cruise_point.Mach'):.2f}")
    print(f"  CL:   {prob.get_val('cruise_point.CL'):.4f} (Target: {TARGET_CL_CRUISE})")
    print(f"  CD:   {prob.get_val('cruise_point.CD'):.5f}")
    print(f"  CM:   {prob.get_val('cruise_point.CM'):.4f}")
    print(f"  Elev: {prob.get_val('cruise_point.Elevator'):.4f} deg")


    print("\n--- Final Values for Loiter Point ---")
    print(f"  Mach: {prob.get_val('loiter_point.Mach'):.2f}")
    print(f"  CL:   {prob.get_val('loiter_point.CL'):.4f} (Target: {TARGET_CL_LOITER})")
    print(f"  CD:   {prob.get_val('loiter_point.CD'):.5f}")
    print(f"  CM:   {prob.get_val('loiter_point.CM'):.4f}")
    print(f"  Elev: {prob.get_val('loiter_point.Elevator'):.4f} deg")
    
    # prob.list_outputs(residuals=True, prom_name=True)
    # prob.list_inputs(prom_name=True)

    # Cleanup, if necessary (e.g., if OVLGroup opens files or external processes)
    # prob.cleanup()
    print("\nScript finished.")
```
