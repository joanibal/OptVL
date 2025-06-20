import numpy as np
from optvl.optvl_class import OVLSolver # Adjust import if OVLSolver is elsewhere

def create_simple_wing_geom_data():
    """Defines a simple rectangular wing geometry."""
    aircraft_data = {
        "name": "ProgrammaticWing",
        "reference_conditions": {
            "mach": 0.1,
            "Sref": 10.0,
            "Cref": 1.0,
            "Bref": 10.0,
            "Xref": 0.25,
            "Yref": 0.0,
            "Zref": 0.0,
            "CD0": 0.01
        },
        "symmetry": {
            "iy_sym": 1, # Symmetric about Y=0 (simulating a half-model)
                         # Or use iy_sym = 0 and define full wing for simplicity initially
        },
        "surfaces": [
            {
                "name": "MainWing",
                "n_chordwise": 10,
                "chord_spacing": "cosine", # This key needs to be mapped to CSPACE in set_geometry_from_data
                "surface_spanwise_spacing": True, # This key needs to be mapped to LSURFSPACING
                "n_spanwise": 16, # Panels on the half-wing if iy_sym=1, or total if iy_sym=0
                "span_spacing": "cosine", # This key needs to be mapped to SSPACE in set_geometry_from_data
                # "y_duplicate_surface": {"y_dup": 0.0}, # Not needed if iy_sym=1 handles it, or if defining full wing.
                                                       # If iy_sym=1, AVL handles duplication.
                                                       # If iy_sym=0, and we want a YDUP, then this is used.
                                                       # For this example, let iy_sym = 1 and omit y_duplicate_surface.
                "angle": 0.0, # Overall incidence for the wing
                "sections": [
                    { # Root section
                        "le_point": [0.0, 0.0, 0.0],
                        "chord": 1.0,
                        "incidence_angle": 2.0, # degrees
                        # n_spanwise_panels and span_spacing_parameter not needed if surface_spanwise_spacing is True
                        "airfoil": {"naca": "0012"},
                    },
                    { # Tip section
                        "le_point": [0.0, 5.0, 0.0], # Half-span is 5.0
                        "chord": 1.0,
                        "incidence_angle": 2.0, # degrees
                        "airfoil": {"naca": "0012"},
                    }
                ]
            }
        ],
        "global_controls": [] # No controls for this simple example
    }
    return aircraft_data

if __name__ == "__main__":
    print("Creating programmatic geometry data...")
    geom_data_input = create_simple_wing_geom_data()

    # Translate the example data to the format expected by set_geometry_from_data
    # This is a temporary step, set_geometry_from_data should ideally handle this structure.
    # For now, we adapt the input to what set_geometry_from_data was built for.

    translated_geom_data = {
        "name": geom_data_input["name"],
        "Sref": geom_data_input["reference_conditions"]["Sref"],
        "Cref": geom_data_input["reference_conditions"]["Cref"],
        "Bref": geom_data_input["reference_conditions"]["Bref"],
        "Xref": geom_data_input["reference_conditions"]["Xref"],
        "Yref": geom_data_input["reference_conditions"]["Yref"],
        "Zref": geom_data_input["reference_conditions"]["Zref"],
        "Mach": geom_data_input["reference_conditions"]["mach"],
        "CD0": geom_data_input["reference_conditions"]["CD0"],
        "IYsym": geom_data_input["symmetry"]["iy_sym"],
        "IZsym": 0, # Assuming not symmetric in XZ plane for this example
        "Zsym": 0.0, # Assuming symmetry plane is at Z=0 if IZsym=1
        "surfaces": []
    }

    for surf_in in geom_data_input["surfaces"]:
        surf_out = {
            "name": surf_in["name"],
            "Nchordwise": surf_in["n_chordwise"],
            # Cspace: 1.0 for uniform, -1.0 for cosine (example values)
            "Cspace": -1.0 if surf_in["chord_spacing"] == "cosine" else 1.0,
            "Lsurfspacing": surf_in["surface_spanwise_spacing"],
            "Nspanwise": surf_in["n_spanwise"], # Corresponds to NVS if Lsurfspacing is true
             # Sspace: 1.0 for uniform, -1.0 for cosine
            "Sspace": -1.0 if surf_in.get("span_spacing", "uniform") == "cosine" else 1.0,
            "angle": surf_in["angle"],
            "sections": [],
            "controls": [] # Assuming no controls per surface from input for now
        }
        # If iy_sym = 1, YDUPLICATE is not set at the surface level in AVL input for the primary surface.
        # AVL handles the duplication. If we were defining a full wing (iy_sym=0) and wanted
        # one surface to be a duplicate of another, then YDUPLICATE would be used.
        # For iy_sym=1, the surface defined is the main one, and AVL mirrors it.

        for sec_in in surf_in["sections"]:
            sec_out = {
                "Xle": sec_in["le_point"][0],
                "Yle": sec_in["le_point"][1],
                "Zle": sec_in["le_point"][2],
                "chord": sec_in["chord"],
                "ainc": sec_in["incidence_angle"],
                # Nspan and Sspace for sections are only used if Lsurfspacing is False
                "Nspan": surf_in.get("n_spanwise_section", 0), # Example, assuming 0 if not defined
                "Sspace": -1.0 if surf_in.get("span_spacing_section", "uniform") == "cosine" else 1.0,
                "afile": f"NACA {sec_in['airfoil']['naca']}" if "naca" in sec_in["airfoil"] else sec_in["airfoil"].get("datfile", ""),
                # Optional CLaf, CDCL can be added here if needed
            }
            surf_out["sections"].append(sec_out)
        translated_geom_data["surfaces"].append(surf_out)

    # Add global controls if any (translating from geom_data_input["global_controls"])
    # For this example, it's empty. If it had entries, they would need to be mapped to
    # the control definition structure expected by set_geometry_from_data,
    # particularly for DNAME setup.

    print("Initializing OVLSolver without a geometry file...")
    try:
        solver = OVLSolver(geo_file=None, debug=True) # Enable debug for more verbose output
    except Exception as e:
        print(f"Error initializing OVLSolver: {e}")
        print("Make sure the OptVL library (including Fortran components) is compiled and accessible.")
        import traceback
        traceback.print_exc()
        exit()

    print("Setting geometry from data...")
    try:
        solver.set_geometry_from_data(translated_geom_data)
        print("Geometry set successfully.")
    except Exception as e:
        print(f"Error in set_geometry_from_data: {e}")
        import traceback
        traceback.print_exc()
        exit()

    # Set up a simple run case
    print("Setting up analysis case...")
    try:
        solver.set_constraint("alpha", 5.0) # Alpha = 5 degrees

        print("Executing run...")
        solver.execute_run()
        print("Run executed.")

        print("Fetching total forces...")
        total_forces = solver.get_total_forces()
        print("Total Forces:")
        if total_forces:
            for k, v in total_forces.items():
                print(f"  {k}: {v:.4f}")
        else:
            print("  No force data retrieved.")

        output_avl_file = "programmatic_wing.avl"
        print(f"Writing generated geometry to {output_avl_file}...")
        solver.write_geom_file(output_avl_file)
        print(f"Geometry written to {output_avl_file}.")

    except Exception as e:
        print(f"An error occurred during analysis: {e}")
        import traceback
        traceback.print_exc()

    print("Example finished.")
