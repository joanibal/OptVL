set -xe

PROJECT_DIR="$1"

cd $PROJECT_DIR/tests

# install tesing dependencies
pip install psutil openmdao!=3.38
# pip install testflo

# testflo -v .

#HACK: if the tests are not split up the CI runs out of memory...
# python -m unittest -v




# test package built and installed correctly
python test_import.py
# python -m unittest -v test_import.py
# python -m unittest -v test_io.py

# # test mem ussage of pyavl and test framework
# python -m unittest -v test_tear_down.py

# # test basic avl functionality
# python -m unittest -v test_parameters.py
# python -m unittest -v test_analysis.py
# python -m unittest -v test_surf_geom.py
# python -m unittest -v test_constraints.py
# python -m unittest -v test_stab_derivs.py
# python -m unittest -v test_body_axis_derivs.py
# # test eigenmode analysis
# python -m unittest -v test_eigen_analysis.py

# # tests for adjoint
# python -m unittest -v test_new_subroutines.py
# #HACK: if the tests are not split up the windows version has
# #      an error loading the shared lib.
# python -m unittest -v test_partial_derivs.TestFunctionPartials
# python -m unittest -v test_partial_derivs.TestResidualPartials
# python -m unittest -v test_consurf_partial_derivs.TestResidualDPartials
# python -m unittest -v test_consurf_partial_derivs.TestConSurfDerivsPartials
# python -m unittest -v test_stab_derivs_partial_derivs.TestResidualUPartials
# python -m unittest -v test_stab_derivs_partial_derivs.TestStabDerivDerivsPartials

# python -m unittest -v test_body_axis_derivs_partial_derivs.py
# python -m unittest -v test_total_derivs.py

# # test openmdao wrapper and basic optimization results
# python -m unittest -v test_om_wrapper.py

