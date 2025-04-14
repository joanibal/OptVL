# Development Docs
to add a variable to avl shared memory (common block)

1. add it to the COMMON line AND declare its type above
2. navigate to `src/includes` and run `gen_ad_inc.py` 

## releasing a new version
1. Bump the version number in pyproject.toml AND in meson.build.
    - I cannot figure out a way to single source this so it will have to be done in two steps for now.
2. Then create a new release on github. Create a new taged version as part of the release
3. Check that that new version has been released on PYPI as it should