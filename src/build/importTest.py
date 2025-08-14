#! /usr/bin/env python
import sys

name = "libavl"
print("Testing if module %s can be imported..." % name)
import_cmd = "import %s" % name
exec(import_cmd)

# try:
#     exec(import_cmd)
# except ImportError:
#     print("Error: libadflow was not imported correctly")
#     sys.exit(1)
# end try

print("\033[1;32mModule %s was passed import test.\033[0m" % name)
