
# The Code
This code is an OO version of the simple Ising Spin model for magnetic domains.
It's spread across files and uses some OO features.
All of the code files are in a ".src" directory. Module and .o files are also placed here.
It can be built using "./build IsingSpin.f90" and should run correctly. The final
executable is called IsingSpin and should be in the root directory (not src)


# The Docs
The documentation is INCOMPLETE to avoid too much unneeded text in the docs, but
demonstrates all of the tricks discussed in this repo.

## Notes on the Fortran Implementation

* Functions on types use C++ style namespacing with "::" i.e. "typename::function"
* Doxygen doesn't extract from _within_ the Program block, but does extract from within functions
* Doxygen doesn't understand Fortran OO, so I have used a custom command to document the Type-bound functions


