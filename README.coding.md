# Versioning

Recommended commits are tagged "x.y.z".

* New z = Bugs fixed.
* New y = Features added.
* New x = Backward compatibility is not guaranteed.

Each tagged commit normally passes QA tests.  Users are
encouraged to always test the latest version.

# Coding style

leansdr uses C++ for namespaces and type-safe polymorphism.
No attempt is made to follow popular object-oriented practices.

* Member variables are not prefixed with "m_".

* Destructors are not implemented and memory management is minimal.
  In practice, after the signal processing flow graph is instantiated,
  no allocation/deallocation is expected until exit.

* There are no unnecessary getter/setter methods.

* Dependencies are kept to a minimum (no STL, no iostream).

# Known limitations

* The code is not intended to be thread-safe.
