leansdr uses C++ for namespaces and type-safe polymorphism.
No attempt is made to follow popular object-oriented practices.

* Member variables are not prefixed with "m_".

* Destructors are not implemented and memory management is minimal.
  In practice, after the signal processing flow graph is instantiated,
  no allocation/deallocation is expected until exit.

* There are no unnecessary getter/setter methods.

* Dependencies are kept to a minimum (no STL, no iostream).

Other notes:

* The code is not thread-safe.

