# Installation

dlufd builds against gstrickd's photon-stream headers and stage ids. CMake first looks for an
installed gstrickd package (`find_package(gstrickd 0.1 CONFIG)`); if none is found and
`gstrickd_DIR` is not set, it asks the Python that runs CMake for `python -m gstrickd --cmake-dir`
(pip install gstrickd). If that also fails, it fetches gstrickd at the pinned commit with
FetchContent and reports which source it uses. Setting `-DDLUFD_REQUIRE_INSTALLED_GSTRICKD=ON`
turns the fetch off, so the build needs an installed package (for offline or sdist builds).
