# Installation

dlufd builds against the stream component (gstrickd headers and stage ids). CMake uses, in this order, the stream
target of the combined obdeect build, the stream sources in `../stream` of the repository, or an
installed gstrickd package (`gstrickd_DIR`, then `python -m gstrickd --cmake-dir`). A source
distribution does not contain the stream sources, so it needs the installed package.
