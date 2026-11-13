# Installation

obdeect propagation builds against the stream component (obdeect stream headers and stage ids). CMake uses, in this order, the stream
target of the combined obdeect build, the stream sources in `../stream` of the repository, or an
installed obdeect stream package (`obdeect_stream_DIR`, then `python -m obdeect.stream --cmake-dir`). A source
distribution does not contain the stream sources, so it needs the installed package.
