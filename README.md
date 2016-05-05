ra4_draw
========

Repository for plotting and table-making utilities combining the flexibility of analysis_code's simple and flexible string-based plotting with the speed of ra4_macros' "looping" scripts.

#### Code documentation
Doxygen-based documentation is available at [Adam's UCSB-HEP webpage](http://hep.ucsb.edu/people/ald77/documentation/doc_ra4_draw/).

#### Setup and compilation
Compilation requires c++11 and ROOT, but not CMSSW. To compile, simply run

    ./compile.sh

#### Making histograms
An example script is available under src/test.cxx. To execute, compile and then run

    ./run/test.exe
