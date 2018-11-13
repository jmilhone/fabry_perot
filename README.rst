**********
Fabry
**********

Package for analyzing Fabry-Perot interference patterns

Installation
**************

Prerequisites
================

* PyMultiNest can be installed via pip or via the github repository.

.. code-block:: 
    
    $ git clone https://github.com/JohannesBuchner/PyMultiNest/
    $ cd PyMultiNest
    $ python setup.py install


Use the "--user" switch if you want to install locally.

* PyMultiNest requires MultiNest to run. The simple instructions for building and installing are

.. code-block:: 
    
    $ git clone https://github.com/JohannesBuchner/MultiNest
    cd MultiNest/build
    cmake ..
    make


More detailed instructions are located `here <http://johannesbuchner.github.io/pymultinest-tutorial/install.html#on-your-own-computer>`_.

* Rawpy requires libraw. If you looking for a specific version, it can be installed from the source repository.
    
.. code-block::
    
    $ git clone https://github.com/LibRaw/LibRaw.git libraw
    $ git clone https://github.com/LibRaw/LibRaw-cmake.git libraw-cmake
    $ cd libraw
    $ git checkout 0.19.0
    $ cp -R ../libraw-cmake/* .
    $ cmake .
    $ sudo make install

Afterwards rawpy can be installed using ``pip install rawpy --no-binary rawpy``. 




