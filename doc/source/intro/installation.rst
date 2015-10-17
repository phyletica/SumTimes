.. _installation:

************
Installation
************


.. _prerequisites:

Prerequisites
=============

|sumtimes|_ has been tested on versions 2.7 and 3.4 of |python|. There are two |python|
dependencies:

*   |dendro|_ :cite:`Dendropy`
*   |pyyaml|_

If you do not have these modules, you will need to install them before you can
use |sumtimes|_. You can do this via the command line::

    $ sudo pip install dendropy
    $ sudo pip install pyyaml

If you do not have admin privileges, or you prefer to install locally, you can
use::

    $ pip install --user dendropy
    $ pip install --user pyyaml


Installing sumtimes.py
======================

You should have |python|_ and :ref:`the required modules<prerequisites>`
installed.

First, if you have |git|_ installed, clone the Git repository::

    $ git clone https://github.com/phyletica/SumTimes.git

If you do not have |git|_, you can download a snapshot of the repository as a
tar or zip archive by clicking the respective folder icons at the top of this
page. Either way, move into the downloaded directory::

    $ cd SumTimes

Next, let's make sure |sumtimes|_ is working as expected on your machine:

.. parsed-literal::

    $ ./|st| --run-tests

This will run a suite of tests.
All of |sumtimes|_ is self-contained within the |lst| script, so to install,
all you need to do is copy (or link) the |lst| file to a location that is in
your PATH. For example:

.. parsed-literal::

    $ sudo cp |st| /usr/local/bin

Or:

.. parsed-literal::

    $ cp |st| ~/bin

Alternatively, you can run |lst| from your current directory by using the
``./sumtimes.py`` like we did above when we ran the tests. However, the
remainder of this documentation assumes the script is in your PATH, so if it is
not, you will have to remember to add the "``./``" bit.

You should be able to call up the help menu of |lst|:

.. parsed-literal::

    $ ./|st| -h

You should see output that looks like::

    usage: sumtimes.py [-h] [--np NP] [--seed SEED] YAML-CONFIG-PATH
    
    sumtimes.py Version 0.1.0
    
    positional arguments:
      YAML-CONFIG-PATH  Path to the YAML-formatted config file.
    
    optional arguments:
      -h, --help        show this help message and exit
      --np NP           The maximum number of processes to run in parallel. The
                        default is the smaller of the number of tree files or the
                        number of CPUs available on the machine.
      --seed SEED       Random number seed to use for the analysis.
