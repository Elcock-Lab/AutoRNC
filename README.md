# AutoRNC
Fortran package for generating ribosome nascent-chain models

# Installation instructions:

compiling the code:
-------------------

cd SOURCE/
for i in *.f ; do gfortran $i -o $i.exe ; done
rename .f.exe .exe *exe
cd ..

note that three ribosome templates are provided: 5UYM, 3JBV, 7OIZ
each template is provided in two forms, e.g. for 5UYM:

folder 5UYM     : this contains an unmodified 5UYM structure
folder 5UYM_L24 : this contains a version with the L24 loop rebuilt
                  so that it doesn't occlude the tunnel exit

download the dipeptide structure libraries from zenodo:
-------------------------------------------------------

make sure you have space! they occupy ~3.4 GB...
you can find the libraries here:

https://zenodo.org/records/8101871

try running the examples:
-------------------------

There are five example folders provided, each of which corresponds
to a panel of Figure 2 from the McDonnell/Elcock paper. Each folder
contains a README which outlines the commands used and all the
necessary arguments. Note that the most challenging example case
(#5) can take a while to produce a conformation

Running AutoRNC basically comes down to executing two commands:

(1) make_comfile_AHEMODEL_nascent_chains_for_release_final_2023.exe

This is a script-generating code that takes a bunch of arguments
provided by the user and generates a single script that can be used
to generate RNC models.

*********************************************************************
note that you will need to edit the command line arguments in the examples
to reflect the paths to the various folders on your system...
*********************************************************************

successful execution of the make_comfile...exe code will generate a
script with the following name:

comfile_build_nascent_chain

(2) run the script thus:

nohup ./comfile_build_nascent_chain &

note that when the script runs it will generate a new folder
called AHEMODEL_001 ; the two .pdb files that you want are:

AHEMODEL_001/P12345.AHEMODEL.pdb  : this contains the nascent chain model(s) - if
                                    multiple models are made they are separated by
                                    MODEL/ENDMDL statements in the .pdb file so can
                                    be viewed as a movie trajectory in VMD
AHEMODEL_001/P12345.HETATMs.pdb   : this contains the ribosome & any blocking atoms

note that the reason that the nascent chain is written to a separate file from the
ribosome is that in typical usage one will be making many nascent chain conformations,
so there is no point continually writing out the ribosome coordinates when they will
be unchanged...

