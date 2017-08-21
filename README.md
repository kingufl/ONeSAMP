ONeSAMP 2.0 version 1.0.1

By David Tallmon, Mark Heim, Christina Boucher, et al.

=============
=============
== License ==
=============
=============

This code is licensed under the GNU Public License, version 3.

---

ONeSAMP 2.0 computes the effective population size of gene data sets.
Copyright (C) 2016 Mark Heim, David Tallmon, Christina Boucher, et al.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---

==================
==================
== Introduction ==
==================
==================

The enclosed software will calculate the nucleotide effective diversity
population size (N_e) on an input file of unlinked, neutral SNPs or
microsatellites specified in the GENEPOP format, version 4.

The software works better for small populations (~100 individuals or less)
and small numbers of loci (~2000 loci or less).

Currently, ONeSAMP 2.0 is only available for Linux and Unix systems. It
may not be too difficult to eventually operate as a single node on a Mac.

=============
=============
== Install ==
=============
=============

To use the software, you must have the following software packages installed.
Alternative versions of these systems may also work.

gcc 4.8.3 with OpenMP                      https://gcc.gnu.org/
make 3.82                                  https://www.gnu.org/software/make/
perl 5.22.1                                https://www.perl.org/
Rscript 3.1.2                              https://www.r-project.org/
g++ 4.8.3                                  https://gcc.gnu.org/
gtest 1.0                                  https://github.com/google/googletest

Some of these programs may already be installed.

Programs that may have needed to have been used in prior releases:

simuPOP         (no longer needed as of this release)
simcoal2        (no longer needed as of this release)
python          (no longer needed as of this release)

====================
====================
== Usage overview ==
====================
====================

Please follow below steps.
1. generate simulalated populations with below command :
	./driver.sh generateSimulatedPopulationsOfKnownNe
	Before running this command make sure to set export NeVals="?" in driver.sh so that simulator can generate population for that Ne Value.
2. Run OneSamp : follow below steps.
==========
= STEP 1 =
==========
BUILD ONESAMP2

Build ONeSAMP2 using the supplied makefile using these commands:

cd ONeSAMP2; make

On the first run, the following lines should appear if everything is installed
correctly. (If not, run "make clean" and run "make".)

gcc -fopenmp -c ./refactor/engine/refactor_main.c -fopenmp -o ./refactor/engine/refactor_main.o
gcc -fopenmp -c ./refactor/engine/refactor_engine.c -fopenmp -o ./refactor/engine/refactor_engine.o
g++ -fopenmp -c ./refactor/engine/refactor_engine_test.cc -fopenmp -o ./refactor/engine/refactor_engine_test.o
gcc -fopenmp -c ./refactor/macro/refactor_macro.c -fopenmp -o ./refactor/macro/refactor_macro.o
g++ -fopenmp -c ./refactor/macro/refactor_macro_test.cc -fopenmp -o ./refactor/macro/refactor_macro_test.o
gcc -fopenmp -c ./refactor/arguments/refactor_arguments.c -fopenmp -o ./refactor/arguments/refactor_arguments.o
g++ -fopenmp -c ./refactor/arguments/refactor_arguments_test.cc -fopenmp -o ./refactor/arguments/refactor_arguments_test.o
gcc -fopenmp -c ./refactor/parser/refactor_parser.c -fopenmp -o ./refactor/parser/refactor_parser.o
g++ -fopenmp -c ./refactor/parser/refactor_parser_test.cc -fopenmp -o ./refactor/parser/refactor_parser_test.o
gcc -fopenmp -c ./refactor/memory/refactor_memory.c -fopenmp -o ./refactor/memory/refactor_memory.o
g++ -fopenmp -c ./refactor/memory/refactor_memory_test.cc -fopenmp -o ./refactor/memory/refactor_memory_test.o
gcc -fopenmp -c ./refactor/stats/refactor_stats.c -fopenmp -o ./refactor/stats/refactor_stats.o
g++ -fopenmp -c ./refactor/stats/refactor_stats_test.cc -fopenmp -o ./refactor/stats/refactor_stats_test.o
gcc  ./refactor/engine/refactor_main.o ./refactor/engine/refactor_engine.o ./refactor/arguments/refactor_arguments.o ./refactor/parser/refactor_parser.o ./refactor/memory/refactor_memory.o ./refactor/macro/refactor_macro.o ./refactor/stats/refactor_stats.o -fopenmp -o ./refactor/release/refactor_main  -lm
gcc ./refactor/engine/refactor_coalescent_simulator.c -fopenmp -o ./refactor/release/refactor_coalescent_simulator -lm
g++ -fopenmp -c ./refactor/engine/refactor_tests_main.cc -fopenmp -o ./refactor/engine/refactor_tests_main.o
g++  ./refactor/engine/refactor_tests_main.o ./refactor/engine/refactor_engine_test.o ./refactor/engine/refactor_engine.o ./refactor/macro/refactor_macro.o ./refactor/memory/refactor_memory.o ./refactor/stats/refactor_stats.o  ./refactor/arguments/refactor_arguments.o ./refactor/parser/refactor_parser.o -fopenmp -o ./refactor/engine/refactor_engine_test  -lgtest
g++  ./refactor/arguments/refactor_arguments.o ./refactor/macro/refactor_macro_test.o ./refactor/macro/refactor_macro.o ./refactor/engine/refactor_tests_main.o -fopenmp -o ./refactor/macro/refactor_macro_test  -lgtest
g++  ./refactor/arguments/refactor_arguments.o ./refactor/engine/refactor_tests_main.o ./refactor/arguments/refactor_arguments_test.o ./refactor/macro/refactor_macro.o -fopenmp -o ./refactor/arguments/refactor_arguments_test  -lgtest
g++  ./refactor/engine/refactor_tests_main.o ./refactor/parser/refactor_parser.o ./refactor/arguments/refactor_arguments.o ./refactor/parser/refactor_parser_test.o ./refactor/memory/refactor_memory.o ./refactor/macro/refactor_macro.o -fopenmp -o ./refactor/parser/refactor_parser_test  -lgtest
g++  ./refactor/engine/refactor_tests_main.o ./refactor/memory/refactor_memory_test.o ./refactor/memory/refactor_memory.o ./refactor/arguments/refactor_arguments.o ./refactor/macro/refactor_macro.o -fopenmp -o ./refactor/memory/refactor_memory_test  -lgtest
g++  ./refactor/engine/refactor_tests_main.o ./refactor/stats/refactor_stats_test.o ./refactor/stats/refactor_stats.o ./refactor/arguments/refactor_arguments.o ./refactor/macro/refactor_macro.o ./refactor/memory/refactor_memory.o -fopenmp -o ./refactor/stats/refactor_stats_test  -lgtest
g++  ./refactor/engine/refactor_engine_test.o ./refactor/macro/refactor_macro_test.o ./refactor/arguments/refactor_arguments_test.o ./refactor/parser/refactor_parser_test.o ./refactor/memory/refactor_memory_test.o ./refactor/stats/refactor_stats_test.o ./refactor/engine/refactor_tests_main.o ./refactor/engine/refactor_engine.o ./refactor/arguments/refactor_arguments.o ./refactor/parser/refactor_parser.o ./refactor/memory/refactor_memory.o ./refactor/macro/refactor_macro.o ./refactor/stats/refactor_stats.o -fopenmp -o ./refactor/release/refactor_all_test  -lgtest

==========
= STEP 2 =
==========
EDIT THE DRIVER TO REFLECT LOCATIONS OF COMPILERS, INTERPRETERS, AND STORAGE

Now, navigate to the subdirectory refactor/release and, with a text editor,
edit driver.sh:

cd refactor/release
gedit driver.sh

There are 4 steps of editing that must take place in the driver.sh file near
the top of the code. The one that is absolutely critical is step 2. A user
MUST enter in a legal path to 3 locations for temporary storage that are NOT
being used by any other program.

==========
= STEP 3 =
==========
COPY THE INPUT FILE IN ONeSAMP2/refactor/release/ directory
and cd to this directory by following command
cd ONeSAMP2/refactor/release/

CHECK IF THERE ARE SYNTAX ERRORS IN INPUT FILE
If your file doesn't have syntax errors then this step in not required

Suppose population1.gen is the input file then Run the following command:

./refactor_main -x < population1.gen

Fix any formatting errors that arise in the input file.

A successful run should include a line with -l and -i flags.

==========
= STEP 4 =
==========
Make sure you have only those gen files which for which you need to run ONeSAMP in release folder.

Then, append a special (.reduced) suffix using the Unix move command, if you have more than one gen file then you need to do it for all:

mv population1.gen population1.gen.reduced

==========
= STEP 5 =
==========
ENTER SIMULATION PARAMETERS

In "driver.sh", scroll down to a section marked:

###########
# ONeSAMP #
###########

There, you can edit several parameters. Make sure there aren't any spaces
between the parameter and the equals sign.

What you must edit are the following:

ONESAMP2COAL_MINALLELEFREQUENCY
(What is the minimum allele frequency necessary in order to keep that allele
in the simulation?)

mutationRate
(What is the microsatellite or SNP mutation rate? Add double
quotes to it, e.g., "0.000000012")

rangeNe
(In what range do we expect N_e to be in? This is the posterior distribution
of N_e. Add a comma but no spaces between the numbers, e.g., 6,1000. Also,
both numbers here MUST be even.)

theta
(In what range do we expect theta, a parameter used in Hoppe's urn model in the 
simulated populations, to be? Add a comma but no spaces between the
numbers, e.g., 0.000048,0.0048.)

numOneSampTrials
(How many populations do we want ONeSAMP2 to simulate? A good rule of thumb
may be at least 50000, but it may be good to go higher. See Tallmon et al. for
a graph of accuracy versus number of trials in regards to the original version
of ONeSAMP.)

duration
(In what range of number of generatiosn should ONeSAMP2 simulate a population
bottleneck? The default is a random number chosen from 2 to 8, inclusive. Add a
comma but no spaces between the numbers, e.g., 2,8.)

processPriority
(How intensively should ONeSAMP2 use machine resources? 19 is low priority, 0
is normal priority.)

microsatsOrSNPs
(Does this file contain microsats or SNPs? Use 'm' for microsatellites or 's'
for SNPs)

blocksize
(How many populations should ONeSAMP2 simulate before saving to disk? To reduce
RAM usage, set this to a low number.)

===========
= STEP 6  =
===========
CALCULATE REFERENCE POINT FOR POPULATION

Run the following command to calculate a reference point for (should take no more than a minute):

./driver.sh analyzeUnknownPopulation

===========
= STEP 7  =
===========
START ONESAMP

Ensure that step 6 is complete. Then, run

./driver.sh runOnesamp start

If network traffic is too high, run

Now, wait awhile. You can check ONeSAMP's progress by running:

./driver.sh runOnesamp progress

in another terminal in the same directory.

===========
= STEP 8  =
===========
STOP ONESAMP

Once the progress reads 100 percent or you want to stop the program anytime in middle, stop the onesamp with the command

./driver.sh runOnesamp stop

Alternatively, we can also check if the program has finished by checking if refactor_main and driver.sh programs
 are running or not by the following command :
 ps -e | grep driver
 ps -e | refactor
===========
= STEP 9 =
===========
RUN STATISTICAL ANALYSIS

Then, run the data analysis script in R to determine N_e. It should show
a message like the following:

Below are the mean, median, and 95 credible limits
for the posterior distribution of the effective
population size from OneSamp

mean        median      lower95CL   upper95CL
24.16       31.27       15.91       34.80

Run below command, it will create result.txt:
./driver.sh estimateNe

==========================
==========================
== Usage on the website ==
==========================
==========================

Currently under development as of writing; this should be up soon. Check back
soon at this URL:

http://plaza.ufl.edu/surajk95/onesamp/
