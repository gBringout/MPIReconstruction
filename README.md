# MPI Reconsutrion #

This repository contains exemplary Matlab codes for simple MPI reconstructions.

Data are provided [here] (http://www.tuhh.de/ibi/research/mpi-data-format.html) or [here] (https://github.com/KsenijaGraefe/SingleSidedData) and have to be placed in the data folder. The file `measurement_5.h5` from the first adresse has to be renamed to `measurement.h5`.

# reco.m #

Running this script, you should obtain after roughly one minutes those three graphs:

<img src="/results/SM.jpg" height="200">

<img src="/results/SpectrumMeasure.jpg" height="200">

<img src="/results/Reco.jpg" height="200">

The last one present the results of reconstruction of the concentration map a tracer using the system matrix/calibration approach. Using the same measurements, the inconsistent system of linear equations is solved using the signal acquired by a single channel of the scanner.

Five algorithms are used to solved it:
 1. A least square approach,
 2. an Algebraic Reconstruction Technique (ART) also known as the Kaczmarz's algorithm,
 3. a modified ART algorithm, forcing a positive and real approximation of the solution at the start of each iteration,
 4. a modified ART algorithm, regularizing and forcing a positive and real approximation of the solution,
 5. pseudoinverse approach.


# Reco_IWMPI2016.m #

Running the script will required quite around 4 GB of RAM and 100 minutes and will produce the results presented at IWMPI 2016.


You can try whatever you want to improve these reconstruction!
