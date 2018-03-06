#!/bin/bash

#Clean
rm -rf output

#Compile files
make

./run

#Read variables to be used
source 'info.data'

#Create output directory
mkdir output

#Realisations if noise
if [ $REALISATIONS -ne 0 ] ;
then
    for j in `seq 1 1 $REALISATIONS`
    do
        echo 'Realisation num:' $j
        #Run executable
        ./run
        cat eigenvalues.data >> eigenvalues_realisations.data
    done
fi

#Graphical outputs (not done if realisations)
if [ $GRAPHICAL = 'T' ] ;
then
    #Plot stuff
    echo ' >> Plot'

    mv src/Plotting/*.pyc ./bin

    python2 bin/probabilities.pyc $VECTORS
    python2 bin/eigenvalues.pyc $VECTORS
    python2 bin/dynamics.pyc $TOTALTIME $INITIALVEC $N
    python2 bin/exmap.pyc

    mkdir output/Plots
    mv ./*.png ./output/Plots
fi

mv ./*.data ./*.out ./output

#Clean everything
make clean
rm -rf ./bin


