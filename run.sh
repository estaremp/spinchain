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

mv src/Plotting/*.pyc ./bin

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

    python2 python_scripts/averages.py eigenvalues_realisations.data $VECTORS

    if [ $GRAPHICAL = 'T' ] ;
    then
    #Plot averaged
    python2 bin/averagedEigenvalues.pyc $VECTORS
    fi

fi

#Graphical outputs (not done if realisations)
if [ $GRAPHICAL = 'T' ] ;
then
    #Plot stuff
    echo ' >> Plot'

    python2 bin/probabilities.pyc $VECTORS
    python2 bin/eigenvalues.pyc $VECTORS
    python2 bin/dynamics.pyc $TOTALTIME $INITIALVEC $N
    python2 bin/exmap.pyc
    if [ $EOF = 'T' ] ;
    then
        python2 bin/eof.pyc $TOTALTIME $INITIALVEC
    fi
fi

#If Plots have been created move them to the required folder
if [ $GRAPHICAL = 'T' ] || [ $REALISATIONS -ne 0 ] ;
then
    mkdir output/Plots
    mv ./*.png ./output/Plots
fi

mv ./*.data ./*.out ./output

#Clean everything
make clean
rm -rf ./bin


