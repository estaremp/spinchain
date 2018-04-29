#!/bin/bash

#Clean
make clean
rm -rf ./bin
rm -rf output

#Compile files
make

./run

#This should probably be read using grep
#Read variables to be used
source 'info.data'

#Create output directory
mkdir output

mv src/Plotting/*.pyc ./bin

if [ $MAX_EOF = 'T' ] && [ $REALISATIONS -eq 1 ];
then
python2 python_scripts/maxEOF.py eof.data $REALISATIONS $TOTALTIME $OFFNOISE $DIAGNOISE
fi

#Realisations if noise
if [ $REALISATIONS -ne 1 ] ;
then
    for j in `seq 1 1 $REALISATIONS`
    do
        echo 'Realisation num:' $j
        #Run executable
        ./run

        cat eigenvalues.data >> eigenvalues_realisations.data

        if [ $SINGLE = 'T' ] ;
        then
        cat eof.data >> eof_realisations_SP.data
        elif [ $MAX_EOF = 'T' ] ;
        then
        cat eof.data >> eof_realisations_FD.data
        fi
    done

    python2 python_scripts/averageEigenvalues.py eigenvalues_realisations.data $VECTORS $REALISATIONS
    if [ $SINGLE = 'T' ] ;
    then
    python2 python_scripts/averageEOF.py eof_realisations_SP.data $TA $OFFNOISE $DIAGNOISE
    elif [ $MAX_EOF = 'T' ] ;
    then
    python2 python_scripts/maxEOF.py eof_realisations_FD.data $REALISATIONS $TOTALTIME $OFFNOISE $DIAGNOISE
    fi

    if [ $GRAPHICAL = 'T' ] ;
    then
    #Plot averaged
    python2 bin/averagedEigenvalues.pyc $VECTORS
    fi

    #If Plots have been created move them to the required folder
    mkdir output/Plots
    mv ./*.png ./output/Plots

fi

#Graphical outputs (not done if single time realisations)
if [ $GRAPHICAL = 'T' ] && [ $SINGLE = 'F' ]  ;
then
    #Plot stuff
    echo ' >> Plot'

    python2 bin/probabilities.pyc $VECTORS
    python2 bin/eigenvalues.pyc $VECTORS
    python2 bin/dynamics.pyc $TOTALTIME $INITIALVEC $N
    python2 bin/exmap.pyc $TOTALTIME
    if [ $EOF = 'T' ] ;
    then
        python2 bin/eof.pyc $TOTALTIME $INITIALVEC
    fi

    #If Plots have been created move them to the required folder
    mkdir output/Plots
    mv ./*.png ./output/Plots
fi


mv ./*.data ./*.out ./output

#Clean everything
make clean
rm -rf ./bin


