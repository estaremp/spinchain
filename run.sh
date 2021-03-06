#!/bin/bash

#Clean
#make clean
rm -rf ./bin

#Move state into main folder for initialisation
mv ./output/final_state.data ./

rm -rf ./output

#Compile files
make

./run

#This should probably be read using grep
#Read variables to be used
source 'info.data'

#Create output directory
mkdir output

cp src/Plotting/*.py ./bin


if [ $MAX_EOF = 'T' ] && [ $REALISATIONS -eq 1 ] ;
then
    if [ $METHOD = 'DIAG' ] ;
    then
    python python_scripts/maxEOF.py eof.data $REALISATIONS $TOTALTIME $OFFNOISE $DIAGNOISE
    fi
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

    python python_scripts/averageEigenvalues.py eigenvalues_realisations.data $VECTORS $REALISATIONS
    if [ $SINGLE = 'T' ] ;
    then
    python python_scripts/averageEOF.py eof_realisations_SP.data $TA $OFFNOISE $DIAGNOISE
    elif [ $MAX_EOF = 'T' ] ;
    then
    python python_scripts/maxEOF.py eof_realisations_FD.data $REALISATIONS $TOTALTIME $OFFNOISE $DIAGNOISE
    fi

    if [ $GRAPHICAL = 'T' ] ;
    then
    #Plot averaged
    python bin/averagedEigenvalues.py $VECTORS
    fi

    #If Plots have been created move them to the required folder
    mkdir output/Plots
    mv ./*.png ./output/Plots

fi

#Graphical outputs (not done if single time realisations or integration method)
if [ $GRAPHICAL = 'T' ] && [ $SINGLE = 'F' ] ;
then
    #Plot stuff
    echo ' >> Plot'

    if [ $METHOD = 'DIAG' ] ;
    then
    python bin/probabilities.py $VECTORS
    python bin/eigenvalues.py $VECTORS
        if [ $EOF = 'T' ] ;
        then
        python bin/eof.py $TOTALTIME $INITIALVEC
        fi
    fi
    python bin/dynamics.py $TOTALTIME $INITIALVEC $N
    python bin/exmap.py $TOTALTIME


    #If Plots have been created move them to the required folder
    mkdir output/Plots
    mv ./*.png ./output/Plots
fi


mv ./*.data ./*.out ./output

#Clean everything
make clean
rm -rf ./bin


