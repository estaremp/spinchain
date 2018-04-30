 #!/bin/bash

#RUN OVER DISORDER WEIGHTED AGAINST $DELT
for e in `seq 0.0 0.1 0.5`;
do
    echo 'Noise:' $e

    #modify the code with level of disorder (diagonal-ENERC) and (offdiagonal-NOISEV)
    sed 's/SED_ENERC/'"$e"' /g' BASECODEb.f90 > BASECODEout.f90

done



