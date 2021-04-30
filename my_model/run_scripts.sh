echo 'You can run auto_complete.sh first to complete the BioM model.
    mc start stop var. Model
    6    x    y    3    ER
    6    x    y    2    GB
    6    x    y    3    Conf.
    6    x    y    4    GC
    6    x    y    5    BC
    6    x    y    6    BGC

Start and stop is up to user, this is just the number of iterations.
Other inputs are for each model.
*Microcconnectome (mc) is only for 6 at the moment, with a view to incuding 0-5 later.'

read -p 'Enter mc Number (0-6): ' number1
read -p 'Start: ' number2
read -p 'Stop: ' number3
read -p 'Enter variant (column 4): ' number4
read -p 'Enter model: ' number5

var1=1
var2=2
var3=3
var4=4

echo

if [ "$number4" -eq $var1 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        echo $i
        python3 Models.py "$number1" "$number5"
        python3 statistics.py "$number5"
        cd ../
    done
fi
if [ "$number4" -eq $var2 ]
then
    for i in $(seq ${number2} ${number3}); do
        echo $PWD
        cd ../reference_models/
        echo $PWD
        python3 neuron_swapping.py "$number1"
        python3 concatenate_py.py "$number1"
        cd ../my_model/output/GB/general_reconstruction/
        rm *
        cd ../../../src/
        echo $PWD
        python3 statistics.py "$number5"

        cd ../../reference_models/
        echo $PWD
        echo
    done
fi
if [ "$number4" -eq $var3 ]
then
    for i in $(seq ${number2} ${number3}); do
        echo $i
        cd ../reference_models/
        python3 ER.py "$number1"
        cd ../my_model/src/
        python3 statistics.py "$number5"
        cd ../
    done
fi
if [ "$number4" -eq $var4 ]
then
    cd ../reference_models/
    python3 biological_model.py "$number1"
    python3 conc.py "$number1"
    cd ../my_model/src/
    python3 split_by_layer.py
    cd ../output/Bio_M/reconstruction/
    rm -r *.csv
    cd ../../../src
    python3 statistics.py "$number5"
fi
