echo 'Running order:      mc start stop group   Model
                    6  1       1      4       7 BioM model (Must be run First)
                    6  1       1      3       1 Erdos-Renyi model
                    6  1       1      2       2 General Biol model
                    6  1       1      1       3 config model
                    6  1       1      1       4 geometric config model
                    6  1       1      5       5 block config model
                    6  1       1      6       6 block geometric config model'

read -p 'Enter mc (0-6): ' number1
read -p 'Start: ' number2
read -p 'Stop: ' number3
read -p 'Enter group: ' number4
read -p 'Enter model: ' number5

var1=1
var2=2
var3=3
var4=4
var5=5
var6=6


if [ "$number4" -eq $var1 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5"
        python3 statistics.py "$number5"
        cd ../
    done
fi
if [ "$number4" -eq $var2 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5"
        julia swap_concat.jl
        cd ../output/GB/general_reconstruction/
        rm *
        cd ../../../src/
        python3 statistics.py "$number5"
        cd ../
    done
fi
if [ "$number4" -eq $var3 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5"
        python3 statistics.py "$number5"
        cd ../
    done
fi
if [ "$number4" -eq $var4 ]
then
    cd src/
    python3 Models.py "$number1" "$number5"
    julia concatenation.jl
    python3 -c "import functions; functions.complete_blocks(6)"
    cd ../output/Bio_M/reconstruction/
    rm *.csv
    cd ../../../src
    python3 statistics.py "$number5"
    cd ../
fi


if [ "$number4" -eq $var5 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5"
        python3 statistics.py "$number5"
        cd ../output/BC/blocks/
        rm *.npy
        cd ../new_blocks/
        rm *.npy
        cd ../../../../
    done
fi

if [ "$number4" -eq $var6 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5"
        python3 statistics.py "$number5"
        cd ../output/BC/blocks/
        rm *.npy
        cd ../../BGC/new_blocks/
        rm *.npy
        cd ../../../
    done
fi
