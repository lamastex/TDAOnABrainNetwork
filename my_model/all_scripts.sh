cat <<EOF

                                    Model Menu
EOF
cat <<EOF

===============================================================================================
||  mc  | start | stop | group |  Model  |  Graphics  |      Model Name                      ||
===============================================================================================
|| 0-6  |   1   |  1   |   4   |    7    |     1/0    | Bio-M MC (Must be run First)         ||
|| 0-6  |   1   |  1   |   3   |    1    |     1/0    | Erdős–Rényi Model                    ||
|| 0-6  |   1   |  1   |   2   |    2    |     1/0    | General Biological Model             ||
|| 0-6  |   1   |  1   |   1   |    3    |     1/0    | Configuration Model                  ||
|| 0-6  |   1   |  1   |   1   |    4    |     1/0    | Geometric Configuration Model        ||
|| 0-6  |   1   |  1   |   5   |    5    |     1/0    | Block Configuration Model            ||
|| 0-6  |   1   |  1   |   6   |    6    |     1/0    | Block Geometric Configuration Model  ||
===============================================================================================
EOF

read -p 'Enter mc (0-6): ' number1
read -p 'Start: ' number2
read -p 'Stop: ' number3
read -p 'Enter group: ' number4
read -p 'Enter model: ' number5
read -p 'Graphics: ' number6

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
        python3 Models.py "$number1" "$number5" "$number6"
        python3 statistics.py "$number1" "$number5" "$number6"
        cd ../
    done
fi
if [ "$number4" -eq $var2 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5" "$number6"
        julia swap_concat.jl
        cd ../output/GB/general_reconstruction/
        rm *
        cd ../../../src/
        python3 statistics.py "$number1" "$number5" "$number6"
        cd ../
    done
fi
if [ "$number4" -eq $var3 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5" "$number6"
        python3 statistics.py "$number1" "$number5" "$number6"
        cd ../
    done
fi
if [ "$number4" -eq $var4 ]
then
    cd src/
    python3 Models.py "$number1" "$number5" "$number6"
    julia concatenation.jl
    python3 -c "import functions"
    cd ../output/Bio_M/reconstruction/
    rm *.csv
    cd ../../../src
    python3 statistics.py "$number1" "$number5" "$number6"
    cd ../
fi


if [ "$number4" -eq $var5 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5" "$number6"
        python3 statistics.py "$number1" "$number5" "$number6"
        cd ../output/BC/new_blocks/
        rm *.npy
        cd ../../../
    done
fi

if [ "$number4" -eq $var6 ]
then
    for i in $(seq ${number2} ${number3}); do
        cd src/
        python3 Models.py "$number1" "$number5" "$number6"
        python3 statistics.py "$number1" "$number5" "$number6"
        cd ../output/BGC/new_blocks/
        rm *.npy
        cd ../../../
    done
fi
