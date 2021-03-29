echo 'Script runtime is around 10 minutes'

read -p 'Enter mc Number (0-6): ' number1
read -p 'Enter bin counts (-1, -2,...): ' number2
read -p 'Enter probability (-1, -2, ...) (same as bins): ' number3
read -p 'Enter start: ' number4
read -p 'Enter stop: ' number5
read -p 'Enter model: ' number6

var1=1
var2=2
var3=3

if [  "$number6" -eq $var1 ]
then
    cd src/
    python3 compute_stats.py "$number1" "$number2" "$number3" "$number4" "$number5" "$number6"
elif [  "$number6" -eq $var2 ]
then
    cd src/
    python3 compute_stats.py "$number1" "$number2" "$number3" "$number4" "$number5" "$number6"
elif [  "$number6" -eq $var3 ]
then
    cd src/
    python3 block_geometric_configuration.py "$number1" "$number4" "$number5"
fi
