echo 'This script runtime can vary considerably with regards to each MC. An optimal setting so far is "2 f 1 -1 -1". The -1s are used in order to not lose accuracy when recomputing the connections.

read -p 'Enter mc Number (0-6): ' number1
read -p 'Enter m_type (f/e/i): ' number2
read -p 'Enter instance (1-5,11-15,21-25): ' number3
read -p 'Enter bin counts (-1, -2,...): ' number4
read -p 'Enter probability (-1, -2, ...) (same as bins): ' number5

if [ -z "$number1" ] || [ -z "$number2" ] || [ -z "$number3" ] || [ -z "$number4" ] || [ -z "$number5" ]
then
    echo 'Inputs cannot be blank please try again'
    exit 0
else
    cd src/
    python3 h5_to_csv.py "$number1" "$number2"
    python3 concatenate_csv.py "$number1" "$number2"
    python3 reconfig_pre.py "$number1" "$number3" "$number4" "$number5" "$number2"
    python3 compute_stats.py "$number1" "$number3" "$number2"
fi
