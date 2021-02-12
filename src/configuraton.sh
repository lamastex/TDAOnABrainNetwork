echo 'This script will take approximately 15 minutes to run, so patience is required! 5 arguments are required; connectome number, m_type, your instance, number of bins to remove and number of probability values to remove (minimum for both is -1)'

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
    python3 h5_to_csv.py "$number1" "$number2"
    python3 concatenate_csv.py "$number1" "$number2"
    python3 config_pre.py "$number1" "$number3" "$number4" "$number5" "$number2"
    python3 rebuild_array_compute_ts.py "$number1" "$number3" "$number2"
fi
