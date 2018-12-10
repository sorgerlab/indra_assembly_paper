#!/bin/bash
export VIRTUALENVWRAPPER_PYTHON=/usr/local/bin/python3
export WORKON_HOME=/pmc/virtualenvs/
source /usr/local/bin/virtualenvwrapper.sh
workon py36
export PYTHONPATH=/pmc/indra
for i in  `seq 1 50`;
do
    nohup python -u run_reading.py $i &>> reader_log$i.txt &
done
