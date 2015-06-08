#!/bin/bash
bash test-shell.sh || exit
python setup.py test || exit