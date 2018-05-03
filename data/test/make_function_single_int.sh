#!/bin/bash

testfile=function_single_int.tab
echo -e "peptide\tgo\tint" > $testfile
echo -e "A\tGO:0008152\t10" >> $testfile
echo -e "B\tGO:0022610\t20" >> $testfile