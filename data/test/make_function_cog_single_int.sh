#!/bin/bash

testfile=function_cog_single_int.tab
echo -e "peptide\tcog\tint" > $testfile
echo -e "A\tO\t10" >> $testfile
echo -e "B\tC\t20" >> $testfile
echo -e "C\tO\t40" >> $testfile