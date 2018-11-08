#!/bin/bash
testfile=fun_tax_interaction.txt
echo -e "peptide\tcog\tgenus\tint1\tint2\tint3\tint4\tint5\tint6" > $testfile
echo -e "A\tM, O\tStreptococcus\t105\t103\t102\t2\t1\t4" >> $testfile
echo -e "B\tJ\tStreptococcus\t1\t5\t3\t105\t102\t103" >> $testfile