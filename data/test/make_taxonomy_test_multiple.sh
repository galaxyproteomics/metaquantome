#!/bin/bash

# echo -e 'peptide\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tint1\tint2\tint3\tint4\n'\
# 'NOTISSREAL\tbacteria\tproteobacteria\tepsilonproteobacteria\tcampylobacterales\thelicobacteraceae\theliobacter\tpylori\t2000\t5000\t9000\t1000\n'\
# 'ISSNOTREALE\tbacteria\tfirmicutes\tclostridia\tclostridiales\tclostridiaceae\tclostridium\tdifficile\t7000\t2000\t1000\t500\n'\
# 'ISSVERYNOTREALE\tbacteria\tfirmicutes\tclostridia\tclostridiales\tclostridiaceae\t\t\t1000\t3000\t0\t0\n'\
# 'AAAADQLSIWADR\tbacteria\tfirmicutes\tnegativicutes\tveillonaellales\tveillonellaceae\tveillonella\t\t0\t0\t100\t0' > taxonomy_test_multiple.tab

# in order:
# helicobacter pylori (only proteobacteria)
# clostridioides difficile
# clostridiaceae
# veillonella

echo -e 'peptide\tlca\tint1\tint2\tint3\tint4\n'\
'NOTISSREAL\t210\t2000\t5000\t9000\t1000\n'\
'ISSNOTREALE\t1496\t7000\t2000\t1000\t500\n'\
'ISSVERYNOTREALE\t31979\t1000\t3000\t0\t0\n'\
'AAAADQLSIWADR\t29465\t0\t0\t100\t0'\
> taxonomy_test_multiple.tab
