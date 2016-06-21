#cd ~/Desktop/FYP_final/data/bcode/bin/min_coverage

# backup
cp ../../lib/exec/background_subtraction/background_subtraction.cpp background_subtraction_backup.cpp

cp background_subtraction005.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll005.txt

cp background_subtraction007.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll007.txt

# cp background_subtraction009.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll009.txt

cp background_subtraction010.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll010.txt

# cp background_subtraction020.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll020.txt

cp background_subtraction030.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll030.txt

# cp background_subtraction040.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll040.txt

cp background_subtraction050.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll050.txt

# cp background_subtraction060.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll060.txt

cp background_subtraction070.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll070.txt

# cp background_subtraction080.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll080.txt

# cp background_subtraction090.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll090.txt

# cp background_subtraction100.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll100.txt

# restore
cp background_subtraction_backup.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp