#cd ~/Desktop/FYP_final/data/bcode/bin/diff_thresh/

# backup
cp ../../lib/exec/background_subtraction/background_subtraction.cpp background_subtraction_backup.cpp

cp background_subtraction50.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll50.txt

# cp background_subtraction52.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll52.txt

# cp background_subtraction54.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll54.txt

# cp background_subtraction56.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll56.txt

# cp background_subtraction58.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll58.txt

cp background_subtraction60.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll60.txt

# cp background_subtraction62.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll62.txt

# cp background_subtraction64.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll64.txt

# cp background_subtraction66.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll66.txt

# cp background_subtraction68.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll68.txt

# cp background_subtraction70.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
# bash ../../../../run.sh
# cp ../outputAll.txt outputAll70.txt

# restore
cp background_subtraction_backup.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp