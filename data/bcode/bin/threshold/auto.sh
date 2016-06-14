cd ~/Desktop/FYP_final/data/bcode/bin/threshold

# backup
cp ../../lib/exec/background_subtraction/background_subtraction.cpp background_subtraction_backup.cpp

cp background_subtraction20.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll20.txt

cp background_subtraction21.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll21.txt

cp background_subtraction22.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll22.txt

cp background_subtraction23.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll23.txt

cp background_subtraction24.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll24.txt

cp background_subtraction25.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll25.txt

cp background_subtraction26.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll26.txt

cp background_subtraction27.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll27.txt

cp background_subtraction28.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll28.txt

cp background_subtraction29.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll29.txt

cp background_subtraction30.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp
bash ../../../../run.sh
cp ../outputAll.txt outputAll30.txt

# restore
cp background_subtraction_backup.cpp ../../lib/exec/background_subtraction/background_subtraction.cpp