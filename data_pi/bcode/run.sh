git add /home/pi/Documents/FYP_final/data_pi/bcode/lib/exec/background_subtraction/background_subtraction.cpp

cd /home/pi/Documents/FYP_final/data_pi/bcode/lib/exec/background_subtraction
make

cd /home/pi/Documents/FYP_final/data_pi/bcode/bin
#./background_subtraction.ln -i testAll_pi.txt -o outputAll_pi.txt
./background_subtraction.ln -i test01_pi.txt -o output01_pi.txt

cd /home/pi/Documents/FYP_final/data_pi/bcode/bin
#git add testAll_pi.txt outputAll_pi.txt
git add test01_pi.txt output01_pi.txt

