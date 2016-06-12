git add /home/karenng/Desktop/FYP_final/data/bcode/lib/exec/background_subtraction/background_subtraction.cpp

cd /home/karenng/Desktop/FYP_final/data/bcode/lib/exec/background_subtraction
make

cd /home/karenng/Desktop/FYP_final/data/bcode/bin
# ./background_subtraction.ln -i testAll.txt -o outputAll.txt
./background_subtraction.ln -i test01.txt -o output01.txt
# ./background_subtraction.ln -i test02.txt -o output02.txt
# ./background_subtraction.ln -i test03.txt -o output03.txt
# ./background_subtraction.ln -i test04.txt -o output04.txt
# ./background_subtraction.ln -i test05.txt -o output05.txt

cd /home/karenng/Desktop/FYP_final/data/bcode/bin
# git add testAll.txt outputAll.txt
git add test01.txt output01.txt
# git add test02.txt output02.txt
# git add test03.txt output03.txt
# git add test04.txt output04.txt
# git add test05.txt output05.txt
