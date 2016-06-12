git add /home/pi/Documents/FYP_final/data_pi/bcode/lib/exec/background_subtraction/background_subtraction.cpp

cd /home/pi/Documents/FYP_final/data_pi/bcode/lib/exec/background_subtraction
make -j 4

cd /home/pi/Documents/FYP_final/data_pi/bcode/bin
./background_subtraction.ln -i ./datasets/K01.txt -o ./datasets/K01_out.txt
./background_subtraction.ln -i ./datasets/K07.txt -o ./datasets/K07_out.txt
./background_subtraction.ln -i ./datasets/KDC2-KD04.txt -o ./datasets/KDC2-KD04_out.txt
./background_subtraction.ln -i ./datasets/Mo02.txt -o ./datasets/Mo02_out.txt
./background_subtraction.ln -i ./datasets/Mo03.txt -o ./datasets/Mo03_out.txt
./background_subtraction.ln -i ./datasets/Mo06.txt -o ./datasets/Mo06_out.txt
./background_subtraction.ln -i ./datasets/Mo12.txt -o ./datasets/Mo12_out.txt
./background_subtraction.ln -i ./datasets/Mo13.txt -o ./datasets/Mo13_out.txt
./background_subtraction.ln -i ./datasets/MT09.txt -o ./datasets/MT09_out.txt
./background_subtraction.ln -i ./datasets/MT10.txt -o ./datasets/MT10_out.txt
./background_subtraction.ln -i ./datasets/MT30.txt -o ./datasets/MT30_out.txt
./background_subtraction.ln -i ./datasets/MTEXT06.txt -o ./datasets/MTEXT06_out.txt
./background_subtraction.ln -i ./datasets/PI05.txt -o ./datasets/PI05_out.txt
./background_subtraction.ln -i ./datasets/PI10.txt -o ./datasets/PI10_out.txt
./background_subtraction.ln -i ./datasets/PI30.txt -o ./datasets/PI30_out.txt
./background_subtraction.ln -i ./datasets/PI35.txt -o ./datasets/PI35_out.txt
./background_subtraction.ln -i ./datasets/SY09.txt -o ./datasets/SY09_out.txt
./background_subtraction.ln -i ./datasets/SY13.txt -o ./datasets/SY13_out.txt
./background_subtraction.ln -i ./datasets/SY16.txt -o ./datasets/SY16_out.txt
./background_subtraction.ln -i ./datasets/SY31.txt -o ./datasets/SY31_out.txt
./background_subtraction.ln -i ./datasets/T02.txt -o ./datasets/T02_out.txt
./background_subtraction.ln -i ./datasets/T05.txt -o ./datasets/T05_out.txt
./background_subtraction.ln -i ./datasets/T08.txt -o ./datasets/T08_out.txt
./background_subtraction.ln -i ./datasets/T11.txt -o ./datasets/T11_out.txt
./background_subtraction.ln -i ./datasets/T12.txt -o ./datasets/T12_out.txt
./background_subtraction.ln -i ./datasets/T15.txt -o ./datasets/T15_out.txt
./background_subtraction.ln -i ./datasets/W09.txt -o ./datasets/W09_out.txt
./background_subtraction.ln -i ./datasets/W11.txt -o ./datasets/W11_out.txt
./background_subtraction.ln -i ./datasets/W11b.txt -o ./datasets/W11b_out.txt
./background_subtraction.ln -i ./datasets/W12.txt -o ./datasets/W12_out.txt
./background_subtraction.ln -i ./datasets/W13.txt -o ./datasets/W13_out.txt
./background_subtraction.ln -i ./datasets/W13all.txt -o ./datasets/W13all_out.txt
./background_subtraction.ln -i ./datasets/W17.txt -o ./datasets/W17_out.txt


cd /home/pi/Documents/FYP_final/data_pi/bcode/bin
git add ./datasets/
