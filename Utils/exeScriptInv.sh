#!/bin/bash
python3 ../Utils/mat_gen.py 512x512
cp ../Utils/mat1.txt ../Serial/mat1_512.txt 
cp ../Utils/mat1.txt ../OpenMP/mat1_512.txt 
python3 ../Utils/mat_gen.py 1024x1024
cp ../Utils/mat1.txt ../Serial/mat1_1024.txt 
cp ../Utils/mat1.txt ../OpenMP/mat1_1024.txt 
python3 ../Utils/mat_gen.py 2048x2048
cp ../Utils/mat1.txt ../Serial/mat1_2048.txt 
cp ../Utils/mat1.txt ../OpenMP/mat1_2048.txt 
python3 ../Utils/mat_gen.py 3072x3072
cp ../Utils/mat1.txt ../Serial/mat1_3072.txt 
cp ../Utils/mat1.txt ../OpenMP/mat1_3072.txt 
python3 ../Utils/mat_gen.py 4096x4096
cp ../Utils/mat1.txt ../Serial/mat1_4096.txt 
cp ../Utils/mat1.txt ../OpenMP/mat1_4096.txt 
gcc -O2 -w ../Serial/mat_inv.c -o matinv -lm
cp ../Utils/matinv ../Serial/
gcc -O2 -w -fopenmp ../OpenMP/mat_inv.c -o matinv -lm
cp ../Utils/matinv ../OpenMP/
chmod +x ../Serial/matinv
chmod +x ../OpenMP/matinv
MAX=3
echo -e "*************** SERIAL - 1 THREAD ***************\n" >> ./resultInversion.txt
cd "../Serial"
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
export OMP_NUM_THREADS=2
echo -e "\n\n*************** OMP - 2 THREADs ***************\n" >> ../Utils/resultInversion.txt
cd "../OpenMP/"
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
export OMP_NUM_THREADS=4
echo -e "\n\n*************** OMP - 4 THREADs ***************\n" >> ../Utils/resultInversion.txt
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
export OMP_NUM_THREADS=8
echo -e "\n\n*************** OMP - 8 THREADs ***************\n" >> ../Utils/resultInversion.txt
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
export OMP_NUM_THREADS=16
echo -e "\n\n*************** OMP - 16 THREADs ***************\n" >> ../Utils/resultInversion.txt
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
export OMP_NUM_THREADS=24
echo -e "\n\n*************** OMP - 24 THREADs ***************\n" >> ../Utils/resultInversion.txt
echo -e "\nN 512\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_512.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_1024.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_2048.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_3072.txt >> ../Utils/resultInversion.txt; done
echo -e "\nN 4096\n" >> ../Utils/resultInversion.txt
for i in {0..$MAX}; do ./matinv mat1_4096.txt >> ../Utils/resultInversion.txt; done
echo "END" >> ./result.txt
