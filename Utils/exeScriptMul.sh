#!/bin/bash
python3 ../Utils/mat_gen.py 256x256 256x256
cp ../Utils/mat1.txt ../Serial/mat1_256.txt 
cp ../Utils/mat2.txt ../Serial/mat2_256.txt
cp ../Utils/mat1.txt ../OpenMP/mat1_256.txt 
cp ../Utils/mat2.txt ../OpenMP/mat2_256.txt  
python3 ../Utils/mat_gen.py 512x512 512x512
cp ../Utils/mat1.txt ../Serial/mat1_512.txt 
cp ../Utils/mat2.txt ../Serial/mat2_512.txt
cp ../Utils/mat1.txt ../OpenMP/mat1_512.txt 
cp ../Utils/mat2.txt ../OpenMP/mat2_512.txt  
python3 ../Utils/mat_gen.py 1024x1024 1024x1024
cp ../Utils/mat1.txt ../Serial/mat1_1024.txt 
cp ../Utils/mat2.txt ../Serial/mat2_1024.txt
cp ../Utils/mat1.txt ../OpenMP/mat1_1024.txt 
cp ../Utils/mat2.txt ../OpenMP/mat2_1024.txt 
python3 ../Utils/mat_gen.py 2048x2048 2048x2048
cp ../Utils/mat1.txt ../Serial/mat1_2048.txt 
cp ../Utils/mat2.txt ../Serial/mat2_2048.txt
cp ../Utils/mat1.txt ../OpenMP/mat1_2048.txt 
cp ../Utils/mat2.txt ../OpenMP/mat2_2048.txt 
python3 ../Utils/mat_gen.py 3072x3072 3072x3072
cp ../Utils/mat1.txt ../Serial/mat1_3072.txt 
cp ../Utils/mat2.txt ../Serial/mat2_3072.txt
cp ../Utils/mat1.txt ../OpenMP/mat1_3072.txt 
cp ../Utils/mat2.txt ../OpenMP/mat2_3072.txt 
gcc -O2 -w ../Serial/mat_mul.c -o matmul
cp ../Utils/matmul ../Serial/
gcc -O2 -w -fopenmp ../OpenMP/mat_mul.c -o matmul -lm
cp ../Utils/matmul ../OpenMP/
gcc -O2 -w -fopenmp ../OpenMP/mat_mul_blocks.c -o matmulOMP45 -lm
cp ../Utils/matmulOMP45 ../OpenMP/
chmod +x ../Serial/matmul
chmod +x ../OpenMP/matmul
chmod +x ../OpenMP/matmulOMP45
MAX=3
echo -e "*************** SERIAL - 1 THREAD ***************\n" >> ./resultMultiplication.txt
cd "../Serial"
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=2
echo -e "\n\n*************** OMP - 2 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
cd "../OpenMP/"
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=4
echo -e "\n\n*************** OMP - 4 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=8
echo -e "\n\n*************** OMP - 8 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=16
echo -e "\n\n*************** OMP - 16 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=24
echo -e "\n\n*************** OMP - 24 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\n\n****************************** OMP WITH TASK - PARALLELIZATION PER BLOCKS ******************************\n" >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=2
echo -e "\n\n*************** OMP - 2 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
cd "../OpenMP/"
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=4
echo -e "\n\n*************** OMP - 4 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=8
echo -e "\n\n*************** OMP - 8 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=16
echo -e "\n\n*************** OMP - 16 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
export OMP_NUM_THREADS=24
echo -e "\n\n*************** OMP - 24 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt; done
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
for i in {0..$MAX}; do ./matmulOMP45 mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt; done
echo "END" >> ./resultMultiplication.txt
