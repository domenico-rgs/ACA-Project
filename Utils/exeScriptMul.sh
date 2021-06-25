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
gcc -O3 -w ../Serial/mat_mul.c -o matmul
cp ../Utils/matmul ../Serial/
gcc -O3 -w -fopenmp ../OpenMP/mat_mul.c -o matmul -lm
cp ../Utils/matmul ../OpenMP/
chmod +x ../Serial/matmul
chmod +x ../OpenMP/matmul
echo -e "*************** SERIAL - 1 THREAD ***************\n" >> ./resultMultiplication.txt
cd "../Serial"
echo -e "\nN 256\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=2
echo -e "\n\n*************** OMP - 2 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
cd "../OpenMP/"
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=4
echo -e "\n\n*************** OMP - 4 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=8
echo -e "\n\n*************** OMP - 8 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=16
echo -e "\n\n*************** OMP - 16 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
export OMP_NUM_THREADS=24
echo -e "\n\n*************** OMP - 24 THREADs ***************\n" >> ../Utils/resultMultiplication.txt
echo -e "\nN 256\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_256.txt mat2_256.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 512\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_512.txt mat2_512.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 1024\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_1024.txt mat2_1024.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 2048\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_2048.txt mat2_2048.txt >> ../Utils/resultMultiplication.txt
echo -e "\nN 3072\n\n" >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
./matmul mat1_3072.txt mat2_3072.txt >> ../Utils/resultMultiplication.txt
echo "END" >> ./resultMultiplication.txt
mail -s 'Computation complete' domenico.ragusa01@universitadipavia.it <<< 'Your computation batch is now finished. Check the results'
