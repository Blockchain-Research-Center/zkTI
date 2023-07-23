#!/bin/sh

# first build the C++ project and get the ./zkTI executable program.

answer_path="../datasets/d_duck_identification/answer.csv"
truth_path="../datasets/d_duck_identification/truth.csv"
task_n=(108 100 100 100 80 100 60 40 50 10 5)
worker_n=(39 38 30 25 25 20 25 25 10 10 2)

n=${#task_n[@]}

for (( i=0; i<${n}; i++ ))
do
  ./build/src/zkTI ${answer_path} ${truth_path} ${task_n[$i]} ${worker_n[$i]}

  echo
  echo
  echo
done