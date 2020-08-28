# Finding-Significant-Items

## Introduction

LTC is a succinct data structure to find top significant or frequent or persistent items in a high speed data stream. It is first published in ICDE 2019, and the Journal version of this paper is under review. This repository is the code for the Journal version paper.

    Yang, Tong, et al. "Finding Significant Items in Data Streams." 2019 IEEE 35th International Conference on Data Engineering (ICDE). IEEE, 2019.

## Source Code.

LTC.cpp contains base LTC algorithm as well as the optimizations. It is cleaned up and the comments and variable names in it should make it understandable. In addition, it has the functions used to run expeirments in the Journal paper. Users can refer to the main() function to run experiments. 

## Dataset

- stack-new.txt: It contains 10,000,000 flows from the stackexchange website. Each flow contains the item ID and the arriving time. We commend setting the number of periods to 1000, and the memory size to more than 20KB.
- ddos.txt: It is converted from the a Kaggle dataset [here](https://www.kaggle.com/devendra416/ddos-datasets).

Since the maximal file size in the Github is 100MB, we have to compress these files.

## Contact

Everyone is welcomed to open issues in this repo!
