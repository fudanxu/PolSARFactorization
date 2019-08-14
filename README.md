# Polarimetric SAR Image Factorization 
The folder stores the implementation process of Polarimetric SAR Image Factorization. 

If you use this code, please refer to below papers:
### Polarimetric SAR Image Factorization
[1]  Feng Xu, Qian Song, and Ya-Qiu Jin,  "Polarimetric SAR Image Factorization",  IEEE Transactions on Geoscience and Remote Sensing, Vol. 55, No. 9, pp. 5026-5041, 2017.  
### NMF optimization using graph regularized non-negative matrix factorization with divergence formulation locality preserving.
[1]  Deng Cai, Xiaofei He, Xuanhui Wang, Hujun Bao, and Jiawei Han. "Locality Preserving Nonnegative Matrix Factorization", Proc. 2009 Int. Joint Conf. on Arti_cial Intelligence (IJCAI-09), Pasadena, CA, July 2009. 
[2]  Deng Cai, Xiaofei He, Jiawei Han, Thomas Huang. "Graph Regularized Non-negative Matrix Factorization for Data Representation", IEEE Transactions on Pattern Analysis and Machine Intelligence, , Vol. 33, No. 8, pp. 1548-1560, 2011.  

## Requirements
- Matlab

## Datatsets
The datasets adopted in our paper "Polarimetric SAR Image Factorization" are collected by UAVSAR, which can be downloaded by the website 
https://uavsar.jpl.nasa.gov.
- C_TestData1.mat  The covariance matrix of a UAVSAR PolSAR image with size 300*300.
- C_TestData1.mat  The covariance matrix of a UAVSAR PolSAR image with size 900*900.

## Acknowledgement 
Implementation of non-negative matrix factorization is borrowed from Cai Deng's code (http://www.cad.zju.edu.cn/home/dengcai/Data/data.html).
