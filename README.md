# DexterousManipulation

This code is implementation of paper : 
Seunghwan Lee, Ri Yu, Jungnam Park, Mridul Aanjaneya, Eftychios Sifakis and Jehee Lee, 
Dexterous Manipulation and Control with Volumetric Muscles, 
ACM Transactions on Graphics (SIGGRAPH 2018), Volume 37, Issue4, August 2018, Article No. 57. 
https://dl.acm.org/citation.cfm?id=3201330

If you are interested in this implementation, please refer to the paper.

# Prerequisite


DART version 6 http://dartsim.github.io/ 

Eigen version 3

IPOPT https://projects.coin-or.org/Ipopt

Boost

As for the Projective Dynamics, I refers to the code of Tiantian Liu(http://tiantianliu.cn/).

# How to build

```console
mkdir build 
cd build
cmake ..
```

If you want to simulate without volumetric muscles, you can disable USE_MUSCLE in CMakeLists.txt(or in ccmake)

Unfortunately, I'm not managing this code anymore, but I'll be happy to answer your questions for compiling, linking, or details of implementation. Feel free to ask. 

Email : lsw9021@gmail.com
