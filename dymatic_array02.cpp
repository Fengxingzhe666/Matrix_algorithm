// dymatic_array02.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// 主程序使用和测试matrix.h

#include <iostream>
#include "matrix.h"
#include <fstream>

int main()
{
    using namespace std;
    Matrix a,b,r;
    double t[4][4] = { {1,0,0,0},{1,2,0,0},{2,1,3,0},{1,2,1,4} };
    double t1[1][1] = { 2 };
    a.set(t);
    a.disp();
    //a.set({ {1,2},{3,4} });
    b.set({ { 2,2,2 }, { 3,3,3 } });
    r = a.accompany();//伴随矩阵
    r.disp();
    r = inv(a);//逆矩阵
    r.disp();
    r = inv(a) * a;//重载*运算符，实现矩阵相乘。逆矩阵*本体=单位矩阵（对角全1）
    r.disp();
    r = a + b;
    r.disp();
}
