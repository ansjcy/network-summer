//
//  demo.cpp
//  coding
//
//  Created by Anakin on 16/8/9.
//  Copyright © 2016年 Anakin. All rights reserved.
//
#include <map>
#include <iostream>
#include <stdio.h>
using namespace std;
int main()
{
    map<int, int> aa;
    aa[0] = 1;
    aa[10] = 2;
    map<int, int> &bb = aa;
    bb[10] = 1;
    cout << aa[10];
}
