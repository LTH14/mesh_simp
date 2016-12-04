//
//  main.cpp
//  mesh_simp
//
//  Created by litianhong on 2016/11/26.
//  Copyright © 2016年 litianhong. All rights reserved.
//

#include <iostream>
#include <cstring>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include "SimpleObject.h"
#include "Vec3f.h"
using namespace cv;
using namespace std;
using namespace SimpleOBJ;

int main() {
    // insert code here...
    CSimpleObject obj;
    float r;
    char fn[100] = "/Users/litianhong/Desktop/C++/mesh_simp/mesh_simp/";
    char input[50];
    cout << "Please type the name of file to contract" << endl;
    scanf("%s", input);
    cout << "Please type the expected ratio of remaining phases" << endl;
    cin >> r;
    strcat(fn, input);
    obj.LoadFromObj(fn);
    int n = obj.m_nowTriangles;
    int i = 0;
    while (obj.m_nowTriangles > r*n)
    {
        i++;
        obj.contract();
    }
    char saveplace[100] = "/Users/litianhong/Desktop/C++/mesh_simp/";
    cout << "max edges" << obj.max_edge_num << endl;
    cout << "max faces" << obj.max_face_num << endl;
    strcat(saveplace, input);
    obj.SaveToObj(saveplace);
    exit(0);
}
