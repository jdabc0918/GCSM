// stdafx.h : 標準のシステム インクルード ファイルのインクルード ファイル、または
// 参照回数が多く、かつあまり変更されない、プロジェクト専用のインクルード ファイル
// を記述します。
//

#pragma once
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include "targetver.h"
#include <stdio.h>
#include <tchar.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <opencv2/opencv.hpp>
//#include <opencv2/nonfree/nonfree.hpp>
using namespace std;
using namespace cv;

#ifdef _DEBUG
#pragma comment(lib,"opencv_world300d.lib")
#else
//Releaseモードの場合
#pragma comment(lib,"opencv_world300.lib")
#endif