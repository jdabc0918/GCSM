// stdafx.h : �W���̃V�X�e�� �C���N���[�h �t�@�C���̃C���N���[�h �t�@�C���A�܂���
// �Q�Ɖ񐔂������A�����܂�ύX����Ȃ��A�v���W�F�N�g��p�̃C���N���[�h �t�@�C��
// ���L�q���܂��B
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
//Release���[�h�̏ꍇ
#pragma comment(lib,"opencv_world300.lib")
#endif