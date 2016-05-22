// GCSaliencyMap.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"

const cv::Vec3b colorPalette[] = {
    cv::Vec3b(0,0,255),
    cv::Vec3b(0, 255,0),
    cv::Vec3b(255,0,0),
    cv::Vec3b(255, 0, 255),
    cv::Vec3b(0, 255, 255),
    cv::Vec3b(255,255,0),
    cv::Vec3b(128, 0, 255),
    cv::Vec3b(0, 255, 128),
    cv::Vec3b(255, 128, 0),
    cv::Vec3b(255, 0, 128),
    cv::Vec3b(0, 128, 255),
    cv::Vec3b(128, 255, 0),
    cv::Vec3b(0, 0, 128),
    cv::Vec3b(0, 128, 0),
    cv::Vec3b(128, 0, 0)
};

void createMask(cv::Mat img, cv::Mat &dst)
{
    cv::Vec3b targetColor(36, 28, 237);
    
    dst = cv::Mat::zeros(img.size(), CV_8UC1);
    for (int i = 0; i < img.cols* img.rows; i++)
    {
        if (img.at<cv::Vec3b>(i) == targetColor)
            dst.at<uchar>(i) = 255;
    }
}

void coloringLabels(cv::Mat labels, cv::Mat &dst)
{
    cv::Mat ret = cv::Mat::zeros(labels.size(), CV_8UC3);

    for (int i = 0; i < ret.cols*ret.rows; i++)
    {
        int lab = labels.at<int>(i);

        cv::Vec3b c = cv::Vec3b(0, 0, 0);
        if (lab >= 0)  c = colorPalette[ lab % 15];

        ret.at<cv::Vec3b>(i) = c;
    }

    dst = ret.clone();
}

int countColor(cv::Mat img, std::vector<std::pair<cv::Vec3b, int>> &colorList, cv::Mat mask = cv::Mat())
{
    if (mask.type() != CV_8UC3)
    {
        cv::cvtColor(mask, mask, cv::COLOR_GRAY2RGB);
    }
    
    cv::Mat Img = img & mask;
    int sz = Img.rows * Img.cols;
    colorList.clear();
    int color_num = 0;

    for (int i = 0; i < sz; i++)
    {
        if (mask.at<cv::Vec3b>(i) == cv::Vec3b(255,255,255))
        {
            cv::Vec3b c = img.at<cv::Vec3b>(i);
            color_num = (int)colorList.size();

            if (color_num == 0)
            {
                std::pair<cv::Vec3b, int> obj = std::make_pair(c, 1);
                colorList.push_back(obj);
            }
            else
            {
                for (int j = 0; j < color_num; j++)
                {
                    cv::Vec3b cj = colorList[j].first;
                    if (c == cj)
                    {
                        colorList[j].second++;
                        break;
                    }
                    else if (j == color_num - 1)
                    {
                        std::pair<cv::Vec3b, int> obj = std::make_pair(c, 1);
                        colorList.push_back(obj);
                    }
                }
            }
        }
    }

    return (int)colorList.size();
}

double colorDist(cv::Vec3b c1, cv::Vec3b c2)
{
    double ret = 0.0;
    for (int i = 0; i < 3; i++) ret += pow((double)c1[i] - (double)c2[i], 2);
    return sqrt(ret);
}

void RegionBaseContrast(cv::Mat img, cv::Mat label_image, int label_num)
{
    //  ！ラベル番号は連番前提
    
    int height = img.rows;
    int width = img.cols;
    int sz = height*width;

    //  ラベルごとにマスクを作る
    std::vector<cv::Mat> Masks(label_num, cv::Mat());
    for (int i = 0; i < label_num; i++)
    {
        cv::Mat msk = cv::Mat::zeros(label_image.size(), CV_8UC3);
        for (int j = 0; j < sz; j++)
        {
            if (label_image.at<int>(j) == i) msk.at<cv::Vec3b>(j) = cv::Vec3b(255,255,255);
        }
        Masks[i] = msk;
    }

    //  重みwr(面積)の計算
    std::vector<int> wr(label_num, 0);
    for (int i = 1; i < label_num; i++)
    {
        cv::Mat msk = Masks[i];
        for (int j = 0; j < sz; j++)
        {
            if (msk.at<cv::Vec3b>(j) == cv::Vec3b(255,255,255)) wr[i]++;
        }
    }

    //  重みwsr(画像中心からの距離)の計算
    std::vector<double> wsr(label_num, 0.0);
    cv::Point2d center(0.5, 0.5);
    for (int i = 1; i < label_num; i++)
    {
        cv::Mat msk = Masks[i];
        double di = 0.0;
        for (int j = 0; j < sz; j++)
        {
            if (msk.at<cv::Vec3b>(j) == cv::Vec3b(255, 255, 255))
            {
                double dx = (double)(j % width) / width - center.x;
                double dy = (double)(j / width) / height - center.y;
                di += sqrt(dx*dx + dy*dy);
            }
        }
        di /= wr[i];
        wsr[i] = exp(-9 * di);
    }
    
    //  Ds
    cv::Mat Ds = cv::Mat::zeros(cv::Size(label_num, label_num), CV_64FC1);
    for (int i = 1; i < label_num; i++) //i=0は真っ黒なので無視
    {
        cv::Mat ri_msk = Masks[i];
        cv::Point2d center_i(0.0, 0.0);
        for (int k = 0; k < sz; k++)
        {
            if(ri_msk.at<cv::Vec3b>(k) == cv::Vec3b(255,255,255))
                center_i += cv::Point2d(k % width, k / width);
        }
        center_i /= wr[i];

        for (int j = i + 1; j < label_num; j++)
        {
            cv::Mat rj_msk = Masks[j];
            cv::Point2d center_j(0.0, 0.0);
            for (int k = 0; k < sz; k++)
            {
                if (rj_msk.at<cv::Vec3b>(k) == cv::Vec3b(255, 255, 255))
                    center_j += cv::Point2d(k % width, k / width);
            }
            center_j /= wr[j];

            double dx = (center_i.x - center_j.x) / width;
            double dy = (center_i.y - center_j.y) / height;
            Ds.at<double>(i, j) = sqrt(dx*dx + dy*dy);
            std::cout << "Ds(r" << i << ", r" << j << ") = " << Ds.at<double>(i, j) << std::endl;
        }
    }
    
    //  Dr
    cv::Mat Dr = cv::Mat::zeros(cv::Size(label_num, label_num), CV_64FC1);
    for (int i = 1; i < label_num; i++) //i=0は真っ黒なので無視
    {
        //  一つ目の領域について
        cv::Mat ri_msk = Masks[i];
        int wri = wr[i];
        cv::Mat ri_img = img & ri_msk;      
        std::vector<std::pair<cv::Vec3b, int>> ri_colors;
        int ni = countColor(img, ri_colors, ri_msk);
        cv::imshow("ri", ri_img);
        //  二つ目の領域について
        for (int j = i+1; j < label_num; j++)
        {
            std::cout << "Dr(r" << i << ", r" << j << ")を計算...";
            cv::Mat rj_msk = Masks[j];
            int wrj = wr[j];
            cv::Mat rj_img = img & rj_msk;
            std::vector<std::pair<cv::Vec3b, int>> rj_colors;
            int nj = countColor(img, rj_colors, rj_msk);
            cv::imshow("rj", rj_img);
            cv::waitKey(1);
            for (int k = 0; k < ni; k++)
            {
                cv::Vec3b ci = ri_colors[k].first;
                double fci = (double)ri_colors[k].second / wri;
                for (int m = 0; m < nj; m++)
                {
                    cv::Vec3b cj = rj_colors[m].first;
                    double fcj = (double)rj_colors[m].second / wrj;
                    double dist = colorDist(ci, cj);
                    Dr.at<double>(i, j) += fci * fcj * dist;
                }
            }
            std::cout << Dr.at<double>(i, j) << std::endl;
        }
    }

    //  Sr
    double varS = 0.4;
    std::vector<double> Sr(label_num, 0.0);
    for (int i = 1; i < label_num; i++) //i=0は真っ黒なので無視
    {
        std::cout << "Sr" << i << "を計算...";
        for (int j = 1; j < label_num; j++)
        {
            if (i != j)
            {
                int I = std::min(i, j);
                int J = std::max(i, j);
                double wrj = wr[j];
                double sw = exp(Ds.at<double>(I, J) / (-varS));
                Sr[i] += sw * wrj * Dr.at<double>(I, J);
            }
        }
        Sr[i] *= wsr[i];
        std::cout << Sr[i] << std::endl;
    }

    //  可視化
    double max_Sr = -1.0;
    for (int i = 1; i < label_num; i++)
    {
        max_Sr = std::max(max_Sr, Sr[i]);
    }
    cv::Mat SrMap = cv::Mat::zeros(img.size(), CV_8UC1);
    for (int i = 1; i < label_num; i++)
    {
        cv::Mat msk = Masks[i];
        double ratio = Sr[i] / max_Sr;
        for (int j = 0; j < sz; j++)
        {
            if(msk.at<cv::Vec3b>(j) == cv::Vec3b(255,255,255)) SrMap.at<uchar>(j) = (uchar)(255 * ratio);
        }
    }
    cv::imshow("SrMap", SrMap);
    cv::waitKey(1);
}

int _tmain(int argc, _TCHAR* argv[])
{
    cv::Mat img = cv::imread("image.png");
    cv::Mat marked = cv::imread("marked.png");
    cv::resize(img, img, cv::Size(), 0.5, 0.5, cv::INTER_NEAREST);
    cv::resize(marked, marked, cv::Size(), 0.5, 0.5, cv::INTER_NEAREST);

    //  マスク作成
    cv::Mat mask;
    createMask(marked, mask);

    //  ラベリング
    cv::Mat labels;
    int nLabels = cv::connectedComponents(mask, labels);

    //  watershed
    cv::watershed(img, labels);

    RegionBaseContrast(img, labels, nLabels);

    //  カラーリング
    coloringLabels(labels, labels);

    cv::imshow("image", img);
    cv::imshow("result", labels);
    cv::waitKey();

	return 0;
}

