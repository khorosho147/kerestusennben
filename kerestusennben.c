#include <stdio.h>
#include <math.h>
// 線の本数を定義する
#define N 4
#define PI 3.141592653589793238

int FINLIN = N - 1;

// 点の構造体を定義
struct Point
{
    double x;
    double y;
};

// 線の構造体を定義
struct Line
{
    struct Point start; // 始点
    struct Point end;   // 终点
};

// ベクトルの長さを計算する関数を定義
float vector_length(struct Point p1, struct Point p2)
{
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// ベクトルの角度を計算する関数を定義
float vector_angle(struct Point p1, struct Point p2, struct Point p3, struct Point p4)
{
    // ベクトルのドット積を計算する
    float dot_product = (p2.x - p1.x) * (p4.x - p3.x) + (p2.y - p1.y) * (p4.y - p3.y);
    // ベクトルの長さを計算する
    float length1 = vector_length(p1, p2);
    float length2 = vector_length(p3, p4);
    // 角度の弧度値を計算する
    float angle_radians = acos(dot_product / (length1 * length2));
    // 弧度を角度に変換する
    float angle_degrees = angle_radians * 180 / PI;
    // printf("angle:%.13f\n", angle_degrees);
    //  角度の値を返す
    return angle_degrees;
}

// 傾きを計算
float slope_calculation(struct Point p1, struct Point p2)
{
    return (p2.y - p1.y) / (p2.x - p1.x);
}

// 傾きの弧度
float slope_atan(float x)
{
    return atan(x);
}

float dispx;
float dispy;

float distance_px(float d, float angle)
{
    if (angle >= 1.57079636 && angle <= 1.57079638)
    {
        return 0;
    }
    else
    {
        dispx = d * cosf(angle);
        dispx = fabs(dispx);
        return dispx;
    }
}

float distance_py(float d, float angle)
{

    dispy = d * sinf(angle);
    dispy = fabs(dispy);
    return dispy;
}

struct Point midpoint(struct Point p1, struct Point p2)
{
    struct Point mid;
    mid.x = (p1.x + p2.x) / 2.0;
    mid.y = (p1.y + p2.y) / 2.0;
    return mid;
}

int main()
{
    // 線分の端点の座標を定義する
    struct Line line[N];
    line[0] = (struct Line){{6, 6}, {4, 13}};
    line[1] = (struct Line){{6, 15}, {12, 19}};
    line[2] = (struct Line){{18, 18}, {24, 14}};
    line[3] = (struct Line){{24, 10}, {21, 5}};
    line[4] = (struct Line){{30, 50}, {80, 0}};
    line[5] = (struct Line){{120, 0}, {170, 50}};
    line[6] = (struct Line){{110, 30}, {160, 90}};
    line[7] = (struct Line){{190, 50}, {190, 150}};
    line[8] = (struct Line){{110, 170}, {160, 110}};
    line[9] = (struct Line){{140, 180}, {170, 150}};

    //  線分上の点Pの座標を定義する
    struct Point P[N];
    for (int i = 0; i < N; i++)
    {
        P[i] = midpoint(line[i].start, line[i].end);
    }
    // printf("%f,%f\n", P.x, P.y);
    //  線分の外の点Aの座標を定義する
    struct Point A = {10, 6};
    // 線分の外の点Bの座標を定義する
    struct Point B = {16, 6};

    // ベクトルXYの長さを計算する
    float XY_length[N];
    for (int i = 0; i < N; i++)
    {
        XY_length[i] = vector_length(line[i].start, line[i].end);
    }

    // 移動距離を定義
    float Distance[N];

    // 移動方向定義
    short Directions = 2;

    // 傾き計算
    float slope[N];
    for (int i = 0; i < N; i++)
    {
        slope[i] = slope_calculation(line[i].start, line[i].end);
    }

    // printf("slope :  %f\n", slope1);
    // printf("slope :  %f\n", slope2);
    float slope_angle[N];
    for (int i = 0; i < N; i++)
    {
        slope_angle[i] = slope_atan(slope[i]);
    }

    // printf("slope_angle :  %f\n", slope1_angle);
    // printf("slope_angle :  %f\n", slope2_angle);

    // 入射角と反射角を定義する
    float theta_enter[N];
    float theta_out[N];

    float m = 1.0e-13;

    for (int k = 0; k < 100; k++)
    {
        //移動距離を初期化する
        for (int i = 0; i < N; i++)
        {
            Distance[i] = XY_length[i] / 10.0;
        }

        // １本線
        for (;;)
        {
            // 線分XYと線分APの角度を計算する
            theta_enter[0] = vector_angle(line[0].start, line[0].end, A, P[0]);
            // 線分XYと線分BPの角度を計算する
            theta_out[0] = vector_angle(line[0].start, line[0].end, P[0], P[1]);

            // 入射角が大きい時右に移動
            if (theta_enter[0] - theta_out[0] > m)
            {
                // 前の方向が左方向の場合、移動距離を半分にする
                if (Directions == 0)
                {
                    Distance[0] = Distance[0] / 2.0;

                    // 移動距離が限界になる時、for文から抜ける
                    if (Distance[0] <= 1.0e-13)
                    {
                        // 入射角と反射角の誤差範囲内の時、for文から抜ける
                        break;
                    }
                }
                // 移動方向の記録 (右 = 1)
                Directions = 1;

                if (line[0].end.x - line[0].start.x >= 0)
                {
                    // １象限 ++
                    if (line[0].end.y - line[0].start.y >= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[0].x = P[0].x + distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y + distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x >= line[0].end.x && P[0].y >= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x <= line[0].start.x && P[0].y <= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }

                    // 2象限 +-
                    else if (line[0].end.y - line[0].start.y <= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[0].x = P[0].x + distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y - distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x >= line[0].end.x && P[0].y <= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x <= line[0].start.x && P[0].y >= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                }
                else if (line[0].end.x - line[0].start.x <= 0)
                {
                    // 3象限 --
                    if (line[0].end.y - line[0].start.y <= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[0].x = P[0].x - distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y - distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x <= line[0].end.x && P[0].y <= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x >= line[0].start.x && P[0].y >= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                    // 4象限 - +
                    else if (line[0].end.y - line[0].start.y >= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[0].x = P[0].x - distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y + distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x <= line[0].end.x && P[0].y >= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x >= line[0].start.x && P[0].y <= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                }
            }

            // 反射角が大きい時左に移動
            else if (theta_out[0] - theta_enter[0] > m)
            {
                // 前の方向が右方向の場合、移動距離を半分にする
                if (Directions == 1)
                {
                    Distance[0] = Distance[0] / 2.0; // ２から変更

                    // 移動距離が限界になる時、for文から抜ける
                    if (Distance[0] <= 1.0e-13)
                    {

                        // 入射角と反射角の誤差範囲内の時、for文から抜ける
                        break;
                    }
                }
                // 移動方向の記録 (左 = 0)
                Directions = 0;

                if (line[0].end.x - line[0].start.x >= 0)
                {
                    // １象限 --
                    if (line[0].end.y - line[0].start.y >= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[0].x = P[0].x - distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y - distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x >= line[0].end.x && P[0].y >= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x <= line[0].start.x && P[0].y <= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                    // 2象限 -+
                    else if (line[0].end.y - line[0].start.y <= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[0].x = P[0].x - distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y + distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x >= line[0].end.x && P[0].y <= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x <= line[0].start.x && P[0].y >= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                }
                else if (line[0].end.x - line[0].start.x <= 0)
                {
                    // 3象限 ++
                    if (line[0].end.y - line[0].start.y <= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[0].x = P[0].x + distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y + distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x <= line[0].end.x && P[0].y <= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x >= line[0].start.x && P[0].y >= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                    // 4象限 + -
                    else if (line[0].end.y - line[0].start.y >= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[0].x = P[0].x + distance_px(Distance[0], slope_angle[0]);
                        P[0].y = P[0].y - distance_py(Distance[0], slope_angle[0]);

                        // 端点
                        if (P[0].x <= line[0].end.x && P[0].y >= line[0].end.y)
                        {
                            P[0].x = line[0].end.x;
                            P[0].y = line[0].end.y;

                            break;
                        }
                        else if (P[0].x >= line[0].start.x && P[0].y <= line[0].start.y)
                        {
                            P[0].x = line[0].start.x;
                            P[0].y = line[0].start.y;

                            break;
                        }
                    }
                }
            }

            // 入射角と反射角の誤差範囲が10^-13内の時、結果を出力
            else if (fabs(theta_enter[0] - theta_out[0]) <= 1.0e-13 || fabs(theta_out[0] - theta_enter[0]) <= 1.0e-13)
            {

                // 入射角と反射角の誤差範囲内の時、for文から抜ける
                break;
            }
            else
            {
                printf("error");
                break;
            }
        }

        // 中の線
        for (int i = 1; i < N - 1; i++)
        {

            // 2本線
            for (;;)
            {
                // 線分XYと線分APの角度を計算する
                theta_enter[i] = vector_angle(line[i].start, line[i].end, P[i - 1], P[i]);
                // theta1 = 180 - theta1;
                // 線分XYと線分BPの角度を計算する
                theta_out[i] = vector_angle(line[i].start, line[i].end, P[i], P[i + 1]);

                // 入射角が大きい時右に移動
                if (theta_enter[i] - theta_out[i] > m)
                {
                    // 前の方向が左方向の場合、移動距離を半分にする
                    if (Directions == 0)
                    {
                        Distance[i] = Distance[i] / 2.0;

                        if (Distance[i] <= 1.0e-13)
                        {
                            // 入射角と反射角の誤差範囲内の時、for文から抜ける
                            break;
                        }
                    }
                    // 移動方向の記録 (右 = 1)
                    Directions = 1;

                    if (line[i].end.x - line[i].start.x >= 0)
                    {
                        // １象限 ++
                        if (line[i].end.y - line[i].start.y >= 0)
                        {
                            // PのX座標を移動距離分右に移動
                            P[i].x = P[i].x + distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y + distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x >= line[i].end.x && P[i].y >= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x <= line[i].start.x && P[i].y <= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }

                        // 2象限 +-
                        else if (line[i].end.y - line[i].start.y <= 0)
                        {
                            // PのX座標を移動距離分右に移動
                            P[i].x = P[i].x + distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y - distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x >= line[i].end.x && P[i].y <= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x <= line[i].start.x && P[i].y >= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                    }
                    else if (line[i].end.x - line[i].start.x <= 0)
                    {
                        // 3象限 --
                        if (line[i].end.y - line[i].start.y <= 0)
                        {
                            // PのX座標を移動距離分右に移動
                            P[i].x = P[i].x - distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y - distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x <= line[i].end.x && P[i].y <= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x >= line[i].start.x && P[i].y >= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                        // 4象限 - +
                        else if (line[i].end.y - line[i].start.y >= 0)
                        {
                            // PのX座標を移動距離分右に移動
                            P[i].x = P[i].x - distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y + distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x <= line[i].end.x && P[i].y >= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x >= line[i].start.x && P[i].y <= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                    }
                }

                // 反射角が大きい時左に移動
                else if (theta_out[i] - theta_enter[i] > m)
                {
                    // 前の方向が右方向の場合、移動距離を半分にする
                    if (Directions == 1)
                    {
                        Distance[i] = Distance[i] / 2.0; // ２から変更
                        if (Distance[i] <= 1.0e-13)
                        {

                            // 入射角と反射角の誤差範囲内の時、for文から抜ける
                            break;
                        }
                    }
                    // 移動方向の記録 (左 = 0)
                    Directions = 0;

                    if (line[i].end.x - line[i].start.x >= 0)
                    {
                        // １象限 --
                        if (line[i].end.y - line[i].start.y >= 0)
                        {
                            // PのX座標を移動距離分左に移動
                            P[i].x = P[i].x - distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y - distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x >= line[i].end.x && P[i].y >= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;
                            }
                            else if (P[i].x <= line[i].start.x && P[i].y <= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                        // 2象限 -+
                        else if (line[i].end.y - line[i].start.y <= 0)
                        {
                            // PのX座標を移動距離分左に移動
                            P[i].x = P[i].x - distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y + distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x >= line[i].end.x && P[i].y <= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x <= line[i].start.x && P[i].y >= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                    }
                    else if (line[i].end.x - line[i].start.x <= 0)
                    {
                        // 3象限 ++
                        if (line[i].end.y - line[i].start.y <= 0)
                        {
                            // PのX座標を移動距離分左に移動
                            P[i].x = P[i].x + distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y + distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x <= line[i].end.x && P[i].y <= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x >= line[i].start.x && P[i].y >= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                        // 4象限 + -
                        else if (line[i].end.y - line[i].start.y >= 0)
                        {
                            // PのX座標を移動距離分左に移動
                            P[i].x = P[i].x + distance_px(Distance[i], slope_angle[i]);
                            P[i].y = P[i].y - distance_py(Distance[i], slope_angle[i]);

                            if (P[i].x <= line[i].end.x && P[i].y >= line[i].end.y)
                            {
                                P[i].x = line[i].end.x;
                                P[i].y = line[i].end.y;

                                break;
                            }
                            else if (P[i].x >= line[i].start.x && P[i].y <= line[i].start.y)
                            {
                                P[i].x = line[i].start.x;
                                P[i].y = line[i].start.y;

                                break;
                            }
                        }
                    }
                }

                // 入射角と反射角の誤差範囲が10^-13内の時、結果を出力
                else if (fabs(theta_enter[i] - theta_out[i]) <= 1.0e-13 || fabs(theta_out[i] - theta_enter[i]) <= 1.0e-13)
                {

                    // 入射角と反射角の誤差範囲内の時、for文から抜ける
                    break;
                }
                else
                {
                    printf("error");
                    break;
                }
            }
        }

        // 最後の線
        for (;;)
        {
            // 線分XYと線分APの角度を計算する
            theta_enter[FINLIN] = vector_angle(line[FINLIN].start, line[FINLIN].end, P[FINLIN - 1], P[FINLIN]);
            // theta1 = 180 - theta1;
            // 線分XYと線分BPの角度を計算する
            theta_out[FINLIN] = vector_angle(line[FINLIN].start, line[FINLIN].end, P[FINLIN], B);

            // 入射角が大きい時右に移動
            if (theta_enter[FINLIN] - theta_out[FINLIN] > m)
            {
                // 前の方向が左方向の場合、移動距離を半分にする
                if (Directions == 0)
                {
                    Distance[FINLIN] = Distance[FINLIN] / 2.0;

                    if (Distance[FINLIN] <= 1.0e-13)
                    {
                        break;
                    }
                }
                // 移動方向の記録 (右 = 1)
                Directions = 1;

                if (line[FINLIN].end.x - line[FINLIN].start.x >= 0)
                {
                    // １象限 ++
                    if (line[FINLIN].end.y - line[FINLIN].start.y >= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[FINLIN].x = P[FINLIN].x + distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y + distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x >= line[FINLIN].end.x && P[FINLIN].y >= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x <= line[FINLIN].start.x && P[FINLIN].y <= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }

                    // 2象限 +-
                    else if (line[FINLIN].end.y - line[FINLIN].start.y <= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[FINLIN].x = P[FINLIN].x + distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y - distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x >= line[FINLIN].end.x && P[FINLIN].y <= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x <= line[FINLIN].start.x && P[FINLIN].y >= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                }
                else if (line[FINLIN].end.x - line[FINLIN].start.x <= 0)
                {
                    // 3象限 --
                    if (line[FINLIN].end.y - line[FINLIN].start.y <= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[FINLIN].x = P[FINLIN].x - distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y - distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x <= line[FINLIN].end.x && P[FINLIN].y <= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x >= line[FINLIN].start.x && P[FINLIN].y >= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                    // 4象限 - +
                    else if (line[FINLIN].end.y - line[FINLIN].start.y >= 0)
                    {
                        // PのX座標を移動距離分右に移動
                        P[FINLIN].x = P[FINLIN].x - distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y + distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x <= line[FINLIN].end.x && P[FINLIN].y >= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x >= line[FINLIN].start.x && P[FINLIN].y <= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                }
            }

            // 反射角が大きい時左に移動
            else if (theta_out[FINLIN] - theta_enter[FINLIN] > m)
            {
                // 前の方向が右方向の場合、移動距離を半分にする
                if (Directions == 1)
                {
                    Distance[FINLIN] = Distance[FINLIN] / 2.0; // ２から変更
                    if (Distance[FINLIN] <= 1.0e-13)
                    {
                        break;
                    }
                }
                // 移動方向の記録 (左 = 0)
                Directions = 0;

                if (line[FINLIN].end.x - line[FINLIN].start.x >= 0)
                {
                    // １象限 --
                    if (line[FINLIN].end.y - line[FINLIN].start.y >= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[FINLIN].x = P[FINLIN].x - distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y - distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x >= line[FINLIN].end.x && P[FINLIN].y >= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x <= line[FINLIN].start.x && P[FINLIN].y <= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                    // 2象限 -+
                    else if (line[FINLIN].end.y - line[FINLIN].start.y <= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[FINLIN].x = P[FINLIN].x - distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y + distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x >= line[FINLIN].end.x && P[FINLIN].y <= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x <= line[FINLIN].start.x && P[FINLIN].y >= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                }
                else if (line[FINLIN].end.x - line[FINLIN].start.x <= 0)
                {
                    // 3象限 ++
                    if (line[FINLIN].end.y - line[FINLIN].start.y <= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[FINLIN].x = P[FINLIN].x + distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y + distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x <= line[FINLIN].end.x && P[FINLIN].y <= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                            break;
                        }
                        else if (P[FINLIN].x >= line[FINLIN].start.x && P[FINLIN].y >= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                    // 4象限 + -
                    else if (line[FINLIN].end.y - line[FINLIN].start.y >= 0)
                    {
                        // PのX座標を移動距離分左に移動
                        P[FINLIN].x = P[FINLIN].x + distance_px(Distance[FINLIN], slope_angle[FINLIN]);
                        P[FINLIN].y = P[FINLIN].y - distance_py(Distance[FINLIN], slope_angle[FINLIN]);

                        if (P[FINLIN].x <= line[FINLIN].end.x && P[FINLIN].y >= line[FINLIN].end.y)
                        {
                            P[FINLIN].x = line[FINLIN].end.x;
                            P[FINLIN].y = line[FINLIN].end.y;
                        }
                        else if (P[FINLIN].x >= line[FINLIN].start.x && P[FINLIN].y <= line[FINLIN].start.y)
                        {
                            P[FINLIN].x = line[FINLIN].start.x;
                            P[FINLIN].y = line[FINLIN].start.y;
                            break;
                        }
                    }
                }
            }

            // 入射角と反射角の誤差範囲が10^-13内の時、結果を出力
            else if (fabs(theta_enter[FINLIN] - theta_out[FINLIN]) <= 1.0e-13 || fabs(theta_out[FINLIN] - theta_enter[FINLIN]) <= 1.0e-13)
            {
                break;
            }
            else
            {
                printf("error");
                break;
            }
        }

        // 線分XYと線分APの角度を計算する
        theta_enter[0] = vector_angle(line[0].start, line[0].end, A, P[0]);
        // theta_AP1 = 180 - theta_AP1;
        //  線分XYと線分BPの角度を計算する
        theta_out[0] = vector_angle(line[0].start, line[0].end, P[0], P[1]);

        // 中の線
        //  線分XYと線分APの角度を計算する
        for (int i = 1; i < N - 1; i++)
        {
            theta_enter[i] = vector_angle(line[i].start, line[i].end, P[i - 1], P[i]);
            //  線分XYと線分BPの角度を計算する
            theta_out[i] = vector_angle(line[i].start, line[i].end, P[i], P[i + 1]);
        }

        // 四本線

        // 線分XYと線分APの角度を計算する
        theta_enter[FINLIN] = vector_angle(line[FINLIN].start, line[FINLIN].end, P[FINLIN - 1], P[FINLIN]);
        // theta1 = 180 - theta1;
        // 線分XYと線分BPの角度を計算する
        theta_out[FINLIN] = vector_angle(line[FINLIN].start, line[FINLIN].end, P[FINLIN], B);

        float comparison_enter[N];
        float comparison_out[N];

        float hikaku = 0; // 前回のデータと比較して、変化がないならbreak；

        // 角度と前回一致する時、結界出力
        for (int i = 0; i < N; i++)
        {
            hikaku = hikaku + (comparison_enter[i] - theta_enter[i]) + (comparison_out[i] - theta_out[i]);
        }

        if (hikaku == 0)
        {
            //  結果を出力する
            for (int i = 0; i < N; i++)
            {
                printf("%d本線\n", i + 1);
                printf("P%d:  (%f,%f)\n", i + 1, P[i].x, P[i].y);
                printf("入射角: %.13f°\n", theta_enter[i]);
                printf("反射角: %.13f°\n", theta_out[i]);
            }

            // X座標を出力
            printf("\nX座標を出力\n");
            printf("%f\n", P[0].x);
            printf("%f\n\n", A.x);
            printf("%f\n", line[0].start.x);
            printf("%f\n\n", line[0].end.x);

            for (int i = 0; i < N - 1; i++)
            {
                printf("%f\n", P[i].x);
                printf("%f\n\n", P[i + 1].x);
                printf("%f\n", line[i + 1].start.x);
                printf("%f\n\n", line[i + 1].end.x);
            }
            printf("%f\n", P[N - 1].x);
            printf("%f\n", B.x);

            // Y座標を出力
            printf("\nY座標を出力\n");
            printf("%f\n", P[0].y);
            printf("%f\n\n", A.y);
            printf("%f\n", line[0].start.y);
            printf("%f\n\n", line[0].end.y);

            for (int i = 0; i < N - 1; i++)
            {
                printf("%f\n", P[i].y);
                printf("%f\n\n", P[i + 1].y);
                printf("%f\n", line[i + 1].start.y);
                printf("%f\n\n", line[i + 1].end.y);
            }
            printf("%f\n", P[N - 1].y);
            printf("%f\n\n", B.y);

            break;
        }
        // 前回のデータを格納する
        for (int i = 0; i < N; i++)
        {
            comparison_enter[i] = theta_enter[i];
            comparison_out[i] = theta_out[i];
        }
    }
    return 0;
}