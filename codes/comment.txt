R = 0.3;%半径
L = 1.0;%矩形边长一半

Point(1) = {L, -L, 0};%右下点
Point(2) = {L, L, 0};%右上点
Point(3) = {-L, L, 0};%左上点
Point(4) = {-L, -L, 0};%左下点
Point(5) = {-L + R, -L, 0};%圆弧起点
Point(6) = {-L, -L + R, 0};%圆弧终点
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};%圆弧45度点

Circle(1) = {5, 4, 7};%圆心4，连接5、7的圆弧
Circle(2) = {7, 4, 6};%圆心4，连接7、6的圆弧

Line(3) = {6, 3};%连接6、3的直线
Line(4) = {3, 2};%连接3、2的直线
Line(5) = {2, 1};%连接2、1的直线
Line(6) = {1, 5};%连接1、5的直线
Line(7) = {2, 7};%连接2、7的直线

Curve Loop(1) = {4, 7, 2, 3};%由圆弧和直线组成的曲线环1
Plane Surface(1) = {1};%定义平面1

Curve Loop(2) = {7, -1, -6, -5};%由圆弧和直线组成的曲线环2
Plane Surface(2) = {2};%定义平面2

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;%将该线段的网格划分为123、345、567
Transfinite Surface{1};%使面1网格均匀分布
Transfinite Surface{2};%使面2网格均匀分布

Recombine Surface{1};%重组面1网格
Recombine Surface{2};%重组面2网格

Mesh.ElementOrder = 1;%网格单元阶数为1
Mesh.Algorithm = 8;%设置网格算法为高质量三角形生成算法

// EOF
