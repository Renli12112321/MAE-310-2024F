convq为人工解误差计算器。
运行convq得到disp.mat, 可运行vis_2d_sol得到可视化云图。
quarter-plate-with-hole-quad1.geo 存储了网格生成和边界设置代码。
quarter_plate_with_hole_quad1.m和quarter_plate_with_hole_quad2.m为两种matlab网格，存在bug，一个无QUADS，一个无lines第三列。
driver1.m为第二三问的尝试代码，尚未完成，只有初始设定、坐标变换、精确解、网格导入、IEN等内容。