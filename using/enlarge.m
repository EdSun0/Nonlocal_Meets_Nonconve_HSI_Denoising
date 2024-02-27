% 打开.fig文件
fig = openfig('our.fig', 'reuse'); % 替换 'our.fig' 为更改的fig文件
ax_original = gca; % 获取当前图形中的坐标轴

% 保持当前轴，以便添加图形
hold(ax_original, 'on');

% 定义放大区域
rect = [45 90 20 20]; % 替换为你想放大的区域 [x y width height]

% 在原图上标记放大区域
rectangle(ax_original, 'Position', rect, 'EdgeColor', 'r', 'LineWidth', 2);

% 提取放大区域
magnifiedArea = imcrop(getimage(ax_original), rect);
magnifiedArea = imresize(magnifiedArea, 2); % 2倍放大

% 创建新的轴显示放大区域
fig.Position = [100 100 800 600];

% 定义新轴的位置（右上角）
newAxesPosition = [0.52, 0.667, 0.2, 0.2]; % [left bottom width height]
newAxes = axes('Position', newAxesPosition);

% 显示放大区域
imshow(magnifiedArea, 'Parent', newAxes);

% 标记新图
rectangle('Position', [0, 0, size(magnifiedArea, 1), size(magnifiedArea, 1)], 'EdgeColor', 'r', 'LineWidth', 2, 'Parent', newAxes);