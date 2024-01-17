% ��.fig�ļ�
fig = openfig('our.fig', 'reuse'); % �滻 'our.fig' Ϊ���ĵ�fig�ļ�
ax_original = gca; % ��ȡ��ǰͼ���е�������

% ���ֵ�ǰ�ᣬ�Ա����ͼ��
hold(ax_original, 'on');

% ����Ŵ�����
rect = [45 90 20 20]; % �滻Ϊ����Ŵ������ [x y width height]

% ��ԭͼ�ϱ�ǷŴ�����
rectangle(ax_original, 'Position', rect, 'EdgeColor', 'r', 'LineWidth', 2);

% ��ȡ�Ŵ�����
magnifiedArea = imcrop(getimage(ax_original), rect);
magnifiedArea = imresize(magnifiedArea, 2); % 2���Ŵ�

% �����µ�����ʾ�Ŵ�����
fig.Position = [100 100 800 600];

% ���������λ�ã����Ͻǣ�
newAxesPosition = [0.52, 0.667, 0.2, 0.2]; % [left bottom width height]
newAxes = axes('Position', newAxesPosition);

% ��ʾ�Ŵ�����
imshow(magnifiedArea, 'Parent', newAxes);

% �����ͼ
rectangle('Position', [0, 0, size(magnifiedArea, 1), size(magnifiedArea, 1)], 'EdgeColor', 'r', 'LineWidth', 2, 'Parent', newAxes);