delta = pi / 180;  % Space between elements of antenna
theta = -pi / 2 : delta : pi / 2;  % Array of theta values
target_DoA = [15 * pi / 180];  % Target Direction of Arrival (radian)
interference_source = [-50 * pi / 180, -10 * pi / 180, 40 * pi / 180];  % Interference sources (radian)

% Chuyển giá trị DoA và interference source sang độ
theta_degrees = theta * 180 / pi;  % Chuyển theta sang độ
target_DoA_degrees = target_DoA * 180 / pi;  % Chuyển target_DoA sang độ
interference_source_degrees = interference_source * 180 / pi;  % Chuyển interference_source sang độ

% Tìm index gần nhất trong mảng theta_degrees với target_DoA_degrees
[~, target_index] = min(abs(theta_degrees - target_DoA_degrees));  % Tìm index gần nhất với target_DoA

% Tìm index gần nhất cho từng interference source
interference_indexes = zeros(size(interference_source_degrees));  % Khởi tạo mảng chứa các index
for i = 1:length(interference_source_degrees)
    [~, interference_indexes(i)] = min(abs(theta_degrees - interference_source_degrees(i)));  % Tìm index gần nhất cho mỗi interference source
end

% Hiển thị kết quả
disp('Index của target_DoA:');
disp(target_index);

disp('Index của interference_source:');
disp(interference_indexes);
