clc; clear; close all;

% 이미지가 저장된 폴더 경로
imageFolder = './Plots';
outputVideoName = 'V_VelocityAnimation.mp4';

% 파일 이름 형식에 맞는 모든 이미지 불러오기
imageFiles = dir(fullfile(imageFolder, 'V_Velocity_*.jpg'));

% 파일 이름을 정렬 (숫자가 포함된 경우에도 올바르게)
imageFileNames = {imageFiles.name}; 
imageFileNames = natsortfiles(imageFileNames); % 자연스럽게 정렬

% 정렬된 순서대로 파일 목록을 다시 불러옴
sortedFiles = cellfun(@(x) fullfile(imageFolder, x), imageFileNames, 'UniformOutput', false);

% 비디오 저장을 위한 객체 생성 (프레임 속도 10fps)
videoObj = VideoWriter(outputVideoName, 'MPEG-4');
videoObj.FrameRate = 10; % 초당 10 프레임
open(videoObj);

% 모든 이미지 파일을 읽어서 비디오에 추가
for i = 1:length(sortedFiles)
    img = imread(sortedFiles{i}); % 이미지 읽기
    writeVideo(videoObj, img); % 비디오에 추가
end

% 비디오 저장 완료 후 닫기
close(videoObj);

fprintf('✅ 비디오 저장 완료: %s\n', outputVideoName);