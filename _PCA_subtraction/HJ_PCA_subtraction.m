%% Video Initialization (removing PCA background)
clear; clc; close all;
Msg = 'Select a video file to analyze';
[FileName, PathName] = uigetfile({'*.avi','Video Files(*.avi)';...
                                  '*.*','All Files(*.*)'}, Msg);
cd(PathName)
vidpath = fullfile(PathName, FileName);
do_savevid = 1;
resize_ratio = 0.25;
vid_pca_removed = PCA_bkgsubtraction(vidpath, resize_ratio, do_savevid);
