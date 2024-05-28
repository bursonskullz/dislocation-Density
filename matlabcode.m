% DATE: 4-10-2021
% TIME: 13:52
% PURPOSE: OBTAIN VISUAL GRAPHICS FOR SHOCHWAVE REPORT

strain = load('strain_data');
alphaData = load('alpha_data');
thetaData = load('theta_data');
DData = load('D_data');
BData = load('B_data');
fData = load('forces_data');
coupleStress = load('coupleStress_data');
stressData = load('stress_data');
positiona = load('positionVector_data');

%{
figure(1)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), alphaData(:,3), 32, alphaData(:,3), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = '\alpha_{13} (m^{-1})';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(2)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), thetaData(:,9), 32, thetaData(:,9), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = '\theta_{33} (m^{-2})';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(3)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), coupleStress(:,7), 42, coupleStress(:,7), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'M_{31} (m^{-1})';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}


figure(4)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), coupleStress(:,3), 42, coupleStress(:,3), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'M_{31} (m^{-1})';
set(colorTitleHandle ,'String',titleString,'FontSize',18);


%{
figure(5)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), stressData(:,2), 42, stressData(:,2), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'T_{12} (\mu)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(6)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), stressData(:,1), 42, stressData(:,1), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'T_{11} (\mu)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(7)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), DData(:,1), 42, DData(:,1), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'D_{1231}  (\mu \vert k)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(8)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), BData(:,1), 42, BData(:,1), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'B_{3111} (\mu \vert k)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}

%{
figure(9)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), fData(:,2), 42, fData(:,2), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'F_{2}^\theta (Mpa/ m)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}


%{
figure(10)
colormap(jet);% Set colour map
scatter3(positiona(:,1), positiona(:,2), fData(:,3), 22, fData(:,3), 'filled');
xlabel('x_{1}','FontSize',22);
ylabel('x_{2}','FontSize',22);
view(0,90);
hcb= colorbar; % Throw in the colour bar.
colorTitleHandle = get(hcb,'Title');
titleString = 'F_{2}^\alpha (Mpa/ m)';
set(colorTitleHandle ,'String',titleString,'FontSize',18);
%}
