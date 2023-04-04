clear;clc; close all;
addpath('../Digit_data/0331/test01')
file1 =  'forwardKine.csv';
file2 =  'imuReading.csv';
file3 =  'jointPosFile.csv';
file4 =  'groundTruthFile.csv';
file5 =  'deltaTime.csv';
file6 =  'imuOrient.csv';
file7 =  'contactFile.csv';
file8 =  'baseOrient.csv';
file9 =  'trueFootFile.csv';
file10 = 'timerFile.csv';

leftFoot = readmatrix(file1,"Range",'A:C');
rightFoot = readmatrix(file1,"Range",'D:F');
world_leftFoot = readmatrix(file9,"Range",'A:C');
world_rightFoot = readmatrix(file9,"Range",'D:F');

iumReading = readmatrix(file2);
linAcc = iumReading(:,1:3);
jointPos = readmatrix(file3);
% jointPos(:,11) = [];
truePos = readmatrix(file4,'Range','A:C');
trueVel =  readmatrix(file4,'Range','D:F');
deltaTime = csvread(file5);
imuOrient = readmatrix(file6);
baseOrient = readmatrix(file8);
leftContact = readmatrix(file7,"Range",'A:A');
rightContact = readmatrix(file7,"Range",'B:B');
time = readmatrix(file10);
% predictState = readmatrix(file10);
torso2IMU = -pi/2;
R_torso2IMU = [cos(torso2IMU) 0 sin(torso2IMU);
               0              1 0;
              -sin(torso2IMU) 0 cos(torso2IMU)];

iumReading = (R_torso2IMU * linAcc')';
Q_contact = [0.01, 0 , 0;
             0, 0.01, 0;
             0, 0, 0.01 ];
Q_noContact = [1e10, 0 , 0;
                0, 1e10, 0;
                0, 0, 1e10 ];
q0 = baseOrient(1,1); q1 = baseOrient(1,2); 
q2 = baseOrient(1,3); q3 = baseOrient(1,4);
simLen = size(deltaTime,1);
store = zeros(simLen,12);


Qrotation = [2*(q0^2+q1^2)-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2); 
            2*(q1*q2+q0*q3), 2*(q0^2+q2^2)-1, 2*(q2*q3-q0*q1);
            2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(q0^2+q3^2)-1];
rotationDegX = 0* pi()/180; 
rotationDegY = -90 * pi()/180;
rotationDegZ = 0*pi()/180;

rotationMatrixX =[ 1, 0, 0;
                   0, cos(rotationDegX), -sin(rotationDegX);
                    0, sin(rotationDegX), cos(rotationDegX)];
rotationMatrixY =[ cos(rotationDegY), 0, sin(rotationDegY);
                    0,               1,         0;
                 -sin(rotationDegY), 0, cos(rotationDegY)];
rotationMatrixZ =[  cos(rotationDegZ), -sin(rotationDegZ), 0;
                   sin(rotationDegZ), cos(rotationDegZ), 0;
                    0, 0, 1];

%% configure KF-
p0 = 100*eye(12); 
x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0].'; % initial guess
R = 1*eye(6);
leftFoot_world = zeros(simLen,3);
rightFoot_world = zeros(simLen,3);
% world_leftFoot = (Qrotation * world_leftFoot')';

imu_world = zeros(simLen,3);
for i = 1:simLen
    Q = 0.01*eye(12);
%     Q = Qk_val;
    dt = deltaTime(i);
    A = eye(12);
    time_matrix = [dt, 0 , 0; 
                    0, dt, 0;
                    0, 0, dt];
    A(1:3,4:6) = time_matrix;
    q_vec = baseOrient(i,:)';
    Qrotation = quaternion2RotM(q_vec);
    imuAcc =  Qrotation * iumReading(i,:)';
    lFoot =  Qrotation * leftFoot(i,:)';
    rFoot =  Qrotation * rightFoot(i,:)';
    leftFoot_world(i,:) = lFoot;
    rightFoot_world(i,:) = rFoot;
    imu_world(i,:) = imuAcc;
    B = [0, 0, 0, imuAcc(1)*dt, imuAcc(2)*dt, (imuAcc(3)-9.81)*dt, 0, 0, 0, 0, 0, 0].';
    z = [lFoot(1), lFoot(2), lFoot(3), rFoot(1), rFoot(2), rFoot(3)].';
    H = zeros(6,12);
    H_M1 = -1*eye(3);
    H_M2 = eye(6);
    H(:,1:3) = [H_M1, H_M1].';
    H(:,7:12) = H_M2;
    if (rightContact(i)== 0)
        Q(10:12,10:12) = Q_noContact;
    elseif (rightContact(i)==1)
        Q(10:12,10:12) = Q_contact;
    end
    if (leftContact(i)==0)
        Q(7:9,7:9) = Q_noContact;
    elseif(leftContact(i)==1)
        Q(7:9,7:9) = Q_contact;
    end
    x = A*x0 + B;
    p = A*p0*transpose(A) + Q;
%     normP(i,1) = norm(p(1:6,1:6));
    y = z - H*x;
    s = H*p*transpose(H)+ R;
    k = p *transpose(H)*s^-1;
    x0 = x + k*y;
    store(i,:)= x0.';
    p0 = (eye(12) - k*H)*p;
    y = z - H*x0;
end

%% debug plot
close all
figure 
plot(imu_world(:,1));
hold on
plot(imu_world(:,2));
plot(imu_world(:,3));
legend('x','y','z')
%% plot
figure
subplot(2,1,1);
plot(time,truePos(:,1));
hold on 
plot(time,truePos(:,2));
plot(time,truePos(:,3));
legend('x','y','z','Location','northwest')
title('base ground truth position')
xlabel('time (s)')
ylabel('distance (meters)')
grid on
subplot(2,1,2);
plot(time, store(:,1));
hold on
plot(time, store(:,2));
plot(time, store(:,3));
legend('x','y','z','Location','northwest')
title('base estimated position')
xlabel('time (s)')
ylabel('distance (meters)')
grid on

figure
subplot(2,2,1)
plot(time, truePos(:,1));
hold on 
plot(time, store(:,1));
legend('true','estimated','Location','northwest')
title('base position')
xlabel('time (s)')
ylabel('x distance (meters)')
grid on
ylim([-0.5 3])

subplot(2,2,2)
plot(time, truePos(:,2));
hold on 
plot(time, store(:,2));
legend('true','estimated','Location','northwest')
title('base position')
xlabel('time (s)')
ylabel('y distance (meters)')
grid on
ylim([-0.5 3])

subplot(2,2,[3, 4])
plot(time, truePos(:,3));
hold on 
plot(time, store(:,3));
legend('true','estimated','Location','northwest')
title('base position')
xlabel('time (s)')
ylabel('z distance (meters)')
grid on
ylim([-1 2])

%% foot plot
figure
subplot(2,2,1)
plot(time, world_leftFoot(:,1), 'Linewidth', 3)
hold on
plot(time, store(:,7), 'Linewidth', 2)
legend('true x', 'estimated x','Location','northwest')
title('foot x position')
xlabel('time (s)')
ylabel('x distance (meters)')
ylim([-0.5 2.5])
grid on

subplot(2,2,2)
plot(time, world_leftFoot(:,2), 'Linewidth', 3)
hold on
plot(time, store(:,8), 'Linewidth', 2)
legend('true y', 'estimated y','Location','northwest')
title('foot y position')
xlabel('time (s)')
ylabel('y distance (meters)')
ylim([-0.5 2.5])
grid on

subplot(2,2,3)
plot(time, world_leftFoot(:,3), 'Linewidth', 3)
hold on
plot(time, store(:,9), 'Linewidth', 2)
legend('true z', 'estimated z','Location','northwest')
title('foot z position')
xlabel('time (s)')
ylabel('z distance (meters)')
grid on
% ylim([-1.5 0])
% legend('ture left foot x', 'estimated left foot x', 'ture left foot y', 'estimated left foot y', 'ture left foot z', 'estimated left foot z') 
%% Observability 
% O = H;
% n = 1;
% while n < 12
%     temp = H*A^(n);
%     O = [O;temp];
%     n = n+1;
% end 
%% velocity plot
figure
subplot(2,2,1)
plot(time, trueVel(:,1))
hold on 
plot(time, store(:,4))
title('base velocity')
xlabel('time (s)')
ylabel('x velocity (m/s)')
legend('true','estimated')
ylim([-1 2])
grid on
subplot(2,2,2)
plot(time, trueVel(:,2))
hold on 
plot(time, store(:,5))
title('base velocity')
xlabel('time (s)')
ylabel('y velocity (m/s)')
legend('true','estimated')
ylim([-1 2])
grid on

subplot(2,2,3)
plot(time, trueVel(:,3))
hold on 
plot(time, store(:,6))
title('base velocity')
xlabel('time (s)')
ylabel('z velocity (m/s)')
legend('true','estimated')
grid on
% hold on 
% plot(leftFoot(:,1))


%% function 
function R = quaternion2RotM(q)
%q = w x y z
%     R = [1-2*(q(3)*q(3)+q(4)*q(4)),2*(q(2)*q(3)-q(4)*q(1)),2*(q(2)*q(4)+q(3)*q(1));
%         2*(q(2)*q(3)+q(4)*q(1)),1-2*(q(2)*q(2)+q(4)*q(4)),2*(q(3)*q(4)-q(2)*q(1));
%         2*(q(2)*q(4)-q(3)*q(1)),2*(q(3)*q(4)+q(2)*q(1)),1-2*(q(2)*q(2)+q(3)*q(3))];
    q0 = q(1,1);
    q1 = q(2,1);
    q2 = q(3,1);
    q3 = q(4,1);
    R = [2*(q0^2+q1^2)-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2); 
        2*(q1*q2+q0*q3), 2*(q0^2+q2^2)-1, 2*(q2*q3-q0*q1);
        2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(q0^2+q3^2)-1];
end