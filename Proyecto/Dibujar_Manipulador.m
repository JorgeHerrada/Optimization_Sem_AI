function Dibujar_Manipulador (q)
    a_1 = 0.35;
    a_2 = 0.35;
    a_3 = 0.25;
    
    theta_1 = q(1);
    theta_2 = q(2);
    theta_3 = q(3);
    
    t01 = [a_1*cos(theta_1); a_1*sin(theta_1)];
	
    t02 = [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1); ...
           a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1)];
    
    t03 = [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
           a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];

    hold on 
    grid on

    line([0 t01(1)],[0 t01(2)],'color',[0 0 0],'LineWidth',4)
    line([t01(1) t02(1)],[t01(2) t02(2)],'color',[0 0 0],'LineWidth',4)
    line([t02(1) t03(1)],[t02(2) t03(2)],'color',[0 0 0],'LineWidth',4)

    plot(0,0,'ko','MarkerSize',15,'LineWidth',2)
    plot(t01(1),t01(2),'ko','MarkerSize',15,'LineWidth',2)
    plot(t02(1),t02(2),'ko','MarkerSize',15,'LineWidth',2)
    
    axis([-1 1 -1 1])
