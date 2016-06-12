x = [2.9 97.1];

figure
subplot(1,2,1);
p = pie(x);
title('% Memory Usage on PC');
p(1).FaceColor = 'r';
p(3).FaceColor = 'k';
