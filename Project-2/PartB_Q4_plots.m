contour(x,y,phi_num(:,:,curr),50)
title(sprintf('At t=%s',num2str(t)));
saveas(gcf, sprintf('Contour_implicit_t=%s_PartB_Q3.png', num2str(t)))