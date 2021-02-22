function [gstop,isterminal,direction] = g_fall(t,y)
gstop = y(1);
isterminal = 1;
direction = [];