function [val,isterminal,direction] = g_pendulum(t,y,y0)

val = norm(y-y0)-0.05;
isterminal = 1;
direction = -1;