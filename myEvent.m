function [value, isterminal, direction] = myEvent(t, y, maxy)
value = maxy - y(1);
isterminal = 1;   % Stop the integration
direction  = 0;