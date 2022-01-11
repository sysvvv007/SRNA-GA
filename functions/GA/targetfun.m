function y = targetfun(IP,QP,SimColle)
y= sum(abs(SimColle.IP-IP))+sum(abs(SimColle.QP-QP));
end