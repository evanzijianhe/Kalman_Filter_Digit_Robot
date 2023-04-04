function Qk = Qk_val
    Qrx = 0.001;
    Qry = 0.001;
    Qrz = 0.001;
    Qr = blkdiag(Qrx, Qry, Qrz);
        
    Qvx = 0.01;
    Qvy = 0.01;
    Qvz = 0.01;
    Qv = blkdiag(Qvx, Qvy, Qvz);
    
%     Qqw = 0.01;
%     Qqx = 0.01;
%     Qqy = 0.01;
%     Qqz = 0.01;
%     Qq = blkdiag(Qqw,Qqx, Qqy, Qqz);
    
    Qfoot = 0.01*eye(6);
    Qk = blkdiag(Qr, Qv,Qfoot); 
    
end