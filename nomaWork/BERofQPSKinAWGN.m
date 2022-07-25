
function[snrdb,ber]=BERofQPSKinAWGN()
% clear all;                                            
% close all;                                            
l=10000;
snrdb=1:1:40;
snrlin=10.^(snrdb/10);
for snrdb=1:1:40
    si=2*(round(rand(1,l))-0.5);                      
    sq=2*(round(rand(1,l))-0.5);                                    
    s=si+1i*sq;                               
    w=awgn(s,snrdb,'measured');
    r=w;                                           
    si_=sign(real(r));                                
    sq_=sign(imag(r));                               
    ber1=(l-sum(si==si_))/l;                          
    ber2=(l-sum(sq==sq_))/l;                          
    ber(snrdb)=mean([ber1 ber2]);                         
end
%semilogy(snrdb, ber,'o-')
snrdb=1:1:15;
snrlin=10.^(snrdb./10);
tber=0.5.*erfc(sqrt(snrlin));
% semilogy(snrdb,ber,'-bo',snrdb,tber,'-mh')
% title('QPSK with awgn');
% xlabel('Signal to noise ratio');
% ylabel('Bit error rate');      
% grid on;
end