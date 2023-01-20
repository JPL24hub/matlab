% 1. change input to binary
% 2. append 1 
% 3. append zeros until the total number of digits + 64 is equal modulo 512 is zero
% 4. write number of input bits (in 1) in binary with the 64 bits (in 3) and append
% 5. split into 512 blocks
% 6. for each block spit into 32 bit words (16 by 32)
% 7. expand to 64 words by computing the next words (17, 18, ... 64) as next-word = sigma_1(w(t-2)) + sigma_0(w(t-15)) + w(t-7) + w(t-16)
% 8. initialize a-to-h with H (a-to-f-initial),calculat T1 and T2, shift a-to-h down (ie
% a->b, b->c, c->d, ...) then a<-T1, e<-T2. repeat this process for K(1)+W(1)...K(64)+W(64) in T1 to get a-to-f-new. 
% T1=epsilon_1(e)+choice(e,f,g)+h+K(1..64)+W(1..64),
% T2=epsilon_0(a)+majority(a,b,c)
% 9. compute a-to-f-final = a-to-f-initial + a-to-f-new
% 10. repeat 6 to 9 for each block in 6 using a-to-f-final as a-to-f-initial
% 11. a-to-f-final is the digest

function[digest]=sha256(meg)
% Data preprocessing
    DataInBinary = meg %reshape(dec2bin(meg, 8).'-'0',1,[]) %convert to binaty 
    inputMeg = DataInBinary;
    DataInBinary(end+1) = 1 %append 1
%% 

    L = length(DataInBinary)+64;
    i=0; %how many zeros to append
    while( mod(L,512) ~= 0)
        i=i+1;
        L=L+1;
    end
    z = zeros(1,i); % zeros to append
    input = de2bi(length(inputMeg),64,'left-msb'); %size of input (DataInBinary) in binary
    newData = [DataInBinary z input]; %append zeros then size of input

    row=numel(newData)/512;
    k=512;
    j=1;
    for i=1:row
        data(i,:) = newData(1,j:k);
        k=k+512;
        j=j+512;  
    end


    aTof = zeros(8,32); %initialize a to f
    for i=1:8
       aTof(i,:) = H(i);
    end
    %%
    for ii=1:row
        newData = data(ii,:);
%split into 32 bit long ie 16 by 32
        i=1;
        j=32;
        k=1;
        while(j <= 512)
            w(k,:) = newData(1,i:j);
            i=i+32;
            j=j+32;
            k=k+1; 
        end
%%

%expand to 64 words
        for i=17 : 64
            a1=w(i-7,:);
            a2=w(i-16,:);
            b1=sigma_1(w(i-2,:));
            c1=sigma_0(w(i-15,:));
    
            a = binaryAdd(a1,a2);
            b = binaryAdd(a,b1);
            w(i,:)=binaryAdd(b,c1); 
        end
        w

%%
%compression
        aTofN = aTof;
        aTof_prim = aTof;
      
        for i=1:64
            T_1 = binaryAdd(epsilon_1(aTof(5,:)), choice(aTof(5,:),aTof(6,:),aTof(7,:)));% epsilon(e) + choice(e,f,g)+
            T_2 = binaryAdd(T_1, aTof(8,:)); % h + 
            T_3 = binaryAdd(T_2, K(i)); % k(i) + 
            T1 = binaryAdd(T_3, w(i,:)); % w(i)
            T2 = binaryAdd(epsilon_0(aTof(1,:)), majority(aTof(1,:),aTof(2,:),aTof(3,:)));
    
            aTof_prim(2:8,:) = aTof(1:7,:); % shift lower by 1
            aTof_prim(1,:) = binaryAdd(T1, T2);
            aTof_prim(5,:) = binaryAdd(aTof_prim(5,:), T1);
            aTof = aTof_prim;
        end

        aTofDp = zeros(8,32);
        for i=1:8
            aTofDp(i,:) = binaryAdd(aTofN(i,:), aTof(i,:));
        end
        aTof = aTofDp;

    end

    dgst = reshape(aTof', [1,numel(aTof)]);

    digest = binaryVectorToHex(dgst);
end

