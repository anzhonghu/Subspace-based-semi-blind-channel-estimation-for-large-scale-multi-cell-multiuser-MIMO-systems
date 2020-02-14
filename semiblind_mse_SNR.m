%Semiblind channel estimation for rapid fading multicell MIMO systems
%This file is for the simulation calculation
clear;
close all;
%%constants
L = 3;%cell number
Nd = 100;%symbol length
K = 3;%user number
tao = K;%pilot length
r = 1000;
rc = r * 0.8;%center to vertex distance(m)
rh = 100;%minimum terminal radius of the cell(m)
ra = rc / rh - 1;
gamma = 3.8;%decay exponent
mu = 0;
sigma = 10^(8*0.1);
i_ant = 100;
Num = 1e1;%iteration
Nt = tao + Nd;
alpha_d = 0.9;
SNR_n = 7;
SNR = [0;5;10;15;20;25;30];
%%position of every base
base(1:7,1) = [0;(1i * 2 * rc);(sqrt(3) * rc + 1i * rc);(sqrt(3) * rc - 1i * rc);(-1i * 2 * rc);(-sqrt(3) * rc - 1i * rc);(-sqrt(3) * rc + 1i * rc);];
D = zeros(K,K*L*L);
H = zeros(i_ant,K*L*L);
G = zeros(i_ant,K*L*L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cova_store = zeros(SNR_n, 5);

%Here the iteration begins
for ii = 1 : SNR_n
    amp = 10 ^ (SNR(ii,1)*0.05) / sqrt(i_ant);
    bound = 0;
    bound1 = 0;
    csi_mse = zeros(1,4);
    for jj = 1 : Num
        shadow_amp = lognrnd(mu,sigma);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%position of every terminal, unifrom distribute
        dis(1:K,1:3) = (rem(rand(K,3) * ra, ra) + 1) * rh;
        ang(1:K,1:3) = rand(K,3) * 2 * pi;
        pos(1:K,1:3) = dis .* (exp(1i * ang));
        pos(:,2) = pos(:,2) + base(2,1);
        pos(:,3) = pos(:,3) + base(3,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %channel matrix, large scale fading and small scale fading
        for l1 = 1 : L%BS
            for l2 = 1 : L%user
                H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = 1 / sqrt(2) * (randn(i_ant,K)+1i*randn(i_ant,K));
%                 H_temp = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K)';
%                 [Hn_temp,~,~] = svd(H_temp);
%                 H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = sqrt(i_ant) * Hn_temp(:,i_ant-K+1:i_ant);
                D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = shadow_amp  * diag(((abs(pos(:,l2)-base(l1,1))*0.01).^(-0.5*gamma)));
                G(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) = H(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K) * D(:,(l1-1)*L*K+(l2-1)*K+1:(l1-1)*L*K+l2*K);
            end
        end
        beta = shadow_amp  * ((rc*0.01).^(-0.5*gamma));
        amp_use = amp / beta;%sqrt(rho),or sqrt(rho_d)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate modulated symbols, every N columns corresponds to symbols in
        %one cell
        symbol = sign(randn(K,L*Nd));%BPSK modulation
        for k = 1 : K
            for l = 1 : L*Nd
                if 0 == symbol(k,l)
                    symbol(k,l) = 1;
                else
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        receive_symbol = zeros(i_ant,Nd*L);
        noise =  (randn(i_ant,Nd*L)+1i*randn(i_ant,Nd*L)) / sqrt(2);
        receive_symbol = receive_symbol + noise;
        for j = 1 : L%BS
            for l = 1 : L
                Gjl = G(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K);
                receive_symbol(:,(j-1)*Nd+1:j*Nd) = receive_symbol(:,(j-1)*Nd+1:j*Nd) + amp_use * Gjl  * symbol(:,(l-1)*Nd+1:l*Nd);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate pilots
        pilots = zeros(tao,K);%Chu sequences
        if 0 == mod(tao,2)
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*ll/tao);
                end
            end
        else
            for k = 1 : K
                for l = 0 : tao-1
                    ll = mod(l+k-1, tao);
                    pilots(l+1,k) = exp(1i*pi*ll*(ll+1)/tao);
                end
            end
        end
        R_p = pilots.' * conj(pilots);%the pilots here are of amplitude 1, so the autocorrelation is not an identity matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %received pilots
        receive_pilots = zeros(i_ant,tao*L);
        noise =  (randn(i_ant,tao*L)+1i*randn(i_ant,tao*L)) / sqrt(2);
        receive_pilots = receive_pilots + noise;
        for j = 1 : L%BS
            for l = 1 : L
                Gjl = G(:,(j-1)*L*K+(l-1)*K+1:(j-1)*L*K+l*K);
                receive_pilots(:,(j-1)*tao+1:j*tao) = receive_pilots(:,(j-1)*tao+1:j*tao) + amp_use * Gjl  * pilots.';
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %evd
        H_e_null = zeros(i_ant,K*L);
        H_e_span = zeros(i_ant,K*L);
        H_e_pilot = zeros(i_ant,K*L);
        cita = zeros(K,L);
        H_e = zeros(i_ant,K*L);
        U = zeros(i_ant,K);
        for j = 1 : L
            r_x = zeros(i_ant,i_ant);
            for k = 1 : Nd
                r_x = r_x + receive_symbol(:,(j-1)*Nd+k) * receive_symbol(:,(j-1)*Nd+k)';
            end
            r_x = r_x / Nd;
            [U_temp,D1,~] = svd(r_x);
            U = U_temp(:,1:K);
            %semiblind estimate
            Hpi = receive_pilots(:,(j-1)*tao+1:j*tao) * conj(pilots) / R_p / amp_use;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H_e_pilot(:,(j-1)*K+1:j*K) = Hpi / D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Uc = U_temp(:,K+1:i_ant);
            R_u = zeros(i_ant,i_ant);
            for l = 1 : i_ant-K
                R_u = R_u + Uc(:,l) * Uc(:,l)';
            end
            for k = 1 : K
                sum_H_e_pilot = zeros(i_ant,1);
                R_temp = R_u;
                for j1 = 1 : K
                    R_temp = R_temp +  (H_e_pilot(:,(j-1)*K+j1) * H_e_pilot(:,(j-1)*K+j1)');
                end
                H_e_null(:,(j-1)*K+k) =  i_ant * (R_temp \ H_e_pilot(:,(j-1)*K+k));
            end
            %H_e_null(:,(j-1)*K+1:j*K) = i_ant * U / (H_e_pilot(:,(j-1)*K+1:j*K)'*U);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            H_e_span(:,(j-1)*K+1:j*K) = U * U' * H_e_pilot(:,(j-1)*K+1:j*K);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %EVD based method
            [U_temp,~] = eig(r_x);
            D_comp = i_ant * diag(D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K).^2);
            [~,index] = sort(D_comp);
            U(:,index) = U_temp(:,i_ant-K+1:i_ant);
            %scalar resolvement
            A_sum = zeros(2*K,2*K);
            Ay = zeros(2*K,1);
            for n = 1 : tao
                yn = [real(receive_pilots(:,(j-1)*tao+n));imag(receive_pilots(:,(j-1)*tao+n))];
                An = U * (D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)) * diag(pilots(n,:)) * amp_use;
                A_bar = [real(An),-imag(An);
                    imag(An),real(An)];
                A_sum = A_sum + A_bar.' * A_bar;
                Ay = Ay + A_bar.' * yn;
            end
            cita_temp = A_sum \ Ay;
            cita(:,j) = cita_temp(1:K,1) + cita_temp(K+1:2*K,1) * 1i;
            H_e(:,(j-1)*K+1:j*K) = U * diag(cita(:,j));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Uc = [(U_temp(:,1:i_ant-K))';zeros(K,i_ant)];
%             [~,~,V] = svd(Uc);
%             Vc = V(:,i_ant-K+1:i_ant);
%             H_e_null(:,(j-1)*K+1:j*K) = Vc * Vc' *  H_e_pilot(:,(j-1)*K+1:j*K);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %calculate error rate
        csi_error = zeros(1,4);
        for j = 1 : L
            sum_e_H = zeros(1,4);
            for k = 1 : i_ant
                sum_e_H(1,1) = sum_e_H(1,1) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e(k,(j-1)*K+1:j*K)) * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e(k,(j-1)*K+1:j*K))';
                sum_e_H(1,2) = sum_e_H(1,2) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_null(k,(j-1)*K+1:j*K)) ...
                    * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_null(k,(j-1)*K+1:j*K))';
                sum_e_H(1,3) = sum_e_H(1,3) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_span(k,(j-1)*K+1:j*K)) ...
                    * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_span(k,(j-1)*K+1:j*K))';
                sum_e_H(1,4) = sum_e_H(1,4) + (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_pilot(k,(j-1)*K+1:j*K)) ...
                    * (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K)-H_e_pilot(k,(j-1)*K+1:j*K))';
            end
            csi_error = csi_error + sum_e_H;
        end
        csi_mse = csi_mse + csi_error / L / K / i_ant;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %cramer-rao bound
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                cova1 = zeros(K*i_ant,K*i_ant);
                R1 = zeros(K*i_ant,i_ant*Nt);
                Q1 = zeros(K*Nd,i_ant*Nt);
                for j = 1 : L
                    Djj = diag(D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K));
                    b = [pilots;symbol(:,(j-1)*Nd+1:j*Nd).'];
                    R_e = (amp_use * b * D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';
                    for k = 1 : i_ant
                        R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) = R_e;
                        Q_e = (H(k,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K) * D(:,(j-1)*L*K+(j-1)*K+1:(j-1)*L*K+j*K))';%amp_use not needed
                        for l = 1 : Nd
                            Q1((l-1)*K+1:l*K,(k-1)*Nt+tao+l) = Q_e;
                        end
                    end
                    noise_va = 0;
                    noise_va = noise_va + 1;
                    Q_QH = Q1*Q1';
                    R_RH = zeros(K*i_ant,K*i_ant);
                    R_QH = zeros(K*i_ant,K*Nd);
                    for k = 1 : i_ant
                        R_QH((k-1)*K+1:k*K,:) = R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) * (Q1(:,(k-1)*Nt+1:k*Nt))';
                        R_RH((k-1)*K+1:k*K,(k-1)*K+1:k*K) = R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt) * (R1((k-1)*K+1:k*K,(k-1)*Nt+1:k*Nt))';
                    end
                    RPR1 = R_RH - R_QH / Q_QH * R_QH';
                    cova1 = cova1 + noise_va * eye(K*i_ant) / RPR1;
                end
                bound1 = bound1 + sum(diag(cova1)) / L / K / i_ant;
    end
    csi_mse = csi_mse / Num;
    %bound = real(bound) / Num;
    bound1 = real(bound1) / Num;
    cova_store(ii,:) = [csi_mse,bound1];
    sprintf('%d',ii)
end



h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
semilogy(SNR,cova_store(:,4),'b-.o','MarkerSize',10)
hold on
semilogy(SNR,cova_store(:,1),'k:*')
semilogy(SNR,cova_store(:,2),'k:s','MarkerSize',10)
semilogy(SNR,cova_store(:,3),'k-d','MarkerSize',10)
semilogy(SNR,cova_store(:,5),'r-+')
le = legend('Pilot-based estimator','Estimator in [7]','Estimator in [8]','Proposed estimator','CRB','Location','Southwest');
set(le,'Fontsize',14,'Fontname','Times')
xlabel('SNR(dB)','Fontsize',16,'Fontname','Times')
ylabel('MSE','Fontsize',20,'Fontname','Times')
%print(h,'-dpdf','evd_mse_SNR')