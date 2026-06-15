s = 1e-3; % m
ds = 50e-6; % m
d_distanza_piastre = 0.1e-3; % m
distanza_piastre = 3e-3; % m
U1 = 620; % V
U2 = 401; % V
dU = 5.5; % V
eta = 1.850e-5; % kg / m*s
rho1 = 1.225; % kg/m^3 densità aria
rho2 = 850; % kg/m^3 densità olio
g = 9.81; % m/s^2
b = 82; % mircometri * hPa
p = 1000; % hPa
A = (b / p) * 1e-6;

% parte simbolica per derivate parziali (per calcolo incertezza)
syms s_sym d_sym U_sym t1_sym t2_sym eta_sym rho1_sym rho2_sym g_sym A_sym

v1_sym = s_sym / (2 * t1_sym);
v2_sym = s_sym / (2 * t2_sym);

r0_sym = sqrt(9/2 * ((eta_sym * v2_sym)/((rho2_sym - rho1_sym) * g_sym)));
q0_sym = r0_sym * (6 * pi * eta_sym * d_sym * (v1_sym + v2_sym)) / U_sym;
r_sym = sqrt(r0_sym^2 + A_sym^2/4) - A_sym/2;
q_sym = q0_sym * (1 + A_sym / r_sym)^(-1.5);

% derivate parziali simboliche
dq_ds = diff(q_sym, s_sym);
dq_dd = diff(q_sym, d_sym);
dq_dU = diff(q_sym, U_sym);

%% U1
% rimossa la prima goccia con t1 = 24.883 e t2 = 13.730
t1_u1 = [2.640 5.785 7.668 6.294 2.389 6.177 3.655 15.582 1.759];
t2_u1 = [6.767 7.945 6.940 4.352 8.291 9.405 13.933 6.066 2.326];

v1 = zeros(1, size(t1_u1, 2));
v2 = zeros(1, size(t1_u1, 2));
r0 = zeros(1, size(v2, 2));
q0 = zeros(1, size(v2, 2));
r = zeros(1, size(v2, 2));
q = zeros(1, size(v2, 2));
delta_q = zeros(1, size(t1_u1, 2));

for i = 1:size(t1_u1, 2)
    v1(i) = s / (2 * t1_u1(i));
    v2(i) = s / (2 * t2_u1(i));

    r0(i) = sqrt(9/2 * ((eta * v2(i))/((rho2 - rho1) * g)));
    q0(i) = r0(i) * (6 * pi * eta * distanza_piastre * (v1(i) + v2(i))) / U1;

    r(i) = sqrt((r0(i))^2 + A^2/4) - A/2;
    q(i) = q0(i) * (1 + A / r(i))^(-1.5);

    valori_correnti = {s, distanza_piastre, U1, t1_u1(i), t2_u1(i), eta, rho1, rho2, g, A};
    variabili_sym   = {s_sym, d_sym, U_sym, t1_sym, t2_sym, eta_sym, rho1_sym, rho2_sym, g_sym, A_sym};

    % sostituzione numerica alle derivate parziali simboliche
    dq_ds_num = double(subs(dq_ds, variabili_sym, valori_correnti));
    dq_dd_num = double(subs(dq_dd, variabili_sym, valori_correnti));
    dq_dU_num = double(subs(dq_dU, variabili_sym, valori_correnti));

    delta_q(i) = sqrt((dq_ds_num * ds)^2 + (dq_dd_num * d_distanza_piastre)^2 + (dq_dU_num * dU)^2);
end

% calcolo di n (max(n) dovrebbe essere 7 da manuale)
[q_ordinato, idx_ordinati] = sort(q, 'descend');
curr = 0;
ratio = 0;
ratios = [];

% j = i+1 così si minimizzano le iterazioni (avendo ordinato q)
for i = 1:size(q_ordinato, 2)
    curr = q_ordinato(i);
    for j = i+1:size(q_ordinato, 2)
        ratio = curr / q_ordinato(j);
        ratios(i, j) = ratio;
    end
end

[n_righe, n_colonne] = size(ratios);
counter = 1;

lista_rapporti = []; % ratios senza zeri (solo triangolo superiore)
coppie_indici = [];

for i = 1:n_righe
    for j = i+1:n_colonne
        lista_rapporti(counter) = ratios(i, j);
        coppie_indici(counter, :) = [i, j];

        counter = counter + 1;
    end
end

passo_tolleranza = 0.05;

% necessario per arrotondamento più preciso
lista_rapporti_regolati = round(lista_rapporti / passo_tolleranza) * passo_tolleranza;

% categorie di rapporti (set da lista_rapporti_regolati)
[valori_unici, ~, idx_gruppo] = unique(lista_rapporti_regolati);

Set_Rapporti = table();
riga_tabella = 1;

% tabella rapporti con relative gocce
for g = 1:length(valori_unici)
    elementi_nel_set = find(idx_gruppo == g);

    if length(elementi_nel_set) >= 2
        for k = 1:length(elementi_nel_set)
            idx_lista = elementi_nel_set(k);

            goccia_a = coppie_indici(idx_lista, 1);
            goccia_b = coppie_indici(idx_lista, 2);
            valore_rapporto_reale = lista_rapporti(idx_lista);

            Set_Rapporti.ID_Gruppo(riga_tabella) = g;
            Set_Rapporti.Rapporto_Arrotondato(riga_tabella) = valori_unici(g);
            Set_Rapporti.Rapporto_Reale(riga_tabella) = valore_rapporto_reale;
            Set_Rapporti.Goccia_Numeratore(riga_tabella) = goccia_a;
            Set_Rapporti.Goccia_Denominatore(riga_tabella) = goccia_b;

            riga_tabella = riga_tabella + 1;
        end
    end
end

Set_Rapporti = sortrows(Set_Rapporti, 'ID_Gruppo');

% calcolo n usando metodo monte carlo (tramite dispersione minima)
% l'idea è di minimizzare la dispersione degli e_test andando a trovare il
% vettore di n più adatto dato che q = n*e => e = q / n
dispersione_min = inf;
best_n = zeros(1, 9); % 9 perchè prima misurazione inaccurata eliminata

% si potrebbe ottimizzare usando q ordinato e scartando n non crescenti (verificare)
% 10 milioni di iterazioni per sicurezza
for i = 1:10000000
    n_causali = randi(7, 1, 9);

    e_test = q ./ n_causali;

    % deviazione standard
    dispersione_attuale = std(e_test);

    % eliminando misura inaccurata forse secondo controllo inutile
    if (dispersione_attuale < dispersione_min) && (mean(e_test) > 1.4e-19) && (mean(e_test) < 1.8e-19)
        dispersione_min = dispersione_attuale;
        best_n = n_causali;
    end
end

% stima di c_e, de
c_e = q ./ best_n;
de = delta_q ./ best_n;

w = 1 ./ (de.^2);

e_finale = sum(w .* c_e) / sum(w);
de_finale = 1 / sqrt(sum(w));

disp(e_finale + "+-" + de_finale);