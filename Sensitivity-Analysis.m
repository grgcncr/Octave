% Φορτώνουμε τα δεδομένα από το αρχείο results.mat
load('results.mat', 'av_alt_w', 'cost_w', 'usability_w', 'functionality_w', 'final_priorities', 'alt_matrix', 'criteria_vector', 'cost_vector', 'usability_vector', 'functionality_vector', 'cr_w', 'alt_w');

fac = [3 3 3];;
experts = 15;
nfac = sum(fac,2);
% Ορίζουμε τη συνάρτηση pertub_matrix
function [Pnew, CR] = pertub_matrix(P, s)
  [nx, ny] = size(P);
  DP = zeros(nx, ny); % Δημιουργία πίνακα για τις παραμορφώσεις
  Pnew = ones(nx, ny); % Αρχικοποίηση του Pnew

  for i = 1:nx
    for j = i+1:ny
      DP(i,j) = s * (rand - 0.5); % Παραμόρφωση Pij εντός [-s/2, s/2]
      P(i,j) = P(i,j) * DP(i,j) + P(i,j); % Εφαρμογή παραμόρφωσης στον P
      Pnew(i,j) = closeToValues(P(i,j)); % Προσαρμογή του Pij στην κοντινότερη ακέραια τιμή της κλίμακας
    end
    for j = i+1:ny
      Pnew(j,i) = 1 / Pnew(i,j); % Υπολογισμός αντίστροφου Pij
    end
  end

  [wnew, CR] = eigenmethod(Pnew); % Υπολογισμός eigenmethod
end

% Ορίζουμε τη συνάρτηση rankinversion
function N = rankinversion(T)
  n = length(T);
  for i = 1:n-1
     if T(i) >= T(i+1)
        N = 1; % Αναστροφή κατάταξης
        return;
     end
  end
  N = 0; % Χωρίς αναστροφή κατάταξης
end

% Ορίζουμε τη συνάρτηση closeToValues
function Pijnew = closeToValues(Pij)
  values = [1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7, 8, 9];
  nv = length(values);

  for i = 1:(nv-1)
    if Pij > values(i) && Pij < values(i+1)
      if abs(Pij - values(i)) > abs(values(i+1) - Pij)
        Pijnew = values(i+1);
      elseif abs(Pij - values(i)) < abs(values(i+1) - Pij)
        Pijnew = values(i);
      elseif abs(Pij - values(i)) == abs(values(i+1) - Pij)
        Pijnew = values(i);
      end
    elseif Pij == values(i)
      Pijnew = values(i);
    elseif Pij <= values(1)
      Pijnew = values(1);
    elseif Pij >= values(nv)
      Pijnew = values(nv);
    end
  end
end

% Ορίζουμε τη συνάρτηση eigenmethod
function [w, CR] = eigenmethod(P)
  [n, ~] = size(P);
  [V, D] = eig(P);
  [~, ind] = max(diag(D));
  w = V(:, ind);
  w = w / sum(w); % Κανονικοποίηση των βαρών

  % Υπολογισμός δείκτη συνέπειας (Consistency Ratio, CR)
  lambda_max = max(diag(D));
  CI = (lambda_max - n) / (n - 1);
  RI = [0 0 0.58 0.9 1.12 1.24 1.32 1.41 1.45 1.49];
  if n <= length(RI)
    CR = CI / RI(n);
  else
    CR = CI / (1.49 + 0.11 * (n - 10)); % Επέκταση του RI για n > 10
  end
end

% Ορίζουμε τις παραμέτρους της προσομοίωσης Monte Carlo
N = 1e4; % Αριθμός επαναλήψεων
s_values = [0.2 0.1 0.6]; % Τιμές δύναμης παραμόρφωσης

% Αρχικοποιούμε τον πίνακα για την αποθήκευση των PRR
PRR = zeros(length(s_values), 1);

% Εφαρμόζουμε την ανάλυση ευαισθησίας για κάθε τιμή του s
for k = 1:length(s_values)
  s = s_values(k);
  PRR_simulations = zeros(N, 1);

  % Αρχικοποίηση των νέων βαρών
  Cr_w = cr_w; % Αρχικοποίηση βαρών κριτηρίων
  c_w = cost_w; % Αρχικοποίηση βαρών κόστους
  us_w = usability_w; % Αρχικοποίηση βαρών χρηστικότητας
  fu_w = functionality_w; % Αρχικοποίηση βαρών λειτουργικότητας
  aw = alt_w; % Αρχικοποίηση βαρών εναλλακτικών
  ScenarioValue_ntimes = zeros(N, 4); % Αρχικοποίηση προτεραιοτήτων εναλλακτικών για κάθε επανάληψη
  for iter = 1:N
    for i = 1:experts
      % Παραμόρφωση των πινάκων
      Cr_w(i,:) = pertub_matrix(criteria_vector(i), s); % Παραμορφωμένος πίνακας κριτηρίων
      c_w(i,:) = pertub_matrix(cost_vector(i), s); % Παραμορφωμένος πίνακας κόστους
      us_w(i,:) = pertub_matrix(usability_vector(i), s); % Παραμορφωμένος πίνακας χρηστικότητας
      fu_w(i,:) = pertub_matrix(functionality_vector(i), s); % Παραμορφωμένος πίνακας λειτουργικότητας


      for j = 1:nfac
        aw(i,j) = pertub_matrix(alt_matrix(i,j), s); % Παραμορφωμένος πίνακας εναλλακτικών
      end
    end


    % Εκτίμηση μέσων βαρών και σχετικών βαθμολογιών για τους ειδικούς
    ncr_w = mean(Cr_w, 1);
    nc_w = mean(c_w, 1);
    nsu_w = mean(us_w, 1);
    nfu_w = mean(fu_w, 1);
    naw = mean(cell2mat(aw), 1);
    nfw = [nc_w, nsu_w, nfu_w];

    % Εκτίμηση των προτεραιοτήτων των εναλλακτικών
    for i = 1:4 % Εναλλακτικές
      ScenarioValue_ntimes(iter,i) = 0;
      Nfcur = 0;
      for k = 1:3 % Κριτήρια
        Nfcur = Nfcur + fac(k); % Nfcur δείχνει το μέγιστο j για κάθε k
        for j=j+1:Nfcur % Παράγοντες
          ScenarioValue_ntimes(iter,i) = ScenarioValue_ntimes(iter,i) + ncr_w(k) * nfw(j) * naw(i,j);
        end
      end
    end

    NPRR(iter) = rankinversion(ScenarioValue_ntimes(iter,:)); % Υπολογισμός αναστροφής κατάταξης για κάθε επανάληψη
  end

  % Υπολογισμός του μέσου PRR για την τρέχουσα τιμή του s
  PRR(k) = mean(NPRR);
end

% Εμφανίζουμε τα αποτελέσματα
disp('Perturbation Strength (s)');
disp(s_values');
disp('Rank Reversal Probability (PRR)');
disp(PRR);

% Απεικονίζουμε το PRR ως συνάρτηση του s
plot(s_values, PRR, '-o');
xlabel('Perturbation Strength (s)');
ylabel('Rank Reversal Probability (PRR)');
title('Rank Reversal Probability vs Perturbation Strength');

