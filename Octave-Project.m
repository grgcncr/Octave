% Πίνακες Κριτηρίων, Παραγόντων και Εναλλακτικών
criteria = {'Cost', 'Usability', 'Functionality'};
costFactors = {'Price', 'Maintainance', 'Upgrade'};
usabilityFactors = {'Effectiveness', 'Ease of Use', 'System Responsiveness'};
functionalityFactors = {'Interoperability', 'Security', 'Compatibility'};
alternatives = {'software A', 'software B', 'software C', 'software D'};
Factors = [costFactors, usabilityFactors, functionalityFactors];
experts = 15;

% Συνάρτηση για τον υπολογισμό του Consistency Ratio (CR)
function CR = consistency_ratio(matrix)
  [n, ~] = size(matrix);
  [V, D] = eig(matrix);
  max_eigenvalue = max(diag(D));
  CI = (max_eigenvalue - n) / (n - 1);
  RI_values = [0 0 0.58 0.9 1.12 1.24 1.32 1.41 1.45];
  if n <= length(RI_values)
    RI = RI_values(n);
  else
    error('Η τιμή του n είναι μεγαλύτερη από το μήκος του RI');
  end
  CR = CI / RI;
  if CR > 0.1
    disp('Ασυνεπής πίνακας');
  end
end


% Συνάρτηση δημιουργίας πινάκων των ειδηκών
function expert_matrix = expert_matrix_factory(s)
  expert_matrix = zeros(s, s);
  nan_count = 0;
  for i = 1:s
    for j = 1:s
      if (nan_count < 2) && (randi([1,100]) <= 3)
        expert_matrix(i,j) = NaN;
        nan_count++;
        else
        if i == j
          expert_matrix(i, j) = 1;
        elseif expert_matrix(j, i) == 0
          if randi([1, 2]) == 1
            expert_matrix(i, j) = randi([1, 9]);
          else
            expert_matrix(i, j) = 1 / randi([1, 9]);
          end
        else
          expert_matrix(i, j) = 1 / expert_matrix(j, i);
        end
      end
    end
  end
end


% Συνάρτηση για την συμπλήρωση των ελλιπών τιμών στον πίνακα
function filledMatrix = fill_missing_values(matrix)
  [nRows, nCols] = size(matrix);
  for i = 1:nRows
    for j = 1:nCols
      if isnan(matrix(i, j))
        if i == j
          matrix(i, j) = 1;
        else
          if isnan(matrix(j, i))
            matrix(i, j) = randi([1, 9]);
            matrix(j, i) = 1 / matrix(i, j);
          else
            matrix(i, j) = 1 / matrix(j, i);
          end
        end
      end
    end
  end
  filledMatrix = matrix;
end



% PWC συνάρτηση
function vector = PWC(matrix,experts)
  disp('PWC begin');
  s = size(matrix);
  if experts ~= 1
    vector = cell(1, experts);
    for i = 1:experts
      consistent = false;
      attempts = 0;
      while ~consistent && attempts < 100
        vector{i} = fill_missing_values(expert_matrix_factory(s(2)));
        %disp(consistency_ratio(cost_vector{i}));
        if consistency_ratio(vector{i}) <= 0.1
        consistent = true;
        end
         attempts = attempts + 1;
      end
      if ~consistent
        disp(['Factors :' ,num2str(i),' and ', num2str(attempts), ' attempts']);
      end
    end
  else
    vector = cell(s(2),s(2));
    consistent = false;
    attempts = 0;
    while ~consistent && attempts < 100
      vector = fill_missing_values(expert_matrix_factory(s(2)));
      %disp(consistency_ratio(cost_vector{i}));
      if consistency_ratio(vector) <= 0.1
       consistent = true;
      end
      attempts = attempts + 1;
      if ~consistent
        disp([num2str(attempts), ' attempts']);
      end
    end
  end
end


%Δημιουργία όλων των πινάκων των ειδηκών
criteria_vector = PWC(criteria,experts);
cost_vector = PWC(costFactors,experts);
usability_vector = PWC(usabilityFactors,experts);
functionality_vector = PWC(functionalityFactors,experts);
S = size(Factors);
alt_matrix = cell(experts, S(2));  % Ορίζουμε τον πίνακα alt_matrix με το σωστό μέγεθος

for i = 1:S(2)
  for j = 1:experts
      alt_matrix{j, i} = PWC(alternatives, 1);  % Δημιουργία ενός πίνακα εναλλακτικών για τον εκάστοτε ειδικό
  end
end


% Εμφάνιση των πινάκων
for i = 1:experts
  disp(['Criteria Matrix: ', num2str(i)]);
  disp(criteria_vector{i});
end

for i = 1:experts
  disp(['Cost Factors Matrix: ', num2str(i)]);
  disp(cost_vector{i});
end

for i = 1:experts
  disp(['Usability Factors Matrix: ', num2str(i)]);
  disp(usability_vector{i});
end

for i = 1:experts
  disp(['Functionality Factors Matrix: ', num2str(i)]);
  disp(functionality_vector{i});
end

for i = 1:S(2)
  disp(['  Alternatives Matrix for Factor ', Factors{i}, ':']);
  for j = 1:experts
    disp(['Expert ', num2str(j), ':']);
    disp(alt_matrix{j, i});
  end
end



function weights = compute_weights(matrix, experts)
  % Αρχικοποίηση πίνακα weights με μέγεθος (experts x n), όπου n είναι ο αριθμός των κριτηρίων
  [n, ~] = size(matrix{1,1});
  weights = zeros(experts, n);

  for i = 1:experts
    % Υπολογισμός ιδιοτιμών και ιδιοδιανυσμάτων
    [V, D] = eig(matrix{i});

    % Εύρεση του μέγιστου ιδιοτιμήματος και του αντίστοιχου ιδιοδιανύσματος
    [max_eigenvalue, idx] = max(diag(D));
    max_eigenvector = V(:, idx);

    % Κανονικοποίηση του ιδιοδιανύσματος ώστε να προκύψουν τα βάρη
    sum_eigenvector = sum(max_eigenvector);
    weights(i, :) = max_eigenvector / sum_eigenvector;
  end

  % Εμφάνιση των βαρών
  %disp('Weights:');
  %for i = 1:experts
   % disp(['Expert ', num2str(i), ': ', mat2str(weights(i, :))]);
  %end
end



% Συνάρτηση υπολογισμού μέσου βάρους
function av_weight = average_w(matrix, experts)
  % Υπολογισμός των βαρών για κάθε ειδικό
  weights = compute_weights(matrix, experts);

  % Υπολογισμός του μέσου όρου των βαρών
  av_weight = mean(weights, 1); % Μέσος όρος κατά τη διάσταση των ειδικών

  % Εμφάνιση του μέσου βάρους
  disp('Average Weights:');
  disp(av_weight);
end

%Υπολογισμός μέσων βαρών
av_cr_w = average_w(criteria_vector,experts);

cost_w = average_w(cost_vector,experts);
usability_w = average_w(usability_vector,experts);
functionality_w = average_w(functionality_vector,experts);

av_fac_w = (cost_w + usability_w + functionality_w)/3;
disp('Average Factor weight');
disp(av_fac_w);
alt_w = cell(experts,S(2));
for l = 1:experts
  for p = 1:S(2)
    alt_w(l,p) = compute_weights(alt_matrix(l,p),1);
  end
end

av_alt_w = cell(S(2),1);
disp(size(alt_matrix{1,1}));
for i = 1:S(2)
  disp(['Factor: ', Factors{i}]);
  av_alt_w{i} = average_w(alt_matrix(:,i), experts);
end

% Υπολογισμός τελικών προτεραιοτήτων
final_priorities = zeros(1, 4);
for i = 1:4
  for j = 1:3  % Για κάθε παράγοντα
    final_priorities(i) = final_priorities(i) + av_cr_w(j) * av_fac_w(j) * av_alt_w{j}(i);
  end
end

% Κανονικοποίηση των τελικών προτεραιοτήτων
disp('Final Priorities');
disp(final_priorities = final_priorities / sum(final_priorities));


maut_matrix = [65,80,50,73; 90,55,65,60; 65,100,60,75];


% Συνάρτηση χρησιμότητας (γραμμική κανονικοποίηση)
function utility = U(value, min_val, max_val)
  utility = (value - min_val) / (max_val - min_val);
end

% Συνάρτηση για τον υπολογισμό των τελικών χρησιμοτήτων με τη μέθοδο MAUT
function final_utilities = maut_utilities(alternatives_matrix, weights)
  % alternatives_matrix: πίνακας με τις εναλλακτικές (n x m)
  % weights: πίνακας με τα βάρη των κριτηρίων (1 x n)

  [n, m] = size(alternatives_matrix); % Διαστάσεις του πίνακα εναλλακτικών
  final_utilities = zeros(1, m); % Αρχικοποίηση πίνακα για τις τελικές χρησιμότητες

  % Βρίσκουμε τις ελάχιστες και μέγιστες τιμές για κάθε κριτήριο
  min_vals = min(alternatives_matrix, [], 2);
  max_vals = max(alternatives_matrix, [], 2);

  % Υπολογισμός χρησιμοτήτων για κάθε εναλλακτική
  for j = 1:m
    utility_sum = 0;
    for i = 1:n
      % Υπολογισμός χρησιμότητας για το κριτήριο i χρησιμοποιώντας τη συνάρτηση U
      utility_value = U(alternatives_matrix(i, j), min_vals(i), max_vals(i));

      % Προσθήκη της σταθμισμένης χρησιμότητας στο άθροισμα
      utility_sum = utility_sum + weights(i) * utility_value;
    end
    final_utilities(j) = utility_sum;
  end

  % Εμφάνιση των τελικών χρησιμοτήτων
  disp('Final Utilities:');
  disp(final_utilities);
end


% Υπολογίζουμε τις τελικές χρησιμότητες
final_utilities = maut_utilities(maut_matrix, av_cr_w);
cr_w = compute_weights(criteria_vector, experts);
disp(cr_w(1,:));
% Αποθήκευση των αποτελεσμάτων σε αρχείο .mat
save('results.mat','av_alt_w','cost_w','usability_w','functionality_w','final_priorities','alt_matrix','criteria_vector','cost_vector','usability_vector','functionality_vector','cr_w','alt_w');

