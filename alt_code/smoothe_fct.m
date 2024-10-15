% Smoothes probabilities from heuristic center detection method
function probs_final = smoothe_fct(probs, ~)
    probs_lefts = zeros(size(probs));
    probs_rights = probs_lefts;
    probs_lefts(1) = probs(1);
    for i = 2:length(probs) % To the left, to the left
        probs_lefts(i) = probs_lefts(i-1)*0.99+probs(i);
    end
    probs_rights(end) = probs(end);
    for i = length(probs)-1:-1:1 % To the right, to the right
        probs_rights(i) = probs_rights(i+1)*0.99+probs(i);
    end
    probs_final = probs_lefts + probs_rights-probs; % We need to subtract 'probs', since we counted the ones of 'probs' twice
    probs_final(probs_final<eps) = 0;
end
