tokamaks = {'mastu', 'nstxu', 'sparc'};

% inputs-only runs of all pulses
for tokamak_ = tokamaks 
  tokamak = tokamak_{:};
  pulses = get_pulses(tokamak);
  for pulse = pulses    
    fprintf('Inputs only run: %s pulse %d\n', tokamak, pulse)
    run_pulse(tokamak, pulse, 1);
  end
end


% full solution runs of all pulses
for tokamak_ = tokamaks 
  tokamak = tokamak_{:};
  pulses = get_pulses(tokamak);
  for pulse = pulses    
    fprintf('Full-pulse run: %s pulse %d\n', tokamak, pulse)
    run_pulse(tokamak, pulse, 0);
  end
end


function pulses = get_pulses(tokamak)
  d = dir([getenv('GSROOT') '/tokamaks/' tokamak, '/pulses']);
  pulses = [];
  for i = 1:length(d)
    pulse = str2num(d(i).name);
    if ~isempty(pulse)
      pulses(end+1) = pulse;
    end
  end
end
