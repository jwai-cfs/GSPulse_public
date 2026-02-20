function ys = measure_ys(psizr, ic, tok, t, optimization_signals)
% =========================================================================
% Description:
%
%
%
% Inputs:
%
%
%
% Outputs:
%
%
%
% Additional info:
%
%
% =========================================================================

cv = define_data_indices(optimization_signals, tok);
N = length(t);
ys = cell(N,1);

for i = 1:N

  ydata = struct;

  for j = 1:length(optimization_signals)

    sig = optimization_signals{j};

    if ismember(sig.name, cv.y_names)

      switch sig.calc_type
      
        case 'flux_relative'
          
          ref = structts2struct(sig, {'r1', 'z1', 'r2', 'z2'}, t(i));
          ydata.(sig.name) = bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r1, ref.z1) - ...
            bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r2, ref.z2);
        
        case 'flux_relative_multipoint'
          ref = structts2struct(sig, {'r', 'z', 'r_multiref', 'z_multiref', 'wt_multiref'}, t(i));
          ydata.(sig.name) = bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r, ref.z) - ...
            ref.wt_multiref * bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r_multiref, ref.z_multiref);

        case 'flux_absolute'
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          ydata.(sig.name) = bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r, ref.z);
      
        case 'flux_absolute_avg'
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          ydata.(sig.name) = mean(bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r, ref.z));
        
        case 'vacuum_flux_absolute_avg'          
          ydata.(sig.name) = 0;

        case 'field_absolute_vertical'
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          [~, psi_r] = bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r, ref.z);
          ydata.(sig.name) = sparse(diag(1 ./ (2 * pi * ref.r(:)))) * psi_r;
          
        case {'vacuum_field_absolute_vertical', 'vacuum_field_absolute_radial'}
          ref = structts2struct(sig, {'r'}, t(i));
          ydata.(sig.name) = zeros(size(ref.r(:)));

        case 'field_absolute_radial'
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          [~, ~, psi_z] = bicubicHermite(tok.rg, tok.zg, psizr(:,i), ref.r, ref.z);
          ydata.(sig.name) = -sparse(diag(1 ./ (2 * pi * ref.r(:)))) * psi_z;
      
        case 'coil_currents'
          ydata.(sig.name) = ic(:,i);  
        
        case 'coil_current_combinations'
          ydata.(sig.name) = sig.coil_combinations_matrix * ic(:,i);
  
        otherwise
          error('Calculation type "%s" not yet implemented.', sig.calc_type)
      end
      
      if any(isnan(ydata.(sig.name)(:)))
        error('Measurement of signal "%s" produced nan', sig.name)
      end
    end
  end

  ys{i} = struct2vec(ydata, cv.y_names);
end
