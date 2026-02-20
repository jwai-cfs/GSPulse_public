function cmats = output_model(dpsizrdx, tok, t, optimization_signals)
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


N = length(t);
cmats = cell(N,1);
cv = define_data_indices(optimization_signals, tok);

for i = 1:N

  cdata = struct;

  for j = 1:length(optimization_signals)

    sig = optimization_signals{j};

    if ismember(sig.name, cv.y_names)
      
      switch sig.calc_type
      
        case 'flux_relative'
          
          ref = structts2struct(sig, {'r1', 'z1', 'r2', 'z2'}, t(i));
          cdata.(sig.name) = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r1, ref.z1) - ...
            multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r2, ref.z2);
        
        case 'flux_relative_multipoint'
          ref = structts2struct(sig, {'r', 'z', 'r_multiref', 'z_multiref', 'wt_multiref'}, t(i));
          cdata.(sig.name) = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r, ref.z) - ...
            ref.wt_multiref * multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r_multiref, ref.z_multiref);
          
        case 'flux_absolute'
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          cdata.(sig.name) = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r, ref.z);
      
        case {'flux_absolute_avg', 'vacuum_flux_absolute_avg'}
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          cdata.(sig.name) = mean(multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r, ref.z));
      
        case {'field_absolute_vertical', 'vacuum_field_absolute_vertical'}
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          [~, psi_r] = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r, ref.z);
          cdata.(sig.name) = sparse(diag(1 ./ (2 * pi * ref.r(:)))) * psi_r;
          
        case {'field_absolute_radial', 'vacuum_field_absolute_radial'}
          ref = structts2struct(sig, {'r', 'z'}, t(i));
          [~, ~, psi_z] = multigrid2pt(tok.rg, tok.zg, dpsizrdx, ref.r, ref.z);
          cdata.(sig.name) = -sparse(diag(1 ./ (2 * pi * ref.r(:)))) * psi_z;
      
        case 'coil_currents'
          cdata.(sig.name) = [eye(tok.nc) zeros(tok.nc, tok.nv)];               
        
        case 'coil_current_combinations'
          cdata.(sig.name) = [sig.coil_combinations_matrix zeros(size(sig.coil_combinations_matrix,1), tok.nv)];

        otherwise
          error('Calculation type "%s" not yet implemented.', sig.calc_type)
      end
      
      if any(isnan(cdata.(sig.name)(:)))
        error('Output response model contains nan for "%s" response', sig.name)
      end
    end


  end

  cmats{i} = struct2vec(cdata, cv.y_names);
end
