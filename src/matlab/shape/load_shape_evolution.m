function shapes = load_shape_evolution(t, fps)

% evo = jsondecode(fileread(shape_evolution_fp));
% t = str2num(evo.times);
% fps = fullfile(getenv('GSROOT'), evo.shape_dir, splitlines(evo.shape_fns));
% fps = strrep(fps, ',', '');  % remove commas

shapes = struct;

for i = 1:length(t)
  shape = jsondecode(fileread(fps{i}));

  % read boundary
  shapes.rbbbs.Data(i,:) = shape.rbbbs;
  shapes.zbbbs.Data(i,:) = shape.zbbbs;
  shapes.cp_r.Data(i,:) = shape.cp_r;
  shapes.cp_z.Data(i,:) = shape.cp_z;
  
  % read x-points
  for j = 1:4
    shapes.(['rx' num2str(j)]).Data(i) = shape.manual_pt_inputs_r.(['x_pt' num2str(j)]);
    shapes.(['zx' num2str(j)]).Data(i) = shape.manual_pt_inputs_z.(['x_pt' num2str(j)]);
  end
  
  % read strike points
  for j = 1:8
    shapes.(['rstrike' num2str(j)]).Data(i) = shape.manual_pt_inputs_r.(['strike_pt' num2str(j)]);
    shapes.(['zstrike' num2str(j)]).Data(i) = shape.manual_pt_inputs_z.(['strike_pt' num2str(j)]);
  end
  
  % read additional ref pts
  for j = 1:3
    shapes.(['r_control_pt_ref' num2str(j)]).Data(i) = shape.manual_pt_inputs_r.(['control_pt_ref' num2str(j)]);
    shapes.(['z_control_pt_ref' num2str(j)]).Data(i) = shape.manual_pt_inputs_z.(['control_pt_ref' num2str(j)]);
  end

end

sigs = fieldnames(shapes);
for i = 1:length(sigs)
  shapes.(sigs{i}).Time = t;
end

shapes = check_structts_dims(shapes);